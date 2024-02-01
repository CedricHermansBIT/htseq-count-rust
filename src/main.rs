use bam::record::cigar::Operation;
use bam::record::tags::TagValue;
use bam::BamReader;
use feature::Feature;
use intervaltree::IntervalTree;
use interval::Interval;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use structopt::StructOpt;

mod feature;
mod intervaltree;
mod interval;
mod node;

use node::Node;

fn main() {
    // Command line arguments
    let args = Args::from_args();

    if args.counts_output.is_some() {
        // check if we have write access to the file, otherwise, crash at the start instead of waiting until the end
        if let Err(e) = File::create(args.counts_output.clone().unwrap()) {
            eprintln!("Could not create file: {}", e);
            std::process::exit(1);
        }
    }

    // Try to open the bam file, if it fails, print an error message
    let bam = BamReader::from_path(args.bam.clone(), args.n).expect("Could not read bam file");
    let header = bam.header().clone();
    let reference_names: Vec<String> = header.reference_names().to_owned();
    // bam fields: https://docs.rs/bam/0.3.0/bam/record/struct.Record.html

    // Read the gtf file
    let gtf = read_gtf(&args.gtf, args.t.as_str());

    eprintln!("{:?}",gtf["1"].overlap(24646121, 24646123));
   
    if args.export_feature_map.is_some() {
        let mut file = File::create(args.export_feature_map.clone().unwrap()).expect("Unable to create file");
        for (chr, tree) in gtf.iter() {
            file.write_all(format!("{}:\n", chr).as_bytes()).expect("Unable to write data");
            let _ = tree.top_node.clone().unwrap().write_structure(&mut file, 0);
        }
    }

    // exit(1) to prevent the rest of the program from running
    std::process::exit(1);

    let mut counts = prepare_count_hashmap(&gtf);

    let mut counter = 0;

    count_reads(bam, &mut counter, &mut counts, &args, reference_names, gtf);

    if args.counts_output.is_some() {
        write_counts(counts, args, counter);
    } else {
        print_output(counts, args, counter);
    }
}


// Use the structopt crate to parse command line arguments
#[derive(StructOpt)]
struct Args {
    // Number of threads
    #[structopt(
        short = "n",
        long = "number",
        default_value = "10",
        help = "Number of threads"
    )]
    n: u16,

    // Mode
    // TODO: implement other modes
    #[structopt(short = "m", long = "mode", default_value = "intersection-strict", possible_values = &["intersection-strict", "intersection-nonempty", "union"], help = "Mode to use for counting reads overlapping features. Possible values: intersection-strict, intersection-nonempty, union (default: intersection-strict).")]
    _m: String,

    // Stranded
    // TODO: implement stranded mode
    #[structopt(long = "stranded")]
    _stranded: bool,

    // Quality filter
    #[structopt(
        short = "a",
        long = "amount",
        default_value = "10",
        help = "Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads."
    )]
    a: u8,

    // Type of feature to be used
    #[structopt(
        short = "t",
        long = "type",
        default_value = "exon",
        help = "Feature type (3rd column in GTF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)"
    )]
    t: String,

    // Feature ID
    // TODO: implement actual logic for this option
    #[structopt(
        short = "i",
        long = "id",
        default_value = "gene_name",
        help = "GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). All feature of the right type (see -t option) within the same GTF
    attribute will be added together. The typical way of using this option is to count all exonic reads from each gene and add the exons but other uses are possible
    as well. You can call this option multiple times: in that case, the combination of all attributes separated by colons (:) will be used as a unique identifier,
    e.g. for exons you might use -i gene_id -i exon_number."
    )]
    _i: Vec<String>,

    // Name and type of the bam file
    #[structopt(name = "bam", default_value = "test.bam")]
    bam: String,

    // Name and type of the gtf file
    #[structopt(name = "gtf", default_value = "test.gtf")]
    gtf: String,

    // Secondary alignment mode
    #[structopt(long = "secondary-alignments", default_value = "score", possible_values = &["score", "ignore"], help = "Whether to score secondary alignments (0x100 flag)")]
    secondary_alignments: String,

    // Supplementary alignment mode
    #[structopt(long = "supplementary-alignments", default_value = "score", possible_values = &["score", "ignore"], help = "Whether to score supplementary alignments (0x800 flag)")]
    supplementary_alignments: String,

    // Option to also output a total count of uniquely mapped reads
    #[structopt(
        long = "extended-output",
        help = "Also output a total count of uniquely mapped reads"
    )]
    counts: bool,

    // Output delimiter
    #[structopt(
        short = "d",
        long = "delimiter",
        default_value = "\t",
        help = "Column delimiter in output (default: TAB)."
    )]
    delimiter: String,

    // Output file
    #[structopt(
        short = "c",
        long = "counts_output",
        help = "Filename to output the counts to instead of stdout."
    )]
    counts_output: Option<String>,

    // Export feature map
    #[structopt(
        short = "f",
        long = "export_feature_map",
        help = "Filename to output the feature map for debugging purposes."
    )]
    export_feature_map: Option<String>,
}

fn prepare_count_hashmap(gtf: &HashMap<String, IntervalTree>) -> HashMap<String, i32> {
    let mut counts = HashMap::with_capacity(gtf.len());
    // add all features to the map
    for features in gtf.values() {
        for feature in features.all_intervals.iter() {
            counts.entry(feature.data.as_ref().unwrap().name.clone()).or_insert(0);
        }
    }

    // Add the special keys
    counts.insert("__no_feature".to_string(), 0);
    counts.insert("__ambiguous".to_string(), 0);
    counts.insert("__not_aligned".to_string(), 0);
    counts.insert("__too_low_aQual".to_string(), 0);
    counts.insert("__alignment_not_unique".to_string(), 0);
    counts
}

fn read_gtf(file_path: &str, feature_type_filter: &str) -> HashMap<String, IntervalTree> {
    let mut map: HashMap<String, Vec<Interval>> = HashMap::new();
    let file = File::open(file_path).expect("Kan het bestand niet openen");
    let reader = BufReader::new(file);
    let mut counter = 0;
    for line in reader.lines() {
        counter += 1;
        if counter % 100000 == 0 {
            eprintln!("{} GFF lines processed.", counter);
        }
        if let Ok(line) = line {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let feature_name = fields[2];
            if feature_name != feature_type_filter {
                continue;
            }
            let attributes = fields[8];

            let attributes: Vec<&str> = attributes.split(';').collect();

            let name = attributes
                .iter()
                .find(|&attr| attr.contains("gene_name"))
                .unwrap_or(&"")
                .trim();

            let name = name.replace("gene_name ", "");
            let name = name.trim_matches('"');

            let start = fields[3].parse::<u32>().unwrap_or(0);
            let end = fields[4].parse::<u32>().unwrap_or(0);
            let chr = fields[0].to_string();

            let feature = Feature {
                name: name.to_string(),
                chr: chr.clone(),
                start: start-1,
                end: end,
                start_sorted_index: 0,
                end_sorted_index: 0
            };

            map.entry(chr).or_default().push(Interval::new((start-1) as i32, end as i32, Some(feature)));
        }
    }
    //let node = Node::from_intervals(intervals);
    // print the center node (withouth the children)
    // construct the tree
    // let tree= IntervalTree::from_tuples(tuples);
    // let top_node = tree.top_node.unwrap();
    // eprintln!("tree: center: {}, depth: {}, balance: {}", top_node.x_center, top_node.depth, top_node.balance);
    // eprintln!("Boundary table: {:?}", tree.boundary_table);

    let mut result: HashMap<String, IntervalTree> = HashMap::new();

    for (chr, intervals) in map {
        let tree = IntervalTree::new(Some(intervals));
        result.insert(chr, tree);
        
    }

    eprintln!("{} GFF lines processed.", counter);
    
    result
}

fn should_skip_record(
    record: &bam::Record,
    counts: &mut HashMap<String, i32>,
    args: &Args,
) -> bool {
    // Skip all reads that are not aligned
    if record.ref_id() < 0 {
        *counts.entry("__not_aligned".to_string()).or_insert(0) += 1;
        return true;
    }
    // Skip all reads that are secondary alignments
    if args.secondary_alignments == "ignore" && record.flag().all_bits(0x100) {
        return true;
    }
    // Skip all reads that are supplementary alignments
    if args.supplementary_alignments == "ignore" && record.flag().all_bits(0x800) {
        return true;
    }
    // Skip all reads with MAPQ alignment quality lower than the given minimum value
    if record.mapq() < args.a {
        *counts.entry("__too_low_aQual".to_string()).or_insert(0) += 1;
        return true;
    }
    // Skip all reads that have an optional field "NH" with value > 1
    if let Some(TagValue::Int(i, _)) = record.tags().get(b"NH") {
        if i > 1 {
            *counts
                .entry("__alignment_not_unique".to_string())
                .or_insert(0) += 1;
            return true;
        }
    }
    false
}

fn print_output(counts: HashMap<String, i32>, args: Args, counter: i32) {
    // Print de HashMap
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort_by(|a, b| a.cmp(&b));
    for key in sorted_keys {
        if key.starts_with("__") {
            continue;
        }
        println!("{}{}{}", key, args.delimiter, counts[key]);
    }
    println!("__no_feature{}{}", args.delimiter, counts["__no_feature"]);
    println!("__ambiguous{}{}", args.delimiter, counts["__ambiguous"]);
    println!("__too_low_aQual{}{}",args.delimiter, counts["__too_low_aQual"]);
    println!("__not_aligned{}{}", args.delimiter, counts["__not_aligned"]);
    println!("__alignment_not_unique{}{}",args.delimiter, counts["__alignment_not_unique"]);

    if args.counts {
        println!(
            "Total number of uniquely mapped reads{}{}",
            args.delimiter,
            counter
                - counts["__not_aligned"]
                - counts["__too_low_aQual"]
                - counts["__alignment_not_unique"]
                - counts["__ambiguous"]
                - counts["__no_feature"]
        );
    }
}

fn write_counts(counts: HashMap<String, i32>, args: Args, counter: i32) {
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort_by(|a, b| a.cmp(&b));
    let mut file = File::create(args.counts_output.unwrap()).expect("Unable to create file");
    for key in sorted_keys {
        if key.starts_with("__") {
            continue;
        }
        file.write_all(format!("{}{}{}\n", key, args.delimiter, counts[key]).as_bytes()).expect("Unable to write data");
    }
    file.write_all(format!("__no_feature{}{}\n", args.delimiter, counts["__no_feature"]).as_bytes(),).expect("Unable to write data");
    file.write_all(format!("__ambiguous{}{}\n", args.delimiter, counts["__ambiguous"]).as_bytes()).expect("Unable to write data");
    file.write_all(format!("__too_low_aQual{}{}\n",args.delimiter, counts["__too_low_aQual"]).as_bytes(),).expect("Unable to write data");
    file.write_all(format!("__not_aligned{}{}\n",args.delimiter, counts["__not_aligned"]).as_bytes(),).expect("Unable to write data");
    file.write_all(format!("__alignment_not_unique{}{}\n",args.delimiter, counts["__alignment_not_unique"]).as_bytes(),).expect("Unable to write data");

    if args.counts {
        file.write_all(format!("Total number of uniquely mapped reads{}{}\n",args.delimiter,
            counter - counts["__not_aligned"] - counts["__too_low_aQual"] - counts["__alignment_not_unique"] - counts["__ambiguous"] - counts["__no_feature"])
            .as_bytes(),
        )
        .expect("Unable to write data");
    }
}


fn count_reads(bam: BamReader<File>, counter: &mut i32, counts: &mut HashMap<String, i32>, args: &Args, reference_names: Vec<String>, gtf_end_sorted: HashMap<String, IntervalTree>) {
    for record in bam {
        *counter += 1;
        if *counter % 100000 == 0 {
            //println!("{} records processed.", counter);
            eprintln!("{} records processed.", counter);
        }
        let record = record.unwrap();
        if String::from_utf8_lossy(record.name())=="SRR001432.281211 USI-EAS21_0008_3445:8:7:657:535 length=25" {
            //println!("{}: {}-{}", String::from_utf8_lossy(record.name()), record.start(), record.calculate_end());
            eprintln!("this one");
        }
        if should_skip_record(&record, counts, args) {
            continue;
        }
        // todo
        unimplemented!("todo");
        
        let ref_id = record.ref_id();

        let reference = &reference_names[ref_id as usize];
        if gtf_end_sorted.contains_key(reference) {
            let start_pos: u32 = record.start().try_into().unwrap();
            let end_pos: u32 = record.calculate_end().try_into().unwrap();
            let features = &gtf_end_sorted[reference];

            let mut ambiguous = false;

            //if String::from_utf8_lossy(record.name())=="SRR001432.153301 USI-EAS21_0008_3445:8:4:393:300 length=25" {
                //println!("{}: {}-{}", String::from_utf8_lossy(record.name()), record.start(), record.calculate_end());
                //eprintln!("startindex: {}, endindex: {}", startindex, endindex);
                //eprintln!("feature: {:?}", features.0[startindex]);
                //std::process::exit(1);
            //}
            
            // case 1: read is fully within feature -> count feature
            //   RRR
            //AAAAAAAAA
            // case 2: read is partially within feature -> don't count feature
            //       RRRRR
            //AAAAAAAAA
            // case 3: read overspans intron -> don't count feature
            //  RRRRRR
            //AAAAIIAAAA
            // case 4: read is split over two exons (of same gene) -> count feature
            //  RR--RR
            //AAAAIIAAAA
            // case 5: read is overlaps with two genes, but only one gene covers the entire read -> count feature for that gene
            //  RRR
            //AAAAAAA
            //    BBBBBBB
            // case 6: read is overlaps with two genes, but both genes cover the entire read -> ambiguous
            //  RRR
            //AAAAAAA
            //BBBBBBB
            
            
            // if cigar has length 1 and is a match, we can use the start and end positions of the read
            let cigar = record.cigar();
            let mut feature_name = String::default();
            if cigar.len() == 1 && cigar.iter().next().unwrap().1 == Operation::AlnMatch {
                
                //eprintln!("startindex: {}, endindex: {}", startindex, endindex);
                feature_name = process_partial_read(features, start_pos, end_pos, &mut ambiguous);
            } else if cigar.len()>1 {
                feature_name= String::default();
                // construct partial reads for each cigar element with AlnMatch
                // keep track of feature_names for each partial read, if they are all the same, we can count the feature
                let mut partial_reads: Vec<(u32, u32)> = Vec::new();
                let mut partial_read_feature_names: Vec<String> = Vec::new();

                let mut start_pos = start_pos;
                for cig in cigar.iter() {
                    if cig.1 != Operation::AlnMatch {
                        // Skip all cigar elements that are not matches, but add the length to the start position
                        start_pos += cig.0;
                        continue;
                    }
                    let partial_end_pos = start_pos + cig.0;
                    partial_reads.push((start_pos, partial_end_pos));
                    let temp_feature_name = process_partial_read(features,  start_pos, partial_end_pos, &mut ambiguous);
                    // if ambiguous flag is set, we can stop here, otherwise we can add the feature name to the list
                    if ambiguous {
                        break;
                    } else {
                        partial_read_feature_names.push(temp_feature_name);
                    }
                    start_pos += cig.0;
                }
                // now check if all feature names are the same
                if !ambiguous {
                    let first_feature_name = &partial_read_feature_names[0];
                    if partial_read_feature_names.iter().all(|f| f == first_feature_name) {
                        feature_name = first_feature_name.clone();
                    } else {
                        ambiguous = true;
                    }
                }
                
                // if cigar has length > 1, we need to check each cigar element
                // for cig in cigar.iter() {
                //     if cig.1 != Operation::AlnMatch {
                //         // Skip all cigar elements that are not matches, but add the length to the start position
                //         start_pos += cig.0;
                //         continue;
                //     }
                //     let partial_end_pos = start_pos + cig.0;
                //     if startindex < features.0.len()
                //         && current_feature.start <= partial_end_pos
                //         && current_feature.end >= start_pos
                //     {
                //         let feature_name = &current_feature.name;
                //         while startindex < features.0.len() {
                //             current_feature = &gtf_end_sorted[reference].0[startindex];
                //             if current_feature.start > partial_end_pos {
                //                 break;
                //             }
                //             // check if same feature name, otherwise ambiguous
                //             if feature_name != &current_feature.name {
                //                 //println!("Ambiguous feature found! {}: {}:{}-{} (position: {}); (alternate: {}) (index: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos, feature_name, index);
                //                 ambiguous = true;
                //                 break;
                //             }
                //             //println!("Feature found! {}: {}:{}-{} (position: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos);
                //             feature_count = 1;
                //             startindex += 1;
                //         }
                //     }
                //     start_pos += cig.0;
                // }
            }

            if ambiguous {
                *counts.entry("__ambiguous".to_string()).or_insert(0) += 1;
            } else if feature_name == String::default() {
                *counts.entry("__no_feature".to_string()).or_insert(0) += 1;
            } else {
                *counts.entry(feature_name.clone()).or_insert(0) += 1;
            }
        } else {
            // No reference found for this read
            // TODO: check if we should add this to __no_feature or we should throw an error
            let count = counts.entry("__no_feature".to_string()).or_insert(0);
            *count += 1;
        }
    }
}

fn process_partial_read(features: &IntervalTree, start_pos: u32, end_pos: u32, ambiguous: &mut bool) -> String {
    let mut feature_name = String::default();
    todo!("process_partial_read");
    
    feature_name
}

fn get_end_index(features: &[Feature], _start_pos: u32, end_pos: u32) -> usize{
    // find index of the first feature that starts after the read end
    let endindex = match features.binary_search_by(|f| {
        if f.start > end_pos {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Less
        }
    }) {
        Ok(index) | Err(index) => {
            //eprintln!("f.start>end_pos: {}, f.start: {}, end_pos: {}, index: {}, {}, feature:{:?}", features.len(), features[index].start, end_pos, index, features[index].start > end_pos, features[index]);
            index + 1      
        }
    };
    endindex
}

fn get_start_index(features: &Vec<Feature>, start_pos: u32, _end_pos: u32) -> usize {
    let startindex = match features.binary_search_by(|f| {
        if f.end < start_pos {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    }) {
        Ok(index) | Err(index) => {
            //eprintln!("f.end<start_pos: {}, f.end: {}, start_pos: {}, index: {}, {}, feature:{:?}", features.len(), features[index-1].end, start_pos, index, features[index-1].end < start_pos, features[index-1]);
            if index > 0 {
                index - 1
            } else {
                index
            }
        }
    };
    features[startindex].end_sorted_index
}
