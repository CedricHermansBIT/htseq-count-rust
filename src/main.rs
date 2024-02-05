use bam::record::cigar::Operation;
use bam::record::tags::TagValue;
use bam::{BamReader, RecordWriter, SamWriter};
use feature::Feature;
use intervaltree::IntervalTree;
use interval::Interval;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
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
    let ref_names_to_id: HashMap<String, i32> = reference_names.iter().enumerate().map(|(i, s)| (s.clone(), i as i32)).collect();
    let mut output_sam: Option<SamWriter<BufWriter<File>>> = None;
    if args.output_sam.is_some() {
        output_sam = Some(SamWriter::from_path(args.output_sam.clone().unwrap(), header).expect("Could not create output sam file"));
    }
    // bam fields: https://docs.rs/bam/0.3.0/bam/record/struct.Record.html

    // Read the gtf file
    let gtf = read_gtf(&args.gtf, args.t.as_str(), &ref_names_to_id);

    // let read= 21940455;
    // eprintln!("Searching for reads overlapping position {}-{}...", read, read+25);
    // for overlap in gtf["1"].overlap(read, read+25) {
    //      eprintln!("overlap: {:?}", overlap);
    // }

    
    if args.export_feature_tree.is_some() {
        eprintln!("Exporting feature trees as dot files...");
        let mut file = File::create(args.export_feature_tree.clone().unwrap()).expect("Unable to create file");
        for (chr, tree) in gtf.iter().enumerate() {
            // top node for each tree is the chromosome
            file.write_all(format!("digraph {} {{\n", reference_names[chr]).as_bytes()).expect("Unable to write data");
            let _ = tree.as_ref().unwrap().top_node.clone().unwrap().write_structure(&mut file, 0);
            file.write_all("}\n".as_bytes()).expect("Unable to write data");
        }
    }
    
    // exit(1) to prevent the rest of the program from running for debugging purposes
    //std::process::exit(1);

    let mut counts = prepare_count_hashmap(&gtf);
    let mut read_to_feature: Vec<FeatureType> = Vec::new();
    let mut counter = 0;

    count_reads(bam, &mut counter, &mut counts, &args, gtf, &mut read_to_feature);

    if args.output_sam.is_some() {
        eprintln!("Writing output sam file...");
        // loop through the bam file again and write the reads to the output sam file
        let bam = BamReader::from_path(args.bam.clone(), args.n).expect("Could not read bam file");
        for (i,record) in bam.enumerate() {
            let mut record = record.unwrap();
            let feature = read_to_feature[i].as_bytes();
            record.tags_mut().push_string(b"XF", &feature);
            output_sam.as_mut().unwrap().write(&record).unwrap();
        }


        let mut output_sam = output_sam.unwrap();
        output_sam.flush().unwrap();
        output_sam.finish().unwrap();
    }
    if args.counts_output.is_some() {
        write_counts(counts, args, counter);
    } else {
        print_output(counts, args, counter);
    }
}

enum FeatureType {
    Name(String),
    __no_feature,
    __ambiguous(String),
    __not_aligned,
    __too_low_aQual,
    __alignment_not_unique,
    none,
}

impl FeatureType {
    fn as_bytes(&self) -> Vec<u8> {
        match self {
            FeatureType::Name(s) => s.as_bytes().to_vec(),
            FeatureType::__no_feature => b"__no_feature".to_vec(),
            FeatureType::__ambiguous(s) => format!("__ambiguous[{}]", s).as_bytes().to_vec(),
            FeatureType::__not_aligned => b"__not_aligned".to_vec(),
            FeatureType::__too_low_aQual => b"__too_low_aQual".to_vec(),
            FeatureType::__alignment_not_unique => b"__alignment_not_unique".to_vec(),
            FeatureType::none => b"".to_vec(),
        }
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
    #[structopt(short = "m", long = "mode", default_value = "union", possible_values = &["intersection-strict", "intersection-nonempty", "union"], help = "Mode to use for counting reads overlapping features. Possible values: intersection-strict, intersection-nonempty, union (default: intersection-strict).")]
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
    export_feature_tree: Option<String>,

    #[structopt(
        short = "o",
        long = "output_sam",
        help = "Create a SAM file with the reads and their features."
    )]
    output_sam: Option<String>,
}

fn prepare_count_hashmap(gtf: &Vec<Option<IntervalTree>>) -> HashMap<String, i32> {
    let mut counts = HashMap::with_capacity(gtf.len());
    // add all features to the map
    for tree in gtf {
        if tree.is_none() {
            continue;
        }
        for feature in tree.as_ref().unwrap().all_intervals.iter() {
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

fn read_gtf(file_path: &str, feature_type_filter: &str, ref_names_to_id: &HashMap<String, i32>) -> Vec<Option<IntervalTree>> {
    let mut map: HashMap<i32, Vec<Interval>> = HashMap::new();
    let file = File::open(file_path).expect("Kan het bestand niet openen");
    let mut reader = BufReader::new(file);
    let mut counter = 0;
    let mut line = String::new();
    while reader.read_line(&mut line).unwrap() > 0 {
        counter += 1;
        if counter % 100000 == 0 {
            eprintln!("{} GFF lines processed.", counter);
        }
    
        if line.starts_with('#') {
            line.clear();
            continue;
        }
        let mut fields = line.split('\t');
        let chr_id = ref_names_to_id.get(fields.next().unwrap()).unwrap();
        //eprintln!("chr_id: {}", chr_id);
        let feature_name = fields.nth(1).unwrap();
        //eprintln!("feature_name: {}", feature_name);


        let start = fields.next().unwrap().parse::<i32>().unwrap();
        let end = fields.next().unwrap().parse::<i32>().unwrap()+1;

        if feature_name != feature_type_filter {
            line.clear();
            continue;
        }
        let attributes = fields.nth(3).unwrap();
        //eprintln!("attributes: {}", attributes);

        let name = attributes.split(';')
            .find(|&attr| attr.contains("gene_name"))
            .unwrap_or("")
            .trim()
            .strip_prefix("gene_name ")
            .unwrap_or("")
            .trim_matches('"');

        let feature = Feature {
            name: name.to_string(),
            chr: chr_id.clone(),
            start: min(start, end),
            end: max(start, end),
        };

        map.entry(*chr_id).or_default().push(Interval::new(start, end, Some(feature)));
        line.clear();
    }
    //let node = Node::from_intervals(intervals);
    // print the center node (withouth the children)
    // construct the tree
    // let tree= IntervalTree::from_tuples(tuples);
    // let top_node = tree.top_node.unwrap();
    // eprintln!("tree: center: {}, depth: {}, balance: {}", top_node.x_center, top_node.depth, top_node.balance);
    // eprintln!("Boundary table: {:?}", tree.boundary_table);

    eprintln!("{} GFF lines processed.", counter);
    eprint!("Creating IntervalTree for each chromosome...");
    // prepopulate the map with empty trees
    let mut result: Vec<Option<IntervalTree>> = Vec::with_capacity(ref_names_to_id.len());
    //eprintln!("{:?} chromosomes found.", ref_names_to_id.keys());
    for _ in 0..ref_names_to_id.len() {
        result.push(None);
    }

    for (chr, intervals) in map {
        let tree = IntervalTree::new(Some(intervals));
        result[chr as usize] = Some(tree);
        
    }
    eprintln!("done.");

    //eprintln!("The amount of not empty trees: {}", result.iter().filter(|x| x.is_some()).count());
    result    
}

fn should_skip_record(
    record: &bam::Record,
    counts: &mut HashMap<String, i32>,
    args: &Args,
    read_to_feature: &mut Vec<FeatureType>,
) -> bool {
    // Skip all reads that are not aligned
    if record.ref_id() < 0 {
        read_to_feature.push(FeatureType::__not_aligned);
        *counts.entry("__not_aligned".to_string()).or_insert(0) += 1;
        return true;
    }
    // Skip all reads that are secondary alignments
    if args.secondary_alignments == "ignore" && record.flag().all_bits(0x100) {
        read_to_feature.push(FeatureType::__not_aligned);
        return true;
    }
    // Skip all reads that are supplementary alignments
    if args.supplementary_alignments == "ignore" && record.flag().all_bits(0x800) {
        read_to_feature.push(FeatureType::none);
        return true;
    }
    // Skip all reads with MAPQ alignment quality lower than the given minimum value
    if record.mapq() < args.a {
        read_to_feature.push(FeatureType::__too_low_aQual);
        *counts.entry("__too_low_aQual".to_string()).or_insert(0) += 1;
        return true;
    }
    // Skip all reads that have an optional field "NH" with value > 1
    if let Some(TagValue::Int(i, _)) = record.tags().get(b"NH") {
        if i > 1 {
            read_to_feature.push(FeatureType::__alignment_not_unique);
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


fn count_reads(bam: BamReader<File>, counter: &mut i32, counts: &mut HashMap<String, i32>, args: &Args, gtf: Vec<Option<IntervalTree>>, read_to_feature: &mut Vec<FeatureType>) {

    let processing_function = match args._m.as_str() {
        "intersection-strict" => process_intersection_strict_read,
        "intersection-nonempty" => process_intersection_nonempty_read,
        "union" => process_union_read,
        _ => panic!("Invalid mode"),
    };

    for record in bam {
        //eprintln!("{} records processed.", counter);
        *counter += 1;
        if *counter % 100000 == 0 {
            //println!("{} records processed.", counter);
            eprintln!("{} records processed.", counter);
        }
        let record = record.unwrap();
        // if String::from_utf8_lossy(record.name())=="SRR001432.281211 USI-EAS21_0008_3445:8:7:657:535 length=25" {
        //     //println!("{}: {}-{}", String::from_utf8_lossy(record.name()), record.start(), record.calculate_end());
        //     eprintln!("this one");
        // }
        if should_skip_record(&record, counts, args, read_to_feature) {
            continue;
        }
        // todo
        //unimplemented!("todo");
        
        let ref_id = record.ref_id() as usize;

        if gtf[ref_id].is_some() {
            let start_pos = record.start();
            let end_pos = record.calculate_end() + 1;
            let features = &gtf[ref_id].as_ref().unwrap();

            let mut ambiguous = false;
            
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
            let mut feature = Feature::default();
            if cigar.len() == 1 && cigar.iter().next().unwrap().1 == Operation::AlnMatch {
                
                //eprintln!("startindex: {}, endindex: {}", startindex, endindex);
                feature = processing_function(&features, start_pos, end_pos, &mut ambiguous, &record, counts, read_to_feature);
            } else {
                //todo!("cigar length > 1");
                // construct partial reads for each cigar element with AlnMatch
                // keep track of feature_names for each partial read, if they are all the same, we can count the feature
                let mut start_pos = start_pos;
                                
                // if cigar has length > 1, we need to check each cigar element
                for cig in cigar.iter() {
                    if cig.1 != Operation::AlnMatch {
                        // Skip all cigar elements that are not matches, but add the length to the start position
                        start_pos += cig.0 as i32;
                        continue;
                    }
                    let partial_end_pos = start_pos + cig.0 as i32;
                    
                    let temp_feature_name = processing_function(features,  start_pos, partial_end_pos, &mut ambiguous, &record, counts, read_to_feature);
                    // if ambiguous flag is set, we can stop here, otherwise we can add the feature name to the list
                    if ambiguous {
                        break;
                    } else if feature != Feature::default() && feature.name != temp_feature_name.name {
                        check_ambiguity_union(&features.overlap(start_pos, partial_end_pos), 
                                            start_pos, partial_end_pos, &mut feature, &mut ambiguous, &record,
                                            counts, read_to_feature);
                        if ambiguous {
                            break;
                        }
                    } 
                    else {
                        feature = temp_feature_name;
                    }
                    start_pos = partial_end_pos+1;
                }
            }

            if ambiguous {
                *counts.entry("__ambiguous".to_string()).or_insert(0) += 1;
            } else if feature.name == String::default() {
                *counts.entry("__no_feature".to_string()).or_insert(0) += 1;
                read_to_feature.push(FeatureType::__no_feature);
            } else {
                *counts.entry(feature.name.clone()).or_insert(0) += 1;
                read_to_feature.push(FeatureType::Name(feature.name));
            }
        } else {
            // No reference found for this read
            // TODO: check if we should add this to __no_feature or we should throw an error
            *counts.entry("__no_feature".to_string()).or_insert(0) +=1;
            read_to_feature.push(FeatureType::__no_feature);
        }
    }

    eprintln!("{} records processed.", counter);
}

fn process_union_read(features: &IntervalTree, start_pos: i32, end_pos: i32, ambiguous: &mut bool, record: &bam::Record, counts: &mut HashMap<String, i32>, read_to_feature: &mut Vec<FeatureType>) -> Feature {
    let mut feature = Feature::default();

    let overlapping_features = features.overlap(start_pos, end_pos);

    if overlapping_features.len() == 0 {
        return feature;
    } else if overlapping_features.len() == 1 {
        feature = overlapping_features.iter().next().unwrap().data.as_ref().unwrap().clone();
    } else {
        check_ambiguity_union(&overlapping_features, start_pos, end_pos, &mut feature, ambiguous, &record, counts, read_to_feature);
    }

    //todo!("process_partial_read");
    
    feature
}

fn check_ambiguity_union(overlapping_features: &HashSet<&Interval>, _start_pos: i32, _end_pos: i32, feature: &mut Feature, ambiguous: &mut bool, record: &bam::Record, counts: &mut HashMap<String, i32>, read_to_feature: &mut Vec<FeatureType>) {
    let feature_names: HashSet<String> = overlapping_features.iter().map(|f| f.data.as_ref().unwrap().name.clone()).collect();
    match feature_names.len() {
        0 => {},
        1 => *feature = overlapping_features.iter().next().unwrap().data.as_ref().unwrap().clone(),
        _ => *ambiguous = {
            //use feature_names as the ambiguous feature
            let mut feature_names: Vec<String> = feature_names.into_iter().collect();
            feature_names.sort();
            // write XF:Z:__ambiguous[feature_names] to the tags of the record and write it to the output sam file (separator: +)
            read_to_feature.push(FeatureType::__ambiguous(feature_names.join("+")));
            for feature_name in feature_names {
                *counts.entry(feature_name).or_insert(0) += 1;
            }


            true
        }

    }
}


fn _check_ambiguity(overlapping_features: &HashSet<&Interval>, start_pos: i32, end_pos: i32, feature: &mut Feature, ambiguous: &mut bool, output_sam: &mut Option<SamWriter<BufWriter<File>>>, record: &bam::Record) {
    let mut contained_by = Vec::new();
    for overlap in overlapping_features {
        if overlap.start <= start_pos && overlap.end >= end_pos {
            contained_by.push(overlap.data.as_ref().unwrap().clone());
        }
    }
    
    let feature_names: HashSet<String> = contained_by.iter().map(|f| f.name.clone()).collect();

    match feature_names.len() {
        0 => *ambiguous = {
            // get overlapping_features names
            let mut names = HashSet::new();
            for overlap in overlapping_features {
                names.insert(overlap.data.as_ref().unwrap().name.clone());
            }
            if names.len() == 1 {
                *feature = overlapping_features.iter().next().unwrap().data.as_ref().unwrap().clone();
                false
            } else {
            // sort alphabetically
            let mut names: Vec<String> = names.into_iter().collect();
            names.sort();
            
            // write XF:Z:__ambiguous[feature_names] to the tags of the record and write it to the output sam file (separator: +)
            if let Some(output_sam) = output_sam {
                let mut record = record.clone();
                record.tags_mut().push_string(b"XF", format!("__ambiguous[{}]", names.join("+")).as_bytes());
                output_sam.write(&record).unwrap();
            }
            true
            }
        },
        1 => *feature = contained_by[0].clone(),
        _ => *ambiguous = {
            //use feature_names as the ambiguous feature
            let mut feature_names: Vec<String> = feature_names.into_iter().collect();
            feature_names.sort();
            // write XF:Z:__ambiguous[feature_names] to the tags of the record and write it to the output sam file (separator: +)
            if let Some(output_sam) = output_sam {
                let mut record = record.clone();
                record.tags_mut().push_string(b"XF", format!("__ambiguous[{}]", feature_names.join("+")).as_bytes());
                output_sam.write(&record).unwrap();
            }
            true
        }
    }
}

fn process_intersection_nonempty_read(features: &IntervalTree, start_pos: i32, end_pos: i32, ambiguous: &mut bool, record: &bam::Record, counts: &mut HashMap<String, i32>, read_to_feature: &mut Vec<FeatureType>) -> Feature {
    todo!("process_partial_read for intersection-nonempty");
}

fn process_intersection_strict_read(features: &IntervalTree, start_pos: i32, end_pos: i32, ambiguous: &mut bool, record: &bam::Record, counts: &mut HashMap<String, i32>, read_to_feature: &mut Vec<FeatureType>) -> Feature {
    todo!("process_partial_read for intersection-strict");
}