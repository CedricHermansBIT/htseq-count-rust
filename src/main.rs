use bam::record::cigar::Operation;
use bam::record::tags::TagValue;
use bam::BamReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use structopt::StructOpt;

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
}

// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash)]
struct Feature {
    //type_: String,
    name: String,
    chr: String,
    start: u32,
    end: u32,
}

fn read_gtf(file_path: &str, feature_type_filter: &str) -> HashMap<String, Vec<Feature>> {
    let mut map: HashMap<String, Vec<Feature>> = HashMap::new();
    let file = File::open(file_path).expect("Kan het bestand niet openen");
    let reader = BufReader::new(file);
    let mut counter = 0;
    for line in reader.lines() {
        counter += 1;
        if counter % 100000 == 0 {
            println!("{} GFF lines processed.", counter);
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
                //type_: feature_name.to_string(),
                name: name.to_string(),
                chr,
                start,
                end,
            };
            map.entry(feature.chr.clone()).or_default().push(feature);
        }
    }
    // Sort the features by start position
    for features in map.values_mut() {
        features.sort_by(|a, b| a.end.cmp(&b.end));
    }
    println!("{} GFF lines processed.", counter);
    map
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
    sorted_keys.sort_by(|a, b| a.to_lowercase().cmp(&b.to_lowercase()));
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
    sorted_keys.sort_by(|a, b| a.to_lowercase().cmp(&b.to_lowercase()));
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
    let map = read_gtf(&args.gtf, args.t.as_str());

    let mut counts = HashMap::with_capacity(map.len());
    // add all features to the map
    for features in map.values() {
        for feature in features {
            counts.entry(feature.name.clone()).or_default();
        }
    }

    // Add the special keys
    counts.insert("__no_feature".to_string(), 0);
    counts.insert("__ambiguous".to_string(), 0);
    counts.insert("__not_aligned".to_string(), 0);
    counts.insert("__too_low_aQual".to_string(), 0);
    counts.insert("__alignment_not_unique".to_string(), 0);

    let mut counter = 0;

    for record in bam {
        counter += 1;
        if counter % 100000 == 0 {
            println!("{} records processed.", counter);
        }
        let record = record.unwrap();
        if should_skip_record(&record, &mut counts, &args) {
            continue;
        }
        let ref_id = record.ref_id();

        let reference = &reference_names[ref_id as usize];
        if map.contains_key(reference) {
            let mut start_pos: u32 = record.start().try_into().unwrap();
            let end_pos: u32 = record.calculate_end().try_into().unwrap();
            let features = &map[reference];

            // Start from the index of the first feature that ends after the read start
            let mut index = match features.binary_search_by(|f| {
                if f.end <= start_pos {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            }) {
                Ok(index) => index,
                Err(index) => index,
            };
            //println!("{}: {}-{}/{}", reference, start_pos, end_pos, &map[reference][index].end);
            let mut feature_name = &String::default();
            let mut feature_count = 0;
            let mut ambiguous = false;
            if index >= features.len() {
                // No feature found for this read
                let count = counts.entry("__no_feature".to_string()).or_insert(0);
                *count += 1;
                continue;
            }
            let mut current_feature = &map[reference][index];

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

            // Feature start should be before (or equal to) the read start and feature end should be after (or equal to) the read end
            // This means the read falls fully within the feature (no exon junction)
            if current_feature.start <= start_pos && current_feature.end >= end_pos {
                feature_name = &current_feature.name;
                while index < features.len() && current_feature.start <= end_pos {
                    current_feature = &map[reference][index];
                    // check if same feature name, otherwise ambiguous
                    if feature_name != &current_feature.name {
                        //println!("Ambiguous feature found! {}: {}:{}-{} (position: {}); (alternate: {}) (index: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos, feature_name, index);
                        ambiguous = true;
                        break;
                    }
                    //println!("Feature found! {}: {}:{}-{} (position: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos);
                    feature_count = 1;
                    index += 1;
                }
            // If the read is not fully within the feature, we might need to check for exon junctions
            } else {
                let cigar = record.cigar();
                if cigar.len() > 1 {
                    // match the operation of the cigar elements

                    for cig in cigar.iter() {
                        if cig.1 != Operation::AlnMatch {
                            // Skip all cigar elements that are not matches, but add the length to the start position
                            start_pos += cig.0;
                            continue;
                        }
                        let partial_end_pos = start_pos + cig.0;
                        if index < features.len()
                            && current_feature.start <= partial_end_pos
                            && current_feature.end >= start_pos
                        {
                            let feature_name = &current_feature.name;
                            while index < features.len() {
                                current_feature = &map[reference][index];
                                if current_feature.start > partial_end_pos {
                                    break;
                                }
                                // check if same feature name, otherwise ambiguous
                                if feature_name != &current_feature.name {
                                    //println!("Ambiguous feature found! {}: {}:{}-{} (position: {}); (alternate: {}) (index: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos, feature_name, index);
                                    ambiguous = true;
                                    break;
                                }
                                //println!("Feature found! {}: {}:{}-{} (position: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos);
                                feature_count = 1;
                                index += 1;
                            }
                        }
                        start_pos += cig.0;
                    }
                }
            }
            if ambiguous {
                *counts.entry("__ambiguous".to_string()).or_insert(0) += 1;
            } else if feature_count == 0 {
                *counts.entry("__no_feature".to_string()).or_insert(0) += 1;
            } else {
                *counts.entry(feature_name.clone()).or_insert(0) += feature_count;
            }
        } else {
            // No reference found for this read
            // TODO: check if we should add this to __no_feature or we should throw an error
            let count = counts.entry("__no_feature".to_string()).or_insert(0);
            *count += 1;
        }
    }

    if args.counts_output.is_some() {
        write_counts(counts, args, counter);
    } else {
        print_output(counts, args, counter);
    }
}
