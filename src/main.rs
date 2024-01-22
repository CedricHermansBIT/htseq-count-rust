use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use structopt::StructOpt;
use bam::BamReader;

// Use the structopt crate to parse command line arguments
#[derive(StructOpt)]
struct Args {
    // Number of threads
    #[structopt(short = "n", long = "number", default_value = "10", help="Number of threads")]
    n: u16,

    // Mode
    #[structopt(short = "m", long = "mode", default_value = "intersection-strict")]
    m: String,

    // Stranded
    #[structopt(long = "stranded")]
    stranded: bool,

    // Quality filter
    #[structopt(short = "a", long = "amount", default_value = "10", help="Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads.")]
    a: i32,

    // Type of feature to be used
    #[structopt(short = "t", long = "type", default_value = "exon", help="Feature type (3rd column in GTF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)")]
    t: String,

    // Feature ID
    #[structopt(short = "i", long = "id", default_value = "gene_name", help="GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). All feature of the right type (see -t option) within the same GTF
    attribute will be added together. The typical way of using this option is to count all exonic reads from each gene and add the exons but other uses are possible
    as well. You can call this option multiple times: in that case, the combination of all attributes separated by colons (:) will be used as a unique identifier,
    e.g. for exons you might use -i gene_id -i exon_number.")]
    i: Vec<String>,

    // Name and type of the bam file
    #[structopt(name = "bam", default_value="test.bam")]
    bam: String,

    // Name and type of the gtf file
    #[structopt(name = "gtf", default_value="test.gtf")]
    gtf: String,
}


// Struct to store the features
#[derive(Debug)]
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
    for features in map.values_mut() {
        features.sort_by(|a, b| a.start.cmp(&b.start));
    }
    println!("{} GFF lines processed.", counter);
    map
}

fn main() {
    // Command line arguments
    let args = Args::from_args();
    println!("n = {}", args.n);
    println!("m = {}", args.m);
    println!("stranded = {}", args.stranded);
    println!("a = {}", args.a);
    println!("t = {}", args.t);
    println!("i = {:?}", args.i);
    println!("bam = {}", args.bam);
    println!("gtf = {}", args.gtf);
    
    // Try to open the bam file, if it fails, print an error message
    let bam = BamReader::from_path(args.bam, args.n).expect("Could not read bam file");
    let header = bam.header().clone();
    let reference_names: Vec<String> = header.reference_names().to_owned();
    // bam fields: https://docs.rs/bam/0.3.0/bam/record/struct.Record.html
    //println!("{:?}", reference_names);

    
    // Read the gtf file
    let map = read_gtf(&args.gtf, args.t.as_str());
    // print last feature of chr1
    //println!("{:?}", map["chr1"][map["chr1"].len()-1]);
    
    
    let mut counts = HashMap::new();
    // add all features to the map
    for features in map.values() {
        for feature in features {
            counts.entry(feature.name.clone()).or_default();
        }
    }
    // add __no_feature to the map
    counts.entry("__no_feature".to_string()).or_default();
    // add __ambiguous to the map
    counts.entry("__ambiguous".to_string()).or_default();
    // add __not_aligned to the map
    counts.entry("__not_aligned".to_string()).or_default();

    let mut counter = 0;

    for record in bam {
        counter += 1;
        if counter % 100000 == 0 {
            println!("{} records processed.", counter);
        }
        let record = record.unwrap();
        let ref_id: usize= record.ref_id().try_into().unwrap();
        let reference=reference_names[ref_id].clone();
        if map.contains_key(&reference) {
            let start_pos: u32 = record.start().try_into().unwrap();
            let end_pos: u32 = record.calculate_end().try_into().unwrap();
            // Start from the index of the first feature that starts before the read
            let mut index = match map[&reference].binary_search_by(|f| start_pos.cmp(&f.start)) {
                Ok(index) => index,
                Err(index) => if index > 0 { 
                    //println!("{}", index);
                    index - 1 } else { 
                    println!("index: {}, reference: {}, start_pos: {}, end_pos: {}", index, reference, start_pos, end_pos);
                    std::process::exit(1)
                 },
            };
            let mut feature_name= "".to_string();
            let mut feature_count = 0;
            let mut ambiguous = false;
            while index < map[&reference].len() && map[&reference][index].start <= end_pos {
                // Feature start should be before (or equal to) the read start and feature end should be after (or equal to) the read end
                if map[&reference][index].start <= start_pos && map[&reference][index].end >= end_pos {
                    if feature_name == "" {
                        feature_name = map[&reference][index].name.clone();
                    }
                    // check if same feature name, otherwise ambiguous
                    if feature_name != map[&reference][index].name {
                        //println!("Ambiguous feature found! {}: {}:{}-{} (position: {}); (alternate: {}) (index: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos, feature_name, index);
                        ambiguous = true;
                        break;
                    }
                    //println!("Feature found! {}: {}:{}-{} (position: {})", map[&reference][index].name,reference, map[&reference][index].start, map[&reference][index].end, start_pos);
                    feature_count += 1;
                }
                // Go to the next feature
                index += 1;
            }
            if ambiguous {
                //println!("Ambiguous feature found! {}: {}:{}-{} (position: {})", map[&reference][i].name,reference, map[&reference][i].start, map[&reference][i].end, start_pos);
                let count = counts.entry("__ambiguous".to_string()).or_insert(0);
                *count += 1;
            } else if feature_count == 0 {
                let count = counts.entry("__no_feature".to_string()).or_insert(0);
                *count += 1;
            } else {
                let count = counts.entry(feature_name).or_insert(0);
                *count += feature_count;
            }
        } else {
            let count = counts.entry("__not_aligned".to_string()).or_insert(0);
            *count += 1;
        }
    }

    // Print de HashMap
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort_by(|a, b| a.to_lowercase().cmp(&b.to_lowercase()));
    for key in sorted_keys {
        if key.starts_with("__") {
            continue;
        }
        println!("{}: {}", key, counts[key]);
    }
    println!("{}: {}", "__no_feature", counts["__no_feature"]);
    println!("{}: {}", "__ambiguous", counts["__ambiguous"]);
    println!("{}: {}", "__not_aligned", counts["__not_aligned"]);
}
