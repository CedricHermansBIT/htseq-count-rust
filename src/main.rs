use bam::record::cigar::Operation;
use bam::record::tags::TagValue;
use bam::{RecordReader,BamReader, RecordWriter, SamReader, SamWriter};
use feature::Feature;
use intervaltree::IntervalTree;
use interval::Interval;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::mpsc;
use std::thread;
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
    let mut reads_reader = ReadsReader::from_path(args.bam.clone(), args.n);
    let header = reads_reader.header().clone();
    // check if the header is valid
    check_header_validity(&header, &args);

    let reference_names: Vec<String> = header.reference_names().to_owned();
    let ref_names_to_id: HashMap<String, i32> = reference_names.iter().enumerate().map(|(i, s)| (s.clone(), i as i32)).collect();
    
    // Spawn a thread to write the output sam file
    let (sender, receiver) = mpsc::channel::<FeatureType>();
    //let mut output_sam: Option<SamWriter<BufWriter<File>>> = None;
    // by default, the writer thread does nothing and discards the input
    let output_sam = args.output_sam.clone();
    let input_reads = args.bam.clone();
    let writer_thread = if let Some(output_sam) = output_sam {
        thread::spawn(move || {
            let mut output_sam = SamWriter::from_path(output_sam, header).expect("Could not create output sam file");
            let mut bam = ReadsReader::from_path(input_reads, args.n);
            // read the
            let mut record= bam::Record::new();
            for type_ in receiver {
                record = match bam.read_into(&mut record) {
                    Ok(true) => record,
                    Ok(false) => break,
                    Err(e) => panic!("{}", e),
                };
                let feature = type_.as_bytes();

                // if 'XS' tag (strandedness) already exists, change the type from A (character) to Z (string) // else, do nothing
                // we change this tag from A to Z because htseqcount does this as well
                // TODO: check if we need to support this or not since it is extra work, for no obvious reason
                if record.tags().get(b"XS").is_some() {
                    //get a copy from all tags
                    let tags = record.tags().clone();
                    // clear the old tags
                    record.tags_mut().clear();
                    // add the old tags back, but change the type of the XS tag to Z
                    for (key, value) in tags.iter() {
                        match value {
                            TagValue::Char(c) => {
                                if key == *b"XS" {
                                    // note: c = + or -, but are represented as 43 and 45 in ASCII
                                    record.tags_mut().push_string(b"XS", &[c]);
                                } else {
                                    record.tags_mut().push_char(&key, c);
                                }
                            }
                            TagValue::Int(i, _t) => record.tags_mut().push_num(&key, i as i32),
                            TagValue::Float(f) => record.tags_mut().push_num(&key, f),
                            TagValue::String(s, _t) => record.tags_mut().push_string(&key, s),
                            TagValue::FloatArray(a) => record.tags_mut().push_array(&key, a.raw()),
                            TagValue::IntArray(a) => record.tags_mut().push_array(&key, a.raw()),
                        }
                    }

                    //eprintln!("tags: {:?}", record.tags().raw());
                }
                
                // push also the feature to the XF tag
                record.tags_mut().push_string(b"XF", &feature);
                output_sam.write(&record).unwrap();
            }
        })
    } else {
        thread::spawn(move || { 
            for _ in receiver {
                // do nothing
            }
        })
    };
    // bam fields: https://docs.rs/bam/0.3.0/bam/record/struct.Record.html

    // Read the gtf file
    let gtf = read_gtf(&args.gtf, args.t.as_str(), &ref_names_to_id, &args);

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
            let _ = tree.as_ref().unwrap().top_node.clone().unwrap().write_structure(&mut file, 0, reference_names[chr].clone());
        }
    }
    
    // exit(1) to prevent the rest of the program from running for debugging purposes
    //std::process::exit(1);

    let mut counts = prepare_count_hashmap(&gtf);
    //let mut read_to_feature: Vec<FeatureType> = Vec::new();
    let mut counter = 0;

    count_reads(&mut reads_reader, &mut counter, &mut counts, &args, gtf, sender);
    
    if args.output_sam.is_some() {
        eprintln!("Waiting for writer thread to finish...");
    }
    writer_thread.join().unwrap();
    
    if args.counts_output.is_some() {
        write_counts(counts, args, counter);
    } else {
        print_output(counts, args, counter);
    }
}

fn check_header_validity(header: &bam::Header, args: &Args) {
    if header.lines().count() == 0 {
        if args.bam.ends_with(".bam") {
            eprintln!("The header of the bam file is empty. This is likely due to an invalid bam file.");
        } else {
            eprintln!("The header of the sam file is empty. This is likely due to an invalid sam file. (If you used samtools to convert a bam file to a sam file, make sure to use the -h option to include the header in the output file.)");
        std::process::exit(1);
        }
    }
}

enum ReadsReader {
    BamReader(BamReader<File>),
    SamReader(SamReader<BufReader<File>>),
}

impl ReadsReader {

    fn from_path(path: String, n: u16) -> ReadsReader {
        if path.ends_with(".bam") {
            let reader = BamReader::from_path(path, n);
            ReadsReader::BamReader(reader.unwrap())
        } else if path.ends_with(".sam") {
            let reader = SamReader::from_path(path);
            ReadsReader::SamReader(reader.unwrap())
        } else {
            panic!("File type not supported");
        }
    }

    fn header(&self) -> &bam::Header {
        match self {
            ReadsReader::BamReader(reader) => reader.header(),
            ReadsReader::SamReader(reader) => reader.header(),
        }
    }

    fn read_into(&mut self, record: &mut bam::Record) -> Result<bool, std::io::Error> {
        match self {
            ReadsReader::BamReader(reader) => reader.read_into(record),
            ReadsReader::SamReader(reader) => reader.read_into(record),
        }
    }
}

impl Iterator for ReadsReader {
    type Item = Result<bam::Record, std::io::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            ReadsReader::BamReader(reader) => reader.next(),
            ReadsReader::SamReader(reader) => reader.next(),
        }
    }
}



enum FeatureType {
    Name(String),
    NoFeature,
    Ambiguous(String),
    NotAligned,
    TooLowaQual,
    AlignmentNotUnique,
    None,
}

impl FeatureType {
    fn as_bytes(&self) -> Vec<u8> {
        match self {
            FeatureType::Name(s) => s.as_bytes().to_vec(),
            FeatureType::NoFeature => b"__no_feature".to_vec(),
            FeatureType::Ambiguous(s) => format!("__ambiguous[{}]", s).as_bytes().to_vec(),
            FeatureType::NotAligned => b"__not_aligned".to_vec(),
            FeatureType::TooLowaQual => b"__too_low_aQual".to_vec(),
            FeatureType::AlignmentNotUnique => b"__alignment_not_unique".to_vec(),
            FeatureType::None => b"".to_vec(),
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
    #[structopt(
        short = "s",
        long = "stranded",
        help = "Whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation",
        default_value = "yes",
        possible_values = &["yes", "no", "reverse"]
    )]
    stranded: String,

    // Quality filter
    #[structopt(
        short = "a",
        long = "minaqual",
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
        long = "idattr",
        default_value = "gene_name",
        help = "GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). All feature of the right type (see -t option) within the same GTF
    attribute will be added together. The typical way of using this option is to count all exonic reads from each gene and add the exons but other uses are possible
    as well. You can call this option multiple times: in that case, the combination of all attributes separated by colons (:) will be used as a unique identifier,
    e.g. for exons you might use -i gene_id -i exon_number."
    )]
    i: Vec<String>,

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

    // non-unique parameter
    #[structopt(
        long = "nonunique",
        default_value = "none",
        possible_values = &["none", "all", "fraction","random"],
        help = "Whether and how to score reads that are not uniquely aligned or ambiguously assigned to features (choices: none, all, fraction, random; default: none)"
    )]
    nonunique: String,

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
        long = "samout",
        help = "Create a SAM file with the reads and their features."
    )]
    output_sam: Option<String>,
}

fn prepare_count_hashmap(gtf: &Vec<Option<IntervalTree>>) -> HashMap<String, f32> {
    let mut counts: HashMap<String, f32> = HashMap::with_capacity(gtf.len());
    // add all features to the map
    for tree in gtf {
        if tree.is_none() {
            continue;
        }
        for feature in tree.as_ref().unwrap().all_intervals.iter() {
            counts.entry(feature.data.as_ref().unwrap().name().to_string()).or_insert(0.0);
        }
    }

    // Add the special keys
    counts.insert("__no_feature".to_string(), 0f32);
    counts.insert("__ambiguous".to_string(), 0f32);
    counts.insert("__not_aligned".to_string(), 0f32);
    counts.insert("__too_low_aQual".to_string(), 0f32);
    counts.insert("__alignment_not_unique".to_string(), 0f32);
    counts
}

fn read_gtf(file_path: &str, feature_type_filter: &str, ref_names_to_id: &HashMap<String, i32>, args: &Args) -> Vec<Option<IntervalTree>> {
    let mut map: HashMap<i32, Vec<Interval>> = HashMap::new();
    let file = File::open(file_path).expect("Could not open this file");
    let mut reader = BufReader::new(file);
    let mut counter = 0;
    let mut line = String::default();
    //TODO: deal with -i option correctly
    // For now we just take the first -i option as the feature name
    let attribute = args.i.first().unwrap();
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
        let end = fields.next().unwrap().parse::<i32>().unwrap();
        let strand = fields.nth(1).unwrap();
        if feature_name != feature_type_filter {
            line.clear();
            continue;
        }
        let attributes = fields.nth(1).unwrap();
        //eprintln!("attributes: {}", attributes);


        let name = attributes.split(';')
            .find(|&attr| attr.contains(attribute))
            .unwrap_or("")
            .trim()
            .strip_prefix(format!("{} ", attribute).as_str())
            .unwrap_or("")
            .trim_matches('"');

        let feature = Feature::new(name.to_string(), *chr_id, min(start, end), max(start, end), strand.chars().next().unwrap());


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
    counts: &mut HashMap<String, f32>,
    args: &Args,
    sender: &mpsc::Sender<FeatureType>,
) -> bool {
    // Skip all reads that are not aligned
    if record.ref_id() < 0 {
        _ = sender.send(FeatureType::NotAligned);
        *counts.entry("__not_aligned".to_string()).or_insert(0.0) += 1.0;
        return true;
    }
    // Skip all reads that are secondary alignments
    if args.secondary_alignments == "ignore" && record.flag().all_bits(0x100) {
        _ = sender.send(FeatureType::None);
        return true;
    }
    // Skip all reads that are supplementary alignments
    if args.supplementary_alignments == "ignore" && record.flag().all_bits(0x800) {
        _ = sender.send(FeatureType::None);
        return true;
    }
    // Skip all reads with MAPQ alignment quality lower than the given minimum value
    if record.mapq() < args.a {
        _ = sender.send( FeatureType::TooLowaQual);
        *counts.entry("__too_low_aQual".to_string()).or_insert(0.0) += 1.0;
        return true;
    }
    // Skip all reads that have an optional field "NH" with value > 1
    if let Some(TagValue::Int(i, _)) = record.tags().get(b"NH") {
        if i > 1 {
            _ = sender.send( FeatureType::AlignmentNotUnique);
            *counts
                .entry("__alignment_not_unique".to_string())
                .or_insert(0.0) += 1.0;
            // TODO: in nonunique mode "all" we should also increment the features for each alignment (process the read anyways), but we should check what happens if the read is also ambiguous
            return true;
        }
    }
    false
}

fn print_output(counts: HashMap<String, f32>, args: Args, counter: i32) {
    // Print de HashMap
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort();
    for key in sorted_keys {
        if key.starts_with("__") {
            continue;
        }
        if args.nonunique != "fraction" {
            println!("{}{}{}", key, args.delimiter, counts[key]);
            continue;
        }
        // print the fraction as a float with 1 decimal if not 0
        if counts[key] == 0.0 {
            println!("{}{}{}", key, args.delimiter, counts[key]);
            continue;
        }
        println!("{}{}{:.1}", key, args.delimiter, counts[key]);
    }
    println!("__no_feature{}{}", args.delimiter, counts["__no_feature"]);
    println!("__ambiguous{}{}", args.delimiter, counts["__ambiguous"]);
    println!("__too_low_aQual{}{}",args.delimiter, counts["__too_low_aQual"]);
    println!("__not_aligned{}{}", args.delimiter, counts["__not_aligned"]);
    println!("__alignment_not_unique{}{}",args.delimiter, counts["__alignment_not_unique"]);

    // TODO: check the correctness, since it might depend on nonunique mode
    if args.counts {
        println!(
            "Total number of uniquely mapped reads{}{}",
            args.delimiter,
            counter as f32
                - counts["__not_aligned"]
                - counts["__too_low_aQual"]
                - counts["__alignment_not_unique"]
                - counts["__ambiguous"]
                - counts["__no_feature"]
        );
    }
}

fn write_counts(counts: HashMap<String, f32>, args: Args, counter: i32) {
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort();
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
            counter as f32 - counts["__not_aligned"] - counts["__too_low_aQual"] - counts["__alignment_not_unique"] - counts["__ambiguous"] - counts["__no_feature"])
            .as_bytes(),
        )
        .expect("Unable to write data");
    }
}


fn count_reads(reads_reader: &mut ReadsReader, counter: &mut i32, counts: &mut HashMap<String, f32>, args: &Args, gtf: Vec<Option<IntervalTree>>, sender: mpsc::Sender<FeatureType>) {
    let processing_function = match args._m.as_str() {
        "intersection-strict" => process_intersection_strict_read,
        "intersection-nonempty" => process_intersection_nonempty_read,
        "union" => process_union_read,
        _ => panic!("Invalid mode"),
    };
    let mut record = bam::Record::new();
    let mut overlapping_features: Vec<Feature> = Vec::with_capacity(50);
    loop {
        overlapping_features.clear();
        match reads_reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
    //for record in bam {
        //eprintln!("{} records processed.", counter);
        *counter += 1;
        if *counter % 100000 == 0 {
            //println!("{} records processed.", counter);
            eprintln!("{} records processed.", counter);
        }
        if should_skip_record(&record, counts, args, &sender) {
            continue;
        }
        
        let ref_id = record.ref_id() as usize;

        if gtf[ref_id].is_some() {
            let mut start_pos = record.start() +1;
            //let end_pos = record.calculate_end() + 1;
            let features = &gtf[ref_id].as_ref().unwrap();
            let cigar = record.cigar();
            for cig in cigar.iter() {
                if cig.1 != Operation::AlnMatch {
                    // Skip all cigar elements that are not matches, but add the length to the start position
                    // Soft clips are not added to the start position
                    /* if record.name() == "SRR5724993.43083906".to_string().as_bytes(){

                        eprintln!("start_pos: {}, cig:{:?}", start_pos, cig);
                    } */
                    match cig.1 {
                        Operation::Soft => {},
                        Operation::Insertion => start_pos += 1,
                        _ => start_pos += cig.0 as i32 + 1,
                    }
                    continue;
                }
                let partial_end_pos = start_pos + cig.0 as i32 -1 ;
                /* if record.name() == "SRR5724993.43083906".to_string().as_bytes(){

                    eprintln!("start_pos: {}, end_pos: {}, cig:{:?}", start_pos, partial_end_pos, cig);
                }
                 */
                processing_function(features,  start_pos, partial_end_pos, record.flag().is_reverse_strand(), &mut overlapping_features, args);
                /* if record.name() == "SRR5724993.43083906".to_string().as_bytes(){
                    eprintln!("overlapping_features: {:?}", overlapping_features);
                } */
                start_pos = partial_end_pos;
            }

            // get the unique feature_names from the overlapping features, also filter out the empty names
            let mut unique_feature_names= filter_ambiguity_union(&overlapping_features);
            
            // match based on the amount of unique feature names
            let feature_name_len = unique_feature_names.len();
            match feature_name_len {
                0 => {
                    // if there are no unique feature names, we have no feature
                    *counts.entry("__no_feature".to_string()).or_insert(0.0) += 1.0;
                    _ = sender.send(FeatureType::NoFeature);
                },
                1 => {
                    // if there is only one unique feature name, we have a non-ambiguous read
                    let feature_name = unique_feature_names.first().unwrap();
                    *counts.entry(feature_name.clone()).or_insert(0.0) += 1.0;
                    _ = sender.send(FeatureType::Name(feature_name.clone()));
                },
                _ => {
                    // if there are multiple unique feature names, we have an ambiguous read
                    *counts.entry("__ambiguous".to_string()).or_insert(0.0) += 1.0;
                    match args.nonunique.as_str() {
                        "all" => {
                        // we also increment each feature name by 1
                            for feature_name in unique_feature_names.clone() {
                                *counts.entry(feature_name.clone()).or_insert(0.0) += 1.0;
                            }
                        },
                        "fraction" => {
                            // we increment each feature name by 1 divided by the amount of unique feature names
                            let fractional_count = 1.0 / feature_name_len as f32;
                            for feature_name in unique_feature_names.clone() {
                                *counts.entry(feature_name.clone()).or_insert(0.0) += fractional_count;
                            }
                        },
                        "random" => {
                            // we increment one of the feature names by 1
                            let random_index = rand::random::<usize>() % feature_name_len;
                            let feature_name = unique_feature_names[random_index].clone();
                            *counts.entry(feature_name).or_insert(0.0) += 1.0;
                        },
                        _ => {// do nothing
                        },
                    }
                    unique_feature_names.sort();
                    _ = sender.send(FeatureType::Ambiguous(unique_feature_names.join("+")));
                }
            }
        } else {
            // No reference found for this read
            // TODO: check if we should add this to __no_feature or we should throw an error
            *counts.entry("__no_feature".to_string()).or_insert(0.0) += 1.0;
            _ = sender.send(FeatureType::NoFeature);
        }
    }

    eprintln!("{} records processed.", counter);
}

fn process_union_read(features: &IntervalTree, start_pos: i32, end_pos: i32, is_reverse_strand: bool, overlapping_features: &mut Vec<Feature>, args: &Args) {
    let new_overlap = features.overlap(start_pos, end_pos);
    let strand = if is_reverse_strand { '-' } else { '+' };
    // add all overlapping features to the list
    add_stranded_features(new_overlap, strand, overlapping_features, args);
    
}

fn process_intersection_nonempty_read(_features: &IntervalTree, _start_pos: i32, _end_pos: i32, _strand: bool, _overlapping_features: &mut Vec<Feature>, _args: &Args) {
    todo!("process_partial_read for intersection-nonempty");
}

fn process_intersection_strict_read(features: &IntervalTree, start_pos: i32, end_pos: i32, strand: bool, overlapping_features: &mut Vec<Feature>, args: &Args) {
    //todo!("process_partial_read for intersection-strict");
    let new_contained = features.contains(start_pos, end_pos);
    // change new_contained to a Vec<&Interval>
    let strand = if strand { '-' } else { '+' };
    // add all contained features to the list
    add_stranded_features(new_contained, strand, overlapping_features, args);
}

fn filter_ambiguity_union(
    overlapping_features: &[Feature],
) -> Vec<String> {
    let unique_feature_names: HashSet<String> = overlapping_features.iter().map(|x| x.name().to_string().clone()).filter(|x| !x.is_empty()).collect();
    unique_feature_names.into_iter().collect()
}

fn add_stranded_features(new_overlap: Vec<&Interval>, strand: char, overlapping_features: &mut Vec<Feature>, args: &Args) {
    for overlap in new_overlap {
        let feature = overlap.data.as_ref().unwrap();
        match args.stranded.as_str() {
            "yes" => {
                //eprintln!("feature: {:?}", feature);
                if feature.strand() == strand {
                    overlapping_features.push(feature.clone());
                }
            },
            "reverse" => {
                if feature.strand() != strand {
                    overlapping_features.push(feature.clone());
                }
            },
            "no" => {
                overlapping_features.push(feature.clone());
            },
            _ => {
                panic!("Invalid strandedness");
            }
        }
    }
}