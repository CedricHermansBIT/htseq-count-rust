use bam::record::tags::TagValue;
use bam::{RecordReader,BamReader, RecordWriter, SamReader, SamWriter};
use feature::Feature;
use intervaltree::IntervalTree;
use interval::Interval;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::mpsc;
use std::thread;
use clap::Parser;

mod feature;
mod intervaltree;
mod interval;
mod node;

use node::Node;

fn main() {
    // Command line arguments
    let args = Args::parse();

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
    
    // Only create the assignment channel when --samout is actually used.
    // The old code sent one message per record to a thread that simply discarded
    // it during normal counting, which added synchronization overhead to the hot path.
    let output_sam = args.output_sam.clone();
    let input_reads = args.bam.clone();
    let reader_threads = args.n;
    let (assignment_sender, writer_thread) = if let Some(output_sam) = output_sam {
        let (sender, receiver) = mpsc::channel::<FeatureType>();
        let writer_thread = thread::spawn(move || {
            let mut output_sam = SamWriter::from_path(output_sam, header)
                .expect("Could not create output sam file");
            let mut bam = ReadsReader::from_path(input_reads, reader_threads);
            let mut record = bam::Record::new();
            for type_ in receiver {
                match bam.read_into(&mut record) {
                    Ok(true) => (),
                    Ok(false) => break,
                    Err(e) => panic!("{}", e),
                };
                let feature = type_.as_bytes();
                record.tags_mut().push_string(b"XF", &feature);
                output_sam.write(&record).unwrap();
            }
        });
        (Some(sender), Some(writer_thread))
    } else {
        (None, None)
    };
    // bam fields: https://docs.rs/bam/0.3.0/bam/record/struct.Record.html

    // Read the gtf file
    let gtf = read_gtf(&args.gtf, &args.t, &ref_names_to_id, &args);

    // let read= 21940455;
    // eprintln!("Searching for reads overlapping position {}-{}...", read, read+25);
    // for overlap in gtf["1"].overlap(read, read+25) {
    //      eprintln!("overlap: {:?}", overlap);
    // }

    
    if args.export_feature_tree.is_some() {
        eprintln!("Exporting feature trees as dot files...");
        let mut file = File::create(args.export_feature_tree.clone().unwrap()).expect("Unable to create file");
        for (chr, tree) in gtf.iter().enumerate().take(reference_names.len()) {
            if let Some(tree) = tree {
                if let Some(top_node) = &tree.top_node {
                    let _ = top_node.clone().write_structure(&mut file, 0, reference_names[chr].clone());
                }
            }
        }
    }
    
    // exit(1) to prevent the rest of the program from running for debugging purposes
    //std::process::exit(1);

    let mut counts = prepare_count_hashmap(&gtf);
    //let mut read_to_feature: Vec<FeatureType> = Vec::new();
    let mut counter = 0;

    count_reads(
        &mut reads_reader,
        &mut counter,
        &mut counts,
        &args,
        gtf,
        assignment_sender.as_ref(),
    );

    // Closing the sender lets the SAM writer terminate its receive loop.
    drop(assignment_sender);
    if let Some(writer_thread) = writer_thread {
        eprintln!("Waiting for writer thread to finish...");
        writer_thread.join().unwrap();
    }
    
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
            let reader = BamReader::from_path(path.clone(), n);
            match reader {
                Ok(reader) => ReadsReader::BamReader(reader),
                Err(e) => panic!("{}, File: {}", e, path),
            }
        } else if path.ends_with(".sam") {
            let reader = SamReader::from_path(path.clone());
            match reader {
                Ok(reader) => ReadsReader::SamReader(reader),
                Err(e) => panic!("{}, File: {}", e, path),
            }
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

// Parse command line arguments with clap's derive API
#[derive(Parser)]
struct Args {
    // Number of threads
    #[arg(
        short = 'n',
        long = "threads",
        default_value = "4",
        help = "Number of threads"
    )]
    n: u16,

    // Mode
    #[arg(short = 'm', long = "mode", default_value = "union", value_parser = ["intersection-strict", "intersection-nonempty", "union"], help = "Mode to use for counting reads overlapping features. Possible values: intersection-strict, intersection-nonempty, union (default: intersection-strict).")]
    _m: String,

    // Stranded
    #[arg(
        short = 's',
        long = "stranded",
        help = "Whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation",
        default_value = "yes",
        value_parser = ["yes", "no", "reverse"]
    )]
    stranded: String,

    // Quality filter
    #[arg(
        short = 'a',
        long = "minaqual",
        default_value = "10",
        help = "Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads."
    )]
    a: u8,

    // Type of feature to be used
    #[arg(
        short = 't',
        long = "type",
        default_value = "exon",
        num_args = 1,
        help = "Feature type (3rd column in GTF file) to be used. May be specified multiple times (default: exon)."
    )]
    t: Vec<String>,

    // Feature ID
    // TODO: implement actual logic for this option
    #[arg(
        short = 'i',
        long = "idattr",
        default_value = "gene_id",
        num_args = 1,
        help = "GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). All feature of the right type (see -t option) within the same GTF
    attribute will be added together. The typical way of using this option is to count all exonic reads from each gene and add the exons but other uses are possible
    as well. You can call this option multiple times: in that case, the combination of all attributes separated by colons (:) will be used as a unique identifier,
    e.g. for exons you might use -i gene_id -i exon_number."
    )]
    i: Vec<String>,

    // Name and type of the bam file
    #[arg(value_name = "bam")]
    bam: String,

    // Name and type of the gtf file
    #[arg(value_name = "gtf")]
    gtf: String,

    // Input ordering for paired-end data
    #[arg(
        short = 'r',
        long = "order",
        default_value = "name",
        value_parser = ["name", "pos"],
        help = "Sorting order for paired-end input: query-name grouped ('name') or coordinate sorted ('pos'). Ignored for single-end data."
    )]
    order: String,

    #[arg(
        long = "max-reads-in-buffer",
        default_value = "30000000",
        help = "Maximum number of unmatched mate keys buffered for coordinate-sorted paired-end input."
    )]
    max_buffer_size: usize,

    // Secondary alignment mode
    #[arg(long = "secondary-alignments", default_value = "ignore", value_parser = ["score", "ignore"], help = "Whether to score secondary alignments (0x100 flag)")]
    secondary_alignments: String,

    // Supplementary alignment mode
    #[arg(long = "supplementary-alignments", default_value = "ignore", value_parser = ["score", "ignore"], help = "Whether to score supplementary alignments (0x800 flag)")]
    supplementary_alignments: String,

    // Option to also output a total count of uniquely mapped reads
    #[arg(
        long = "extended-output",
        help = "Also output a total count of uniquely mapped reads"
    )]
    counts: bool,

    // Output delimiter
    #[arg(
        short = 'd',
        long = "delimiter",
        default_value = "\t",
        help = "Column delimiter in output (default: TAB)."
    )]
    delimiter: String,

    // non-unique parameter
    #[arg(
        long = "nonunique",
        default_value = "none",
        value_parser = ["none", "all", "fraction", "random"],
        help = "Whether and how to score reads that are not uniquely aligned or ambiguously assigned to features (choices: none, all, fraction, random; default: none)"
    )]
    nonunique: String,

    // Output file
    #[arg(
        short = 'c',
        long = "counts_output",
        help = "Filename to output the counts to instead of stdout."
    )]
    counts_output: Option<String>,

    // Export feature map
    #[arg(
        short = 'f',
        long = "export_feature_map",
        help = "Filename to output the feature map for debugging purposes."
    )]
    export_feature_tree: Option<String>,

    #[arg(
        short = 'o',
        long = "samout",
        help = "Create a SAM file with the reads and their features."
    )]
    output_sam: Option<String>,
}

fn prepare_count_hashmap(gtf: &Vec<Option<IntervalTree>>) -> HashMap<String, f64> {
    let mut counts: HashMap<String, f64> = HashMap::with_capacity(gtf.len());
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
    counts.insert("__no_feature".to_string(), 0f64);
    counts.insert("__ambiguous".to_string(), 0f64);
    counts.insert("__not_aligned".to_string(), 0f64);
    counts.insert("__too_low_aQual".to_string(), 0f64);
    counts.insert("__alignment_not_unique".to_string(), 0f64);
    counts
}

fn parse_feature_id(raw: &str, id_attributes: &[String], line_number: usize) -> String {
    // The common case is one ID attribute. Scan the attribute field directly
    // instead of allocating a HashMap<String, String> for every GTF row.
    if id_attributes.len() == 1 {
        let wanted = id_attributes[0].as_str();
        let mut found: Option<&str> = None;

        for raw_attr in raw.split(';') {
            let attr = raw_attr.trim();
            if attr.is_empty() {
                continue;
            }

            let mut parts = attr.splitn(2, |c: char| c.is_whitespace() || c == '=');
            let key = parts.next().unwrap_or("").trim();
            if key == wanted {
                found = Some(parts.next().unwrap_or("").trim().trim_matches('"'));
            }
        }

        return found
            .unwrap_or_else(|| {
                panic!(
                    "Feature on line {} does not contain a '{}' attribute",
                    line_number,
                    wanted
                )
            })
            .to_string();
    }

    let mut values: Vec<Option<&str>> = vec![None; id_attributes.len()];
    for raw_attr in raw.split(';') {
        let attr = raw_attr.trim();
        if attr.is_empty() {
            continue;
        }

        let mut parts = attr.splitn(2, |c: char| c.is_whitespace() || c == '=');
        let key = parts.next().unwrap_or("").trim();
        let value = parts.next().unwrap_or("").trim().trim_matches('"');

        for (index, wanted) in id_attributes.iter().enumerate() {
            if key == wanted {
                values[index] = Some(value);
            }
        }
    }

    let mut joined = String::new();
    for (index, wanted) in id_attributes.iter().enumerate() {
        let value = values[index].unwrap_or_else(|| {
            panic!(
                "Feature on line {} does not contain a '{}' attribute",
                line_number,
                wanted
            )
        });
        if index != 0 {
            joined.push(':');
        }
        joined.push_str(value);
    }
    joined
}

fn read_gtf(file_path: &str, feature_type_filter: &[String], ref_names_to_id: &HashMap<String, i32>, args: &Args) -> Vec<Option<IntervalTree>> {
    let mut map: HashMap<i32, Vec<Interval>> = HashMap::new();
    let file = File::open(file_path).expect("Could not open this file");
    let mut reader = BufReader::new(file);
    let mut counter = 0;
    let mut line = String::default();

    let mut chromosome_ids = ref_names_to_id.clone();
    let mut next_chr_id = chromosome_ids.len() as i32;

    while reader.read_line(&mut line).unwrap() > 0 {
        counter += 1;
        if counter % 100000 == 0 {
            eprintln!("{} GFF lines processed.", counter);
        }

        if line.starts_with('#') || line.trim().is_empty() {
            line.clear();
            continue;
        }

        let mut fields = line.trim_end_matches(['\r', '\n']).split('\t');
        let chr_name = fields.next();
        let _source = fields.next();
        let feature_type = fields.next();
        let start_field = fields.next();
        let end_field = fields.next();
        let _score = fields.next();
        let strand_field = fields.next();
        let _frame = fields.next();
        let attributes = fields.next();

        if chr_name.is_none()
            || feature_type.is_none()
            || start_field.is_none()
            || end_field.is_none()
            || strand_field.is_none()
            || attributes.is_none()
            || fields.next().is_some()
        {
            panic!(
                "Invalid GTF/GFF line {}: expected 9 tab-separated fields",
                counter
            );
        }

        let chr_name = chr_name.unwrap();
        let feature_type = feature_type.unwrap();
        if !feature_type_filter
            .iter()
            .any(|wanted_type| wanted_type == feature_type)
        {
            line.clear();
            continue;
        }

        let chr_id = match chromosome_ids.get(chr_name) {
            Some(id) => *id,
            None => {
                let id = next_chr_id;
                chromosome_ids.insert(chr_name.to_string(), id);
                next_chr_id += 1;
                id
            }
        };

        let start = start_field.unwrap().parse::<i32>().unwrap();
        let end = end_field.unwrap().parse::<i32>().unwrap();
        let strand = strand_field.unwrap().chars().next().unwrap_or('.');

        if args.stranded != "no" && strand != '+' && strand != '-' {
            panic!(
                "Feature on line {} has strand '{}', but stranded counting requires '+' or '-'",
                counter,
                strand
            );
        }

        let name = parse_feature_id(attributes.unwrap(), &args.i, counter);

        let feature = Feature::new(
            name,
            chr_id,
            min(start, end),
            max(start, end),
            strand,
        );
        map.entry(chr_id)
            .or_default()
            .push(Interval::new(min(start, end), max(start, end), Some(feature)));
        line.clear();
    }

    eprintln!("{} GFF lines processed.", counter);
    eprint!("Creating IntervalTree for each chromosome...");

    let mut result: Vec<Option<IntervalTree>> = Vec::with_capacity(chromosome_ids.len());
    for _ in 0..chromosome_ids.len() {
        result.push(None);
    }

    for (chr, intervals) in map {
        result[chr as usize] = Some(IntervalTree::new(Some(intervals)));
    }
    eprintln!("done.");
    result
}

fn should_skip_record(
    record: &bam::Record,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) -> bool {
    if !record.flag().is_mapped() || record.ref_id() < 0 {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::NotAligned); }
        *counts.get_mut("__not_aligned").unwrap() += 1.0;
        return true;
    }

    if args.secondary_alignments == "ignore" && record.flag().all_bits(0x100) {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::None); }
        return true;
    }

    if args.supplementary_alignments == "ignore" && record.flag().all_bits(0x800) {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::None); }
        return true;
    }

    if let Some(TagValue::Int(i, _)) = record.tags().get(b"NH") {
        if i > 1 {
            *counts.get_mut("__alignment_not_unique").unwrap() += 1.0;
            if args.nonunique == "none" {
                if let Some(sender) = sender { let _ = sender.send(FeatureType::AlignmentNotUnique); }
                return true;
            }
        }
    }

    if record.mapq() < args.a {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::TooLowaQual); }
        *counts.get_mut("__too_low_aQual").unwrap() += 1.0;
        return true;
    }

    false
}

fn print_output(counts: HashMap<String, f64>, args: Args, counter: i32) {
    // Print de HashMap
    let mut sorted_keys: Vec<_> = counts.keys().collect();
    // Sort the keys case-insensitively
    sorted_keys.sort();
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

    // TODO: check the correctness, since it might depend on nonunique mode
    if args.counts {
        println!(
            "Total number of uniquely mapped reads{}{}",
            args.delimiter,
            counter as f64
                - counts["__not_aligned"]
                - counts["__too_low_aQual"]
                - counts["__alignment_not_unique"]
                - counts["__ambiguous"]
                - counts["__no_feature"]
        );
    }
}

fn write_counts(counts: HashMap<String, f64>, args: Args, counter: i32) {
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
            counter as f64 - counts["__not_aligned"] - counts["__too_low_aQual"] - counts["__alignment_not_unique"] - counts["__ambiguous"] - counts["__no_feature"])
            .as_bytes(),
        )
        .expect("Unable to write data");
    }
}


fn processing_function_for_mode(
    mode: &str,
) -> fn(&IntervalTree, i32, i32, bool, &mut Vec<Vec<Feature>>, &Args) {
    match mode {
        "intersection-strict" => process_intersection_strict_read,
        "intersection-nonempty" => process_intersection_nonempty_read,
        "union" => process_union_read,
        _ => panic!("Invalid mode"),
    }
}

fn ambiguity_function_for_mode(mode: &str) -> fn(&[Vec<Feature>]) -> Vec<String> {
    match mode {
        "intersection-strict" => filter_ambiguity_intersection_strict,
        "intersection-nonempty" => filter_ambiguity_intersection_nonempty,
        "union" => filter_ambiguity_union,
        _ => panic!("Invalid mode"),
    }
}

fn add_record_blocks(
    record: &bam::Record,
    invert_pair_strand: bool,
    gtf: &[Option<IntervalTree>],
    overlapping_features: &mut Vec<Vec<Feature>>,
    args: &Args,
) -> bool {
    if !record.flag().is_mapped() || record.ref_id() < 0 {
        return true;
    }

    let ref_id = record.ref_id() as usize;
    let features = match gtf.get(ref_id).and_then(|tree| tree.as_ref()) {
        Some(features) => features,
        None => return false,
    };

    let processing_function = processing_function_for_mode(args._m.as_str());
    let mut reference_pos = record.start() + 1;
    let record_reverse = record.flag().is_reverse_strand();
    let effective_reverse = if invert_pair_strand {
        !record_reverse
    } else {
        record_reverse
    };

    for cig in record.cigar().iter() {
        let length = cig.0 as i32;
        let operation = cig.1;

        if operation.is_match() && length > 0 {
            let block_end = reference_pos + length - 1;
            processing_function(
                features,
                reference_pos,
                block_end,
                effective_reverse,
                overlapping_features,
                args,
            );
        }

        if operation.consumes_ref() {
            reference_pos += length;
        }
    }

    true
}

fn assign_overlaps(
    overlapping_features: &[Vec<Feature>],
    counts: &mut HashMap<String, f64>,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    let ambiguity_function = ambiguity_function_for_mode(args._m.as_str());
    let mut unique_feature_names = ambiguity_function(overlapping_features);
    let feature_name_len = unique_feature_names.len();

    match feature_name_len {
        0 => {
            *counts.get_mut("__no_feature").unwrap() += 1.0;
            if let Some(sender) = sender { let _ = sender.send(FeatureType::NoFeature); }
        }
        1 => {
            let feature_name = unique_feature_names.first().unwrap();
            *counts.get_mut(feature_name.as_str()).unwrap() += 1.0;
            if let Some(sender) = sender { let _ = sender.send(FeatureType::Name(feature_name.clone())); }
        }
        _ => {
            *counts.get_mut("__ambiguous").unwrap() += 1.0;
            match args.nonunique.as_str() {
                "all" => {
                    for feature_name in &unique_feature_names {
                        *counts.get_mut(feature_name.as_str()).unwrap() += 1.0;
                    }
                }
                "fraction" => {
                    let fractional_count = 1.0 / feature_name_len as f64;
                    for feature_name in &unique_feature_names {
                        *counts.entry(feature_name.clone()).or_insert(0.0) += fractional_count;
                    }
                }
                "random" => {
                    let random_index = rand::random_range(0..feature_name_len);
                    let feature_name = unique_feature_names[random_index].clone();
                    *counts.get_mut(feature_name.as_str()).unwrap() += 1.0;
                }
                _ => {}
            }

            unique_feature_names.sort();
            if let Some(sender) = sender { let _ = sender.send(FeatureType::Ambiguous(unique_feature_names.join("+"))); }
        }
    }
}

fn count_single_record(
    record: &bam::Record,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: &[Option<IntervalTree>],
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    *counter += 1;
    if *counter % 100000 == 0 {
        eprintln!("{} records processed.", counter);
    }

    if should_skip_record(record, counts, args, sender) {
        return;
    }

    let mut overlapping_features = Vec::with_capacity(3);
    if !add_record_blocks(record, false, gtf, &mut overlapping_features, args) {
        *counts.get_mut("__no_feature").unwrap() += 1.0;
        if let Some(sender) = sender { let _ = sender.send(FeatureType::NoFeature); }
        return;
    }

    assign_overlaps(&overlapping_features, counts, args, sender);
}

fn pair_side(record: &bam::Record) -> u8 {
    match (record.flag().first_in_pair(), record.flag().last_in_pair()) {
        (true, false) => 1,
        (false, true) => 2,
        _ => panic!(
            "Paired-end alignment '{}' must have exactly one of the READ1/READ2 flags set",
            String::from_utf8_lossy(record.name())
        ),
    }
}

fn primary_only_pairing(args: &Args) -> bool {
    args.secondary_alignments == "ignore" && args.supplementary_alignments == "ignore"
}

fn should_pre_filter_pair_record(record: &bam::Record, args: &Args) -> bool {
    primary_only_pairing(args)
        && (record.flag().is_secondary() || record.flag().is_supplementary())
}

fn nh_multimap_status(record: &bam::Record) -> Option<bool> {
    match record.tags().get(b"NH") {
        Some(TagValue::Int(value, _)) => Some(value > 1),
        _ => None,
    }
}

fn pair_is_multimapped_htseq_compatible(
    first: Option<&bam::Record>,
    second: Option<&bam::Record>,
) -> bool {
    // HTSeq 2.1.2 checks both NH tags inside one try/except block. If mate 1
    // exists but lacks NH, its KeyError exits the block before mate 2 is
    // inspected. Preserve that behavior for count compatibility.
    if let Some(record) = first {
        match nh_multimap_status(record) {
            Some(true) => return true,
            Some(false) => {}
            None => return false,
        }
    }

    second
        .and_then(nh_multimap_status)
        .unwrap_or(false)
}

fn should_skip_pair(
    first: Option<&bam::Record>,
    second: Option<&bam::Record>,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) -> bool {
    let first_mapped = first.map(|r| r.flag().is_mapped()).unwrap_or(false);
    let second_mapped = second.map(|r| r.flag().is_mapped()).unwrap_or(false);

    if !first_mapped && !second_mapped {
        *counts.get_mut("__not_aligned").unwrap() += 1.0;
        if let Some(sender) = sender { let _ = sender.send(FeatureType::NotAligned); }
        return true;
    }

    if args.secondary_alignments == "ignore"
        && (first.map(|r| r.flag().is_secondary()).unwrap_or(false)
            || second.map(|r| r.flag().is_secondary()).unwrap_or(false))
    {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::None); }
        return true;
    }

    if args.supplementary_alignments == "ignore"
        && (first.map(|r| r.flag().is_supplementary()).unwrap_or(false)
            || second.map(|r| r.flag().is_supplementary()).unwrap_or(false))
    {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::None); }
        return true;
    }

    let multimapped = pair_is_multimapped_htseq_compatible(first, second);
    if multimapped {
        *counts
            .entry("__alignment_not_unique".to_string())
            .or_insert(0.0) += 1.0;
        if args.nonunique == "none" {
            if let Some(sender) = sender { let _ = sender.send(FeatureType::AlignmentNotUnique); }
            return true;
        }
    }

    let low_quality = first.map(|r| r.mapq() < args.a).unwrap_or(false)
        || second.map(|r| r.mapq() < args.a).unwrap_or(false);
    if low_quality {
        *counts.get_mut("__too_low_aQual").unwrap() += 1.0;
        if let Some(sender) = sender { let _ = sender.send(FeatureType::TooLowaQual); }
        return true;
    }

    false
}

fn count_pair(
    first: Option<&bam::Record>,
    second: Option<&bam::Record>,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: &[Option<IntervalTree>],
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    *counter += 1;
    if *counter % 100000 == 0 {
        eprintln!("{} read pairs processed.", counter);
    }

    if should_skip_pair(first, second, counts, args, sender) {
        return;
    }

    let mut overlapping_features = Vec::with_capacity(6);
    let mut all_chromosomes_known = true;

    if let Some(record) = first {
        all_chromosomes_known &= add_record_blocks(
            record,
            false,
            gtf,
            &mut overlapping_features,
            args,
        );
    }
    if let Some(record) = second {
        // HTSeq inverts mate 2 before applying strandedness. With
        // --stranded reverse the later strand comparison inverts it back.
        all_chromosomes_known &= add_record_blocks(
            record,
            true,
            gtf,
            &mut overlapping_features,
            args,
        );
    }

    if !all_chromosomes_known {
        *counts.get_mut("__no_feature").unwrap() += 1.0;
        if let Some(sender) = sender { let _ = sender.send(FeatureType::NoFeature); }
        return;
    }

    assign_overlaps(&overlapping_features, counts, args, sender);
}

fn records_are_mates_name_sorted(first: &bam::Record, second: &bam::Record) -> bool {
    if pair_side(first) == pair_side(second) {
        return false;
    }

    let first_aligned = first.flag().is_mapped();
    let second_aligned = second.flag().is_mapped();
    let first_mate_aligned = first.flag().mate_is_mapped();
    let second_mate_aligned = second.flag().mate_is_mapped();

    if first_aligned != second_mate_aligned || first_mate_aligned != second_aligned {
        return false;
    }

    if !(first_aligned && second_aligned) {
        return true;
    }

    first.ref_id() == second.mate_ref_id()
        && first.start() == second.mate_start()
        && second.ref_id() == first.mate_ref_id()
        && second.start() == first.mate_start()
}

fn process_name_group(
    mut group: VecDeque<bam::Record>,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: &[Option<IntervalTree>],
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    while let Some(record) = group.pop_front() {
        let mate_index = group
            .iter()
            .position(|candidate| records_are_mates_name_sorted(&record, candidate));
        let mate = mate_index.and_then(|index| group.remove(index));

        if pair_side(&record) == 1 {
            count_pair(
                Some(&record),
                mate.as_ref(),
                counter,
                counts,
                args,
                gtf,
                sender,
            );
        } else {
            count_pair(
                mate.as_ref(),
                Some(&record),
                counter,
                counts,
                args,
                gtf,
                sender,
            );
        }
    }
}

fn count_paired_name_sorted(
    reads_reader: &mut ReadsReader,
    first_record: bam::Record,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: &[Option<IntervalTree>],
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    let mut current_name: Option<Vec<u8>> = None;
    let mut group = VecDeque::new();

    let records = std::iter::once(Ok(first_record)).chain(reads_reader.by_ref());
    for result in records {
        let record = result.unwrap_or_else(|e| panic!("{}", e));
        if !record.flag().is_paired() {
            panic!("Mixed single-end and paired-end records are not supported");
        }
        pair_side(&record);

        if should_pre_filter_pair_record(&record, args) {
            continue;
        }

        let name_changed = current_name
            .as_ref()
            .map(|current| current.as_slice() != record.name())
            .unwrap_or(false);

        if name_changed {
            process_name_group(group, counter, counts, args, gtf, sender);
            group = VecDeque::new();
            current_name = Some(record.name().to_vec());
        } else if current_name.is_none() {
            current_name = Some(record.name().to_vec());
        }

        group.push_back(record);
    }

    if !group.is_empty() {
        process_name_group(group, counter, counts, args, gtf, sender);
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
struct MateKey {
    name: Vec<u8>,
    which: u8,
    ref_id: Option<i32>,
    start: Option<i32>,
    mate_ref_id: Option<i32>,
    mate_start: Option<i32>,
    template_len: Option<i32>,
}

fn own_key_from_expected_mate_key(mut key: MateKey) -> MateKey {
    key.which = if key.which == 1 { 2 } else { 1 };
    std::mem::swap(&mut key.ref_id, &mut key.mate_ref_id);
    std::mem::swap(&mut key.start, &mut key.mate_start);
    key.template_len = key.template_len.map(|value| -value);
    key
}

fn expected_mate_key(record: &bam::Record) -> MateKey {
    let aligned = record.flag().is_mapped();
    let mate_aligned = record.flag().mate_is_mapped();
    MateKey {
        name: record.name().to_vec(),
        which: if pair_side(record) == 1 { 2 } else { 1 },
        ref_id: mate_aligned.then_some(record.mate_ref_id()),
        start: mate_aligned.then_some(record.mate_start()),
        mate_ref_id: aligned.then_some(record.ref_id()),
        mate_start: aligned.then_some(record.start()),
        template_len: (aligned && mate_aligned).then_some(-record.template_len()),
    }
}

fn count_paired_position_sorted(
    reads_reader: &mut ReadsReader,
    first_record: bam::Record,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: &[Option<IntervalTree>],
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    let mut buffer: HashMap<MateKey, VecDeque<bam::Record>> = HashMap::new();

    let records = std::iter::once(Ok(first_record)).chain(reads_reader.by_ref());
    for result in records {
        let record = result.unwrap_or_else(|e| panic!("{}", e));
        if !record.flag().is_paired() {
            panic!("Mixed single-end and paired-end records are not supported");
        }
        let side = pair_side(&record);

        if should_pre_filter_pair_record(&record, args) {
            continue;
        }

        let mate_key = expected_mate_key(&record);
        let mate = if let Some(queue) = buffer.get_mut(&mate_key) {
            let mate = queue.pop_front();
            if queue.is_empty() {
                buffer.remove(&mate_key);
            }
            mate
        } else {
            None
        };

        if let Some(mate) = mate {
            if side == 1 {
                count_pair(
                    Some(&record),
                    Some(&mate),
                    counter,
                    counts,
                    args,
                    gtf,
                    sender,
                );
            } else {
                count_pair(
                    Some(&mate),
                    Some(&record),
                    counter,
                    counts,
                    args,
                    gtf,
                    sender,
                );
            }
        } else {
            // Reuse the QNAME allocation from the failed mate lookup instead of
            // allocating it a second time for the record's own buffer key.
            let record_key = own_key_from_expected_mate_key(mate_key);
            buffer
                .entry(record_key)
                .or_default()
                .push_back(record);
            if buffer.len() > args.max_buffer_size {
                panic!(
                    "Maximum paired-end alignment buffer size exceeded ({} mate keys). \
                     Check that --order pos matches the input, or increase --max-reads-in-buffer.",
                    args.max_buffer_size
                );
            }
        }
    }

    for (_, mut queue) in buffer {
        while let Some(record) = queue.pop_front() {
            if pair_side(&record) == 1 {
                count_pair(
                    Some(&record),
                    None,
                    counter,
                    counts,
                    args,
                    gtf,
                    sender,
                );
            } else {
                count_pair(
                    None,
                    Some(&record),
                    counter,
                    counts,
                    args,
                    gtf,
                    sender,
                );
            }
        }
    }
}

fn count_reads(
    reads_reader: &mut ReadsReader,
    counter: &mut i32,
    counts: &mut HashMap<String, f64>,
    args: &Args,
    gtf: Vec<Option<IntervalTree>>,
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    // RecordReader::read_into reuses the record's internal buffers. The bam
    // crate specifically exposes this path to avoid allocating a new Record
    // for every alignment.
    let mut record = bam::Record::new();
    match reads_reader.read_into(&mut record) {
        Ok(true) => {}
        Ok(false) => {
            eprintln!("0 records processed.");
            return;
        }
        Err(e) => panic!("{}", e),
    }

    if record.flag().is_paired() {
        if args.output_sam.is_some() {
            panic!(
                "--samout is not yet supported for paired-end input. \
                 Counting is supported; SAM annotation needs an order-aware writer."
            );
        }

        match args.order.as_str() {
            "name" => count_paired_name_sorted(
                reads_reader,
                record,
                counter,
                counts,
                args,
                &gtf,
                sender,
            ),
            "pos" => count_paired_position_sorted(
                reads_reader,
                record,
                counter,
                counts,
                args,
                &gtf,
                sender,
            ),
            _ => unreachable!(),
        }
        eprintln!("{} read pairs processed.", counter);
    } else {
        loop {
            if record.flag().is_paired() {
                panic!("Mixed single-end and paired-end records are not supported");
            }

            count_single_record(
                &record,
                counter,
                counts,
                args,
                &gtf,
                sender,
            );

            match reads_reader.read_into(&mut record) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => panic!("{}", e),
            }
        }
        eprintln!("{} records processed.", counter);
    }
}

fn process_union_read(features: &IntervalTree, start_pos: i32, end_pos: i32, strand: bool, overlapping_features: &mut Vec<Vec<Feature>>, args: &Args) {
    let new_overlap = features.overlap(start_pos, end_pos);
    let strand = if strand { '-' } else { '+' };
    // add all overlapping features to the list
    add_stranded_features(new_overlap, strand, overlapping_features, args);
    
}

fn process_intersection_nonempty_read(
    features: &IntervalTree,
    start_pos: i32,
    end_pos: i32,
    strand: bool,
    overlapping_features: &mut Vec<Vec<Feature>>,
    args: &Args,
) {
    let overlaps = features.overlap(start_pos, end_pos);
    let read_strand = if strand { '-' } else { '+' };

    let relevant: Vec<&Interval> = overlaps
        .into_iter()
        .filter(|interval| {
            let feature = interval.data.as_ref().unwrap();
            feature_matches_strand(feature, read_strand, args)
        })
        .collect();

    if relevant.is_empty() {
        return;
    }

    let mut boundaries = vec![start_pos, end_pos + 1];
    for interval in &relevant {
        boundaries.push(max(start_pos, interval.start));
        boundaries.push(min(end_pos, interval.end) + 1);
    }
    boundaries.sort_unstable();
    boundaries.dedup();

    for window in boundaries.windows(2) {
        let step_start = window[0];
        let step_end_exclusive = window[1];
        if step_start >= step_end_exclusive || step_start > end_pos {
            continue;
        }

        let mut step_features = Vec::new();
        for interval in &relevant {
            if interval.start <= step_start && interval.end >= step_start {
                step_features.push(interval.data.as_ref().unwrap().clone());
            }
        }

        if !step_features.is_empty() {
            overlapping_features.push(step_features);
        }
    }
}

fn process_intersection_strict_read(features: &IntervalTree, start_pos: i32, end_pos: i32, strand: bool, overlapping_features: &mut Vec<Vec<Feature>>, args: &Args) {
    //todo!("process_partial_read for intersection-strict");
    // Problem now: if we have partial reads: each part must overlap with the same feature, if different parts overlap with different features, we should not count the read as feature but as no_feature
    let new_contained = features.contains(start_pos, end_pos);
    // change new_contained to a Vec<&Interval>
    let strand = if strand { '-' } else { '+' };
    // add all contained features to the list
    add_stranded_features(new_contained, strand, overlapping_features, args);
}

fn filter_ambiguity_union(
    overlapping_features: &[Vec<Feature>],
) -> Vec<String> {
    // flatten the Vec of Vecs to a single Vec
    let overlapping_features: Vec<&Feature> = overlapping_features.iter().flatten().collect();
    let unique_feature_names: HashSet<String> = overlapping_features.iter().map(|x| x.name().to_string().clone()).filter(|x| !x.is_empty()).collect();
    unique_feature_names.into_iter().collect()
}

fn filter_ambiguity_intersection_strict(
    overlapping_features: &[Vec<Feature>],
) -> Vec<String> {
    // if any of the results is empty, we have no feature, so we return an empty Vec
    if overlapping_features.iter().any(|x| x.is_empty()) {
        return Vec::new()
    
    // otherwise, we return the unique feature names
    } else {
        let mut feature_counts: HashMap<String, usize> = HashMap::new();
        let total = overlapping_features.len();

        for features in overlapping_features.iter() {
            let feature_names: HashSet<String> = features.iter().map(|x| x.name().to_string()).collect();
            for name in feature_names {
                *feature_counts.entry(name).or_insert(0) += 1;
            }
        }

        let unique_feature_names: Vec<String> = feature_counts.into_iter()
            .filter(|&(_, count)| count == total)
            .map(|(name, _)| name)
            .filter(|name| !name.is_empty())
            .collect();

        unique_feature_names
    }
}


fn filter_ambiguity_intersection_nonempty(
    overlapping_features: &[Vec<Feature>],
) -> Vec<String> {
    let mut nonempty_steps = overlapping_features.iter().filter(|features| !features.is_empty());

    let first = match nonempty_steps.next() {
        Some(features) => features,
        None => return Vec::new(),
    };

    let mut feature_names: HashSet<String> =
        first.iter().map(|feature| feature.name().to_string()).collect();

    for features in nonempty_steps {
        let current: HashSet<String> =
            features.iter().map(|feature| feature.name().to_string()).collect();
        feature_names = feature_names.intersection(&current).cloned().collect();

        if feature_names.is_empty() {
            return Vec::new();
        }
    }

    feature_names
        .into_iter()
        .filter(|name| !name.is_empty())
        .collect()
}

fn feature_matches_strand(feature: &Feature, strand: char, args: &Args) -> bool {
    match args.stranded.as_str() {
        "yes" => feature.strand() == strand,
        "reverse" => feature.strand() != strand,
        "no" => true,
        _ => panic!("Invalid strandedness"),
    }
}


fn add_stranded_features(new_overlap: Vec<&Interval>, strand: char, overlapping_features: &mut Vec<Vec<Feature>>, args: &Args) {
    // add empty Vec to the overlapping_features list
    overlapping_features.push(Vec::new());
    let index = overlapping_features.len() - 1;
    if new_overlap.len() > 0 {
        for overlap in new_overlap {
            let feature = overlap.data.as_ref().unwrap();
            if feature_matches_strand(feature, strand, args) {
                overlapping_features[index].push(feature.clone());
            }
        }
    }
}
