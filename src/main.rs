use bam::record::tags::TagValue;
use bam::{RecordReader,BamReader, RecordWriter, SamReader, SamWriter};
use feature::Feature;
use intervaltree::IntervalTree;
use interval::Interval;
use std::borrow::Cow;
use std::cmp::{max, min};
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::Instant;
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

    // Read and index the annotation.
    let annotation_started = Instant::now();
    let (gtf, feature_names) = read_gtf(&args.gtf, &args.t, &ref_names_to_id, &args);
    if std::env::var_os("HTSEQ_COUNT_RUST_TIMINGS").is_some() {
        eprintln!(
            "__timing_annotation_total_seconds\t{:.6}",
            annotation_started.elapsed().as_secs_f64()
        );
    }

    // let read= 21940455;
    // eprintln!("Searching for reads overlapping position {}-{}...", read, read+25);
    // for overlap in gtf["1"].overlap(read, read+25) {
    //      eprintln!("overlap: {:?}", overlap);
    // }

    
    // DOT export is deliberately opt-in. Normal counting does no export I/O
    // and does not clone the interval tree.
    if let Some(export_path) = args.export_feature_tree.as_deref() {
        eprintln!("Exporting feature tree to {}...", export_path);
        let mut file = File::create(export_path).expect("Unable to create feature-tree DOT file");
        for (chr, tree) in gtf.iter().enumerate().take(reference_names.len()) {
            if let Some(tree) = tree {
                if let Some(top_node) = &tree.top_node {
                    let _ = top_node.write_structure(&mut file, 0, &reference_names[chr]);
                }
            }
        }
    }
    
    // exit(1) to prevent the rest of the program from running for debugging purposes
    //std::process::exit(1);

    let mut counts = Counts::new(feature_names);
    //let mut read_to_feature: Vec<FeatureType> = Vec::new();
    let mut counter = 0;

    let counting_started = Instant::now();
    count_reads(
        &mut reads_reader,
        &mut counter,
        &mut counts,
        &args,
        gtf,
        assignment_sender.as_ref(),
    );
    if std::env::var_os("HTSEQ_COUNT_RUST_TIMINGS").is_some() {
        eprintln!(
            "__timing_counting_seconds\t{:.6}",
            counting_started.elapsed().as_secs_f64()
        );
    }

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
        long = "export-feature-tree",
        visible_alias = "export_feature_map",
        value_name = "DOT",
        help = "Optional: export the annotation feature tree as a Graphviz DOT file. Disabled by default because export can be slow on large annotations."
    )]
    export_feature_tree: Option<String>,

    #[arg(
        short = 'o',
        long = "samout",
        help = "Create a SAM file with the reads and their features."
    )]
    output_sam: Option<String>,
}

struct Counts {
    feature_names: Vec<Arc<str>>,
    feature_counts: Vec<f64>,
    no_feature: f64,
    ambiguous: f64,
    not_aligned: f64,
    too_low_aqual: f64,
    alignment_not_unique: f64,
}

impl Counts {
    fn new(feature_names: Vec<Arc<str>>) -> Self {
        let feature_counts = vec![0.0; feature_names.len()];
        Counts {
            feature_names,
            feature_counts,
            no_feature: 0.0,
            ambiguous: 0.0,
            not_aligned: 0.0,
            too_low_aqual: 0.0,
            alignment_not_unique: 0.0,
        }
    }

    #[inline]
    fn add_feature(&mut self, feature_id: usize, value: f64) {
        self.feature_counts[feature_id] += value;
    }

    fn special_total(&self) -> f64 {
        self.no_feature
            + self.ambiguous
            + self.not_aligned
            + self.too_low_aqual
            + self.alignment_not_unique
    }
}

fn trim_ascii(mut value: &[u8]) -> &[u8] {
    while value.first().map(|b| b.is_ascii_whitespace()).unwrap_or(false) {
        value = &value[1..];
    }
    while value.last().map(|b| b.is_ascii_whitespace()).unwrap_or(false) {
        value = &value[..value.len() - 1];
    }
    value
}

#[inline]
fn parse_i32_ascii(value: &[u8]) -> i32 {
    let mut result = 0i32;
    for &byte in value {
        if !byte.is_ascii_digit() {
            panic!("Invalid integer field in GTF/GFF");
        }
        result = result
            .checked_mul(10)
            .and_then(|v| v.checked_add((byte - b'0') as i32))
            .expect("GTF/GFF coordinate exceeds i32 range");
    }
    result
}

fn parse_feature_id_bytes<'a>(
    raw: &'a [u8],
    id_attributes: &[String],
    line_number: usize,
) -> Cow<'a, str> {
    if id_attributes.len() == 1 {
        let wanted = id_attributes[0].as_bytes();

        for raw_attr in raw.split(|byte| *byte == b';') {
            let attr = trim_ascii(raw_attr);
            if attr.is_empty() {
                continue;
            }

            let separator = attr
                .iter()
                .position(|byte| byte.is_ascii_whitespace() || *byte == b'=');
            let Some(separator) = separator else {
                continue;
            };

            if &attr[..separator] == wanted {
                let mut value = trim_ascii(&attr[separator + 1..]);
                if value.len() >= 2 && value[0] == b'"' && value[value.len() - 1] == b'"' {
                    value = &value[1..value.len() - 1];
                }
                return Cow::Borrowed(
                    std::str::from_utf8(value)
                        .expect("Feature ID is not valid UTF-8")
                );
            }
        }

        panic!(
            "Feature on line {} does not contain a '{}' attribute",
            line_number,
            id_attributes[0]
        );
    }

    let mut values: Vec<Option<&[u8]>> = vec![None; id_attributes.len()];
    for raw_attr in raw.split(|byte| *byte == b';') {
        let attr = trim_ascii(raw_attr);
        if attr.is_empty() {
            continue;
        }

        let separator = attr
            .iter()
            .position(|byte| byte.is_ascii_whitespace() || *byte == b'=');
        let Some(separator) = separator else {
            continue;
        };
        let key = &attr[..separator];
        let mut value = trim_ascii(&attr[separator + 1..]);
        if value.len() >= 2 && value[0] == b'"' && value[value.len() - 1] == b'"' {
            value = &value[1..value.len() - 1];
        }

        for (index, wanted) in id_attributes.iter().enumerate() {
            if key == wanted.as_bytes() {
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
        joined.push_str(
            std::str::from_utf8(value).expect("Feature ID is not valid UTF-8")
        );
    }
    Cow::Owned(joined)
}

fn read_gtf(
    file_path: &str,
    feature_type_filter: &[String],
    ref_names_to_id: &HashMap<String, i32>,
    args: &Args,
) -> (Vec<Option<IntervalTree>>, Vec<Arc<str>>) {
    let profile_timings = std::env::var_os("HTSEQ_COUNT_RUST_TIMINGS").is_some();
    let parse_started = Instant::now();
    let mut map: HashMap<i32, Vec<Interval>> = HashMap::new();
    let file = File::open(file_path).expect("Could not open this file");
    // GTF/GFF files are large sequential text streams. A larger buffer reduces
    // read syscalls, while reusing a preallocated line avoids early growth.
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut counter = 0;
    let mut line: Vec<u8> = Vec::with_capacity(512);

    let mut chromosome_ids = ref_names_to_id.clone();
    let mut next_chr_id = chromosome_ids.len() as i32;
    let mut feature_ids: HashMap<String, usize> = HashMap::new();
    let mut feature_names: Vec<Arc<str>> = Vec::new();

    while reader.read_until(b'\n', &mut line).unwrap() > 0 {
        counter += 1;
        if counter % 100000 == 0 {
            eprintln!("{} GFF lines processed.", counter);
        }

        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }

        if line.is_empty() || line[0] == b'#' {
            line.clear();
            continue;
        }

        let mut fields = line.split(|byte| *byte == b'\t');
        let chr_field = fields.next();
        let _source = fields.next();
        let feature_type = fields.next();
        let start_field = fields.next();
        let end_field = fields.next();
        let _score = fields.next();
        let strand_field = fields.next();
        let _frame = fields.next();
        let attributes = fields.next();

        if chr_field.is_none()
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

        let feature_type = feature_type.unwrap();
        if !feature_type_filter
            .iter()
            .any(|wanted_type| wanted_type.as_bytes() == feature_type)
        {
            line.clear();
            continue;
        }

        let chr_name = std::str::from_utf8(chr_field.unwrap())
            .expect("GTF/GFF chromosome name is not valid UTF-8");
        let chr_id = match chromosome_ids.get(chr_name) {
            Some(id) => *id,
            None => {
                let id = next_chr_id;
                chromosome_ids.insert(chr_name.to_string(), id);
                next_chr_id += 1;
                id
            }
        };

        let start = parse_i32_ascii(start_field.unwrap());
        let end = parse_i32_ascii(end_field.unwrap());
        let strand = strand_field
            .unwrap()
            .first()
            .copied()
            .unwrap_or(b'.') as char;

        if args.stranded != "no" && strand != '+' && strand != '-' {
            panic!(
                "Feature on line {} has strand '{}', but stranded counting requires '+' or '-'",
                counter,
                strand
            );
        }

        let parsed_name =
            parse_feature_id_bytes(attributes.unwrap(), &args.i, counter);

        let (feature_id, feature_name) =
            if let Some(id) = feature_ids.get(parsed_name.as_ref()) {
                (*id, feature_names[*id].clone())
            } else {
                let owned_name = parsed_name.into_owned();
                let id = feature_names.len();
                let shared_name: Arc<str> = Arc::from(owned_name.as_str());
                feature_ids.insert(owned_name, id);
                feature_names.push(shared_name.clone());
                (id, shared_name)
            };

        let feature = Feature::new(
            feature_id,
            feature_name,
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
    if profile_timings {
        eprintln!(
            "__timing_gtf_parse_seconds\t{:.6}",
            parse_started.elapsed().as_secs_f64()
        );
    }

    eprint!("Creating IntervalTree for each chromosome...");
    let index_started = Instant::now();

    let mut result: Vec<Option<IntervalTree>> = Vec::with_capacity(chromosome_ids.len());
    for _ in 0..chromosome_ids.len() {
        result.push(None);
    }

    for (chr, intervals) in map {
        result[chr as usize] = Some(IntervalTree::new(Some(intervals)));
    }
    eprintln!("done.");
    if profile_timings {
        eprintln!(
            "__timing_index_build_seconds\t{:.6}",
            index_started.elapsed().as_secs_f64()
        );
    }
    (result, feature_names)
}

fn should_skip_record(
    record: &bam::Record,
    counts: &mut Counts,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) -> bool {
    if !record.flag().is_mapped() || record.ref_id() < 0 {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::NotAligned); }
        counts.not_aligned += 1.0;
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
            counts.alignment_not_unique += 1.0;
            if args.nonunique == "none" {
                if let Some(sender) = sender { let _ = sender.send(FeatureType::AlignmentNotUnique); }
                return true;
            }
        }
    }

    if record.mapq() < args.a {
        if let Some(sender) = sender { let _ = sender.send(FeatureType::TooLowaQual); }
        counts.too_low_aqual += 1.0;
        return true;
    }

    false
}

fn sorted_feature_indices(counts: &Counts) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..counts.feature_names.len()).collect();
    indices.sort_unstable_by(|a, b| {
        counts.feature_names[*a]
            .as_ref()
            .cmp(counts.feature_names[*b].as_ref())
    });
    indices
}

fn print_output(counts: Counts, args: Args, counter: i32) {
    for feature_id in sorted_feature_indices(&counts) {
        println!(
            "{}{}{}",
            counts.feature_names[feature_id],
            args.delimiter,
            counts.feature_counts[feature_id]
        );
    }

    println!("__no_feature{}{}", args.delimiter, counts.no_feature);
    println!("__ambiguous{}{}", args.delimiter, counts.ambiguous);
    println!("__too_low_aQual{}{}", args.delimiter, counts.too_low_aqual);
    println!("__not_aligned{}{}", args.delimiter, counts.not_aligned);
    println!(
        "__alignment_not_unique{}{}",
        args.delimiter,
        counts.alignment_not_unique
    );

    if args.counts {
        println!(
            "Total number of uniquely mapped reads{}{}",
            args.delimiter,
            counter as f64 - counts.special_total()
        );
    }
}

fn write_counts(counts: Counts, args: Args, counter: i32) {
    let mut file = File::create(args.counts_output.unwrap())
        .expect("Unable to create file");

    for feature_id in sorted_feature_indices(&counts) {
        writeln!(
            file,
            "{}{}{}",
            counts.feature_names[feature_id],
            args.delimiter,
            counts.feature_counts[feature_id]
        )
        .expect("Unable to write data");
    }

    writeln!(file, "__no_feature{}{}", args.delimiter, counts.no_feature)
        .expect("Unable to write data");
    writeln!(file, "__ambiguous{}{}", args.delimiter, counts.ambiguous)
        .expect("Unable to write data");
    writeln!(
        file,
        "__too_low_aQual{}{}",
        args.delimiter,
        counts.too_low_aqual
    )
    .expect("Unable to write data");
    writeln!(file, "__not_aligned{}{}", args.delimiter, counts.not_aligned)
        .expect("Unable to write data");
    writeln!(
        file,
        "__alignment_not_unique{}{}",
        args.delimiter,
        counts.alignment_not_unique
    )
    .expect("Unable to write data");

    if args.counts {
        writeln!(
            file,
            "Total number of uniquely mapped reads{}{}",
            args.delimiter,
            counter as f64 - counts.special_total()
        )
        .expect("Unable to write data");
    }
}

fn add_record_blocks<'a>(
    record: &bam::Record,
    invert_pair_strand: bool,
    gtf: &'a [Option<IntervalTree>],
    overlapping_features: &mut Vec<Vec<&'a Feature>>,
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
            match args._m.as_str() {
                "intersection-strict" => process_intersection_strict_read(
                    features,
                    reference_pos,
                    block_end,
                    effective_reverse,
                    overlapping_features,
                    args,
                ),
                "intersection-nonempty" => process_intersection_nonempty_read(
                    features,
                    reference_pos,
                    block_end,
                    effective_reverse,
                    overlapping_features,
                    args,
                ),
                "union" => process_union_read(
                    features,
                    reference_pos,
                    block_end,
                    effective_reverse,
                    overlapping_features,
                    args,
                ),
                _ => unreachable!(),
            }
        }

        if operation.consumes_ref() {
            reference_pos += length;
        }
    }

    true
}

fn assign_overlaps(
    overlapping_features: &[Vec<&Feature>],
    counts: &mut Counts,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) {
    // Keep references into the annotation instead of cloning Feature/String
    // values for every read.
    let mut unique_features = match args._m.as_str() {
        "intersection-strict" => filter_ambiguity_intersection_strict(overlapping_features),
        "intersection-nonempty" => filter_ambiguity_intersection_nonempty(overlapping_features),
        "union" => filter_ambiguity_union(overlapping_features),
        _ => unreachable!(),
    };
    let feature_count = unique_features.len();

    match feature_count {
        0 => {
            counts.no_feature += 1.0;
            if let Some(sender) = sender {
                let _ = sender.send(FeatureType::NoFeature);
            }
        }
        1 => {
            let feature = unique_features[0];
            counts.add_feature(feature.id(), 1.0);
            if let Some(sender) = sender {
                let _ = sender.send(FeatureType::Name(
                    counts.feature_names[feature.id()].to_string()
                ));
            }
        }
        _ => {
            counts.ambiguous += 1.0;
            match args.nonunique.as_str() {
                "all" => {
                    for feature in &unique_features {
                        counts.add_feature(feature.id(), 1.0);
                    }
                }
                "fraction" => {
                    let fractional_count = 1.0 / feature_count as f64;
                    for feature in &unique_features {
                        counts.add_feature(feature.id(), fractional_count);
                    }
                }
                "random" => {
                    let random_index = rand::random_range(0..feature_count);
                    let feature = unique_features[random_index];
                    counts.add_feature(feature.id(), 1.0);
                }
                _ => {}
            }

            if let Some(sender) = sender {
                let mut names: Vec<&str> = unique_features
                    .iter()
                    .map(|feature| counts.feature_names[feature.id()].as_ref())
                    .collect();
                names.sort_unstable();
                let _ = sender.send(FeatureType::Ambiguous(names.join("+")));
            }
        }
    }
}

fn count_single_record(
    record: &bam::Record,
    counter: &mut i32,
    counts: &mut Counts,
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
        counts.no_feature += 1.0;
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
    counts: &mut Counts,
    args: &Args,
    sender: Option<&mpsc::Sender<FeatureType>>,
) -> bool {
    let first_mapped = first.map(|r| r.flag().is_mapped()).unwrap_or(false);
    let second_mapped = second.map(|r| r.flag().is_mapped()).unwrap_or(false);

    if !first_mapped && !second_mapped {
        counts.not_aligned += 1.0;
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
        counts.alignment_not_unique += 1.0;
        if args.nonunique == "none" {
            if let Some(sender) = sender { let _ = sender.send(FeatureType::AlignmentNotUnique); }
            return true;
        }
    }

    let low_quality = first.map(|r| r.mapq() < args.a).unwrap_or(false)
        || second.map(|r| r.mapq() < args.a).unwrap_or(false);
    if low_quality {
        counts.too_low_aqual += 1.0;
        if let Some(sender) = sender { let _ = sender.send(FeatureType::TooLowaQual); }
        return true;
    }

    false
}

fn count_pair(
    first: Option<&bam::Record>,
    second: Option<&bam::Record>,
    counter: &mut i32,
    counts: &mut Counts,
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
        counts.no_feature += 1.0;
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
    counts: &mut Counts,
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
    counts: &mut Counts,
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
    counts: &mut Counts,
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
    counts: &mut Counts,
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

fn process_union_read<'a>(
    features: &'a IntervalTree,
    start_pos: i32,
    end_pos: i32,
    strand: bool,
    overlapping_features: &mut Vec<Vec<&'a Feature>>,
    args: &Args,
) {
    let new_overlap = features.overlap(start_pos, end_pos);
    let strand = if strand { '-' } else { '+' };
    add_stranded_features(new_overlap, strand, overlapping_features, args);
}

fn process_intersection_nonempty_read<'a>(
    features: &'a IntervalTree,
    start_pos: i32,
    end_pos: i32,
    strand: bool,
    overlapping_features: &mut Vec<Vec<&'a Feature>>,
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
                step_features.push(interval.data.as_ref().unwrap());
            }
        }

        if !step_features.is_empty() {
            overlapping_features.push(step_features);
        }
    }
}

fn process_intersection_strict_read<'a>(
    features: &'a IntervalTree,
    start_pos: i32,
    end_pos: i32,
    strand: bool,
    overlapping_features: &mut Vec<Vec<&'a Feature>>,
    args: &Args,
) {
    let new_contained = features.contains(start_pos, end_pos);
    let strand = if strand { '-' } else { '+' };
    add_stranded_features(new_contained, strand, overlapping_features, args);
}

fn sorted_unique_features<'a>(features: &[&'a Feature]) -> Vec<&'a Feature> {
    let mut unique: Vec<&'a Feature> = features
        .iter()
        .copied()
        .filter(|feature| !feature.name().is_empty())
        .collect();
    unique.sort_unstable_by_key(|feature| feature.id());
    unique.dedup_by_key(|feature| feature.id());
    unique
}

fn filter_ambiguity_union<'a>(
    overlapping_features: &[Vec<&'a Feature>],
) -> Vec<&'a Feature> {
    let total_features: usize = overlapping_features.iter().map(Vec::len).sum();
    let mut unique = Vec::with_capacity(total_features);
    for feature in overlapping_features.iter().flatten() {
        if !feature.name().is_empty() {
            unique.push(*feature);
        }
    }
    unique.sort_unstable_by_key(|feature| feature.id());
    unique.dedup_by_key(|feature| feature.id());
    unique
}

fn intersect_sorted_features<'a>(
    candidates: &mut Vec<&'a Feature>,
    features: &[&'a Feature],
) {
    if candidates.is_empty() {
        return;
    }

    let current = sorted_unique_features(features);
    if current.is_empty() {
        candidates.clear();
        return;
    }

    let mut write = 0;
    let mut i = 0;
    let mut j = 0;
    while i < candidates.len() && j < current.len() {
        match candidates[i].id().cmp(&current[j].id()) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                candidates[write] = candidates[i];
                write += 1;
                i += 1;
                j += 1;
            }
        }
    }
    candidates.truncate(write);
}

fn filter_ambiguity_intersection_strict<'a>(
    overlapping_features: &[Vec<&'a Feature>],
) -> Vec<&'a Feature> {
    let mut steps = overlapping_features.iter();
    let first = match steps.next() {
        Some(features) if !features.is_empty() => features,
        _ => return Vec::new(),
    };

    let mut candidates = sorted_unique_features(first);
    for features in steps {
        if features.is_empty() {
            return Vec::new();
        }
        intersect_sorted_features(&mut candidates, features);
        if candidates.is_empty() {
            return candidates;
        }
    }
    candidates
}

fn filter_ambiguity_intersection_nonempty<'a>(
    overlapping_features: &[Vec<&'a Feature>],
) -> Vec<&'a Feature> {
    let mut nonempty_steps = overlapping_features
        .iter()
        .filter(|features| !features.is_empty());

    let first = match nonempty_steps.next() {
        Some(features) => features,
        None => return Vec::new(),
    };

    let mut candidates = sorted_unique_features(first);
    for features in nonempty_steps {
        intersect_sorted_features(&mut candidates, features);
        if candidates.is_empty() {
            return candidates;
        }
    }
    candidates
}

fn feature_matches_strand(feature: &Feature, strand: char, args: &Args) -> bool {
    match args.stranded.as_str() {
        "yes" => feature.strand() == strand,
        "reverse" => feature.strand() != strand,
        "no" => true,
        _ => panic!("Invalid strandedness"),
    }
}

fn add_stranded_features<'a>(
    new_overlap: Vec<&'a Interval>,
    strand: char,
    overlapping_features: &mut Vec<Vec<&'a Feature>>,
    args: &Args,
) {
    let mut matching = Vec::new();
    for overlap in new_overlap {
        let feature = overlap.data.as_ref().unwrap();
        if feature_matches_strand(feature, strand, args) {
            matching.push(feature);
        }
    }
    overlapping_features.push(matching);
}

