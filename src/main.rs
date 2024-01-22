// Importeer de benodigde modules
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use structopt::StructOpt;
use bam::BamReader;

// Definieer een structuur om de argumenten op te slaan
#[derive(StructOpt)]
struct Args {
    // Geef de naam en het type van de parameter n
    #[structopt(short = "n", long = "number", default_value = "10", help="Number of threads")]
    n: u16,

    // Geef de naam en het type van de parameter m
    #[structopt(short = "m", long = "mode", default_value = "intersection-strict")]
    m: String,

    // Geef de naam en het type van de parameter --stranded
    #[structopt(long = "stranded")]
    stranded: bool,

    // Geef de naam en het type van de parameter a
    #[structopt(short = "a", long = "amount", default_value = "10", help="Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads.")]
    a: i32,

    // Geef de naam en het type van de parameter t
    #[structopt(short = "t", long = "type", default_value = "exon", help="Feature type (3rd column in GTF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)")]
    t: String,

    // Geef de naam en het type van de parameter i
    #[structopt(short = "i", long = "id", default_value = "gene_name", help="GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). All feature of the right type (see -t option) within the same GTF
    attribute will be added together. The typical way of using this option is to count all exonic reads from each gene and add the exons but other uses are possible
    as well. You can call this option multiple times: in that case, the combination of all attributes separated by colons (:) will be used as a unique identifier,
    e.g. for exons you might use -i gene_id -i exon_number.")]
    i: String,

    // Geef de naam en het type van de parameter bam
    #[structopt(name = "bam", default_value="test.bam")]
    bam: String,

    // Geef de naam en het type van de parameter gtf
    #[structopt(name = "gtf", default_value="test.gtf")]
    gtf: String,
}


// Definieer een structuur om de feature-informatie op te slaan
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
    let mut counter=0;
    for line in reader.lines() {
        counter+=1;
        if counter%100000==0 {
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


            let name = name.trim_matches('"');
            let name = name.replace("gene_name ","");

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
    println!("i = {}", args.i);
    println!("bam = {}", args.bam);
    println!("gtf = {}", args.gtf);
    
    // Lees bam file (eerst, indien probleem, sneller crashen)
    let bam = BamReader::from_path(args.bam, args.n).expect("Could not read bam file");
    let header = bam.header().clone();
    let reference_names: &[String] = header.reference_names();
    //println!("{:?}", reference_names);

    
    // Lees het bestand en sla het resultaat op in een variabele
    let map = read_gtf("test.gtf", args.t.as_str());

    let mut counts: HashMap<String, i32> =HashMap::new();
    
    let mut counter= 0;
    for record in bam {
        counter+=1;
        if counter%100000 == 0 {
            println!("{counter} alignment records processed");
        }
        let record = record.unwrap();
        let ref_id: usize= record.ref_id().try_into().unwrap();
        let reference=reference_names[ref_id].clone();
        if map.contains_key(&reference) {
            let start_pos: u32 = record.start().try_into().unwrap();
            let end_pos: u32 = record.calculate_end().try_into().unwrap();
            if let Ok(index) = map[&reference].binary_search_by(|f| f.start.cmp(&start_pos)) {
                // Begin met deze index en ga door totdat je een feature vindt die groter is dan de eindpositie van de read
                let mut i = index;
                while i < map[&reference].len() && map[&reference][i].start <= end_pos {
                    // Controleer of de feature overlapt met de read
                    if (map[&reference][i].end >= start_pos) && (map[&reference][i].end <= end_pos) {
                        println!("Feature found! {}: {}:{}-{} (position: {})", map[&reference][i].name,reference, map[&reference][i].start, map[&reference][i].end, start_pos);
                        let count = counts.entry(map[&reference][i].name.clone()).or_insert(0);
                        *count += 1;
                    }
                    // Ga naar de volgende feature
                    i += 1;
                }
            }
        }
    }

    
    
    // Print de HashMap
    println!("{:?}", counts);
}
