extern crate clap;
extern crate shrust;
pub mod utils;
use clap::{arg, Command};
use itertools::Itertools;
use ndarray::Array2;
use ndarray::prelude::*;
use ndarray_npy::write_npy;
use utils::*;
use std::f64;
use std::{collections::HashMap, fs::File, io::{BufReader, Write}, sync::Mutex};
use bio::io::fastq;
use indicatif::{ProgressStyle, ProgressIterator};
use rayon::prelude::*;
use std::path::Path;
use anyhow::Result;
use libsufr::{
    suffix_array::SuffixArray, types::{ExtractResult, ExtractOptions, SufrBuilderArgs}, util::read_sequence_file
};
use flate2::read::GzDecoder;


fn kmer_vec_for_record(record: &fastq::Record, percent_mismatch: &f32)->Result<Vec<String>>{
    let seq_len = record.seq().len();
    let kmer_size = kmer_length(seq_len, *percent_mismatch);
    Ok(record.seq().windows(kmer_size).map(|x| String::from_utf8(x.to_vec()).unwrap()).collect())
}

/// Finds all kmer matches of a single read 
fn match_read_kmers(suffarr: &mut SuffixArray, record: &fastq::Record, percent_mismatch: &f32)->Result<Vec<ExtractResult>>{

    let kmer_vec = kmer_vec_for_record(record, percent_mismatch)?;
    let opts = ExtractOptions {
        queries: kmer_vec,
        max_query_len: None,
        low_memory: true,
        prefix_len: Some(100),
        suffix_len: None,
    };
    suffarr.extract(opts)
}

/// Returns a pair of Hashmaps. The first maps the read to its best alignment to each reference (reference_name, (alignment_start_pos, likelihood of alignment)).
/// The second returns the sum of likelihoods of all alignments to each reference.(reference_name, (sum of likelihood of all alignments)). 
fn query_read(suffarr: &mut SuffixArray, refs: &HashMap<String, String>, record: &fastq::Record, percent_mismatch: &f32)->Result<(HashMap<String, (usize, f64)>, HashMap<String, f64>)>{

    let read_len = record.seq().len();
    let read_seq = std::str::from_utf8(record.seq())?;
    let read_qual = record.qual();
    let max_num_mismatches: usize = (read_len as f32 * (percent_mismatch/100_f32)).floor() as usize;

    let best_match: Mutex<HashMap<String, (usize, f64)>> = Mutex::new(HashMap::new());
    let match_likelihood: Mutex<HashMap<String, f64>> = Mutex::new(HashMap::new());


    let matches = match_read_kmers(suffarr, record, percent_mismatch)?;

    matches.par_iter().for_each(|hit| {
        let kmer_start_read = hit.query_num;
        
        for ref_entry in hit.sequences.iter(){
            let seq_name = &ref_entry.sequence_name;
            let seq_start = ref_entry.sequence_start;
            let seq_end = ref_entry.sequence_range.end;
            let ref_len = seq_end-seq_start-1;

            let kmer_start_ref = ref_entry.suffix_offset;
            // if the read alignment starts within the reference and ends within the reference
            if kmer_start_ref>=kmer_start_read && kmer_start_ref+read_len<=ref_len+kmer_start_read{
                // retrieve the full reference
                let align_start = kmer_start_ref-kmer_start_read;
                let seq = refs.get(seq_name).unwrap();
                let ref_match_seg = &seq[align_start..align_start+read_len];
                if num_mismatches(read_seq, ref_match_seg)<=max_num_mismatches{
                    let match_log_prob = compute_match_log_prob(read_seq, read_qual, ref_match_seg);

                    // Update best match found
                    best_match.lock().unwrap().entry(seq_name.to_string())
                        .and_modify(|e| {
                            // e.push((align_start, match_log_prob))
                            if e.1<=match_log_prob{
                                *e = (align_start, match_log_prob);
                            }
                        })
                        .or_insert((align_start, match_log_prob));

                    // Update match score
                    match_likelihood.lock().unwrap().entry(seq_name.to_string())
                        .and_modify(|e| {
                            // e.push((align_start, match_log_prob))
                            *e += match_log_prob;
                        })
                        .or_insert(match_log_prob);

                }
            }
        }
    });

    Ok((best_match.into_inner().unwrap(), match_likelihood.into_inner().unwrap()))
}

fn process_fastq_file(suffarr: &mut SuffixArray, 
                    refs: &HashMap<String, String>, 
                    fastq_file: &Path, 
                    percent_mismatch: &f32, 
                    // outpath: Option<&Mutex<File>>, 
                    ref_ids: &HashMap<String, usize>, 
                    read_ids: &HashMap<String, usize>, 
                    ll_array: &mut Mutex<Array2<f64>>)-> Result<HashMap<(usize, usize), usize>>
{
    let pb = ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})").unwrap();

    let fastq_records = match get_extension_from_filename(fastq_file.to_str().expect("Invalid reads file!")){
       Some("gz") => {
            let f = File::open(fastq_file)?;
            let decoder = GzDecoder::new(f);

            fastq::Reader::from_bufread(BufReader::new(decoder)).records().map(|x| x.unwrap()).collect_vec()
            // fastq::Reader::from_bufread(GzDecoder::new(reader)).records()
       },
       Some("fastq") => {
            let f = File::open(fastq_file)?;
            let reader = BufReader::new(f);
            fastq::Reader::from_bufread(reader).records().map(|x| x.unwrap()).collect_vec()
        },
       _ => {panic!("Invalid file type for reads!")}
    };

    let mut out_aligns: HashMap<(usize, usize), usize> = HashMap::new();

    // let fastq_records = reader
    //     .records()
    //     .map(|x| x.unwrap())
    //     .collect_vec();

    fastq_records
        .iter()
        .progress_with_style(pb)
        .map(|record| {
            let (best_hits, match_likelihoods) = query_read(suffarr, refs, record, percent_mismatch).unwrap();
            // dbg!(&hits);
            (record.id().to_string(), best_hits, match_likelihoods)
        })
        // .par_bridge()
        .for_each(|(read_id, best_hits, match_likelihoods)| {
            match_likelihoods.iter()
                // .map(|(ref_id, positions)| (ref_id.clone(), positions))
                .for_each(|x| {
                    // let outstr = format!("{}\t{}\t{}\t{}\n", read_id, x.0, x.1.0, x.1.1);
                    // dbg!(read_ids, &read_id);
                    let read_idx = read_ids.get(&read_id).unwrap();
                    let ref_idx = ref_ids.get(x.0).unwrap();
                    let best_align = best_hits.get(x.0).unwrap();

                    out_aligns.insert((*read_idx, *ref_idx), best_align.0);
                    // match outpath{
                    //     Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                    //     None => {print!("{}", outstr)},
                    // };

                    ll_array.lock().unwrap()[[read_idx.clone(), ref_idx.clone()]] = *x.1;
                });
        });

    Ok(out_aligns)
}

fn get_proportions(ll_array: &Array2<f64>, num_iter: usize)->Vec<(usize, usize, f64)>{
    let num_reads = ll_array.shape()[0];
    let num_srcs = ll_array.shape()[1];

    let mut props = Array::<f64, _>::zeros((num_iter+1, num_srcs).f());
    let mut w = Array::<f64, _>::zeros((num_reads, num_srcs).f());

    props.slice_mut(s![0, ..num_srcs]).fill(1_f64/num_srcs as f64);

    for i in 0..num_iter{
        // let tmp_w = Array::<f64, _>::ones((num_reads, num_srcs).f());
        let tmp_prop = props.slice_mut(s![i, ..num_srcs]);

        w = ll_array.clone()*tmp_prop;
        // clik = np.sum(w, axis=1)[:, np.newaxis]
        let clik = w.sum_axis(Axis(1)).insert_axis(Axis(1));
        

        // dbg!(w.shape(), clik.shape());

        w = w/clik;

        // props[i+1] = np.mean(w, axis=0)
        props.slice_mut(s![i+1, ..num_srcs]).assign(w.mean_axis(Axis(0)).as_ref().unwrap());

    }

    let mut em_props = Array::<f64, _>::zeros((num_srcs).f());

    em_props.assign(&props.slice(s![num_iter, ..]));

    // w = em_props*np.exp(ll_array)
    w = ll_array.clone()*em_props;
    // clik = np.sum(w, axis=1)[:, np.newaxis]
    let clik = w.sum_axis(Axis(1)).insert_axis(Axis(1));
    // w /= clik
    w = w/clik;

    let mut results = vec![];

    for (i, row) in w.axis_iter(Axis(0)).enumerate() {
        // dbg!(row.shape());
        let (max_idx, _) =
            row.iter()
                .enumerate()
                .fold((0, row[0]), |(idx_max, val_max), (idx, val)| {
                    if &val_max > val {
                        (idx_max, val_max)
                    } else {
                        (idx, *val)
                    }
                });
        // dbg!(i, row.shape());
        // dbg!();
        results.push((max_idx, i, ll_array.ln()[[i, max_idx]]));
    }    


    results
}




fn main() -> Result<()>{
    let matches = Command::new("Maximum Likelihood Metagenomic Classification")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(
            Command::new("build")
                .about("Build suffix tree index from reference fasta file")
                .arg(arg!(-s --source <SRC_FILE> "Source file with sequences(fasta)")
                    .required(true)
                )
                .arg(arg!(-o --out <OUTFILE> "Output index file name")
                    .default_value("")
                    .value_parser(clap::value_parser!(String))
                )

        )
        .subcommand(
            Command::new("fasta")
               .about("Generate Fasta file containing sequences present in index")
               .arg(arg!(-i --index <INDEX_FILE> "Source index file of reference sequences(.sufr)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
        )
        .subcommand(
            Command::new("inspect")
               .about("Inspect a pre-built index")
               .arg(arg!(-i --index <INDEX_FILE> "Source index file of reference sequences(.sufr)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
        )
        .subcommand(
            Command::new("query")
                .arg(arg!(-s --source <SRC_FILE> "Source index file with reference sequences(.sufr)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-p --percent_mismatch <PERCENT_MISMATCH>"Percent mismatch to reference sequences")
                    .required(true)
                    .value_parser(clap::value_parser!(f32))
                    )
                .arg(arg!(-r --reads <READS>"Source file with read sequences(fasta)")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-o --out <OUT_FILE>"Output file")
                    .default_value("out.matches")
                    .value_parser(clap::value_parser!(String))
                    )
                .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2)")
                    .default_value("2")
                    .value_parser(clap::value_parser!(usize))
                    )
        )
        .about("Maximum Likelihood Metagenomic classifier using Suffix trees")
        .get_matches();
  
    match matches.subcommand(){
        Some(("build",  sub_m)) => {
            let src_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();
            
            match outfile{
                "" => {
                    let sequence_delimiter = b'$';
                    let seq_data = read_sequence_file(Path::new(src_file), sequence_delimiter)?;

                    let builder_args = SufrBuilderArgs {
                        text: seq_data.seq,
                        path: Some(format!("{}.sufr", src_file)),
                        low_memory: false,
                        max_query_len: None,
                        is_dna: true,
                        allow_ambiguity: false,
                        ignore_softmask: true,
                        sequence_starts: seq_data.start_positions.into_iter().collect(),
                        sequence_names: seq_data.sequence_names,
                        num_partitions: 16,
                        seed_mask: None,
                        random_seed: 42,
                    };

                    let _ = SuffixArray::new(builder_args)?;

                },
                _ => {
                    let sequence_delimiter = b'$';
                    let seq_data = read_sequence_file(Path::new(src_file), sequence_delimiter)?;

                    let builder_args = SufrBuilderArgs {
                        text: seq_data.seq,
                        path: Some(outfile.to_string()),
                        low_memory: false,
                        max_query_len: None,
                        is_dna: true,
                        allow_ambiguity: false,
                        ignore_softmask: true,
                        sequence_starts: seq_data.start_positions.into_iter().collect(),
                        sequence_names: seq_data.sequence_names,
                        num_partitions: 16,
                        seed_mask: None,
                        random_seed: 42,
                    };

                    let _ = SuffixArray::new(builder_args)?;
                }
            };

        },
        Some(("inspect",  sub_m)) => {
            let src_file = sub_m.get_one::<String>("index").expect("required").as_str();
            summarize_index(src_file)?;

        },
        Some(("fasta",  sub_m)) => {
            let src_file = sub_m.get_one::<String>("index").expect("required").as_str();
            let mut suffarr = SuffixArray::read(src_file, false)?;

            let ref_string = get_refs_from_sa(&mut suffarr)?.into_iter()
                .map(|(id, seq)| format!("\n>{}\n{}\n", id, seq))
                .collect::<String>();

            let mut outfile = File::create(format!("{}.fasta", src_file))?;
            outfile.write_all(ref_string.as_bytes())?;

            // Ok(())

        },
        Some(("query",  sub_m)) => {
            let ref_file = sub_m.get_one::<String>("source").expect("required").as_str();
            let reads_file = sub_m.get_one::<String>("reads").expect("required").as_str();
            let num_threads = sub_m.get_one::<usize>("threads").expect("required");
            let percent_mismatch = sub_m.get_one::<f32>("percent_mismatch").expect("required");
            let outfile = sub_m.get_one::<String>("out").unwrap().as_str();

            let outpath = Some(Mutex::new(File::create(outfile).unwrap()));
            let outstr = format!("ReadID\tRefID\tStart_pos\tPosterior\n");
            match outpath.as_ref(){
                Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                None => {print!("{}", outstr)},
            };


            rayon::ThreadPoolBuilder::new().num_threads(*num_threads).build_global()?;

            let fastq_records = match get_extension_from_filename(reads_file){
                Some("gz") => {
                    let f = File::open(reads_file)?;
                    let decoder = GzDecoder::new(f);
                    fastq::Reader::from_bufread(BufReader::new(decoder)).records().map(|x| x.unwrap()).collect_vec()
                },
                Some("fastq") => {
                    let f = File::open(reads_file)?;
                    let reader = BufReader::new(f);
                    fastq::Reader::from_bufread(reader).records().map(|x| x.unwrap()).collect_vec()
                },
                _ => panic!("Invalid file type for reads!")
            };

            let mut read_ids: HashMap<String, usize> = HashMap::new();
            let mut read_ids_rev: HashMap<usize, String> = HashMap::new();
            for result in fastq_records.iter().enumerate() {
                let record_id = result.1.id().to_string();
                read_ids.insert(record_id.clone(), result.0);
                read_ids_rev.insert(result.0, record_id);
            }


            let mut suffarr: SuffixArray = SuffixArray::read(ref_file, false)?;

            let refs = get_refs_from_sa(&mut suffarr)?;

            let mut ref_ids: HashMap<String, usize> = HashMap::new();
            let mut ref_ids_rev: HashMap<usize, String> = HashMap::new();
            for (idx, result) in refs.keys().enumerate() {
                ref_ids.insert(result.to_string(), idx);
                ref_ids_rev.insert(idx,result.to_string());
            }

            let mut ll_array = Mutex::new(Array2::<f64>::zeros((read_ids.len(), ref_ids.len())));
            ll_array.lock().unwrap().fill(f64::MIN);

            let read_idxs = serde_json::to_string(&read_ids).unwrap();
            let ref_idxs = serde_json::to_string(&ref_ids).unwrap();

            File::create("read_idxs.json").unwrap().write_all(read_idxs.as_bytes())?;
            File::create("ref_idxs.json").unwrap().write_all(ref_idxs.as_bytes())?;

            let out_alignments = process_fastq_file(&mut suffarr, &refs, Path::new(reads_file), percent_mismatch, &ref_ids, &read_ids, &mut ll_array)?;

            let props = get_proportions(&ll_array.lock().unwrap().exp(), 100);

            for (ref_id, read_id, val) in props{
                let outstr = format!("{}\t{}\t{}\t{:.5}\n", read_ids_rev.get(&read_id).unwrap(), ref_ids_rev.get(&ref_id).unwrap(), out_alignments.get(&(read_id, ref_id)).unwrap(), val);
                match outpath.as_ref(){
                    Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes()).unwrap();},
                    None => {print!("{}", outstr)},
                };
            }

            write_npy("ll_array.npy", &ll_array.into_inner().unwrap())?;
        },
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    }
    
    Ok(())
}
