extern crate clap;
extern crate shrust;
pub mod utils;
use clap::{arg, Command};
use itertools::Itertools;
use ndarray::{Array1, Array2, Axis};
use ndarray_npy::write_npy;
use utils::*;
use std::{collections::HashMap, fs::File, io::Write, sync::Mutex};
use bio::io::{fasta,  fastq};
use generalized_suffix_tree::suffix_tree::KGST;
use generalized_suffix_tree::data::tree_item::TreeItem;
use indicatif::{ProgressBar, ProgressStyle, ParallelProgressIterator};
use rayon::prelude::*;
use std::path::Path;

fn build_tree(file: &Path, max_depth: &usize)->KGST<char, String>{
    // println!("Reference Sequence index: {}", file.display());
    let reader = fasta::Reader::from_file(file).unwrap();

    let total_size = reader.records().count();
    

    let pb = ProgressBar::new(total_size as u64).with_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})").unwrap_or(ProgressStyle::default_bar()));
    
    let mut tree: KGST<char, String> = KGST::new('$');

    let reader = fasta::Reader::from_file(file).unwrap();
    
    for result in reader.records() {

        let result_data = result.unwrap();

        let seq: Vec<char> = result_data.seq()
        .to_vec()
        .iter()
        .map(|x| *x as char)
        .collect();
            
        tree.insert(result_data.id().to_string(), seq, &max_depth);

        pb.inc(1);   
    }
    tree
}

/// Returns alignments of one read with log probability of match.
fn query_tree(tree: &KGST<char, String>, q_seq: &[char], num_mismatches: &usize, qual_seq: &[u8])->Vec<(String, Vec<(usize, f64)>)>{
    let mut match_set: Vec<(String, Vec<(usize, f64)>)> = Vec::new();
    let string_len: usize = q_seq.len();
    let chunk_size: usize = string_len/(num_mismatches+1);
    let refs: HashMap<String, Vec<char>> = tree.get_strings()
                                            .values()
                                            .map(|x| (x.0.get_id().clone(), x.0.get_string().iter().map(|a| a.into_inner().cloned().unwrap()).collect_vec()))
                                            .collect::<HashMap<String, Vec<char>>>();
    if string_len>=chunk_size{
        for depth in 0..string_len+1-chunk_size{
            let mut temp_matches: Vec<(String, Vec<(usize, f64)>)> = Vec::new();
            let chunk = &q_seq[depth..depth+(chunk_size)].to_vec();
            for (partial_match_label, partial_match_pos) in tree.substring_match(chunk).into_iter(){
                temp_matches.push((partial_match_label.clone(), partial_match_pos.into_iter()
                                                            .filter(|start_pos| {
                                                                if !(start_pos>=&depth && start_pos+&q_seq.len()<refs.get(&partial_match_label).unwrap().len()){
                                                                    return false;
                                                                }   
                                                                    let mut num_mis: usize = 0;
                                                                    let mut match_string: Vec<char> = vec![];
                                                                    for (a,b) in q_seq.iter().zip(refs.get(&partial_match_label).unwrap()[start_pos-depth..start_pos+q_seq.len()-depth].iter()){
                                                                        match a==b{
                                                                            true => {match_string.push('0')},
                                                                            false => {num_mis+=1;match_string.push('1');}
                                                                        }
                                                                    }
                                                                    match &num_mis<num_mismatches{
                                                                        true => return true,
                                                                        false => return false,
                                                                    };
                                                            })
                                                            .map(|x| {
                                                                (x-depth, match_log_prob(q_seq, qual_seq, &refs.get(&partial_match_label).unwrap()[x-depth..x+q_seq.len()-depth]))
                                                            })
                                                            .collect::<Vec<(usize, f64)>>()));
            }
            match_set.extend(temp_matches);
        }
    }
    match_set.into_iter()
        .filter(|(_id, positions)| !positions.is_empty())
        .collect()
}

/// Returns alignments of all reads with corresponding log probability
fn process_fastq_file(tree:&KGST<char, String>, fastq_file: &Path, percent_mismatch: &f32, outpath: Option<&Mutex<File>>, ref_ids: &HashMap<String, usize>, read_ids: &HashMap<String, usize>, ll_array: &mut Mutex<Array2<f64>>)
{

    let pb = ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})").unwrap();

    let reader = fastq::Reader::from_file(fastq_file).unwrap();

    let fastq_records = reader.records().map(|x| x.unwrap()).collect_vec();

    return fastq_records.par_iter().progress_with_style(pb).map(|read_data| {
            let read_id: String = read_data.id().to_string();
            let read_qual: Vec<u8> = read_data.qual().iter().map(|x| x-33).collect();
            let num_mismatches: usize = (read_data.seq().len() as f32 * (percent_mismatch/100_f32)).floor() as usize;

            let seq: Vec<char> = read_data.seq()
                .to_vec()
                .iter()
                .map(|x| *x as char)
                .collect();
                
            // (ref_id, vec<(start_positions, log_prob)>)
            let hits: Vec<(String, Vec<(usize, f64)>)> = query_tree(tree, &seq, &num_mismatches, &read_qual);
            
            return (read_id, hits);
        })
        .for_each(|(read_id, read_matches)| {
            if read_matches.len()>0{
                read_matches.iter().map(|(ref_id, positions)| (ref_id.clone(), positions.iter().max_by(|x, y| x.1.total_cmp(&y.1)).cloned().unwrap()))
                        .for_each(|x| {
                            // let ref_pos_likelihood = (x.0, x.1.0, x.1.1);
                            let outstr = format!("{}\t{}\t{}\t{}\n", read_id, x.0, x.1.0, x.1.1);
                            let read_idx = read_ids.get(&read_id).unwrap();
                            let ref_idx = ref_ids.get(&x.0).unwrap();

                            ll_array.lock().unwrap()[[read_idx.clone(), ref_idx.clone()]] = x.1.1;

                            // match outpath{
                            //     Some(file) => {file.lock().unwrap().write_all(outstr.as_bytes());},
                            //     None => {print!("{}", outstr)},
                            // };
                        });
                        // .collect_vec();
                        // .max_by(|x, y| x.2.total_cmp(&y.2))
                        // .unwrap();
            }
            
        })
    
    // all_hits
    
}

fn main() {
    let matches = Command::new("Maximum Likelihood Metagenomic Classification")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .arg(arg!(-s --source <SRC_FILE> "Source file with reference sequences(fasta)")
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
            .value_parser(clap::value_parser!(String))
            )
        .arg(arg!(-t --threads <THREADS>"Number of threads (defaults to 2)")
            .value_parser(clap::value_parser!(usize))
            )
        .about("Maximum Likelihood Metagenomic classifier using Suffix trees")
        .get_matches();
    
    let ref_path = Path::new(matches.get_one::<String>("source").expect("required").as_str());
    let read_path = Path::new(matches.get_one::<String>("reads").expect("required").as_str());
    let out_path = matches.get_one::<String>("out");
    let percent_mismatch = matches.get_one::<f32>("percent_mismatch").expect("required");
    let num_threads: Option<&usize> = matches.get_one::<usize>("threads");
    match num_threads{
        Some(threads) => rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap(),
        None => rayon::ThreadPoolBuilder::new().num_threads(2).build_global().unwrap(),
    };

    let outfile = match out_path{
        Some(file_path) => Some(Mutex::new(File::create(file_path).unwrap())),
        _ => None,
    };
    let fastq_reader = fastq::Reader::from_file(read_path).unwrap();

    let max_read_len = fastq_reader.records().par_bridge().map(|x| x.unwrap().seq().len()).max().unwrap();
    let num_mismatches = (max_read_len as f32 * (percent_mismatch/100_f32)).floor() as usize;
    let tree_depth = (max_read_len/(num_mismatches+1))+2;

    let fastq_reader = fastq::Reader::from_file(read_path).unwrap();
    let mut read_ids: HashMap<String, usize> = HashMap::new();
    for result in fastq_reader.records().enumerate() {
        let record_id = result.1.expect("Error during fastq record parsing").id().to_string();
        read_ids.insert(record_id, result.0);
    }

    let ref_reader = fasta::Reader::from_file(ref_path).unwrap();
    let mut ref_ids: HashMap<String, usize> = HashMap::new();
    for result in ref_reader.records().enumerate() {
        let record_id = result.1.expect("Error during fasta record parsing").id().to_string();
        ref_ids.insert(record_id, result.0);
    }

    // likelihood arrays indexed by (read_id,ref_id)
    let mut ll_array = Mutex::new(Array2::<f64>::zeros((read_ids.len(), ref_ids.len())));
    
    let read_idxs = serde_json::to_string(&read_ids).unwrap();
    let ref_idxs = serde_json::to_string(&ref_ids).unwrap();

    File::create("read_idxs.json").unwrap().write_all(read_idxs.as_bytes());
    File::create("ref_idxs.json").unwrap().write_all(ref_idxs.as_bytes());
    let tree: KGST<char, String> = build_tree(ref_path, &tree_depth);
    process_fastq_file(&tree, read_path, percent_mismatch, outfile.as_ref(), &ref_ids, &read_ids, &mut ll_array);
    dbg!(&ll_array);
    // let array = &ll_array.into_inner().unwrap();
    write_npy("ll_array.npy", &ll_array.into_inner().unwrap());
    // let _ = write_matches(out_path, &alignments);
}
