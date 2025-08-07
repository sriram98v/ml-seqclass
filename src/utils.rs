use itertools::izip;
use std::path::Path;
use std::{io::Write, fs};
use std::collections::HashMap;

/// Compute probability of match given ref is true source
pub fn compute_match_log_prob(q_seq: &str, quality_score_vec: &[u8], aligned_ref_seq: &str) -> f64{
    let mut match_log_likelihood = 0_f64;
    for (read_char,reference_char,quality_score) in izip!(q_seq.chars(), aligned_ref_seq.chars(), quality_score_vec){
        match read_char==reference_char{
            true => match_log_likelihood += (1_f64-error_prob(quality_score.clone())).log10(),
            false => match_log_likelihood += (error_prob(quality_score.clone())).log10(),
        }
    }
    return match_log_likelihood
}

/// Compute minimum length of match required for a partial match.
pub fn kmer_length(seq_len: usize, percent_mismatch: f32)->usize{
    let num_mismatches: usize = (seq_len as f32 * (percent_mismatch/100_f32)).floor() as usize;
    seq_len/(num_mismatches+1)
}

pub fn num_mismatches(read_seq: &str, ref_seq: &str)->usize{
    read_seq.chars().zip(ref_seq.chars()).map(|(a,b)| (a!=b) as usize).sum()
}

pub fn write_matches(outpath: Option<&String>, matches: &HashMap<String, Vec<(String, usize, f64)>>)->std::io::Result<()>
{
    match outpath{
        Some(path) => {
            if Path::new(path).exists(){
                fs::remove_file(path)?;
            }
            let mut file_ref = fs::OpenOptions::new().create_new(true).append(true).open(path).expect("Invalid path for result file!");
            let result_header:String = "Read_ID\tRef_ID\thit position\tLog-Likelihood\n".to_string();
            file_ref.write_all(result_header.as_bytes()).expect("write failed");
            matches.into_iter()
                .filter(|(_seq_id, all_matches)| {
                    all_matches[0].0.is_empty()
                })
                .for_each(|(seq_id, all_matches)| {
                    all_matches.iter()
                        .for_each(|(read_id, hit_pos, match_score)| {
                            let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                            file_ref.write_all(out_string.as_bytes()).expect("write failed");
                        });
                    
            });
        },
        None => {
            let result_header:String = "Read_ID\tRef_ID\tposition\tLog-Likelihood\n".to_string();
            print!("{}", result_header);
            matches.into_iter()
                .for_each(|(seq_id, all_matches)| {
                    all_matches.iter()
                        .for_each(|(read_id, hit_pos, match_score)| {
                            let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                            println!("{}", out_string);
                        });
        });
        }
    }
    Ok(())
}

pub fn error_prob(q: u8)->f64{
    10_f64.powf(-(q as f64/10_f64))
}

pub fn error_log_prob(q: u8)->f64{
    -(q as f64/10_f64)
}

pub fn complement(q_seq: Vec<char>)->Vec<char>{
    q_seq.iter()
        .map(|x|{
            match x{
                'T' => 'A',
                'C' => 'G',
                'A' => 'T',
                'G' => 'C',
                _ => 'E',
            }
        })
        .collect()
}

// fn build_sa(ref_file: &Path, outfile: &String)->Result<()>{
//     let sequence_delimiter = b'$';
//     let seq_data = read_sequence_file(ref_file, sequence_delimiter)?;

//     let builder_args = SufrBuilderArgs {
//         text: seq_data.seq,
//         path: Some(outfile.to_string()),
//         low_memory: false,
//         max_query_len: None,
//         is_dna: true,
//         allow_ambiguity: false,
//         ignore_softmask: true,
//         sequence_starts: seq_data.start_positions.into_iter().collect(),
//         sequence_names: seq_data.sequence_names,
//         num_partitions: 16,
//         seed_mask: None,
//         random_seed: 42,
//     };
//     let suffix_array = SuffixArray::new(builder_args)?;

//     Ok(())
// }

// fn load_sufr(fname: &Path)->Result<SuffixArray>{
//     SuffixArray::read(fname.to_str().expect(".sufr file cannot be found!"), true)
// }
