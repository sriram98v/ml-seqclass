use itertools::izip;
use std::path::Path;
use std::{io::Write, fs};
use std::collections::HashMap;

/// Compute probability of match given ref is true source
pub fn match_log_prob(q_seq: &[char], quality_score_vec: &[u8], aligned_ref_seq: &[char]) -> f64{
    let mut match_log_likelihood = 0_f64;
    for (read_char,reference_char,quality_score) in izip!(q_seq, aligned_ref_seq, quality_score_vec){
        match read_char==reference_char{
            true => match_log_likelihood += (1_f64-error_prob(quality_score.clone())).log10(),
            false => match_log_likelihood += (error_prob(quality_score.clone())).log10(),
        }
    }
    return match_log_likelihood
}

pub fn write_matches(outpath: Option<&String>, matches: &HashMap<String, (String, usize, f64)>)->std::io::Result<()>
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
                .filter(|(_seq_id, (read_id, _hit_pos, _match_score))| {
                    read_id.is_empty()
                })
                .for_each(|(seq_id, (read_id, hit_pos, match_score))| {
                let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                file_ref.write_all(out_string.as_bytes()).expect("write failed");
            });
        },
        None => {
            let result_header:String = "Read_ID\tRef_ID\tposition\tLog-Likelihood\n".to_string();
            print!("{}", result_header);
            matches.into_iter().for_each(|(seq_id, (read_id, hit_pos, match_score))| {
                let out_string:String = format!("{}\t{}\t{}\t{}\n", seq_id, read_id, hit_pos, match_score);
                print!("{}", out_string);
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