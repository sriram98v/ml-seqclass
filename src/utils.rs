use anyhow::Result;
use format_num::NumberFormat;
use itertools::izip;
use libsufr::suffix_array::SuffixArray;
use libsufr::types::SuffixSortType;
use std::path::Path;
use std::{io::Write, fs};
use std::collections::HashMap;
use tabled::Table;
use std::ffi::OsStr;

/// Compute probability of match given ref is true source
pub fn compute_match_log_prob(q_seq: &str, quality_score_vec: &[u8], aligned_ref_seq: &str) -> f32{
    let mut match_log_likelihood = 0_f32;
    for (read_char,reference_char,quality_score) in izip!(q_seq.chars(), aligned_ref_seq.chars(), quality_score_vec){
        match read_char==reference_char{
            true => match_log_likelihood += (1_f32-error_prob(quality_score.clone())).log10(),
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

pub fn write_matches(outpath: Option<&String>, matches: &HashMap<String, Vec<(String, usize, f32)>>)->std::io::Result<()>
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

pub fn error_prob(q: u8)->f32{
    10_f32.powf(-(q as f32/10_f32))
}

pub fn error_log_prob(q: u8)->f32{
    -(q as f32/10_f32)
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

pub fn get_refs_from_sa(suffarr: &mut SuffixArray)->Result<HashMap<String, String>>{
    let seq_starts = suffarr.metadata()?.sequence_starts;
    let seq_names = suffarr.metadata()?.sequence_names;
    let text_len = suffarr.metadata()?.text_len;
    let mut last_stop = text_len;
    let refs = seq_names.iter().rev().zip(seq_starts.iter().rev())
        .map(|(label, start)| {
            let seq = suffarr.string_at(*start, Some(last_stop-start-1)).unwrap();
            last_stop = *start;
            (label.to_string(), seq)
        }).collect::<HashMap<String, String>>();
    
    Ok(refs)
}

pub fn get_ref_starts_from_sa(suffarr: &mut SuffixArray)->Result<HashMap<usize, String>>{
    let seq_starts = suffarr.metadata()?.sequence_starts;
    let seq_names = suffarr.metadata()?.sequence_names;
    let text_len = suffarr.metadata()?.text_len;
    let mut last_stop = text_len;
    let refs = seq_names.iter().rev().zip(seq_starts.iter().rev())
        .map(|(_, start)| {
            let seq = suffarr.string_at(*start, Some(last_stop-start-1)).unwrap();
            last_stop = *start;
            (*start, seq)
        }).collect::<HashMap<usize, String>>();
    
    Ok(refs)
}

pub fn summarize_index(sufr_file: &str) -> Result<()> {
    let suffix_array = SuffixArray::read(sufr_file, false)?;
    let meta = suffix_array.metadata()?;
    let num_fmt = NumberFormat::new();
    let mut rows = vec![vec!["Filename".to_string(), meta.filename.clone()]];
    rows.push(vec![
        "Modified".to_string(),
        meta.modified.format("%Y-%m-%d %H:%M").to_string(),
    ]);
    rows.push(vec![
        "File Size".to_string(),
        format!("{} bytes", num_fmt.format(",.0", meta.file_size as f32)),
    ]);
    rows.push(vec![
        "File Version".to_string(),
        meta.file_version.to_string(),
    ]);
    rows.push(vec!["DNA".to_string(), meta.is_dna.to_string()]);
    rows.push(vec![
        "Allow Ambiguity".to_string(),
        meta.allow_ambiguity.to_string(),
    ]);
    rows.push(vec![
        "Ignore Softmask".to_string(),
        meta.ignore_softmask.to_string(),
    ]);
    rows.push(vec![
        "Text Length".to_string(),
        num_fmt.format(",.0", meta.text_len as f32),
    ]);
    rows.push(vec![
        "Len Suffixes".to_string(),
        num_fmt.format(",.0", meta.len_suffixes as f32),
    ]);

    match meta.sort_type {
        SuffixSortType::Mask(seed_mask) => {
            rows.push(vec!["Seed mask".to_string(), seed_mask.mask])
        }
        SuffixSortType::MaxQueryLen(max_query_len) => rows.push(vec![
            "Max query len".to_string(),
            num_fmt.format(",.0", max_query_len as f32),
        ]),
    };

    rows.push(vec![
        "Num sequences".to_string(),
        num_fmt.format(",.0", meta.num_sequences as f32),
    ]);
    let seq_starts = meta
        .sequence_starts
        .into_iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    rows.push(vec![
        "Sequence starts".to_string(),
        textwrap::wrap(&seq_starts, 40).join("\n"),
    ]);
    let seq_names = meta.sequence_names.join(", ");
    rows.push(vec![
        "Sequence names".to_string(),
        textwrap::wrap(&seq_names, 40).join("\n"),
    ]);

    let table = Table::from_iter(rows);
    println!("{table}");

    Ok(())
}

pub fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename)
        .extension()
        .and_then(OsStr::to_str)
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
