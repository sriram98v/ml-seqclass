# Maximum Likelihood Metagenomic Classification

This repository contains a suffix tree based metagenomic classifier using maximum likelihood. 

## Installation

To install you must first have cargo and rustup installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After installing the above command you can run the following to install seq_class:
```bash
cargo install --git=https://github.com/sriram98v/ml-seqclass
```

Alternatively, you can install seq_class by cloning this repository and building it locally:
```bash
git clone https://github.com/sriram98v/ml-seqclass
cd ml-seqclass
cargo install --path=./
```

## Usage
In order to start classifying a set of reads, you need to build an index from you set of reference sequences. You can build an index using the following command.
```bash
seq_class -s <reference sequence file path> -r <path to fastq file> -p <percent mismatch> -o <output file> -t <num threads (optional)>
```
This will produce two JSON files and one ```.npy``` file.

Then to find maximum likelihood matches of each read to a reference, use the python script provided
```bash
usage: get-proportions.py [-h] [-o OUTPUT_DIR] [-n NUM_ITER]

EM proportion estimations

options:
  -h, --help            show this help message and exit
  -o, --output_dir OUTPUT_DIR
                        Output directory
  -n, --num_iter NUM_ITER
                        Output directory
```

You can refer the man page for seq_class for more details by running
```bash
seq_class -h
```
