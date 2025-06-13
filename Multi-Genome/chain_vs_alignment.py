import os
import random
from collections import defaultdict
from pathlib import Path
from Bio import Align
from matplotlib import pyplot as plt


MAX_NUMS = 100


def get_alignment(seq1, seq2):
    aligner = Align.PairwiseAligner(scoring="blastn", mode="global")
    return aligner.score(seq1, seq2)


#returns a dictionary mapping each ACR with its sequence
def create_seq_dict(seq_input_dir):
    seq_dict = {}
    search_dir = Path(seq_input_dir)
    for file in search_dir.glob("*.group.fa"):
        with open(file) as seq_file:
            line = seq_file.readline()
            while line and line[0] == ">":
                sequence = seq_file.readline().strip()
                seq_name = line[1:].strip()
                seq_dict[seq_name] = sequence
                line = seq_file.readline()
    return seq_dict


def get_data(input_dir, seq_input_dir, output_dir):
    
    seq_dict = create_seq_dict(seq_input_dir)
    
    ACRs = os.listdir(input_dir)
    random.shuffle(ACRs)
    
    count = 0
    for file in ACRs:
        with open(f"{output_dir}/{file}", "w") as output_file:
            with open(f"{input_dir}/{file}") as ACR_file:
                output_file.write("Genome1\tGenome2\tChain_len\tAlignment_score")
                for line in ACR_file:
                    line_arr = line.split("\t")
                    if len(line_arr) < 2:
                        continue
                    seq_name_1 = line_arr[0]
                    seq_name_2 = line_arr[1]
                    chain_len = int(line_arr[2])
                    seq_1 = seq_dict[seq_name_1].upper()
                    seq_2 = seq_dict[seq_name_2].upper()
                    alignment = get_alignment(seq_1, seq_2)
                    output_file.write(f"\n{seq_name_1}\t{seq_name_2}\t{chain_len}\t{alignment}")
        count += 1
        print(f"Processed {count} files")

def create_graph(input_dir, output_dir):
    alignment_scores = []
    chain_lengths = []

    search_dir = Path(input_dir)
    for file_name in search_dir.glob("*.tsv"):
        with open(file_name) as file:
            next(file)
            for line in file:
                line_arr = line.split("\t")
                chain_lengths.append(float(line_arr[2]))
                alignment_scores.append(float(line_arr[3]))

    plt.scatter(chain_lengths, alignment_scores, s=.1)
    plt.xlabel("Chain Lengths")
    plt.ylabel("Alignment Scores")
    plt.title("Chain Length vs Alignment Score -- Intra-ACR")
    plt.savefig(f"{output_dir}/chain_vs_align_intra_min1.png")



# INPUT_DIR = "/home/mwarr/Data/Chaining_min1_intra"  
# SEQ_INPUT_DIR = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes"     
# OUTPUT_DIR = "/home/mwarr/Data/Align_vs_Chain_min1_intra"

# get_data(INPUT_DIR, SEQ_INPUT_DIR, OUTPUT_DIR)

INPUT_DIR = "/home/mwarr/Data/Align_vs_Chain_min1_intra"
OUTPUT_DIR = "../Visualizations"

create_graph(INPUT_DIR, OUTPUT_DIR)