import os
import random
from collections import defaultdict
from pathlib import Path
from Bio import Align


MAX_NUMS = 100


def get_alignment(seq1, seq2):
    aligner = Align.PairwiseAligner(scoring="blastn", mode="global")
    return aligner.score(seq1, seq2)


def create_seq_dict(seq_input_dir):
    seq_dict = {}
    search_dir = Path(seq_input_dir)
    for file in search_dir.glob("*"):
        with open(file) as seq_file:
            line = seq_file.readline()
            while line and line[0] == ">":
                sequence = seq_file.readline().strip()
                seq_name = line[1:].strip()
                seq_dict[seq_name] = sequence
                line = seq_file.readline()
    return seq_dict


def get_data(input_dir, seq_input_dir):
    
    seq_dict = create_seq_dict(seq_input_dir)
    chain_len_dict = defaultdict(list)
    
    ACRs = os.listdir(input_dir)
    random.shuffle(ACRs)
    
    count = 0
    for file in ACRs:
        with open(f"{input_dir}/{file}") as ACR_file:
            line_count = 0
            for line in ACR_file:
                line_arr = line.split("\t")
                if len(line_arr) < 2:
                    continue
                seq_name_1 = line_arr[0]
                seq_name_2 = line_arr[1]
                chain_len = int(line_arr[2])
                if len(chain_len_dict[chain_len]) < MAX_NUMS:
                    seq_1 = seq_dict[seq_name_1]
                    seq_2 = seq_dict[seq_name_2]
                    chain_len_dict[chain_len].append(get_alignment(seq_1, seq_2))
                line_count += 1
        count += 1
        print(f"Processed {count} files")
    print(chain_len_dict)
        

get_data("../Chaining_min1_intra", "../ACR_sequences")