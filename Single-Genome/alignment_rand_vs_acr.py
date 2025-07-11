from Bio import Align #type: ignore
import random
import time
import math
import sys
from pathlib import Path
from useful_functions import *


def time_test():
    bases = ['A', 'C', 'T', 'G']
    target = ""
    query = ""
    for i in range(1000):
        target += bases[random.randint(0, 3)]
        query += bases[random.randint(0, 3)]
    
    print("Beginning alignment")
    start_time = time.time()
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"
    score = aligner.score(target, query)
    print(score, flush=True)
    print(time.time() - start_time)


'''
Takes in a file with a list of regions and outputs a fasta file of those regions.
Regions should be formatted chr#_####to#### or chr#_####to####_rand
'''
def make_acr_fasta(input_file, genome_dict, output_file):
    with open(output_file, "w") as output:
        with open(input_file) as input:
            for line in input:
                if line != "\n":
                    chr = line[ : line.index("_")]
                    start = int(line[line.index("_") + 1 : line.index("to")])
                    if "_rand" in line:
                        stop = int(line[line.index("to") + 2 : line.index("_rand")])
                    stop = int(line[line.index("to") + 2 : ])
                    seq = genome_dict[chr][start - 1 : stop]
                output.write(f">{line.strip()}\n{seq.strip()}\n")

'''
Takes in a fasta file (input_file) and outputs <num_files> files to <output_dir>.
The output files together will contain the same sequences as the original fasta file.
<num_sequences> is the number of sequences in <input_file>.
'''
def split_fasta(input_file, output_dir, num_files, num_sequences):
    seq_per_file = math.ceil(num_sequences / num_files)
    
    with open(input_file) as input:
        id = input.readline()
        for i in range(num_files):
            with open(f"{output_dir}/temp_{i}", "w") as output:
                for _ in range(seq_per_file):
                    output.write(id)
                    line = input.readline()
                    while line and line != "\n" and line[0] != ">":
                        output.write(line)
                        line = input.readline()
                    if line == "\n" or not line:
                        break
                    id = line

'''
Takes in two fasta files (input_file and ref_file) and outputs the alignment score between
every sequence in <input_file> with every sequence in <ref_file> to <output_file>.
'''
def output_align(input_file, ref_file, output_file):
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"
    
    #read in files and create dictionaries
    ref_dict = create_genome_dict(ref_file)
    input_dict = create_genome_dict(input_file)
    
    output = []
    
    #pairwise align
    for id2, seq2 in input_dict.items():
        for id1, seq1 in ref_dict.items():
            score = aligner.score(seq1, seq2)
            output.append((id1, id2, score))

    #output to file
    with open(output_file, "w") as out_file:
        for item in output:
            out_file.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")

def get_lowest_align_scores(input_dir):
    min = 0
    for file_path in Path(input_dir).glob("*freq.tsv"):
        with open(file_path) as file:
            for line in file:
                line_arr = line.split("\t")
                if len(line_arr) < 2:
                    continue
                score = float(line_arr[0])
                if score < min:
                    min = score
    print(min)

def adjust_scores(adjustment, input_dir):
    for file_path in Path(input_dir).glob("*freq.tsv"):
        with open(f"{file_path.parent}/{file_path.stem}_2.tsv", "w") as output:
            print(f"{file_path.parent}/{file_path.stem}_2.tsv")
            with open(file_path) as file:
                for line in file:
                    line_arr = line.split("\t")
                    if len(line_arr) < 2:
                        continue
                    output.write(f"{float(line_arr[0]) + adjustment}\t{line_arr[1]}")


if __name__ == "__main__":
    output_align(sys.argv[1], sys.argv[2], sys.argv[3])
    #split_fasta("/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/seta_90.fa", "/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/seta_temp", 50, 28165)
