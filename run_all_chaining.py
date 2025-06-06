from pathlib import Path
from collections import defaultdict
import numpy as np
import random
from chaining import chain
import os
import time

def run_chaining_all(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    #loop through all files
    for acr_combo in search_dir.glob('*'):
        genome_pairs = defaultdict(list)
        with open(acr_combo) as acr_combo_file:
            for line in acr_combo_file:
                line_arr = line.split("##")
                motif_match_1 = line_arr[0].strip().split("\t")
                motif_match_2 = line_arr[1].strip().split("\t")

                try:
                    assert(len(motif_match_1) % 2 == 0 and len(motif_match_2) % 2 == 0)
                except:
                    print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(acr_combo)}")
                for i in range(0, len(motif_match_1), 2):
                    for j in range(0, len(motif_match_2), 2):
                        #get the pair and sort it
                        genome_pair = [motif_match_1[i], motif_match_2[j]]
                        genome_pair_sorted = tuple(sorted(genome_pair))

                        pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()))
                        #see if the indices need to be swapped to match up with how the 
                        #pair was sorted
                        pair_argsort = np.argsort(genome_pair)
                        if pair_argsort[0] == 1:
                            #swap the order
                            pair_indices = (pair_indices[1], pair_indices[0])
                            
                        genome_pairs[genome_pair_sorted].append(pair_indices)
        with open(f"{output_dir}/chr3_chr_4_output.tsv", "w") as output_file:
            for pair in genome_pairs.keys():
                chain_len = chain(genome_pairs[pair])
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\n")




def create_test_file(num_lines, num_per_line):
    with open("test_folder/chr3_chr4.txt", "w") as file: 
        for k in range(num_lines):
            for i in range(num_per_line):
                rand_str = ""
                for j in range(5):
                    rand_str += chr(random.randint(97, 122))
                file.write(f"{rand_str}\t{random.randint(1, 100)}\t")
            file.write("##")
            for i in range(num_per_line):
                rand_str = ""
                for j in range(5):
                    rand_str += chr(random.randint(97, 122))
                file.write(f"{rand_str}\t{random.randint(1, 100)}\t")
            file.write("\n")

create_test_file(4, 79000)
start_time = time.time()
run_chaining_all("test_folder", "test_output_folder")
print(f"Time: {time.time() - start_time}")
