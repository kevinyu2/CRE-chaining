import sys
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
sys.path.append("../")
from pathlib import Path
from collections import defaultdict
import numpy as np
import random
from chaining import chain_driver, chain_local_driver
import os
import time
import itertools

'''
Step 3 in the pipeline.

Uses the anchor data to chain all 79,000 regions (79 genomes with 1,000 ACRs each). Outputs a 
file for each ACR pair. Each file contains every genome pair for the given ACR with the chain length
and the number of anchors.
'''


# Batch helper
def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


####################################################################
# Inter

# Worker function for parallel chaining
def single_pair_chain(args) :
    input_file, output_dir = args
    # {(genome1, genome2) : [anchors]}
    genome_pairs = defaultdict(list)
    with open(input_file) as acr_combo_file:
        #file is formatted: genome1 \t location \t ... ## ... genome4 \t location \t
        for line1,line2 in itertools.zip_longest(*[acr_combo_file]*2):
            line_arr = line1.split("##")
            if len(line_arr) == 1:
                continue
            motif_match_1 = line_arr[0].strip().split("\t")
            motif_match_2 = line_arr[1].strip().split("\t")

            try:
                assert(len(motif_match_1) % 2 == 0 and len(motif_match_2) % 2 == 0)
            except:
                print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(input_file)}")
            
            #consider all pairs of motifs from 1 and 2
            for i in range(0, len(motif_match_1), 2):
                for j in range(0, len(motif_match_2), 2):
                    genome_pair = (motif_match_1[i], motif_match_2[j])

                    pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()))
                        
                    genome_pairs[genome_pair].append(pair_indices)
    
    #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
    with open(f"{output_dir}/{input_file.stem}.tsv", "w") as output_file:
        # Keep track of seen items
        seen_dict = {}
        for pair, anchors in genome_pairs.items():
            key = tuple(anchors)
            if key in seen_dict :
                output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(anchors)}\n")
            else :
                chain_len = chain_driver(anchors, False)
                seen_dict[key] = chain_len
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(anchors)}\n")


# driver for parallel chaining
def run_chaining_all(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    acr_combos = list(search_dir.glob('*'))

    # Do this to include the output director as an arg
    args = [(file, output_dir) for file in acr_combos]

    with Pool(cpu_count()) as pool:
        for _ in tqdm(pool.imap_unordered(single_pair_chain, args), total=len(acr_combos), miniters = 1000):
            pass

##############################################################
# Weighted

# Worker for the weighted version
def single_pair_chain_weighted(args) :
    input_file, output_dir = args

    genome_pairs = defaultdict(list)
        
    with open(input_file) as acr_combo_file:
        #file is formatted: genome1 \t location \t ... ## ... genome4 \t location \t
        for line1,line2 in itertools.zip_longest(*[acr_combo_file]*2):
            line_arr = line1.split("##")
            # if len(line_arr) == 1:
            #     continue
            motif_match_1 = line_arr[0].strip().split("\t")
            motif_match_2 = line_arr[1].strip().split("\t")
            score = float(line2.strip())

            try:
                assert(len(motif_match_1) % 2 == 0 and len(motif_match_2) % 2 == 0)
            except:
                print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(input_file)}")
            
            #consider all pairs of motifs from 1 and 2
            for i in range(0, len(motif_match_1), 2):
                for j in range(0, len(motif_match_2), 2):
                    genome_pair = (motif_match_1[i], motif_match_2[j])

                    pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()), score)
                        
                    genome_pairs[genome_pair].append(pair_indices)
    
    #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
    with open(f"{output_dir}/{input_file.stem}.tsv", "w") as output_file:
        # Keep track of seen items
        seen_dict = {}
        for pair, anchors in genome_pairs.items():
            key = tuple(anchors)
            if key in seen_dict :
                output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(anchors)}\n")
            else :
                chain_len = chain_driver(anchors, True)
                seen_dict[key] = chain_len
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(anchors)}\n")

# Driver for the weighted version
def run_chaining_all_weighted(input_root_dir, output_dir):
    search_dir = Path(input_root_dir)
    acr_combos = list(search_dir.glob('*'))

    # Do this to include the output director as an arg
    args = [(file, output_dir) for file in acr_combos]

    with Pool(cpu_count()) as pool:
        for _ in tqdm(pool.imap_unordered(single_pair_chain_weighted, args), total=len(acr_combos)):
            pass
   
##################################################################
# Local version

def single_pair_chain_local(args) :
    input_file, output_dir, match, mismatch, gap = args
    # {(genome1, genome2) : [anchors]}
    genome_pairs = defaultdict(list)
    with open(input_file) as acr_combo_file:
        #file is formatted: genome1 \t location \t ... ## ... genome4 \t location \t
        for line1,line2 in itertools.zip_longest(*[acr_combo_file]*2):
            line_arr = line1.split("##")
            if len(line_arr) == 1:
                continue
            motif_match_1 = line_arr[0].strip().split("\t")
            motif_match_2 = line_arr[1].strip().split("\t")

            try:
                assert(len(motif_match_1) % 2 == 0 and len(motif_match_2) % 2 == 0)
            except:
                print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(input_file)}")
            
            #consider all pairs of motifs from 1 and 2
            for i in range(0, len(motif_match_1), 2):
                for j in range(0, len(motif_match_2), 2):
                    genome_pair = (motif_match_1[i], motif_match_2[j])

                    pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()))
                        
                    genome_pairs[genome_pair].append(pair_indices)
    
    #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
    with open(f"{output_dir}/{input_file.stem}.tsv", "w") as output_file:
        # Keep track of seen items
        seen_dict = {}
        for pair, anchors in genome_pairs.items():
            key = tuple(anchors)
            if key in seen_dict :
                output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key][0]}\t{seen_dict[key][1]}\n")
            else :
                chain_len = chain_local_driver(anchors, match, mismatch, gap, False)
                seen_dict[key] = chain_len
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len[0]}\t{chain_len[1]}\n")


# Driver for the local version
def run_chaining_all_local(input_root_dir, match, mismatch, gap, output_dir):
    search_dir = Path(input_root_dir)
    acr_combos = list(search_dir.glob('*'))

    # Do this to include the output director as an arg
    args = [(file, output_dir, match, mismatch, gap,) for file in acr_combos]

    with Pool(cpu_count()) as pool:
        for _ in tqdm(pool.imap_unordered(single_pair_chain_local, args), total=len(acr_combos)):
            pass

########################################################

# run_chaining_all('/home/mwarr/Data/Anchors_min1', '/home/mwarr/Data/Chaining_min1_par')

run_chaining_all_local("/home/mwarr/Data/Anchors_min1_local", 5, -2, -1, "/home/mwarr/Data/Chaining_min1_local_par")




















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
            file.write(f"\n{random.random()*2}\n")