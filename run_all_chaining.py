from pathlib import Path
from collections import defaultdict
import numpy as np
import random
from chaining import chain, chain_weighted
import os
import time
import itertools

'''
Step 3 in the pipeline.

Uses the anchor data to chain all 79,000 regions (79 genomes with 1,000 ACRs each). Outputs a 
file for each ACR pair. Each file contains every genome pair for the given ACR with the chain length
and the number of anchors.
'''

def run_chaining_all(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    #loop through all files (combinations of ACRs)
    count = 0
    start_time = time.time()
    for acr_combo in search_dir.glob('*'):
        count += 1
        if count % 100 == 0:
            print(f"Finished {count} files after {time.time() - start_time} seconds", flush=True)
        
        genome_pairs = defaultdict(list)
        with open(acr_combo) as acr_combo_file:
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
                    print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(acr_combo)}")
                
                #consider all pairs of motifs from 1 and 2
                for i in range(0, len(motif_match_1), 2):
                    for j in range(0, len(motif_match_2), 2):
                        genome_pair = (motif_match_1[i], motif_match_2[j])

                        pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()))
                            
                        genome_pairs[genome_pair].append(pair_indices)
        
        #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
        with open(f"{output_dir}/{acr_combo.stem}.tsv", "w") as output_file:
            # Keep track of seen items
            seen_dict = {}
            for pair in genome_pairs.keys():
                key = tuple(genome_pairs[pair])
                if key in seen_dict :
                    output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(genome_pairs[pair])}\n")
                else :
                    chain_len = chain(genome_pairs[pair])
                    seen_dict[key] = chain_len
                    output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(genome_pairs[pair])}\n")


            # for pair in genome_pairs.keys():
            #     chain_len = chain(genome_pairs[pair])
            #     output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(genome_pairs[pair])}\n")



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

def run_chaining_all_weighted(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    #loop through all files (combinations of ACRs)
    count = 0
    start_time = time.time()
    for acr_combo in search_dir.glob('*'):
        count += 1
        if count % 100 == 0:
            print(f"Finished {count} files after {time.time() - start_time} seconds", flush=True)
        genome_pairs = defaultdict(list)
        
        with open(acr_combo) as acr_combo_file:
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
                    print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(acr_combo)}")
                
                #consider all pairs of motifs from 1 and 2
                for i in range(0, len(motif_match_1), 2):
                    for j in range(0, len(motif_match_2), 2):
                        genome_pair = (motif_match_1[i], motif_match_2[j])

                        pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()), score)
                            
                        genome_pairs[genome_pair].append(pair_indices)
        
        #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
        with open(f"{output_dir}/{acr_combo.stem}.tsv", "w") as output_file:
            # Keep track of seen items
            seen_dict = {}
            for pair in genome_pairs.keys():
                key = tuple(genome_pairs[pair])
                if key in seen_dict :
                    output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(genome_pairs[pair])}\n")
                else :
                    chain_len = chain_weighted(genome_pairs[pair])
                    seen_dict[key] = chain_len
                    output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(genome_pairs[pair])}\n")



# Different because input file is slightly different
def run_chaining_all_intra(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    #loop through all files (combinations of ACRs)
    count = 0
    start_time = time.time()
    for acr_combo in search_dir.glob('*'):
        count += 1
        if count % 10 == 0:
            print(f"Finished {count} files after {time.time() - start_time} seconds", flush=True)
        # print(count, flush=True)
        # print(acr_combo, flush=True)
        inner_loop(input_root_dir, output_dir, acr_combo)


def inner_loop(input_root_dir, output_dir, acr_combo) :
    genome_pairs = defaultdict(list)
    with open(acr_combo) as acr_combo_file:
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
                print(f"ERROR! There are not the same number of genomes as locations in an entry in {str(acr_combo)}")
            

            #error if there are too many locations. The processing would take too long.
            if len(motif_match_1) > 30000 :
                with open(f"{output_dir}/{acr_combo.stem}.tsv", "w") as file:
                    file.write(f"File skipped due to runtime\nThis motif is in {len(motif_match_1)} locations in the genomes (too many!)")
                    print(f"{len(motif_match_1)} locations in the genomes. {acr_combo.stem} skipped due to runtime.")
                    return
                    
            #consider all pairs of motifs from 1 and 2
            for i in range(0, len(motif_match_1), 2):
                for j in range(0, len(motif_match_2), 2):
                    if motif_match_1[i] < motif_match_2[j] :
                        genome_pair = (motif_match_1[i], motif_match_2[j])
                        pair_indices = (int(motif_match_1[i+1]), int(motif_match_2[j+1].strip()))
                    else :
                        genome_pair = (motif_match_2[j], motif_match_1[i])
                        pair_indices = (int(motif_match_2[j+1].strip()), int(motif_match_1[i+1]))

                    genome_pairs[genome_pair].append(pair_indices)
                    # print(pair_indices)

    #output in this format: genome1 \t genome2 \t chain_length \t #_of_anchors
    with open(f"{output_dir}/{acr_combo.stem}.tsv", "w") as output_file:
        # Keep track of seen items
        seen_dict = {}
        for pair in genome_pairs.keys():
            key = tuple(genome_pairs[pair])
            if key in seen_dict :
                output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(genome_pairs[pair])}\n")
            else :
                chain_len = chain(genome_pairs[pair])
                seen_dict[key] = chain_len
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(genome_pairs[pair])}\n")



run_chaining_all_intra("../Anchors_min1_intra", "../Chaining_min1_intra")