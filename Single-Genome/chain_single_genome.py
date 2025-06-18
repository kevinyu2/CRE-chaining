from pathlib import Path
import sys
sys.path.append("../")
from chaining import chain_driver, chain_local_driver

from collections import defaultdict
import itertools
import numpy as np
import glob
import os
import math
import time

'''
The one genome pipeline.

Uses preprocessed data (motifs for each genome for each ACR) and the motif pair data 
to output the anchors for each pair of ACRs across all genomes
'''
# Takes in an xstreme folder     
# Gets motif locations
def get_motif_loc_dict(data_dir) :
    # Holds where each motif is located {MOTIF: [(acr, loc), (acr, loc)], MOTIF2 ...}
    motif_loc_dict = defaultdict(list)

    print("Filling Motif Dict")

    pattern = os.path.join(data_dir, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob.glob(pattern)

    for fimo_file in fimo_files:
        
        # Fill in the dict
        with open(fimo_file, "r") as f:
            # Ignore first line
            next(f)
            for line in f: 
                line_arr = line.rstrip().split('\t')
                # The file ends tsv info early
                if len(line_arr) < 5: 
                    break
                motif_loc_dict[line_arr[0]].append((line_arr[2], int(line_arr[3])))

    return motif_loc_dict


def get_motif_loc_dict_local(data_dir) :
    # Holds set of positions for each ACR. These will be converted to local start positions
    # {acr: {loc, loc, loc...}}
    acr_dict = defaultdict(set)

    pattern = os.path.join(data_dir, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob.glob(pattern)

    for fimo_file in fimo_files:
        
        with open(fimo_file, "r") as f:
            # Ignore first line
            next(f)
            for line in f: 
                line_arr = line.rstrip().split('\t')
                if len(line_arr) < 5: 
                    break
                acr_dict[line_arr[2]].add(int(line_arr[3]))

    # {acr: {start_pos : motif_no}}
    acr_no_dict = defaultdict(dict)
    for acr, locs in acr_dict.items() :
        acr_no_dict[acr] = {v: i for i, v in enumerate(sorted(locs))}


    # Holds where each motif is located {MOTIF: [(acr, loc), (acr, loc)], MOTIF2 ...}
    motif_loc_dict = defaultdict(list)

    print("Filling Motif Dict")

    for fimo_file in fimo_files:
        
        # Fill in the dict
        with open(fimo_file, "r") as f:
            # Ignore first line
            next(f)
            for line in f: 
                line_arr = line.rstrip().split('\t')
                if len(line_arr) < 5: 
                    break
                motif_loc_dict[line_arr[0]].append((line_arr[2], acr_no_dict[line_arr[2]][int(line_arr[3])]))

    return motif_loc_dict

# Finds pairwise chain length (global)
# Takes in a dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors
def chain_all_pairs(motif_loc_dict, out_file) :

    print("Finding Anchors")

    anchor_dict = defaultdict(list)
    start_time = time.time()
    count = 0

    for motif_name, loc_list in motif_loc_dict.items() :
        for (acr1, loc1), (acr2, loc2) in itertools.combinations(loc_list, 2): 
            if acr1 < acr2 :
                anchor_dict[(acr1, acr2)].append((loc1, loc2))
            elif acr1 > acr2 :
                anchor_dict[(acr2, acr1)].append((loc2, loc1))
        count += 1
        print(f"Finished {count} motifs after {time.time() - start_time} seconds", flush=True)

    print("Chaining")
    start_time = time.time()

    num_chained = 0
    with open(out_file, 'w') as out:
        for pair, anchors in anchor_dict.items() :
            num_chained += 1
            if num_chained % 1000000 == 0 :
                print(f"Finished {num_chained} chains after {time.time() - start_time} seconds", flush=True)

            out.write(f"{pair[0]}\t{pair[1]}\t{chain_driver(anchors, False)}\t{len(anchors)}\n")

# Local version
# Printed file in out_file in the format acr1\tacr2\tchain_score\tchain_len\tno_anchors
def chain_all_pairs_local(motif_loc_dict, match, mismatch, gap, out_file) :

    print("Finding Anchors")

    anchor_dict = defaultdict(list)
    start_time = time.time()
    count = 0

    for motif_name, loc_list in motif_loc_dict.items() :
        for (acr1, loc1), (acr2, loc2) in itertools.combinations(loc_list, 2): 
            if acr1 < acr2 :
                anchor_dict[(acr1, acr2)].append((loc1, loc2))
            elif acr1 > acr2 :
                anchor_dict[(acr2, acr1)].append((loc2, loc1))
        count += 1
        print(f"Finished {count} motifs after {time.time() - start_time} seconds", flush=True)

    print("Chaining")
    start_time = time.time()

    num_chained = 0
    with open(out_file, 'w') as out:
        for pair, anchors in anchor_dict.items() :
            num_chained += 1
            if num_chained % 1000000 == 0 :
                print(f"Finished {num_chained} chains after {time.time() - start_time} seconds", flush=True)

            chain_scores = chain_local_driver(anchors, match, mismatch, gap, False)
            out.write(f"{pair[0]}\t{pair[1]}\t{chain_scores[0]}\t{chain_scores[1]}\t{len(anchors)}\n")



def chain_global_driver(input_dir, out_file) :
    motif_l_dict = get_motif_loc_dict(input_dir)
    chain_all_pairs(motif_l_dict, out_file)


def chain_local_driver(input_dir, out_file, match, mismatch, gap) :
    motif_l_dict = get_motif_loc_dict_local(input_dir)
    chain_all_pairs_local(motif_l_dict, match, mismatch, gap, out_file)


                             
# chain_global_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/mwarr/Data/Chaining_one.tsv")

chain_local_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/mwarr/Data/Chaining_one_local.tsv", 5, -2, -1)


# motif_l_dict = {'ABC' : [("a", 4), ("b", 6), ("a", 3), ("c", 10)], 
#                 'BCD' : [("a", 5), ("c", 6), ("d", 6)],
#                 'CDE' : [("a", 11), ("d", 4), ("b", 11), ("b", 12)]}  

