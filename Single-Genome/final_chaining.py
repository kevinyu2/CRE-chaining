# Parallelizes chaining but NOT finding anchors

from pathlib import Path
from glob import glob
import sys
sys.path.append("../")
from chaining import chain_driver, chain_local_driver
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import defaultdict
from itertools import combinations, product
import numpy as np
import os
import math
import time
import pickle
import shutil



'''
The one genome pipeline, parallelized.

Uses preprocessed data (motifs for each genome for each ACR) and the motif pair data 
to output the anchors for each pair of ACRs across all genomes
'''
# Takes in an xstreme folder     
# Gets motif locations
def get_motif_loc_dict(data_dir) :
    # Holds where each motif is located {MOTIF: {acr: [loc, ..., loc] acr:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))

    print("Filling Motif Dict")

    pattern = os.path.join(data_dir, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob(pattern)

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
                motif_loc_dict[line_arr[0].split('-')[1]][line_arr[2]].append(int(line_arr[3]))
    
    # Remove duplicates from repeated sequences (basically remove overlaps)
    for motif, single_motif_dict in motif_loc_dict.items() :
        motif_len = len(motif)
        for acr, loc_list in single_motif_dict.items() :

            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_motif_dict[acr] = [x for x in loc_list if x not in to_remove]

    return motif_loc_dict


def get_motif_loc_dict_local(data_dir) :
    # This calls the original one, but we must change the indices 
    start_dict = get_motif_loc_dict(data_dir)

    # Stores all locations of each acr as a set for numbering
    # {acr: {loc, loc, loc...}}
    acr_dict = defaultdict(set)

    for motif, single_motif_dict in start_dict.items() :
        for acr, loc_list in single_motif_dict.items() :
            acr_dict[acr].update(loc_list)

    # Matches positions to motif number
    # {acr: {start_pos : motif_no}}
    acr_no_dict = defaultdict(dict)
    for acr, locs in acr_dict.items() :
        acr_no_dict[acr] = {v: i for i, v in enumerate(sorted(locs))}
    
    # The localized motif dict to return
    motif_loc_dict = defaultdict(lambda: defaultdict(list))
    
    for motif, single_motif_dict in start_dict.items() :
        for acr, loc_list in single_motif_dict.items() :
            for loc in loc_list :
                motif_loc_dict[motif][acr].append(acr_no_dict[acr][loc])


    return motif_loc_dict




# Find anchors given a loc dict. For local and global
def find_anchors(motif_loc_dict) :
    print("Finding Anchors")

    anchor_dict = defaultdict(list)
    start_time = time.time()

    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :
        tqdm.write(f"Processing: {motif_name}")

        for acr1, acr2 in combinations(single_motif_dict.keys(), 2):
            # Get the Cartesian product of the two value lists
            if acr1 < acr2 :
                anchor_dict[(acr1, acr2)].extend(list(product(single_motif_dict[acr1], single_motif_dict[acr2])))
            else :
                anchor_dict[(acr2, acr1)].extend(list(product(single_motif_dict[acr2], single_motif_dict[acr1])))
    
    with open(f"{start_time}_anchor_dict.pkl", "wb") as f:
        pickle.dump(anchor_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    return anchor_dict



# Parallelized chaining worker
def single_pair_chaining(args):
    pair, anchors = args
    result = chain_driver(anchors, False)
    return f"{pair[0]}\t{pair[1]}\t{result}\t{len(anchors)}\n"



# Finds pairwise chain length (global)
# Takes in a dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors
def chain_all_pairs(full_anchor_dict, out_file) :

    print("Chaining")

    # Parallelized chaining
    pair_items = list(full_anchor_dict.items())
    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        for line in tqdm(pool.imap_unordered(single_pair_chaining, pair_items), total=len(pair_items)):
            out.write(line)
           
# Parallelized local chaining worker
def single_pair_local_chaining(args):
    pair, anchors, match, mismatch, gap = args
    result = chain_local_driver(anchors, match, mismatch, gap, False)
    return f"{pair[0]}\t{pair[1]}\t{result[0]}\t{result[1]}\n"



# Local version
# Printed file in out_file in the format acr1\tacr2\tchain_score\tchain_len\tno_anchors
def chain_all_pairs_local(full_anchor_dict, match, mismatch, gap, out_file) :

    print("Chaining")

    parallelized_inputs = [(pair, anchors, match, mismatch, gap) for pair, anchors in full_anchor_dict.items()]
    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        for line in tqdm(pool.imap_unordered(single_pair_local_chaining, parallelized_inputs), total=len(parallelized_inputs)):
            out.write(line)

# If you have intermediate data (pickled anchor dict director or pickled full dict), can pass the directory name / file location in and it will skip that step
def chain_global_driver(input_dir, out_file, full_anchor_dict_pkl = None) :
    if full_anchor_dict_pkl == None :
        motif_l_dict = get_motif_loc_dict(input_dir)
        full_anchor_dict = find_anchors(motif_l_dict)
    else :
        full_anchor_dict = pickle.load(full_anchor_dict_pkl)
    chain_all_pairs(full_anchor_dict, out_file)


def chain_local_driver(input_dir, out_file, match, mismatch, gap, full_anchor_dict_pkl = None) :
    if full_anchor_dict_pkl == None :
        motif_l_dict = get_motif_loc_dict_local(input_dir)
        full_anchor_dict = find_anchors(motif_l_dict)
    else :
        full_anchor_dict = pickle.load(full_anchor_dict_pkl)
    chain_all_pairs_local(full_anchor_dict, match, mismatch, gap, out_file)


                             
chain_global_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "./Chaining_one_par.txt")

# chain_local_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/mwarr/Data/Chaining_one_local.tsv", 5, -2, -1)

# motif_l_dict = {'ABC' : {"a" : [3, 4, 12], "b" : [6], "c" : [10]},
#                 'BCD' : {"a" : [5], "c" : [6], "d" : [12]},
#                 'CDE' : {"b" : [12, 17], "a" : [11], "d" : [4]}}


# chain_all_pairs(motif_l_dict, "blah")


# anchor_dict = {('A', 'B') : [(4, 5), (6, 7), (8, 9), (1, 10)], 
#                ('B', 'C') : [(1, 15), (2, 2)], 
#                ('A', 'C') : [(6, 7), (8, 9)]}

# chain_all_pairs(anchor_dict, 'Test_chain.txt')


# merge_anchors('anchor_pkl_dir')