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


###########################################################
# Global dict

full_anchor_dict = defaultdict(list)


###########################################################
# Helpers for parallelization

def init_worker(d):
    global full_anchor_dict
    full_anchor_dict = d

def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

#######################################################
# Motif dict creation

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


######################################################################
# Anchors

# Find anchors given a loc dict. For local and global
def find_anchors(motif_loc_dict) :
    print("Finding Anchors")

    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :
        tqdm.write(f"Processing: {motif_name}")

        for acr1, acr2 in combinations(single_motif_dict.keys(), 2):
            # Get the Cartesian product of the two value lists
            if acr1 < acr2 :
                full_anchor_dict[(acr1, acr2)].extend(list(product(single_motif_dict[acr1], single_motif_dict[acr2])))
            else :
                full_anchor_dict[(acr2, acr1)].extend(list(product(single_motif_dict[acr2], single_motif_dict[acr1])))
    
    # with open(f"{start_time}_anchor_dict.pkl", "wb") as f:
        # pickle.dump(anchor_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


####################################################################
# Global chaining


# Parallelized chaining worker
def batch_pair_chaining(pairs):
    results = []
    for pair in pairs:
        chain_len = chain_driver(full_anchor_dict[pair], False)
        results.append(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(full_anchor_dict[pair])}\n")

    return results

# Finds pairwise chain length (global)
# Uses global anchor dict dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors
def chain_all_pairs(out_file) :

    print("Chaining")

    # Parallelized chaining
    pair_items = list(full_anchor_dict.keys())
    batches = list(chunkify(pair_items, 10000000))
    with Pool(cpu_count(), initializer=init_worker, initargs=(full_anchor_dict,)) as pool, open(out_file, 'w') as out:
        for lines in tqdm(pool.imap_unordered(batch_pair_chaining, batches), total=len(batches)):
            out.writelines(lines)



#######################################################################
# Local chaining


# Parallelized chaining worker
def batch_pair_local_chaining(args):
    pairs, match, mismatch, gap = args
    results = []
    for pair in pairs:
        chain_results = chain_local_driver(full_anchor_dict[pair], match, mismatch, gap, False)
        results.append(f"{pair[0]}\t{pair[1]}\t{chain_results[0]}\t{chain_results[1]}\n")

    return results


# Local version
# Printed file in out_file in the format acr1\tacr2\tchain_score\tchain_len\tno_anchors
def chain_all_pairs_local(match, mismatch, gap, out_file) :

    print("Chaining")
    pair_items = list(full_anchor_dict.keys())
    batches = list(chunkify(pair_items, 10000000))

    # Add match, mismatch, and gap
    parallelized_inputs = [(batch, match, mismatch, gap) for batch in batches]
    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        for lines in tqdm(pool.imap_unordered(batch_pair_local_chaining, parallelized_inputs), total=len(parallelized_inputs)):
            out.writelines(lines)


##############################################################
# Drivers

# If you have intermediate data (pickled anchor dict director or pickled full dict), can pass the directory name / file location in and it will skip that step
def chain_global_driver(input_dir, out_file, full_anchor_dict_pkl = None) :
    if full_anchor_dict_pkl == None :
        motif_l_dict = get_motif_loc_dict(input_dir)
        find_anchors(motif_l_dict)
    else :
        full_anchor_dict = pickle.load(full_anchor_dict_pkl)
    chain_all_pairs(out_file)


def chain_local_driver(input_dir, out_file, match, mismatch, gap, full_anchor_dict_pkl = None) :
    if full_anchor_dict_pkl == None :
        motif_l_dict = get_motif_loc_dict_local(input_dir)
        find_anchors(motif_l_dict)
    else :
        full_anchor_dict = pickle.load(full_anchor_dict_pkl)
    chain_all_pairs_local(match, mismatch, gap, out_file)

##########################################################################

                             
chain_global_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "./Chaining_one_par.txt")

# chain_local_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/mwarr/Data/Chaining_one_local.tsv", 5, -2, -1)

# motif_l_dict = {'ABC' : {"a" : [3, 4, 12], "b" : [6], "c" : [10]},
#                 'BCD' : {"a" : [5], "c" : [6], "d" : [12]},
#                 'CDE' : {"b" : [12, 17], "a" : [11], "d" : [4]}}
# find_anchors(motif_l_dict)



# full_anchor_dict = {('A', 'B') : [(4, 5), (6, 7), (8, 9), (1, 10)], 
#                ('B', 'C') : [(1, 15), (2, 2)], 
#                ('A', 'C') : [(6, 7), (8, 9)]}
# chain_all_pairs('./test_chain_out_par.tsv')


