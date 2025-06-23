# Parallelizes chaining but NOT finding anchors

from pathlib import Path
from glob import glob
import sys
sys.path.append("../")
from chaining import chain_driver_np, chain_local_driver_np
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import defaultdict
from itertools import combinations, product
import os
from multiprocessing import Manager
from multiprocessing import shared_memory
import numpy as np


'''
The one genome pipeline, parallelized.

Uses preprocessed data (motifs for each genome for each ACR) and the motif pair data 
to output the anchors for each pair of ACRs across all genomes
'''


###########################################################
# Global memory variables

global_anchor_array = None
global_anchor_shm = None
global_anchor_name = 'anchors'
BATCH_SIZE = 100000

###########################################################
# Helpers for parallelization

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
                if int(line_arr[3]) not in motif_loc_dict [line_arr[0].split('-')[1]][line_arr[2]] :
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

#####################################################################
# Calculate number of anchors


def init_anchor_array(total_anchors):
    global global_anchor_array, global_anchor_shm
    global_anchor_shm = shared_memory.SharedMemory(create=True, size=total_anchors * 2 * 4, name = global_anchor_name)  # 2 ints per anchor, 4 bytes each
    global_anchor_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=global_anchor_shm.buf)

def cleanup_anchor_array():
    global global_anchor_shm
    global_anchor_shm.close()
    global_anchor_shm.unlink()

# Calculate the total number of anchors
# Returns an int, the number of anchors, and a dict, which maps each acr to its start and stop points in the anchor space
def calculate_no_anchors(motif_loc_dict) :
    print("Calculating Anchor Dict Allocation Space")

    # Running total needed room
    total = 0
    # To track individual space needed
    # {(acr1, acr2) : space}
    acr_space_dict = defaultdict(int)
    # To give start and stop indices
    # {(acr1, acr2) : (start, stop)}
    acr_start_stop_dict = defaultdict(lambda: (0, 0))
    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :
        sizes = {acr: len(single_motif_dict[acr]) for acr in list(single_motif_dict.keys())}

        for acr1, acr2 in combinations(single_motif_dict.keys(), 2) :
            key = (acr1, acr2) if acr1 < acr2 else (acr2, acr1)
            pair_count = sizes[acr1] * sizes[acr2]
            acr_space_dict[key] += pair_count
            total += pair_count
    
    curr_start = 0
    for key, space in acr_space_dict.items() :
        acr_start_stop_dict[key] = (curr_start, curr_start + space)
        curr_start += space
    
    return total, acr_start_stop_dict


######################################################################
# Anchors

# Find anchors given a loc dict. For local and global
def find_anchors(anchor_loc_dict, motif_loc_dict) :
    print("Finding Anchors")

    curr_filled_per_pair = defaultdict(int)

    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :

        for acr1, acr2 in combinations(single_motif_dict.keys(), 2):
            key = (acr1, acr2) if acr1 < acr2 else (acr2, acr1)

            anchors1 = single_motif_dict[acr1 if acr1 < acr2 else acr2]
            anchors2 = single_motif_dict[acr2 if acr1 < acr2 else acr1]

            # Write to the correct slice of the global array
            n1, n2 = len(anchors1), len(anchors2)
            n_total = n1 * n2

            anchors = np.empty((n_total, 2), dtype=np.int32)
            start_idx = anchor_loc_dict[key][0] + curr_filled_per_pair[key]
            end_idx = start_idx + n_total

            i = 0
            for a1 in anchors1:
                for a2 in anchors2:
                    anchors[i, 0] = a1
                    anchors[i, 1] = a2
                    i += 1

            # One write
            global_anchor_array[start_idx:end_idx] = anchors
            curr_filled_per_pair[key] += n_total


####################################################################
# Global chaining


# Parallelized chaining worker
def batch_pair_chaining(args):
    pairs, total_anchors = args
    shm = shared_memory.SharedMemory(name='anchors')

    # Create NumPy array wrapper (no copy)
    shared_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=shm.buf)
    results = []
    for acr1, acr2, start, stop in pairs:        
        chain_len = chain_driver_np(shared_array[start:stop], False)
        results.append(f"{acr1}\t{acr2}\t{chain_len}\t{stop - start}\n")

    return results

# Finds pairwise chain length (global)
# Uses global anchor dict dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors
def chain_all_pairs(anchor_loc_dict, total_anchors, out_file) :

    print("Chaining")

    # Parallelized chaining
    anchor_items = [(a, b, start, stop) for (a, b), (start, stop) in anchor_loc_dict.items()]

    batches = list(chunkify(anchor_items, BATCH_SIZE))
    parallelized_inputs =  [(batch, total_anchors) for batch in batches]

    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        try :
            for lines in tqdm(pool.imap_unordered(batch_pair_chaining, parallelized_inputs), total=len(batches)):
                out.writelines(lines)
        finally:
            pool.close()
            pool.join()           


#######################################################################
# Local chaining


# Parallelized chaining worker
def batch_pair_local_chaining(args):
    pairs, total_anchors, match, mismatch, gap = args
    shm = shared_memory.SharedMemory(name='anchors')
    shared_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=shm.buf)

    results = []
    for acr1, acr2, start, stop in pairs:
        chain_results = chain_local_driver_np(shared_array[start:stop], match, mismatch, gap, False)
        results.append(f"{acr1}\t{acr2}\t{chain_results[0]}\t{chain_results[1]}\n")

    return results


# Local version
# Printed file in out_file in the format acr1\tacr2\tchain_score\tchain_len\tno_anchors
def chain_all_pairs_local(anchor_loc_dict, total_anchors, match, mismatch, gap, out_file) :

    print("Chaining")
    anchor_items = [(a, b, start, stop) for (a, b), (start, stop) in anchor_loc_dict.items()]
    batches = list(chunkify(anchor_items, BATCH_SIZE))

    # Add match, mismatch, and gap
    parallelized_inputs =  [(batch, total_anchors, match, mismatch, gap) for batch in batches]
    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        try :
            for lines in tqdm(pool.imap_unordered(batch_pair_local_chaining, parallelized_inputs), total=len(parallelized_inputs)):
                out.writelines(lines)
        finally :
            pool.close()
            pool.join()


##############################################################
# Drivers

# If you have intermediate data (pickled anchor dict director or pickled full dict), can pass the directory name / file location in and it will skip that step
# def chain_global_driver(input_dir, out_file) :
#     motif_l_dict = get_motif_loc_dict(input_dir)
#     find_anchors(motif_l_dict)
#     chain_all_pairs(out_file)



def chain_global_driver(input_dir, out_file) :
    motif_l_dict = get_motif_loc_dict(input_dir)
    space, anchor_loc_dict = calculate_no_anchors(motif_l_dict)
    init_anchor_array(space)
    find_anchors(anchor_loc_dict, motif_l_dict)
    chain_all_pairs(anchor_loc_dict, space, out_file)
    cleanup_anchor_array()

def chain_local_driver(input_dir, out_file, match, mismatch, gap) :
    motif_l_dict = get_motif_loc_dict_local(input_dir)
    space, anchor_loc_dict = calculate_no_anchors(motif_l_dict)
    init_anchor_array(space)
    find_anchors(anchor_loc_dict, motif_l_dict)
    chain_all_pairs_local(anchor_loc_dict, space, match, mismatch, gap, out_file)
    cleanup_anchor_array()

##########################################################################

# Test code
                             
chain_global_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "./Chaining_one_par_f.tsv")

# chain_local_driver("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/mwarr/Data/Chaining_one_local.tsv", 5, -2, -1)

# motif_l_dict = {'ABC' : {"a" : [3, 4, 12], "b" : [6], "c" : [10]},
#                 'BCD' : {"a" : [5], "c" : [6], "d" : [12]},
#                 'CDE' : {"b" : [12, 17], "a" : [11], "d" : [4]}}

# motif_l_dict = get_motif_loc_dict('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/')

# space, anchor_loc_dict = calculate_no_anchors(motif_l_dict)
# init_anchor_array(space)

# find_anchors(anchor_loc_dict, motif_l_dict)
# print(global_anchor_array)
# chain_all_pairs(anchor_loc_dict, space, 'test_chain_out_par.tsv')
# cleanup_anchor_array()

# print(full_anchor_dict)

# full_anchor_dict = {('A', 'B') : [(4, 5), (6, 7), (8, 9), (1, 10)], 
#                ('B', 'C') : [(1, 15), (2, 2)], 
#                ('A', 'C') : [(6, 7), (8, 9)]}
# chain_all_pairs('./test_chain_out_par.tsv')


