# Finds sequence of motifs to use in BLAST style search

from pathlib import Path
from glob import glob
import sys
import os
from collections import defaultdict
import numpy as np
import heapq


# How many base pairs to be considered a break
BREAK_LENGTH = 3000

def find_dict_len(dict) :
    sum = 0
    for k, v in dict.items() :
        for vk, vv in v.items() :
            sum += len(vv)
    print(sum)

def get_motif_loc_dict(data_dir) :
    # Holds where each motif is located {chr: {MOTIF: [loc, loc, ...], MOTIF:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(set))

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

                motif_loc_dict[line_arr[2]][line_arr[0]].add(int(line_arr[3]))

    # Turn them to lists
    for chr, single_chr_dict in motif_loc_dict.items() :
        for motif, loc_list in single_chr_dict.items() :
            single_chr_dict[motif] = list(loc_list)

    # find_dict_len(motif_loc_dict)
    # Remove duplicates from repeated sequences (basically remove overlaps)
    for chr, single_chr_dict in motif_loc_dict.items() :
        for motif, loc_list in single_chr_dict.items() :
            motif_len = len(motif)

            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_chr_dict[motif] = sorted([x for x in loc_list if x not in to_remove])

    # print(motif_loc_dict['Chr1']['AGAGAGAGAGAGAGA'])

    # find_dict_len(motif_loc_dict)

    # print(motif_loc_dict.keys())
    
    return motif_loc_dict

def merge_dict_lists(list_dict):
    heap = []
    result = []

    # Initialize the heap with the first element of each list
    for key, lst in list_dict.items():
        if lst:  # Skip empty lists
            heapq.heappush(heap, (lst[0], key, 0))  # (value, source_key, index_in_list)

    while heap:
        value, key, idx = heapq.heappop(heap)
        result.append((value, key))

        # If there are more items in this list, push the next one
        if idx + 1 < len(list_dict[key]):
            next_value = list_dict[key][idx + 1]
            heapq.heappush(heap, (next_value, key, idx + 1))

    return result

def print_sequences(motif_loc_dict, output_dir) :
    for chr, single_chr_dict in motif_loc_dict.items() :
        with open(f"{output_dir}/{chr}", 'w') as out :
            result_list = merge_dict_lists(single_chr_dict)
            previous = 10000000
            for motif_loc_pair in result_list :
                # IF we need to add a break
                if motif_loc_pair[0] - previous >= BREAK_LENGTH :
                    out.write("BREAK\n")
                out.write(f'{motif_loc_pair[1]}\t{motif_loc_pair[0]}\n')
                previous = motif_loc_pair[0]


motif_loc_dict = get_motif_loc_dict('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_ArabidopsisDAPv1/')
print_sequences(motif_loc_dict, '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/ArabidopsisDAPv1_fimo/')