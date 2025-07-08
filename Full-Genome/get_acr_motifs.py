from collections import defaultdict
import os
from glob import glob
import heapq

def get_motif_loc_dict(data_dir) :
    # Holds where each motif is located {acr: {motif: [loc, ..., loc], motif:[loc, ..., loc]}}
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
                if int(line_arr[3]) not in motif_loc_dict[line_arr[2]][line_arr[0].split('-')[1]] :
                    motif_loc_dict[line_arr[2]][line_arr[0].split('-')[1]].append(int(line_arr[3]))
    
    # Remove duplicates from repeated sequences (basically remove overlaps)
    for acr, single_acr_dict in motif_loc_dict.items() :
        for motif, loc_list in single_acr_dict.items() :
            motif_len = len(motif)
            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_acr_dict[motif] = [x for x in loc_list if x not in to_remove]

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

def print_sequences(motif_loc_dict, outfile) :
    with open(outfile, 'w') as out :
        for acr, single_acr_dict in motif_loc_dict.items() :
            out.write(f"ACR: {acr}\n")
            result_list = merge_dict_lists(single_acr_dict)
            for motif_loc_pair in result_list :
                out.write(f'{motif_loc_pair[1]}\n')


mld = get_motif_loc_dict("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/")
print_sequences(mld, '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/acr_motif_sequences.txt')