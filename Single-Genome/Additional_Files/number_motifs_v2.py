from pathlib import Path
from glob import glob
import sys
from tqdm import tqdm
from collections import defaultdict
import os
import numpy as np

def get_motif_loc_dict(data_dir) :
    # Holds where each motif is located {MOTIF: {acr: [loc, ..., loc] acr:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))


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

    acr_num_motifs_dict = defaultdict(int)

    for motif, single_motif_dict in motif_loc_dict.items() :
        for acr, loc_list in single_motif_dict.items() :
            acr_num_motifs_dict[acr] += len(loc_list)

    with open('/home/mwarr/Data/single_genome_motif_counts_v2.tsv', 'w') as f: 
        for key, value in acr_num_motifs_dict.items() :
            f.write(f'{key}\t{value}\n')

get_motif_loc_dict('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/')