# Finds the motifs that were ACTUALLY used

from pathlib import Path
from glob import glob
import sys
from collections import defaultdict
from itertools import combinations, product
import os
import numpy as np

def get_motif_loc_dict(data_dir, out_file) :
    motif_loc_dict = defaultdict(set)

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
                motif_loc_dict[line_arr[0]].add(line_arr[1])
                if len(motif_loc_dict[line_arr[0]]) > 1 :
                    print(line)
                    
    
    with open(out_file, 'w') as out :
        for motif, locs in motif_loc_dict.items() :
            if len(locs) > 1 :
                print(f"Multiple sources: {motif}")
            else :
                item = list(locs)[0]
                if 'STREME' in item :
                    to_write = 'streme_out/streme.xml'
                else :
                    to_write = 'meme_out/meme.xml'
                out.write(f"{motif}\t{to_write}\t{item}\n")

get_motif_loc_dict("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/", "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motifs.txt")