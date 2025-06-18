import sys
sys.path.append("../")
from pathlib import Path
from collections import defaultdict
import numpy as np
import random
from chaining import chain_driver, chain_local_driver
import os
import time
import itertools
from concurrent.futures import ProcessPoolExecutor

def run_one(args):
    input_root_dir, output_dir, acr_combo = args
    inner_loop(input_root_dir, output_dir, acr_combo)


# Different because input file is slightly different
def run_chaining_all_intra(input_root_dir, output_dir):

    search_dir = Path(input_root_dir)
    acr_combos = list(search_dir.glob('*'))

    jobs = [(input_root_dir, output_dir, acr) for acr in acr_combos]


    #loop through all files in parallel (combinations of ACRs)
    count = 0
    start_time = time.time()


    with ProcessPoolExecutor(max_workers = 120) as executor:
        for count, _ in enumerate(executor.map(run_one, jobs), start=1):
            if count % 100 == 0:
                elapsed = time.time() - start_time
                print(f"Finished {count} files after {elapsed:.2f} seconds", flush=True)

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
        for pair, anchors in genome_pairs.items():
            key = tuple(anchors)
            if key in seen_dict :
                output_file.write(f"{pair[0]}\t{pair[1]}\t{seen_dict[key]}\t{len(anchors)}\n")
            else :
                chain_len = chain_driver(anchors, False)
                seen_dict[key] = chain_len
                output_file.write(f"{pair[0]}\t{pair[1]}\t{chain_len}\t{len(anchors)}\n")


run_chaining_all_intra('/home/mwarr/Data/Anchors_min1_intra/', '/home/mwarr/Data/parallel_test')