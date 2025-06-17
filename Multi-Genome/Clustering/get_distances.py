import pandas as pd
from pathlib import Path
import numpy as np
import time


# alpha ((longest - chain)/longest) + (1 - alpha) ((shortest - chain)/shortest) 
def get_distances_intra(alpha, max_length_dict, acr_file, outfile) :
    with open(outfile, 'w') as outfile :
        with open(acr_file, 'r') as infile :
            for line in infile:
                line_arr = line.rstrip().split('\t')

                if len(line_arr) < 2 :
                    outfile.write('File skipped due to runtime\n')
                    return
                if line_arr[0] != line_arr[1] :
                    chain_len = int(line_arr[2])
                    longest = max(max_length_dict[line_arr[0]], max_length_dict[line_arr[1]])
                    shortest = min(max_length_dict[line_arr[0]], max_length_dict[line_arr[1]])
                    outfile.write(f"{line_arr[0]}\t{line_arr[1]}\t{alpha * ((longest - chain_len)/longest) + (1 - alpha) * ((shortest - chain_len)/shortest)}\n")

# Driver for getting all distances in an intra folder
def get_all_distances_intra(alpha, input_folder, output_folder, max_length_dict) :
    for chain_file in Path(input_folder).rglob("*"):
        get_distances_intra(alpha, max_length_dict, chain_file, Path(output_folder) / (chain_file.stem + '.tsv'))


def get_distance_inter(alpha, max_length_dict, acr_file) :
    total = 0.0
    with open(acr_file, 'r') as f :
        lines = f.readlines()
    for line in lines :
        genome1, genome2, chain, anchors = line.split('\t')
        max1 = max_length_dict[genome1]
        max2 = max_length_dict[genome2]
        if max1 > max2:
            longest = max1
            shortest = max2
        else:
            longest = max2
            shortest = max1
        total += alpha * ((longest - int(chain))/longest) + (1 - alpha) * ((shortest - int(chain))/shortest)
    return total/ len(lines)

# Driver for pairwise distance file inter
def get_all_distances_inter(alpha, chaining_folder, outfile, max_length_dict) :

    with open(outfile, 'w') as out:
        start_time = time.time()
        count = 0
        for chain_file in Path(chaining_folder).rglob("*"):
            count += 1
            if count % 1000 == 0 :
                print(f"Finished {count} files after {time.time() - start_time} seconds", flush=True)
            acr_1 = chain_file.stem.split('_Chr')[0]
            acr_2 = 'Chr' + chain_file.stem.split('_Chr')[1]
            out.write(f"{acr_1}\t{acr_2}\t{get_distance_inter(alpha, max_length_dict, chain_file)}\n")


def make_max_length_dict(dict_file) :
    motif_df = pd.read_csv(dict_file, sep = '\t')
    motifs_dict = {}

    for i, row in motif_df.iterrows() :
        motifs_dict[row['genome']] = int(row['tomtom_motif_count_unique_start'])

    return motifs_dict

max_length_dict = make_max_length_dict('/home/mwarr/motif_counts.tsv')
get_all_distances_intra(0, '/home/mwarr/Data/Chaining_min1_intra/', '/home/mwarr/Data/Clustering/Distances_min1_intra_alpha0', max_length_dict)
#get_all_distances_inter(0, '/home/mwarr/Data/Chaining_min1/', '/home/mwarr/Data/Clustering/Distances_min1_alpha0.tsv', max_length_dict)