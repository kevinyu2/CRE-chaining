import pandas as pd
from pathlib import Path
import numpy as np

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
    distances = []
    with open(acr_file, 'r') as f :
        for line in f :
            line_arr = line.split('\t')
            max1 = max_length_dict[line_arr[0]]
            max2 = max_length_dict[line_arr[1]]
            longest = max(max1, max2)
            shortest = min(max1, max2)
            distances.append(alpha * ((longest - int(line_arr[2]))/longest) + (1 - alpha) * ((shortest - int(line_arr[2]))/shortest))
    return np.mean(distances)

# Driver for pairwise distance file inter
def get_all_distances_inter(alpha, chaining_folder, outfile, max_length_dict) :
    with open(outfile, 'w') as out:
        # out.write('ACR1\tACR2\tDistance\n')
        for chain_file in Path(chaining_folder).rglob("*"):
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
# get_all_distances_intra(0.5, '/home/mwarr/Data/Chaining_min1_intra/', '/home/mwarr/Data/Clustering/Distances_min1_intra_alpha50', max_length_dict)
get_all_distances_inter(0.5, '/home/mwarr/Data/Chaining_min1_mini/', '/home/mwarr/Data/Clustering/Distances_min1_mini_alpha50.tsv', max_length_dict)