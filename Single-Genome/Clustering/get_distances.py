import pandas as pd
from pathlib import Path
import numpy as np
from collections import defaultdict
import time
from tqdm import tqdm


def make_max_length_dict(dict_file) :
    motifs_dict = {}

    with open(dict_file, 'r') as f:
        for line in f:
            line_arr = line.rstrip().split('\t')
    
            motifs_dict[line_arr[0]] = int(line_arr[1])

    return motifs_dict


def get_all_distances(alpha, chaining_file, outfile, max_length_dict) :

    with open(outfile, 'w') as out:
        with open(chaining_file, 'r') as input :
            for line in tqdm(input, miniters= 1000000) :
                genome1, genome2, chain, _ = line.rstrip().split('\t')
                max1 = max_length_dict[genome1]
                max2 = max_length_dict[genome2]
                if max1 > max2:
                    longest = max1
                    shortest = max2
                else:
                    longest = max2
                    shortest = max1
            
                out.write(f"{genome1}\t{genome2}\t{alpha * ((longest - int(chain))/longest) + (1 - alpha) * ((shortest - int(chain))/shortest)}\n")

# Here the beta is a linear multiplier: we are doing beta/chain score 
def get_all_distances_local(beta, chaining_file, outfile) :
    with open(outfile, 'w') as out:
        with open(chaining_file, 'r') as input :
            for line in tqdm(input, miniters= 1000000) :
                genome1, genome2, chain, _ = line.rstrip().split('\t')
                out.write(f"{genome1}\t{genome2}\t{beta / float(chain)}\n")

#########################################################

# max_length_dict = make_max_length_dict('/home/mwarr/Data/single_genome_motif_counts_v2.tsv')
# get_all_distances(0, '/home/kyu/CRE-chaining/Single-Genome/Chaining_one_par_f.tsv', '/home/mwarr/Data/distances_alpha0_one.tsv', max_length_dict)

get_all_distances_local(5, '/home/mwarr/Data/Chaining_one_local.tsv', '/home/mwarr/Data/distances_beta5_local_one.tsv')