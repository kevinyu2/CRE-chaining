
from collections import defaultdict
from pathlib import Path
from glob import glob
from itertools import islice
import numpy as np
import time

##############################################################################
# Settings

# Initial kmer size to set as a seed
KMER_SIZE = 4

# ACR set to consider. If all, set to None
ACR_SET = '/home/mwarr/Data/seta.txt'
# ACR_SET = None

# Ignore BREAK motifs
IGNORE_BREAK = True

# How much larger the sequence window is to the ACR window
WINDOW_SIZE_MULTIPLIER = 2

# Scoring
MATCH = 5
MISMATCH = -2
GAP = -1

# Score cutoff
SCORE_THRESH = 100

#############################################################################

# Returns: dict for kmer locations, dict to match indexes to locations, list of ACRs for each chromosome
# kmer dict's index has the index of the first element
def parse_input(input_folder) :
    # {(KMER) : [(chromosome, index), (chromosome, index)]}
    kmer_dict = defaultdict(list)

    # {(chromosome, index) : location}
    location_dict = {}

    # {chromosome: [motif, motif, ..., motif]}
    motif_lists = defaultdict(list)

    print("Reading Input Sequence")


    # Create the motif lists and location_dict
    for file_path in Path(input_folder).iterdir():
        chromosome = str(file_path).rstrip().split('/')[-1]
        with open(file_path, 'r') as motif_file :

            motif_count = 0
            for line in motif_file :
                if IGNORE_BREAK :
                    if 'BREAK\n' not in line:
                        # should be [motif, location]
                        line_arr = line.rstrip().split('\t')

                        # Fill locaiton dict
                        location_dict[(chromosome, motif_count)] = int(line_arr[1])

                        # fill in motif_list
                        motif_lists[chromosome].append(line_arr[0])

                        motif_count += 1

    # Make them numpy instead
    for chr, motif_list in motif_lists.items() :
        motif_lists[chr] = np.array(motif_list)
                        
    for chr, motif_list in motif_lists.items() :
        list_len = len(motif_list)
        iterators = [islice(motif_list, i, None) for i in range(KMER_SIZE)]
        for idx, kmer in enumerate(zip(*iterators)) :
            kmer_dict[kmer].append((chr, idx))

    return kmer_dict, location_dict, motif_lists


def find_matches(acr_file, outfile, kmer_dict, location_dict, motif_lists) :
    print("Finding Matches")
    kmer_count = 0

    # If restricting by set
    if ACR_SET != None :
        acrs_to_use = set()
        with open(ACR_SET, 'r') as acr_set :
            for line in acr_set :
                acrs_to_use.add(line.rstrip())

    # {ACR: [motif, motif...]}
    acr_dict = defaultdict(list)
    with open(acr_file, 'r') as acrf :
        curr_acr = ""
        for line in acrf :
            # End of one
            if 'ACR: ' in line:
                # RESET
                curr_acr = line.rstrip().split('ACR: ')[1]

            else :
                if ACR_SET != None :
                    if curr_acr in acrs_to_use :
                        acr_dict[curr_acr].append(line.rstrip())
                else :
                    acr_dict[curr_acr].append(line.rstrip())

    # Get reverse:
    for acr in list(acr_dict.keys()) :
        acr_dict[f"{acr}_rev"] = acr_dict[acr][::-1]

    # Calculate number of seeds
    for acr, acr_motif_list in acr_dict.items() :
        iterators = [islice(acr_motif_list, i, None) for i in range(KMER_SIZE)]
        for idx, kmer in enumerate(zip(*iterators)) :

            kmer_count += len(kmer_dict[kmer])


    print(f"No. seeds: {kmer_count}")

    start_time = time.time()
    # Store lengths of each chromosome's motif sequences
    sequence_lengths = {}
    for chr, m_list in motif_lists.items() :
        sequence_lengths[chr] = len(m_list)
            
    done_seeds = 0
    hits = 0
    # Actually match the seeds and extend
    with open(outfile, 'w') as out :
        # Print out info
        out.write(f"kmer_length: {KMER_SIZE}\twindow_size: {WINDOW_SIZE_MULTIPLIER}\tscore_thresh: {SCORE_THRESH}\tmatch: {MATCH}\tmismatch: {MISMATCH}\tgap: {GAP}\tignore_break: {IGNORE_BREAK}\tacr_set: {ACR_SET}\n")
        out.write("ACR\tScore\tChr\tStart\tStop\n")
        for acr, acr_motif_list in acr_dict.items() :
            iterators = [islice(acr_motif_list, i, None) for i in range(KMER_SIZE)]
            len_acr_motif = len(acr_motif_list)
            for acr_idx, kmer in enumerate(zip(*iterators)) :

                for chr, sequence_idx in kmer_dict[kmer] :
                    # Get start and stop locations for query. Don't let bounds extend beyond the bounds of the array
                    # Note: start is inclusive, stop is not (as is python custom)
                    start = max(0, sequence_idx - (acr_idx * WINDOW_SIZE_MULTIPLIER))
                    stop = min(sequence_lengths[chr], sequence_idx + KMER_SIZE + ((len_acr_motif - (acr_idx + KMER_SIZE)) * WINDOW_SIZE_MULTIPLIER))

                    print(f"{chr}, {sequence_idx} | {acr}, {acr_idx} : {kmer} | {start} to {stop}")

                    local_score, local_start, local_stop = local_align(motif_lists[chr][start:stop], acr_motif_list)

                    if local_score >= SCORE_THRESH :
                        hits += 1
                        out.write(f"{acr}\t{local_score}\t{chr}\t{location_dict[(chr, start + local_start)]}\t{location_dict[(chr, start + local_stop)]}\n")

                    done_seeds += 1
                    if done_seeds % 1000 == 0 :
                        print(f"Done: {done_seeds} in {time.time() - start_time} seconds ({hits} hits)")


# SW alignment
def local_align(list_seq, list_acr) :
    l1 = len(list_seq)
    l2 = len(list_acr)
    score_matrix = np.zeros((l1 + 1, l2 + 1), dtype=int)
    
    max_score = 0
    max_pos = (0, 0)

    # Fill the scoring matrix
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            if list_seq[i - 1] == list_acr[j - 1]:
                score = MATCH
            else:
                score = MISMATCH

            diag = score_matrix[i - 1][j - 1] + score
            up   = score_matrix[i - 1][j] + GAP
            left = score_matrix[i][j - 1] + GAP

            score_matrix[i][j] = max(0, diag, up, left)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback to find start position
    i, j = max_pos
    end_i, end_j = i, j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current = score_matrix[i][j]
        if list_seq[i - 1] == list_acr[j - 1]:
            score = MATCH
        else:
            score = MISMATCH

        if current == score_matrix[i - 1][j - 1] + score:
            i -= 1
            j -= 1
        elif current == score_matrix[i - 1][j] + GAP:
            i -= 1
        else:
            j -= 1

    start_i, start_j = i, j

    return max_score, start_i, end_i - 1

##################################################################################
# Drivers

kmer_dict, location_dict, motif_lists = parse_input('/home/kyu/test_blast/test_chrs')
find_matches('/home/kyu/test_blast/test_acr.txt', '/home/kyu/test_blast/test_out_nonpar.tsv', kmer_dict, location_dict, motif_lists)


# kmer_dict, location_dict, motif_lists = parse_input('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motif_sequences')

# find_matches('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/acr_motif_sequences.txt', 
#              '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits.tsv', 
#              kmer_dict, location_dict, motif_lists)



# print(kmer_dict[('AAACCCTAAWT', 'AAACCCTAAWT', 'AAACCCTAAWT')])