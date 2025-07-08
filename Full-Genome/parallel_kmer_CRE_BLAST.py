
from collections import defaultdict
from pathlib import Path
from glob import glob
from itertools import islice
import numpy as np
import time
from multiprocessing import shared_memory
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

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


###########################################################################
# Global memory variables

# [[chromosome, index], [chromosome, index]] (all are ints)
global_kmer_array = None
global_kmer_shm = None
global_kmer_name = 'kmer'

# [motif, motif...]
global_motif_array = None
global_motif_shm = None
global_motif_name = 'motif'

BATCH_SIZE = 200

# Fills the global motif array and kmer array, returns the location dict
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

# Turn into global array
# Returns two things needed to reconstruct the kmer dict: a match to the start and end locations of each kmer and a dict to convert chromosomes to ints
def init_global_kmer(kmer_dict, motif_to_int) :
    print("Initializing Global Kmer")

    # First calculate the size
    kmer_array_size = 0
    for locs in kmer_dict.values() :
        kmer_array_size += len(locs)
    
    global global_kmer_array, global_kmer_shm
    global_kmer_shm = shared_memory.SharedMemory(create=True, size=kmer_array_size * 2 * 4, name = global_kmer_name)  # 4 bytes each
    global_kmer_array = np.ndarray((kmer_array_size, 2), dtype=np.int32, buffer=global_kmer_shm.buf)

        
    # Map locations, fill array
    kmer_array_locs = {}
    # Map chromosomes to ints for the global array
    chr_to_int = {}
    int_to_chr = {}

    curr_chr_num = 0
    curr_start = 0
    for kmer_motif, locs in kmer_dict.items() :
        kmer = tuple(motif_to_int[m] for m in kmer_motif)
        for i,loc in enumerate(locs):
            if loc[0] not in chr_to_int :
                chr_to_int[loc[0]] = curr_chr_num
                int_to_chr[curr_chr_num] = loc[0]
                curr_chr_num += 1

            global_kmer_array[curr_start + i][0] = chr_to_int[loc[0]]
            global_kmer_array[curr_start + i][1] = loc[1]

        kmer_array_locs[kmer] = (curr_start, curr_start + len(locs))
        curr_start += len(locs)


    return kmer_array_size, chr_to_int, int_to_chr, kmer_array_locs


# Turn into global array
# Returns mappings from motif to int, int to motif, and a place where all the motifs are
def init_global_motif(motif_lists) :
    print("Initializing Global Motif")

    # First calculate the size
    motif_array_size = 0
    for motif_list in motif_lists.values() :
        motif_array_size += len(motif_list)

    global global_motif_array, global_motif_shm
    global_motif_shm = shared_memory.SharedMemory(create=True, size=motif_array_size * 4, name = global_motif_name)  # 4 bytes each
    global_motif_array = np.ndarray(motif_array_size, dtype=np.int32, buffer=global_motif_shm.buf)

    # Maps motif location to shared array location
    motif_array_locs = {}
    # Convert to ints
    motif_to_int = {}
    int_to_motif = {}
    # Tracks motif number
    curr_motif_idx = 0
    curr_array_idx = 0
    for chr, motif_list in motif_lists.items() :
        motif_array_locs[chr] = (curr_array_idx, curr_array_idx + len(motif_list))
        for motif in motif_list :
            # Found a new motif
            if motif not in motif_to_int :
                motif_to_int[motif] = curr_motif_idx
                int_to_motif[curr_motif_idx] = motif
                curr_motif_idx += 1
            
            global_motif_array[curr_array_idx] = motif_to_int[motif]
            curr_array_idx += 1

    return motif_array_size, motif_to_int, int_to_motif, motif_array_locs
    

# Final cleanup step for shared memory
def cleanup_global_arrays():
    global global_kmer_shm, global_motif_shm
    global_kmer_shm.close()
    global_kmer_shm.unlink()

    global_motif_shm.close()
    global_motif_shm.unlink()



def batched_worker(args) :
    # acr_motif_lists is an array of pairs of (acr, [motif list])
    acr_motif_lists, kmer_array_size, kmer_array_locs, sequence_lengths, int_to_chr, motif_array_locs, motif_array_size, KMER_SIZE, WINDOW_SIZE_MULTIPLIER, SCORE_THRESH = args
    results = []

    # Open shared memory
    kmer_shm = shared_memory.SharedMemory(name='kmer')
    shared_kmer_array = np.ndarray((kmer_array_size, 2), dtype=np.int32, buffer=kmer_shm.buf)

    motif_shm = shared_memory.SharedMemory(name='motif')
    shared_motif_array = np.ndarray(motif_array_size, dtype=np.int32, buffer=motif_shm.buf)

    # Iterate thorugh all in batch
    for pair in acr_motif_lists :
        acr, motif_list = pair

        iterators = [islice(motif_list, i, None) for i in range(KMER_SIZE)]
        len_acr_motif = len(motif_list)
        for acr_idx, kmer in enumerate(zip(*iterators)) :
            if kmer in kmer_array_locs :
                for kmer_array_loc in range(kmer_array_locs[kmer][0], kmer_array_locs[kmer][1]) :
                    chr_int = shared_kmer_array[kmer_array_loc][0]
                    chr = int_to_chr[chr_int]
                    sequence_idx = shared_kmer_array[kmer_array_loc][1]

                    start = max(0, sequence_idx - (acr_idx * WINDOW_SIZE_MULTIPLIER))
                    stop = min(sequence_lengths[chr], sequence_idx + KMER_SIZE + ((len_acr_motif - (acr_idx + KMER_SIZE)) * WINDOW_SIZE_MULTIPLIER))

                    # print(motif_list)
                    # print(shared_motif_array[motif_array_locs[chr][0] + start : motif_array_locs[chr][0] + stop])
                    # print()

                    local_score, local_start, local_stop = local_align(shared_motif_array[motif_array_locs[chr][0] + start : motif_array_locs[chr][0] + stop], motif_list)
                    if local_score >= SCORE_THRESH :
                        results.append((acr, local_score, chr, start, stop, local_start, local_stop))
    return results

def batch_dict_items(d):
    items = list(d.items())
    for i in range(0, len(items), BATCH_SIZE):
        yield items[i:i + BATCH_SIZE]

def find_matches(motif_lists, acr_file, outfile, kmer_array_size, kmer_array_locs, int_to_chr, motif_array_locs, motif_array_size, location_dict, motif_to_int) :
    print("Finding Matches")

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
                        # Check if motif is found in the dict, if not, safe to just make it -1
                        if line.rstrip() in motif_to_int :
                            acr_dict[curr_acr].append(motif_to_int[line.rstrip()])
                        else :
                            acr_dict[curr_acr].append(-1)

                else :
                    if line.rstrip() in motif_to_int :
                        acr_dict[curr_acr].append(motif_to_int[line.rstrip()])
                    else :
                        acr_dict[curr_acr].append(-1)

    # Get reverse:
    for acr in list(acr_dict.keys()) :
        acr_dict[f"{acr}_rev"] = acr_dict[acr][::-1]

    # Store lengths of each chromosome's motif sequences
    sequence_lengths = {}
    for chr, m_list in motif_lists.items() :
        sequence_lengths[chr] = len(m_list)
            
    # Actually match the seeds and extend
    with open(outfile, 'w') as out :
        # Print out info
        out.write(f"kmer_length: {KMER_SIZE}\twindow_size: {WINDOW_SIZE_MULTIPLIER}\tscore_thresh: {SCORE_THRESH}\tmatch: {MATCH}\tmismatch: {MISMATCH}\tgap: {GAP}\tignore_break: {IGNORE_BREAK}\tacr_set: {ACR_SET}\n")
        out.write("ACR\tScore\tChr\tStart\tStop\n")            


        parallelized_inputs =  [(batch, kmer_array_size, kmer_array_locs, sequence_lengths, int_to_chr, motif_array_locs, motif_array_size, KMER_SIZE, WINDOW_SIZE_MULTIPLIER, SCORE_THRESH) for batch in batch_dict_items(acr_dict)]

        with Pool(cpu_count()) as pool:
            try :
                for lines in tqdm(pool.imap_unordered(batched_worker, parallelized_inputs), total=len(parallelized_inputs)):
                    # print(lines)
                    for line in lines :
                        acr, local_score, chr, start, stop, local_start, local_stop = line
                        out.write(f"{acr}\t{local_score}\t{chr}\t{location_dict[(chr, start + local_start)]}\t{location_dict[(chr, start + local_stop)]}\n")
            finally:
                pool.close()
                pool.join()




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



def run_all(input_folder, acr_file, outfile) :
    kmer_dict, location_dict, motif_lists = parse_input(input_folder)
    motif_array_size, motif_to_int, int_to_motif, motif_array_locs = init_global_motif(motif_lists)
    kmer_array_size, chr_to_int, int_to_chr, kmer_array_locs = init_global_kmer(kmer_dict, motif_to_int)
    find_matches(motif_lists, acr_file, outfile, kmer_array_size, kmer_array_locs, int_to_chr, motif_array_locs, motif_array_size, location_dict, motif_to_int)
    cleanup_global_arrays()

################################################################################

run_all('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motif_sequences',
        '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/acr_motif_sequences.txt',
        '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_par.tsv')
# find_matches('/home/kyu/test_blast/test_acr.txt', '/home/kyu/test_blast/test_out.tsv', kmer_dict, location_dict, motif_lists)


# kmer_dict, location_dict, motif_lists = parse_input('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motif_sequences')

# find_matches('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/acr_motif_sequences.txt', 
#              '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits.tsv', 
#              kmer_dict, location_dict, motif_lists)



# print(kmer_dict[('AAACCCTAAWT', 'AAACCCTAAWT', 'AAACCCTAAWT')])