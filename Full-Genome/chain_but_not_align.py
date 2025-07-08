
from collections import defaultdict
from pathlib import Path
from glob import glob
from itertools import islice
import numpy as np
import time
from multiprocessing import shared_memory
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from Bio import Align


##############################################################################
# Settings

# Additional chian scoring threshold for cutoff
SCORE_THRESH = 100

# How similar to be to ignore it
SIMILARITY_CUTOFF = 0.9

# SEQUENCE_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa'
SEQUENCE_FILE = '/home/kyu/Test_scripts/test_full_genome.fa'

# Hits
HITS_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_ArabidopsisPBM.tsv'

chrs_to_use = {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}

##############################################################################
# [nt, nt...]
global_seq_array = None
global_seq_shm = None
global_seq_name = 'seq'

BATCH_SIZE = 2000

def parse_sequence() :
    sequence_dict = defaultdict(list)
    sequence_sizes = defaultdict(int)
    with open(SEQUENCE_FILE, 'r') as seq_file :
        curr_chr = ""
        to_use = True
        for line in seq_file :
            if '>' in line :
                curr_chr = line.rstrip().split('>')[1]
                if curr_chr in chrs_to_use :
                    to_use = True
                else :
                    to_use = False
            elif to_use :
                sequence_dict[curr_chr].append(line.rstrip())
                sequence_sizes[curr_chr] += len(line.rstrip())
    
    total_size = 0
    for single_size in sequence_sizes.values() :
        total_size += single_size


    try: 
        shm = shared_memory.SharedMemory(name='seq')
        shm.unlink()  # This removes the shared memory segment
        shm.close()
    except:
        _ = 5

    # Create the shared memory
    global global_seq_array, global_seq_shm
    global_seq_shm = shared_memory.SharedMemory(create=True, size=total_size * 1, name = global_seq_name)  # 1 byte each
    global_seq_array = np.ndarray(total_size, dtype='S1', buffer=global_seq_shm.buf)

    # {Chr: (start, end)}
    seq_loc_dict = {}
    curr_start = 0
    for chr, size in sequence_sizes.items() :
        full_seq = ''.join(sequence_dict[chr])
        global_seq_array[curr_start:curr_start + size] = np.frombuffer(bytes(full_seq, encoding='ascii'), dtype='S1')
        seq_loc_dict[chr] = (curr_start, curr_start + size)
        curr_start += size

    return seq_loc_dict
    


# Final cleanup step for shared memory
def cleanup_global_arrays():
    global_seq_shm.close()
    global_seq_shm.unlink()

def parse_input() :
    # To ensure no repeats
    lines_dict = set()

    # (unbatched) inputs: [(qchr, qstart, qend, rchr, rstart, rend)]
    queries = []
    with open(HITS_FILE, 'r') as hits :
        next(hits)
        next(hits)
        for line in hits :
            if line not in lines_dict :
                lines_dict.add(line)
                line_arr = line.rstrip().split('\t')
                qchr = line_arr[0].split('_')[0]
                qstart = int(line_arr[0].split('_')[1].split('to')[0])
                qend = int(line_arr[0].rstrip().split('to')[1])
                queries.append((qchr, qstart, qend, line_arr[2], int(line_arr[3]), int(line_arr[4])))

def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]



seq_loc_dict = parse_sequence()
parse_input()
cleanup_global_arrays()