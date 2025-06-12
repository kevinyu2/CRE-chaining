from collections import defaultdict
import time
from chaining import chain
import pandas as pd

# This algorithm finds matching k-mer start locations then discounts the ones that get exended by tracking with a 2D array
def find_mems(k, seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    # Add all kmers of seq 1 into a dict
    k_dict = defaultdict(list)
    for i in range(0, n - k + 1):
        kmer = seq1[i:i+k]
        k_dict[kmer].append(i)

    # Set to track k-mers
    precense_set = set()

    # List of mems
    mem_starts = []

    for j in range(0, m - k + 1):
        kmer = seq2[j:j+k]
        matches = k_dict[kmer]
        for i in matches:
            if (i-1, j-1) not in precense_set:
                mem_starts.append((i, j))

            precense_set.add((i, j))

    return mem_starts

#with open("seq1.txt") as seq1_file:
    #with open("seq2.txt") as seq2_file:
        #seq1_str = seq1_file.read().strip().replace('\n', "")
        #seq2_str = seq2_file.read().strip().replace('\n', "")
        #mems = find_mems(9, seq1_str,seq2_str)
        #print(chain(mems))

with open("a_thaliana_chr1.fa") as thal_file:
    thal_file.readline()
    seq = thal_file.read()
    df = pd.read_csv("accessible_regions.bed", sep="\t")
    
    start_1 = df["start"][0] - 2000
    end_1 = df["end"][0] + 2000
    region_1 = seq[start_1:end_1]
    
    start_2 = df["start"][1] - 2000
    end_2 = df["end"][1] + 2000
    region_2 = seq[start_2:end_2]

    mems = find_mems(9, region_1, region_2)
    print(len(mems))
    print(chain(mems))




