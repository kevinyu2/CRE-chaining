from collections import defaultdict
import random 

# Size and chromosomes to model
HITS_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_DAPv1.tsv'

SEQUENCE_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa'

OUT_FILE ='/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_RANDOM_DAPv1.tsv'

chrs_to_use = {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}

# Get the sizes of all the sequences to not break
sequence_sizes = defaultdict(int)
with open(SEQUENCE_FILE, 'r') as seq_file :
    curr_chr = ""
    for line in seq_file :
        if '>' in line :
            curr_chr = line.rstrip().split('>')[1]
        else :
            sequence_sizes[curr_chr] += len(line.rstrip())

line_dict = set()

with open(HITS_FILE, 'r') as hits, open(OUT_FILE, 'w') as out :
    out.write(f"RANDOM-{hits.readline()}")
    out.write(f"{hits.readline()}")

    for line in hits :
        if line not in line_dict :
        
            line_arr = line.rstrip().split('\t')
            hit_len = int(line_arr[4]) - int(line_arr[3])
            chr = line_arr[2]
            if chr in chrs_to_use :

                # print(f"{chr}\t{sequence_sizes[chr]}\t{sequence_sizes[chr] - hit_len}")

                random_number = random.randint(0, sequence_sizes[chr] - hit_len - 1)

                out.write(f"{line_arr[0]}\t{line_arr[1]}\t{chr}\t{random_number}\t{random_number + hit_len}\n")
            line_dict.add(line)

