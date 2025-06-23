import glob
import os
from collections import defaultdict
from itertools import product
#Note: This hasn't been thoroughly tested because I changed my approach

#Takes in up to two xstreme folders (input an empty string to use only one xstreme folder)
#Takes in two sets. Every element in set1 will be chained with every element in set2.
def get_motif_loc_dict(data_dir1, data_dir2, set1, set2) :
    # Holds where each motif is located {MOTIF: {acr: [loc, ..., loc] acr:[loc, ..., loc]}}
    motif_loc_dict_set1 = defaultdict(lambda: defaultdict(list)) #holds set 1 ACRs
    motif_loc_dict_set2 = defaultdict(lambda: defaultdict(list)) #holds set 2 ACRs

    print("Filling Motif Dict")

    pattern1 = os.path.join(data_dir1, 'fimo_out_*/', 'fimo.tsv')
    pattern2 = os.path.join(data_dir2, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob.glob(pattern1)
    fimo_files.extend(glob.glob(pattern2))

    
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
                motif = line_arr[0]
                if "-" in motif:
                    motif = motif[motif.index("-") + 1:]
                acr = line_arr[2]
                location = int(line_arr[3])
                if acr in set1:
                    motif_loc_dict_set1[motif][acr].append(location)
                if acr in set2:
                    motif_loc_dict_set2[motif][acr].append(location)


    # Remove duplicates from repeated sequences (basically remove overlaps)
    #dictionary 1
    for motif, single_motif_dict in motif_loc_dict_set1.items() :
        motif_len = len(motif)
        for acr, loc_list in single_motif_dict.items() :

            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_motif_dict[acr] = [x for x in loc_list if x not in to_remove]
    #dictionary 2
    for motif, single_motif_dict in motif_loc_dict_set2.items() :
        motif_len = len(motif)
        for acr, loc_list in single_motif_dict.items() :

            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_motif_dict[acr] = [x for x in loc_list if x not in to_remove]

    return (motif_loc_dict_set1, motif_loc_dict_set2)


# Finds chain length between two sets (global)
# Takes in the dictionaries from previous function and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors where acr1 is from set 1
# and acr2 is from set 2
def chain_two_sets(motif_loc_dict1, motif_loc_dict2, out_file) :

    print("Finding Anchors")

    anchor_dict = defaultdict(list)
    start_time = time.time()
    count = 0

    for motif_name, single_motif_dict in motif_loc_dict1.items():
        if motif_name in motif_loc_dict2:
            single_motif_dict2 = motif_loc_dict_2[motif_name]
            # Get the Cartesian product of the two value lists
            for acr1, acr2 in product(single_motif_dict.keys(), single_motif_dict2.keys()):
                # Get the Cartesian product of the two value lists
                anchor_dict[(acr1, acr2)].extend(list(product(single_motif_dict[acr1], single_motif_dict2[acr2])))
        count += 1
        print(f"Finished {count}, motif: {motif_name} after {time.time() - start_time} seconds", flush=True)

    
    print("Chaining")
    start_time = time.time()

    num_chained = 0
    with open(out_file, 'w') as out:
        for pair, anchors in anchor_dict.items() :
            num_chained += 1
            if num_chained % 1000000 == 0 :
                print(f"Finished {num_chained} chains after {time.time() - start_time} seconds", flush=True)

            out.write(f"{pair[0]}\t{pair[1]}\t{chain_driver(anchors, False)}\t{len(anchors)}\n")


#for testing
data_dir1 = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr3_3759930to3761697._xstreme"
data_dir2 = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr3_3825364to3826552._xstreme"
set1 = set()
set2 = set()
with open("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr3_3759930to3761697.group.txt") as file:
    for line in file:
        set1.add(line.strip())
with open("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr3_3825364to3826552.group.txt") as file:
    for line in file:
        set2.add(line.strip())

dicts = get_motif_loc_dict(data_dir1, data_dir2, set1, set2)


#get_motif_loc_dict("/home/mwarr/one_genome_xstreme", "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr3_3825364to3826552._xstreme", set(), set())