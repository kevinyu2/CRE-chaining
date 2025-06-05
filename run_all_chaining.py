from pathlib import Path
from collections import defaultdict
import numpy as np

def run_chaining_all(root_dir):

    genome_pairs = defaultdict(list)

    search_dir = Path(root_dir)
    #loop through all files
    for acr_combo in search_dir.glob('*'):
        with open(acr_combo) as acr_combo_file:
            for line in acr_combo_file:
                line_arr = line.split("\t")
                #get the pair and sort it
                pair = [line_arr[0], line_arr[2]]
                pair_sorted = tuple(sorted(pair))

                pair_indices = (line_arr[1], line_arr[3].strip())
                #see if the indices need to be swapped to match up with how the 
                #pair was sorted
                pair_argsort = np.argsort(pair)
                pair_order = (line_arr[1], line_arr[2])
                
                if pair_order[0] == 1:
                    #swap the order
                    temp = pair_indices[0]
                    pair_indices[0] = pair_indices[1]
                    pair_indices[1] = temp 
                    
                
                genome_pairs[pair_sorted].append(pair_indices)
    print(genome_pairs)

#with open("chr1_chr2.txt", "w") as file:
    #file.write("abc\t5\tdef\t6\nabc\t4\thij\t3\nklm\t3\tabc\t7\nhij\t2\tdef\t7")
run_chaining_all("test_folder")
