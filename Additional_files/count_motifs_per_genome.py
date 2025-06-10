from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import Align
import sys


def create_tomtom_set() :

    # Get all _xtreme folders
    search_dir = Path('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes')
    xstreme_dirs = [p for p in search_dir.glob('*') if p.is_dir() and p.name.endswith('._xstreme')]

    tomtom_set = set()

    for x_dir in xstreme_dirs:
        tomtom_path = x_dir / 'meme_tomtom_out' / 'tomtom.tsv'
        if tomtom_path.exists():
            try : 
                tomtom_df = pd.read_csv(tomtom_path, sep = '\t')

                tomtom_df = tomtom_df.dropna(subset=['Target_consensus'])
                for i, row in tomtom_df.iterrows() :
                    tomtom_set.add(row['Query_ID'])

            # Usually means the tsv file failed and is empty
            except :
                print(f"File failed: {tomtom_path}")

    return tomtom_set

def find_motif_counts(output_file, tomtom_set) :
    with open(output_file, "w") as f:

        search_dir = Path("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes")
        x_streme_folders = [folder for folder in search_dir.glob("*._xstreme") if folder.is_dir()]

        # {(acr, genome): (total_count, tomtom_count)}
        motif_count_dict = defaultdict(lambda: [0, 0])

        for x_treme_folder in x_streme_folders:
            # Find the acr name
            start_index = str(x_treme_folder).index("Chr")
            end_index = str(x_treme_folder).index("._xstreme")
            acr_name = str(x_treme_folder)[start_index:end_index]
            
            #Find all the motifs for each genome in this ACR and record
            for fimo_folder in x_treme_folder.glob("fimo*"):
                with open(f"{fimo_folder}/fimo.tsv", 'r') as read_file:
                    #read header
                    try:
                        read_file.readline()
                    except:
                        print(f"{fimo_folder}/fimmo.tsv is empty")

                    #iterate through each row in the table
                    for line in read_file:
                        line_arr = line.split("\t")
                        if len(line_arr) < 5:
                            break
                        #get only the sequence (the sequence is sometimes preceeded by a number)
                        index = 0
                        for char in line_arr[0]:
                            if char.isalpha():
                                break
                            else:
                                index += 1
                        sequence = line_arr[0][index:]

                        motif_count_dict[(acr_name, line_arr[2])][0] += 1
                        if str(sequence) in tomtom_set :
                            motif_count_dict[(acr_name, line_arr[2])][1] += 1


        f.write(f"acr\tgenome\ttotal_motif_count\ttomtom_motif_count\n")

        for pair in motif_count_dict.keys() :
            f.write(f"{pair[0]}\t{pair[1]}\t{motif_count_dict[pair][0]}\t{motif_count_dict[pair][1]}\n")


tomtom_set = create_tomtom_set()
print("Found tomtom set")
find_motif_counts('../motif_counts.tsv', tomtom_set)