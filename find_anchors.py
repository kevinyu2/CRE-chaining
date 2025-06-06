from pathlib import Path
from collections import defaultdict
import itertools
import numpy as np


# Should be tsv of motif1, motif2, start1, start2, score. Motif1 and motif2 are sorted in alphabetical order
# Returns a dict that pairs (motif1, motif2) : (start1, start2, score)
def read_motif_pairs(min_score):
    print("Reading Scores File")
    SCORES_FILE_NAME = "scores.tsv"

    pairwise_score_dict = {}

    with open(SCORES_FILE_NAME) as pair_scores:
        next(pair_scores)
        for line in pair_scores :
            line_arr = line.rstrip().split('\t')
            if float(line_arr[4]) > min_score:
                pairwise_score_dict[(line_arr[0], line_arr[1])] = (int(line_arr[2]), int(line_arr[3]), float(line_arr[4]))
    return pairwise_score_dict


# def find_anchors_single_genome_pair(data_dir, sequence1, genome1, sequence2, genome2, score_dict) :
#     file_path_1 = Path(data_dir) / sequence1 / (genome1 + ".tsv")
#     file_path_2 = Path(data_dir) / sequence2 / (genome2 + ".tsv")

#     anchors = []

#     # Has tuples of (motif, start)
#     genome_1_motif_list = []
#     with file_path_1.open('r') as f :
#         for line in f:
#             line_arr = line.rstrip().split('\t')
            
#             # OK to be inefficient in the small case, will just brute force
#             genome_1_motif_list.append((line_arr[0], line_arr[1]))
    
#     # Scan second for anchors
#     with file_path_2.open('r') as f :
#         for line in f:
#             line_arr = line.rstrip().split('\t')

#             for motif_loc_pair in genome_1_motif_list :
#                 if tuple(sorted([line_arr[0], motif_loc_pair[0]])) in score_dict :
#                     anchors.append((motif_loc_pair[1], line_arr[1], score_dict[tuple(sorted([line_arr[0], motif_loc_pair[0]]))][0]))

#     return anchors


# Finds pairwise for all chromosomes            
# Make sure to clear the output folder first since we are appending to files 
# Printed output is a file for each chr chr pair, each line has genome\tloc\tgenome\tloc\t...genome\tloc##genome\tloc\t....genome\tloc\nscore\n       
def find_anchors_all_pairs(data_dir, out_dir, score_dict) :
    # Holds where each motif is located {MOTIF: {chr : [(genome, loc), (genome, loc)]), chr2 : [...]}, MOTIF2 ...}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))

    print("Filling Motif Dict")

    for file_path in Path(data_dir).rglob("*"):
        if file_path.is_file() and file_path.suffix == ".tsv":
            # Save folder name and file name to give to the dict
            # Folder name is the acr (chr), file name is the genome
            folder_name = file_path.parent.name
            file_name = file_path.stem

            # Fill in the dict
            with open(file_path, "r") as f:
                # Ignore first line
                next(f)
                for line in f:
                    # Each line is motif, location
                    line_arr = line.rstrip().split('\t')
                    motif_loc_dict[line_arr[0]][folder_name].append((file_name, int(line_arr[1])))
    

    print("Finding Anchors")
    print(f"Number of Pairs Total: {len(score_dict.keys())}")
    progress_counter = 0

    # Faster to iterate through the keys of the match dict
    for pair in score_dict.keys() :
        if progress_counter % 10000 == 0 :
            if progress_counter != 0 :
                print(f"Done: {progress_counter}")
        progress_counter += 1
        
        ## REMINDER: {MOTIF: {chr : [(genome, loc), (genome, loc)]), chr2 : [...]}, MOTIF2 ...}
        # If pair is the same, iterate through combinations
        # We don't need to add start points
        if pair[0] == pair[1] :
            # Iterate through all pairs of acrs
            for acr1, acr2 in itertools.combinations(motif_loc_dict[pair[0]].keys(), 2):
                
                # No need to worry about same acrs because combinations
                sorted_pair = sorted([acr1, acr2])
                out_file = Path(out_dir) / (sorted_pair[0] + '_' + sorted_pair[1] + '.tsv')
                
                # Don't need to worry about offsets
                with open(out_file, "a") as f:
                    # Write all genome loc pairs 
                    for i, genome_and_loc in enumerate(motif_loc_dict[pair[0]][sorted_pair[0]]) :
                        # This guarantees the tabs are correct
                        if i == 0 :
                            f.write(f"{genome_and_loc[0]}\t{genome_and_loc[1]}")
                        else :
                            f.write(f"\t{genome_and_loc[0]}\t{genome_and_loc[1]}")
                    # separator
                    f.write("##")
                    for i, genome_and_loc in enumerate(motif_loc_dict[pair[1]][sorted_pair[1]]) :
                        if i == 0 :
                            f.write(f"{genome_and_loc[0]}\t{genome_and_loc[1]}")
                        else :
                            f.write(f"\t{genome_and_loc[0]}\t{genome_and_loc[1]}")
                    f.write(f"\n{score_dict[pair][2]}\n")

        # Case where they are different
        else :
            for acr1 in motif_loc_dict[pair[0]].keys() :
                for acr2 in motif_loc_dict[pair[1]].keys() :
                    if acr1 != acr2 :
                        # Sort this, we then need to add the starts to the correct one (so get the argsort)
                        sorted_order = np.argsort([acr1, acr2])
                        sorted_pair = sorted([acr1, acr2])
                        out_file = Path(out_dir) / (sorted_pair[0] + '_' + sorted_pair[1] + '.tsv')

                        
                        # Worry about offsets
                        with open(out_file, "a") as f:
                            # Write all genome loc pairs 
                            for i, genome_and_loc in enumerate(motif_loc_dict[pair[sorted_order[0]]][sorted_pair[0]]) :
                                # This guarantees the tabs are correct
                                if i == 0 :
                                    f.write(f"{genome_and_loc[0]}\t{genome_and_loc[1] + score_dict[pair][sorted_order[0]]}")

                                else :
                                    f.write(f"\t{genome_and_loc[0]}\t{genome_and_loc[1] + score_dict[pair][sorted_order[0]]}")

                            # separator
                            f.write("##")
                            for i, genome_and_loc in enumerate(motif_loc_dict[pair[sorted_order[1]]][sorted_pair[1]]) :

                                if i == 0 :
                                    f.write(f"{genome_and_loc[0]}\t{genome_and_loc[1] + score_dict[pair][sorted_order[1]]}")
                                else :
                                    f.write(f"\t{genome_and_loc[0]}\t{genome_and_loc[1] + score_dict[pair][sorted_order[1]]}")
                            f.write(f"\n{score_dict[pair][2]}\n")
                    
        # If pair different
        # We need to add start points
        # else :
        #     for location_1 in motif_loc_dict[pair[0]] :
        #         for location_2 in motif_loc_dict[pair[1]] :
        #             if location_1[0] != location_2[0] :
        #                 # Sort this, we then need to add the starts to the correct one (so get the argsort)
        #                 sorted_order = np.argsort([location_1[0], location_2[0]])
        #                 sorted_pair = sorted([location_1, location_2])
        #                 out_file = Path(out_dir) / (sorted_pair[0][0] + '_' + sorted_pair[1][0] + '.tsv')
                        
        #                 with open(out_file, "a") as f:
        #                     f.write(f"{sorted_pair[0][1]}\t{sorted_pair[0][2] + score_dict[pair][sorted_order[0]]}\t{sorted_pair[1][1]}\t{sorted_pair[1][2]  + score_dict[pair][sorted_order[1]]}\t{score_dict[pair][2]}\n")
        

DATA_DIR = '../Preprocessed_Data'
OUT_DIR = '../Anchors_min1'
MIN_SCORE = 1.0



# score_dict = {("ACCG", "ACCG") : (0, 0, 8), 
#               ("GGAC", "GGAC") : (0, 0, 8), 
#               ("AAAA", "AAAA") : (0, 0, 8), 
#               ("AAAA", "ACCG") : (3, 0, 2),
#                 ("AAAAACC", "TAACC") : (3, 1, 5) 
#               }


score_dict = read_motif_pairs(MIN_SCORE)
find_anchors_all_pairs(DATA_DIR, OUT_DIR, score_dict)
