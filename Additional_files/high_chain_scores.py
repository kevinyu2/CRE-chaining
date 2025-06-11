from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np

# Prints files with high score and ratio based on provided thresholds
def get_high_scores(input_dir, outfile, thresh_score, thresh_ratio):

    with open(outfile, 'w') as f:

        f.write("File_name\tMean_ratio\tMean_score\tFile_len")
        input_dir = Path(input_dir)
        count = 0
        for file_name in input_dir.glob("*"):
            count += 1
            if count % 1000 == 0 :
                print(f"Completed: {count}")
            with open(file_name) as acr_file:
                ratio_list = []
                score_list = []
                line_count = 0
                for line in acr_file:
                    line_count += 1
                    line_arr = line.split("\t")

                    if len(line_arr) > 1 :

                        ratio_list.append(float(line_arr[2])/float(line_arr[3]))
                        score_list.append(float(line_arr[2]))
                
                # print()
                # print(ratio_list)
                # print(score_list)
                if np.mean(ratio_list) > thresh_ratio and np.mean(score_list) > thresh_score :
                    f.write(f"{file_name}\t{np.mean(ratio_list)}\t{np.mean(score_list)}\t{line_count}\n")

# Gets the number you need to be above
def get_stdev_threshold(num_stdevs, anchors) :
    base = np.sqrt(anchors) * 2 - 1.771088 * (anchors ** (1/6))
    return base + (num_stdevs * (0.5222718 * (anchors ** (0.2470597))))
 

# Prints files based on number of standard deviations from expected chain length
def get_high_scores_statistically(input_dir, outfile, thresh_score, thresh_stdevs):

    with open(outfile, 'w') as f:

        f.write("File_name\tMean_ratio\tMean_score\tFile_len")
        input_dir = Path(input_dir)
        count = 0
        for file_name in input_dir.glob("*"):
            count += 1
            if count % 1000 == 0 :
                print(f"Completed: {count}")
            with open(file_name) as acr_file:
                anchors_list = []
                score_list = []
                line_count = 0
                for line in acr_file:
                    line_count += 1
                    line_arr = line.split("\t")

                    if len(line_arr) > 1 :

                        anchors_list.append(float(line_arr[3]))
                        score_list.append(float(line_arr[2]))
                
                # print()
                # print(ratio_list)
                # print(score_list)
                if np.mean(score_list) > get_stdev_threshold(thresh_stdevs, np.mean(anchors_list)) and np.mean(score_list) > thresh_score :
                    f.write(f"{file_name}\t{np.mean(score_list)/np.mean(anchors_list)}\t{np.mean(score_list)}\t{line_count}\n")


# get_high_scores("../../Chaining_min1", "./high_chain_scores_min1.txt", 5, 0.5)
get_high_scores_statistically("../../Chaining_min1", "/home/kyu/high_chain_scores_min1_stdev1.txt", 5, 1)
