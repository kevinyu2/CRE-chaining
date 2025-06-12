from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# Prints files with high score and ratio based on provided thresholds
def get_high_scores(input_tsv, outfile, thresh_score, thresh_ratio):

    with open(outfile, 'w') as f:
        f.write("File_name\tMean_ratio\tMean_score\tFile_len\n")

        df = pd.read_csv(input_tsv, sep = '\t')
                
        for i, row in df.iterrows() :
            if row['Mean_ratio'] > thresh_ratio and row["Mean_score"] > thresh_score :
                f.write(f"{row['File_name']}\t{row['Mean_ratio']}\t{np.mean(row["Mean_score"])}\t{row['File_len']}\n")

# Gets the number you need to be above
def get_stdev_threshold(num_stdevs, anchors) :
    base = np.sqrt(anchors) * 2 - 1.771088 * (anchors ** (1/6)) + 0.6
    return base + (num_stdevs * (0.5222718 * (anchors ** (0.2470597))))
 

# Prints files based on number of standard deviations from expected chain length
def get_high_scores_statistically(input_tsv, outfile, thresh_score, thresh_stdevs):

    with open(outfile, 'w') as f:

        f.write("File_name\tMean_ratio\tMean_score\tFile_len\n")
        
        
        df = pd.read_csv(input_tsv, sep = '\t')
                
        for i, row in df.iterrows() :      
            if row['Mean_anchors'] > get_stdev_threshold(thresh_stdevs, row['Mean_anchors']) and row['Mean_score'] > thresh_score :
                f.write(f"{row['File_name']}\t{row['Mean_ratio']}\t{np.mean(row["Mean_score"])}\t{row['File_len']}\n")


get_high_scores_statistically('/home/mwarr/average_min1.tsv', "/home/kyu/high_chain_scores_min1_stdev1.txt", 5, 1)
