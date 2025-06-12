from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import math
import pandas as pd

'''
Makes a file that has the averages across all inter gene comparisons
'''

# Gets the expected value and standard deviation for the anchor number
# Returns a n by 2 array, where arr[i][0] means for i anchors, the expected value and arr[i][1] means for i anchors, the stdev
def get_exp_stdev(n) :
    table = np.zeros((n, 2))
    assert (n > 1)

    # Don't fill these with formula
    table[0][0] = 0
    table[0][1] = 1
    table[1][0] = 1
    table[1][1] = 1

    for i in range(2, n) :
        table[i][0] = np.sqrt(i) * 2 - 1.771088 * (i ** (1/6)) + 0.6
        table[i][1] = 0.5222718 * (i ** (0.2470597))
    
    return table

# Gets a dict of genome/acr to number of motifs 
def get_motif_count_dict() :
    motif_counts = {}
    df = pd.read_csv('/home/mwarr/motif_counts.tsv', sep = '\t')
    for i, row in df.iterrows() :
        motif_counts[row['genome']] = int(row['tomtom_motif_count'])

    return motif_counts


def get_means(input_dir, outfile):

    expected_stdev_table = get_exp_stdev(100000)
    motif_counts = get_motif_count_dict()

    with open(outfile, 'w') as f:

        f.write("File_name\tMean_ratio\tMean_anchors\tMean_score\tMean_zscore\tMean_small_percent\tMean_large_percent\nFile_len\n")
        input_dir = Path(input_dir)
        count = 0
        for file_name in input_dir.glob("*"):
            count += 1
            if count % 1000 == 0 :
                print(f"Completed: {count}")
            with open(file_name) as acr_file:
                ratio_list = []
                score_list = []
                anchors_list = []
                stdevs_list = []
                small_percent = []
                large_percent = []
                line_count = 0
                for line in acr_file:
                    line_count += 1
                    line_arr = line.split("\t")

                    if len(line_arr) > 1 :
                        chain_len = int(line_arr[2])
                        anchor_num = int(line_arr[3].strip())
                        ratio_list.append(chain_len/anchor_num)
                        score_list.append(chain_len)
                        anchors_list.append(anchor_num)
                        # (actual - expected) / stdev
                        stdevs_list.append(float(chain_len - expected_stdev_table[anchor_num][0]) / 
                                           float(expected_stdev_table[anchor_num][1]))
                        pair_motif_counts = sorted([motif_counts[line_arr[0]], motif_counts[line_arr[1]]])
                        small_percent.append(chain_len/pair_motif_counts[0])
                        large_percent.append(chain_len/pair_motif_counts[1])                        


                f.write(f"{file_name}\t{np.mean(ratio_list)}\t{np.mean(anchors_list)}\t{np.mean(score_list)}\t{np.mean(stdevs_list)}\t{np.mean(small_percent)}\t{np.mean(large_percent)}\t{line_count}\n")


get_means('/home/mwarr/Chaining_min1', '/home/mwarr/average_min1.tsv')
