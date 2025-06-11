from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import math
import pandas as pd

'''
Makes graph for anchor count vs chain length
'''

# input_tsv should be the mean file, outfile is the graph you want
def plot_anchor_chain(input_tsv, outfile):

    df = pd.read_csv(input_tsv, sep = '\t')
    

    x_max = math.ceil(max(df['Mean_anchors']))
    x_range = range(1, x_max + 1)

    y_exp = [get_exp(x) for x in x_range]
    y_1stdev = [get_stdev_threshold(1, x) for x in x_range]
    y_2stdev = [get_stdev_threshold(2, x) for x in x_range]

    plt.figure()
    plt.scatter(df['Mean_anchors'], df['Mean_score'], s = 0.1, alpha = 0.1, c = 'blue')
    plt.xscale('log')
    plt.plot(x_range, y_exp, label = 'Expected Chain Length', c = 'red', alpha = 0.5)
    # plt.plot(x_range, y_1stdev, label = '1 Stdev', c = 'red', alpha = 0.5)
    # plt.plot(x_range, y_2stdev, label = '2 Stdev', c = 'red', alpha = 0.5)
    plt.legend()
    plt.xlabel('Number of Anchors')
    plt.ylabel('Avg Chain Length')
    plt.title("Average Chain Length vs Anchor Count for All Inter Gene Pairs")
    plt.savefig(outfile)

def get_exp(anchors) :
    return (np.sqrt(anchors) * 2 - 1.771088 * (anchors ** (1/6)) + 0.6)



# Gets the number you need to be above
def get_stdev_threshold(num_stdevs, anchors) :
    base = np.sqrt(anchors) * 2 - 1.771088 * (anchors ** (1/6)) + 0.6
    return base + (num_stdevs * (0.5222718 * (anchors ** (0.2470597))))
 
plot_anchor_chain('/home/mwarr/average_min1.tsv', '/home/kyu/anchor_count_to_chain_len.png')