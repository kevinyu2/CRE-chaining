import numpy as np
import pandas as pd
import random
random.seed(10)

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
import seaborn as sns #type:ignore
import matplotlib.pyplot as plt

def gen_matrix(distance_file):
    # -----------------------------
    # Load distance matrix
    # -----------------------------
    
    # Maps each label to an index
    index_dict = {}
    curr_idx = 0

    with open(distance_file, 'r') as distances :
        for line in distances:
            line_arr = line.split('\t')
            if line_arr[0] not in index_dict :
                index_dict[line_arr[0]] = curr_idx
                curr_idx += 1
            if line_arr[1] not in index_dict :
                index_dict[line_arr[1]] = curr_idx
                curr_idx += 1
    # Maps indexes to labels
    label_dict = {v: k for k, v in index_dict.items()}

    labels = index_dict.keys()

    dist_matrix = np.full((curr_idx, curr_idx), 1.0)
    np.fill_diagonal(dist_matrix, 0)

    with open(distance_file, 'r') as distances :
        for line in distances:
            line_arr = line.split('\t')
            if float(line_arr[2]) < 0:
                line_arr[2] = 0
            dist_matrix[index_dict[line_arr[0]]][index_dict[line_arr[1]]] = line_arr[2]
            dist_matrix[index_dict[line_arr[1]]][index_dict[line_arr[0]]] = line_arr[2]
            
    return dist_matrix

"""
    Parameters
    ----------
    distancespath : str
        Path to the distances
    sample_size : int
        Number of points to randomly sample.
    method : str
        Linkage method (e.g., 'average', 'single', 'complete', 'ward').
"""
def main(dist_matrix, sample_size=1000, method="average"):

    n = dist_matrix.shape[0]
    assert dist_matrix.shape[0] == dist_matrix.shape[1], "Matrix must be square"
    # -----------------------------
    # Randomly sample indices
    # -----------------------------
    if sample_size == None:
        sampled_dist_matrix = dist_matrix
    else:
        sampled_indices = random.sample(range(n), sample_size)

        # Subset the distance matrix
        sampled_dist_matrix = dist_matrix[np.ix_(sampled_indices, sampled_indices)]

    # -----------------------------
    # Convert to condensed form
    # -----------------------------
    # scipy expects the upper triangle (excluding diagonal)
    condensed_dist = squareform(sampled_dist_matrix, checks=False)

    # -----------------------------
    # Compute linkage
    # -----------------------------
    Z = linkage(condensed_dist, method=method)

    # -----------------------------
    # Plot clustermap
    # -----------------------------
    g = sns.clustermap(
        sampled_dist_matrix,
        row_linkage=Z,
        col_linkage=Z,
        cmap="viridis",
        xticklabels=False,
        yticklabels=False
    )

    g.ax_row_dendrogram.set_visible(False) # Hide row dendrogram
    g.ax_col_dendrogram.set_visible(False) # Hide column dendrogram

    if sample_size == None:
        plt.savefig(f"clustermap_all.png")
    else:
        plt.savefig(f"clustermap_{sample_size}.png")


if __name__ == "__main__":
    distance_file = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90-10/ref_set_dist_glob_weighted_alpha0.tsv'
    dist_matrix = gen_matrix(distance_file)
    main(dist_matrix, sample_size=None)
