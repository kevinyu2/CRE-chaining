from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from pathlib import Path
import time
import numpy as np
from matplotlib import pyplot as plt

def create_index_dict(genome_dir, file_name):
    genome_ind = {}
    #label each ACR in the genomes with an index
    with open(f"{genome_dir}/{file_name.stem}.group.txt") as file:
        count = 0
        for line in file:
            genome_ind[line.strip()] = count
            count += 1
    return genome_ind

def get_score_matrix(file_name, genome_ind):
    size = len(genome_ind.keys())
    scores = np.zeros((size, size))
    with open(file_name) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            genome1 = line_arr[0]
            genome2 = line_arr[1]
            score = float(line_arr[2])
            genome1_ind = genome_ind[genome1]
            genome2_ind = genome_ind[genome2]
            scores[genome1_ind][genome2_ind] = score
            scores[genome2_ind][genome1_ind] = score
    return scores

def cluster_intra_eval(k, input_dir, genome_dir):
    search_dir = Path(input_dir)
    start_time = time.time()
    count = 0
    for file_name in search_dir.glob("*"):
        genome_ind_dict = create_index_dict(genome_dir, file_name)]
            scores = get_score_matrix(file_name, genome_ind_dict)
        if len(scores[0]) < 20:
            continue
        kmeans = KMeans(n_clusters = i, n_init='auto', random_state = 0)
        labels = kmeans.fit_predict(scores)
        sil_score = silhouette_score(scores, labels)
        return sil_score


                

# #cluster_intra_evaluate(2, "../test_folder", "../test_output_folder")
# k = []
# sil_scores = []
# for i in range(2, 20):
cluster_intra_eval(2, "/home/mwarr/Data/Clustering/Distances_min1_intra_alpha50", "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes" )
#     k.append(i)
#     sil_scores.append(sil_score)
#     plt.plot(k, sil_scores)
#     plt.savefig("../Visuals/silhouette.png")