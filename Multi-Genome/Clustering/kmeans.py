from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from pathlib import Path
import time
import numpy as np
from matplotlib import pyplot as plt
import warnings
import umap
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

#returns a dictionary where the keys are the genomes and the value is an index
#associated with that genome
def create_index_dict(genome_dir, file_path):
    genome_ind = {}
    #label each ACR in the genomes with an index
    with open(f"{genome_dir}/{file_path.stem}.group.txt") as file:
        count = 0
        for line in file:
            genome_ind[line.strip()] = count
            count += 1
    return genome_ind

#Creates a matrix of the distance scores between each genome pair
def get_score_matrix(file_path, genome_ind):
    size = len(genome_ind.keys())
    scores = np.zeros((size, size))
    with open(file_path) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            genome1 = line_arr[0]
            genome2 = line_arr[1]
            score = round(float(line_arr[2]), 2)
            if score < 0:
                score = 0
            genome1_ind = genome_ind[genome1]
            genome2_ind = genome_ind[genome2]
            scores[genome1_ind][genome2_ind] = score
            scores[genome2_ind][genome1_ind] = score
    return scores

#Gives the silhouette score (-1 to 1, with 1 being the best) of the clustering within a  
#single ACR, using distances as the features
def cluster_intra_eval(k, file_path, genome_dir):
    genome_ind_dict = create_index_dict(genome_dir, file_path)
    scores = get_score_matrix(file_path, genome_ind_dict)
    #Don't cluster if there aren't enough unique scores
    kmeans = KMeans(n_clusters = k, n_init='auto', random_state = 0)
    try:
        labels = kmeans.fit_predict(scores)
    except:
        return (-2, 0)
    sil_score = silhouette_score(scores, labels, metric="precomputed")
    return (sil_score, labels)

#Cluster the distances within a single ACR with distances as features
def cluster_intra_output(k, file_path, genome_dir, output_dir, silhouette):
    genome_ind_dict = create_index_dict(genome_dir, file_path)
    scores = get_score_matrix(file_path, genome_ind_dict)
    kmeans = KMeans(n_clusters = k, n_init='auto', random_state = 0)
    labels = kmeans.fit_predict(scores)
    with open(f"{output_dir}/{file_path.stem}_k{k}_sil{round(silhouette, 2)}.tsv", "w") as out_file:
        for genome, index in genome_ind_dict.items():
            label = labels[index]
            out_file.write(f"{genome}\t{label}\n")

def visualize_clusters(n_clusters, distance_matrix, output_dir, labels):
    reducer = umap.UMAP(metric='precomputed', random_state=42)
    coords = reducer.fit_transform(distance_matrix)

    cmap = plt.get_cmap('tab10', n_clusters)  # categorical colormap with n colors
    bounds = np.arange(-0.5, n_clusters, 1)   # boundaries between colors
    norm = BoundaryNorm(bounds, cmap.N)

    scatter = plt.scatter(coords[:, 0], coords[:, 1], c=labels, cmap=cmap, norm=norm)

    cbar = plt.colorbar(scatter, ticks=np.arange(n_clusters))
    cbar.set_label('Cluster Number')
    cbar.ax.set_yticklabels(np.arange(n_clusters))

    plt.savefig(output_dir)

#Tries different k values, finds the best k value based on max silhouette score, 
#and then outputs the cluster labels to a file 
def cluster_all_driver(input_dir, genome_dir, output_dir):
    search_dir = Path(input_dir)
    start_time = time.time()
    count = 0
    for file_path in search_dir.glob("*"):
        max_score = -2
        k_max = 0
        for k in range(2, 17):
            sil_score = cluster_intra_eval(k, file_path, genome_dir)
            #this means that a warning was caught -- k is greater than the number of unique
            #clusters. This means we want to stop increasing k
            if sil_score == -2:
                break
            #update if the silhouette score is greater than current max
            if sil_score > max_score:
                max_score = sil_score
                k_max = k
        if max_score != -2:
            cluster_intra_output(k_max, file_path, genome_dir, output_dir, max_score)
        else:
            #this means there is really only one cluster (all have very similar distances)
            print(f"{file_path} could not be clustered")
        count += 1
        print(f"Finished {count} in {time.time() - start_time} seconds")
    



INPUT_DIR = "/home/mwarr/Data/Clustering/Distances_min1_intra_alpha50/Chr2_5345616to5346848.tsv"
OUTPUT_DIR = "/home/mwarr/Data/Clustering/Labels_kmeans_min1_intra_alpha50"
GENOME_DIR = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes"


genome_ind = create_index_dict(GENOME_DIR, Path(INPUT_DIR))
dist_matrix = get_score_matrix(Path(INPUT_DIR), genome_ind)
kmeans = KMeans(n_clusters = 16, n_init='auto', random_state = 0)
labels = kmeans.fit_predict(dist_matrix)
visualize_clusters(16, dist_matrix, "/home/mwarr/clusters.png", labels)

#cluster_all_driver(INPUT_DIR, GENOME_DIR, OUTPUT_DIR)

# sil_scores = []
# for i in range(2, 20):
#     sil_score = cluster_intra_eval(i, Path(INPUT_DIR), GENOME_DIR)
#     if sil_score == -2:
#         break
#     sil_scores.append(sil_score)

# plt.plot(sil_scores)
# plt.savefig("/home/mwarr/silhouette.png")



