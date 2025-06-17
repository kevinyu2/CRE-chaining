from hierarchical import *
from dbscan import *
from kmeans import cluster_intra_eval
from pathlib import Path
import time
from matplotlib import pyplot as plt
import numpy as np
import sys

#Note, a sil score of -2 indicates that the either all were put into one cluster, or all were put into
#their own cluster every time we tried to cluster


def count_unique_values(list):
    values = set()
    count = 0
    for item in list:
        if item not in values:
            values.add(item)
            count += 1
    return count

distance_dir = Path("/home/mwarr/Data/Clustering/Distances_min1_intra_alpha0")
GENOME_DIR = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes"
best_hier = []
hier_clusters = []
best_db = []
db_clusters = []
best_kmeans = []
kmeans_clusters = []

start_time = time.time()
count = 0

for file in distance_dir.glob("*"):
    hier_scores = []
    dbscan_scores = []
    kmeans_scores = []
    for i in range(1, 410, 10):
        threshold = i * .002
        #hierarchial clustering
        try:
            cluster_output = hier_cluster(distance_file, threshold)
            sil_score = cluster_output[2]
            num_clusters = count_unique_values(cluster_output[0])
            hier_scores.append((sil_score, num_clusters))
        except ValueError:
            hier_scores.append((-2, 0))
            print("Value Error", file = sys.stderr)

        #db scan
        # try:
        #     cluster_output = db_cluster(distance_file, threshold)
        #     sil_score = cluster_output[2]
        #     num_clusters = count_unique_values(cluster_output[0])
        #     dbscan_scores.append((sil_score, num_clusters))
        # except ValueError:
        #     print("Value Error", file = sys.stderr)


    #kmeans clustering
    for i in range(2, 20):
        try:
            cluster_output = cluster_intra_eval(i, file, GENOME_DIR)
            sil_score = cluster_output[0]
            if sil_score == -2:
                kmeans_scores.append(sil_score)
                break
            else:
                num_clusters = count_unique_values(cluster_output[1])
                kmeans_scores.append((sil_score, num_clusters))
        except ValueError:
            kmeans_scores.append((-2,0))
            print("Value Error", file = sys.stderr)
    
    max_hier = max(hier_scores)
    # max_db = max(dbscan_scores)
    max_kmeans = max(kmeans_scores)
    
    best_hier.append(max_hier[0])
    hier_clusters.append(max_hier[1])
    # best_db.append(max_db[0])
    # db_clusters.append(max_db[1])
    best_kmeans.append(max_kmeans[0])
    kmeans_clusters.append(max_kmeans[1])

    count += 1

    print(f"Finished {count} files in {time.time() - start_time} seconds")

plt.clf()
plt.hist(best_hier)
plt.savefig("/home/mwarr/hier_sil_min1_intra_alpha0")
plt.xlabel("Silhouette Score")
plt.ylabel("Frequency")
plt.yscale("log")
plt.clf()

plt.scatter(hier_clusters, best_hier)
plt.savefig("/home/mwarr/hier_clus_min1_intra_alpha0")
plt.xlabel("Number of clusters")
plt.ylabel("Silhouette Score")
plt.clf()

# plt.hist(best_db)
# plt.savefig("/home/mwarr/db_sil_min1_intra_alpha0")
# plt.clf()

# plt.scatter(db_clusters, best_db)
# plt.savefig("/home/mwarr/db_clus_min1_intra_alpha0")
# plt.clf()

plt.hist(best_kmeans)
plt.savefig("/home/mwarr/kmeans_sil_min1_intra_alpha0")
plt.xlabel("Silhouette Score")
plt.ylabel("Frequency")
plt.yscale("log")
plt.clf()

plt.scatter(kmeans_clusters, best_kmeans)
plt.savefig("/home/mwarr/kmeans_clus_min1_intra_alpha0")
plt.xlabel("Number of clusters")
plt.ylabel("Silhouette Score")
plt.clf()


