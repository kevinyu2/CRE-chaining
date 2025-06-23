import numpy as np
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
from sklearn.metrics import silhouette_score
import umap
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

distance_threshold = 0.2
distance_file = '/home/mwarr/Data/distances_alpha0_one.tsv'

print(f"Using Distance Threshold {distance_threshold}")
print()

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

distance_matrix = np.full((curr_idx, curr_idx), 1.0)
np.fill_diagonal(distance_matrix, 0)

with open(distance_file, 'r') as distances :
    for line in distances:
        line_arr = line.split('\t')
        if float(line_arr[2]) < 0:
            line_arr[2] = 0
        distance_matrix[index_dict[line_arr[0]]][index_dict[line_arr[1]]] = line_arr[2]
        distance_matrix[index_dict[line_arr[1]]][index_dict[line_arr[0]]] = line_arr[2]

# Cluster
clustering = AgglomerativeClustering(
    metric='precomputed', 
    linkage='average',
    distance_threshold=distance_threshold,
    n_clusters=None
)

labels = clustering.fit_predict(distance_matrix)

# for i,lab in enumerate(labels) :
#     print(f"{label_dict[i]}\t{lab}\n", end = '')

# {cluster_num: (idx1,idx2,...)}
cluster_dict = defaultdict(set)
for i, lab in enumerate(labels) :
    cluster_dict[lab].add(i)

for cluster_no, cluster in cluster_dict.items() :
    if len(cluster) > 5 :
        intra_distances = []
        inter_distances = []
        for acr1 in cluster :
            for acr2 in cluster :
                if acr1 != acr2 :
                    intra_distances.append(distance_matrix[acr1][acr2])

        for acr1 in range(curr_idx) :
            for acr2 in range(curr_idx) :
                if acr1 != acr2 and acr1 not in cluster and acr2 in cluster :
                    inter_distances.append(distance_matrix[acr1][acr2])
        
        print(f"Cluster no: {cluster_no}")
        for i, c in enumerate(cluster) :
            if i != 0 :
                print("\t", end = '')
            print(f"{label_dict[c]}", end = '')
        print()
        print(f"Intra Cluster Average Distance: {np.mean(intra_distances)}\nInter Cluster Average Distance: {np.mean(inter_distances)}\n")

if max(labels) != 0:
    score = silhouette_score(distance_matrix, labels, metric='precomputed')
    print(f"Silhouette score: {score:.3f}")
else:
    score = -2   

