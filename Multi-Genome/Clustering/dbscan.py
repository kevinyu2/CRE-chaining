import numpy as np
from sklearn.cluster import DBSCAN
from collections import defaultdict
from sklearn.metrics import silhouette_score
import umap
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


distance_file = '/home/mwarr/Data/Clustering/Distances_min1_mini_alpha50.tsv'
distance_threshold = 0.7 # What the distance must be to make it into a cluster

def db_cluster(distance_file, distance_threshold):
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
    labels = DBSCAN(eps=distance_threshold, min_samples=1, metric='precomputed').fit_predict(distance_matrix)


    # for i,lab in enumerate(labels) :
    #     print(f"{label_dict[i]}\t{lab}\n", end = '')

    # {cluster_num: (idx1,idx2,...)}
    cluster_dict = defaultdict(set)
    for i, lab in enumerate(labels) :
        cluster_dict[lab].add(i)

    for cluster_no, cluster in cluster_dict.items() :
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
        
        # print(f"Cluster no: {cluster_no}")
        # for i, c in enumerate(cluster) :
        #     if i != 0 :
        #         print("\t", end = '')
        #     print(f"{label_dict[c]}", end = '')
        # print()
        # print(f"Intra Cluster Average Distance: {np.mean(intra_distances)}\nInter Cluster Average Distance: {np.mean(inter_distances)}\n")

    if max(labels) != 0:
        score = silhouette_score(distance_matrix, labels, metric='precomputed')
        print(f"Silhouette score: {score:.3f}")  
    else:
        score = -2  

    return (labels, distance_matrix, score)   



   

def visualize(labels, distance_matrix, distance_threshold):
    n_clusters = len(np.unique(labels))

    reducer = umap.UMAP(metric='precomputed', random_state=42)
    coords = reducer.fit_transform(distance_matrix)

    cmap = plt.get_cmap('tab10', n_clusters)  # categorical colormap with n colors
    bounds = np.arange(-0.5, n_clusters, 1)   # boundaries between colors
    norm = BoundaryNorm(bounds, cmap.N)

    scatter = plt.scatter(coords[:, 0], coords[:, 1], c=labels, cmap=cmap, norm=norm, s = 0.1)

    cbar = plt.colorbar(scatter, ticks=np.arange(n_clusters))
    cbar.set_label('Cluster Number')
    cbar.ax.set_yticklabels(np.arange(n_clusters))
    plt.title(f"UMAP Projection of Heirarchical, Distance = {distance_threshold}")

    plt.savefig('/home/mwarr/clusters_dbscan.png')
    
