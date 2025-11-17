import numpy as np
from sklearn.cluster import KMeans
from collections import defaultdict
from sklearn.metrics import silhouette_score
from tqdm import tqdm


#####################################################################

# Number of clusters. If it's a list, the data points will be clustered len(num_clusters) times,
# once for each distance threshold
num_clusters = [50, 70, 90, 110, 130, 150]
distance_file = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90_10/ref_set_dist_glob_weighted_alpha0.tsv'
OUT_DIR = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/cluster_exp/'


# Only output clusters of a certain size
CLUSTER_MIN_SIZE = 4

#If True, the evaluation will be output instead of the cluster file
evaluate = True

#######################################################################

for num_cluster in num_clusters:
    if evaluate:
        OUTFILE = f"{OUT_DIR}kmeans_clusters_{num_cluster}.tsv"
        with open(OUTFILE, 'w') as out:
            out.write(f"Clustering with Kmeans with {num_cluster} Clusters")

    print(f"Clustering with Kmeans with {num_cluster} Clusters")

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
    clustering = KMeans(n_clusters = num_cluster, n_init="auto")

    labels = clustering.fit_predict(distance_matrix)

    cluster_dict = defaultdict(set)
    for i, lab in enumerate(labels) :
        cluster_dict[lab].add(i)

    if evaluate:
        try:
            print(f"Silhouette Score: {silhouette_score(distance_matrix, labels)}", flush=True)
        except Exception as e:
            print(f"Exception: {e}", flush=True)
        
        count = 0
        for cluster in cluster_dict.values():
            if len(cluster) > CLUSTER_MIN_SIZE:
                count += 1
        print(f"Number of significant clusters: {count}\n", flush=True)

    else:
        print("Outputting...")
        with open(OUTFILE, 'w') as out:

            for cluster_no, cluster in tqdm(cluster_dict.items()) :
                if len(cluster) > CLUSTER_MIN_SIZE :

                    if len(cluster) > 1:
                        # Get the representative
                        representative = -1
                        min_distances = float('inf')
                        for c in cluster :
                            curr_sum = 0
                            # calculate average distances
                            for d in cluster :
                                if c != d :
                                    curr_sum += distance_matrix[c][d] 

                            if curr_sum < min_distances :
                                representative = c
                                min_distances = curr_sum
                    # if the cluster size is just 1
                    elif len(cluster) == 1 :
                        representative = cluster[0]
                    else :
                        representative = -1
                        label_dict[-1] = "None"

                    out.write(f"Cluster no: {cluster_no}, Representative: {label_dict[representative]}\n")
                    for i, c in enumerate(cluster) :
                        if i != 0 :
                            out.write("\t")
                        out.write(f"{label_dict[c]}")
                    out.write("\n")