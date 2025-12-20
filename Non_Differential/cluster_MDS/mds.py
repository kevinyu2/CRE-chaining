import numpy as np
from numpy.linalg import eigh, pinv
from sklearn.manifold import MDS
import time
from matplotlib import pyplot as plt
import seaborn as sns #type: ignore
import csv


def mds(distance_file, out_dir, evaluate, dimensions):
    # Create distance matrix

    start_time = time.time()

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

    print(f"Created distance matrix after {round(time.time() - start_time, 2)} s", flush=True)

    #####################################################################
    # Perform MDS

    if evaluate:
        stress = []
        for num in range(200, 5100, 100):
            start_time = time.time()
            mds = MDS(n_components=num, metric=False, normalized_stress=True, dissimilarity="precomputed")
            mds.fit(distance_matrix)
            print(f"Scikit non-metric MDS with {num} dimensions finished in {time.time() - start_time} s", flush=True)
            print(f"Normalized Stress: {mds.stress_}", flush=True)
            stress.append(mds.stress_)

        plt.plot(range(200, 5100, 100), stress)
        plt.title("Non-Metric MDS Stress")
        plt.xlabel("Dimensions")
        plt.ylabel("Stress")
        plt.savefig("mds_stress.png", dpi=200)
    else:
        #get MDS embedding
        start_time = time.time()
        mds = MDS(n_components=dimensions, metric=False, normalized_stress=True, dissimilarity="precomputed")
        mds.fit(distance_matrix)
        print(f"Scikit non-metric MDS with {dimensions} dimensions finished in {time.time() - start_time} s", flush=True)
        print(f"Normalized Stress: {mds.stress_}", flush=True)
        embedding = mds.embedding_

        #Save embedding
        row_headers = np.array([label_dict[ind] for ind in range(len(embedding))]).reshape(-1, 1)
        embedding = np.hstack([row_headers, embedding])
        np.savetxt(f"{out_dir}/mds_embedding_{dimensions}.csv", embedding, delimiter=",", fmt="%s")
        return mds.embedding_

def read_embedding(in_file):
    data = []
    with open(in_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            # Skip row header (first column)
            data.append([float(x) for x in row[1:]])

    return np.array(data)

def create_heatmap(embedding, dimensions):
        #Use embedding to create heatmap
        sns.clustermap(embedding)
        plt.savefig(f"cluster_map_{dimensions}.png")


#####################################################################
#Drivers and Inputs

distance_file = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90-10/ref_set_dist_glob_weighted_alpha0.tsv'
embedding_file = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/cluster_mds_random_50-50/mds_embedding_1500.csv"
out_dir = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/cluster_mds_random_50-50"

def evaluate_mds(dimensions):
    mds(distance_file, out_dir, evaluate=True, dimensions=dimensions)

def read_mds_create_heatmap(dimensions):
    embedding = read_embedding(embedding_file)
    print("read embedding", flush=True)
    create_heatmap(embedding, dimensions)

def perform_mds_create_heatmap(dimensions):
    embedding = mds(distance_file, out_dir, evaluate=False, dimensions=dimensions)
    create_heatmap(embedding, dimensions)


if __name__ == "__main__":
    for dimensions in [50, 100, 200]:
        perform_mds_create_heatmap(dimensions)
#####################################################################