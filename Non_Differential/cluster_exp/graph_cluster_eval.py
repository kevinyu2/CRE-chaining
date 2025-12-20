from matplotlib import pyplot as plt

num_clusters = range(200, 2200, 200)

filepath = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/cluster_exp/cluster_eval.txt"

def plot(eval_method, hierarchical, kmeans, spectral):
    plt.clf()
    plt.plot(num_clusters, hierarchical, label="Hierarchical")
    plt.plot(num_clusters, kmeans, label="K-Means")
    plt.plot(num_clusters, spectral, label="Spectral")
    plt.legend()
    if eval_method[-1] == ":":
        eval_method = eval_method[:-1]
    plt.title(f"Clustering Evaluation: {eval_method}")
    plt.xlabel("Number of Clusters")
    plt.ylabel("Score")
    plt.savefig(f"{eval_method[:3].lower()}_eval.png", dpi=300)

for eval_method in ["Silhouette", "Davies-Bouldin:", "Calinski-Harabasz:"]:
    list_index = -1
    lists = [[], [], []] #holds three lists: hierarchical scores, kmeans scores, and spectral scores
    with open(filepath) as in_file:
        for line in in_file:
            if line == "\n":
                continue
            line_arr = line.rstrip().split(" ")
            if len(line_arr) == 1: #new clustering method
                list_index += 1
                continue
            if line_arr[0] == eval_method:
                lists[list_index].append(float(line_arr[-1]))
    plot(eval_method, lists[0], lists[1], lists[2])
            
            



