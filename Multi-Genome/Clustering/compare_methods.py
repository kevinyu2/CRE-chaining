from hierarchical import *
from dbscan import *
from kmeans import *
from spectral import *
from pathlib import Path
import time
from matplotlib import pyplot as plt
import numpy as np
import sys

#Note: a sil score of -2 indicates that the either all were put into one cluster, or all were put into
#their own cluster every time we tried to cluster (essentially, no clustering)


def count_unique_values(list):
    values = set()
    count = 0
    for item in list:
        if item not in values:
            values.add(item)
            count += 1
    return count

def output_to_file(output_dir, method, best_scores):
    with open(f"{output_dir}/{method}_sil_clus_min1_intra_alpha0.tsv", "w") as file:
            file.write("silhouette_score\tnum_clusters")
            for pair in best_scores:
                #pair is (sil_score, #clusters)
                file.write(f"\n{pair[0]}\t{pair[1]}")


def output_best_clusters(distance_dir, genome_dir, output_dir):
    #List of tuples -- one for each file
    #[(sil_score, #clusters), (...), ...]
    best_hier = []
    best_kmeans = []
    best_db = []
    best_spec = []

    start_time = time.time()
    count = 0

    for file in distance_dir.glob("*"):
        #lists of tuples with the sil scores and number of clusters
        hier_scores = []
        db_scores = []
        kmeans_scores = []
        spec_scores = []
        ind_dict = create_index_dict(genome_dir, file)
        distance_matrix = get_dist_matrix(file, ind_dict)

        for i in range(1, 410, 10):
            threshold = i * .002
            print(f"threshold {threshold}")
            #hierarchial clustering
            try:
                cluster_output = hier_cluster_get_sil(distance_matrix, threshold)
                sil_score = cluster_output[0]
                num_clusters = count_unique_values(cluster_output[1])
                hier_scores.append((sil_score, num_clusters))
                print(f"Hier sil score: {sil_score}")
            except ValueError: 
            #raised by silhouette score; means either all were in same cluster or
            #all were in their own cluster
                hier_scores.append((-2, 0))
                print("Value Error")

            #dbscan clustering
            try:
                cluster_output = db_cluster_get_sil(distance_matrix, threshold)
                sil_score = cluster_output[0]
                num_clusters = count_unique_values(cluster_output[1])
                db_scores.append((sil_score, num_clusters))
                print(f"dbscan sil score: {sil_score}")
            except ValueError: 
            #raised by silhouette score; means either all were in same cluster or
            #all were in their own cluster
                db_scores.append((-2, 0))
                print("Value Error")


        
        for i in range(2, 20):
            #kmeans clustering
            try:
                cluster_output = kmeans_cluster_get_sil(distance_matrix, i)
                sil_score = cluster_output[0]
                num_clusters = count_unique_values(cluster_output[1])
                kmeans_scores.append((sil_score, num_clusters))
            except ValueError:
                #raised by silhouette score; means all were in one cluster
                kmeans_scores.append((-2,0))
                print("Value Error", file = sys.stderr)
            
            #spectral clustering
            try:
                cluster_output = spec_cluster_get_sil(distance_matrix, i)
                sil_score = cluster_output[0]
                num_clusters = count_unique_values(cluster_output[1])
                spec_scores.append((sil_score, num_clusters))
            except ValueError:
                #raised by silhouette score; means all were in one cluster
                spec_scores.append((-2,0))
                print("Value Error", file = sys.stderr)
        
        best_hier.append(max(hier_scores))
        best_kmeans.append(max(kmeans_scores))
        best_db.append(max(db_scores))
        best_spec.append(max(spec_scores))

        count += 1

        print(f"Finished {count} files in {time.time() - start_time} seconds")

    output_to_file(output_dir, "kmeans", best_kmeans)
    output_to_file(output_dir, "hier", best_hier)
    output_to_file(output_dir, "db", best_db)
    output_to_file(output_dir, "spec", best_spec)
    
    return (best_hier, best_kmeans, best_db, best_spec)


def create_graphs(best_hier, best_kmeans, best_db, best_spec):
    best_hier, hier_clusters = zip(*best_hier)
    best_db, db_clusters = zip(*best_db)
    best_kmeans, kmeans_clusters = zip(*best_kmeans)
    best_spec, spec_clusters = zip(*best_spec)

    plt.clf()
    plt.hist(best_hier, bins=100)
    plt.xlabel("Silhouette Score")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.savefig("/home/mwarr/hier_sil_min1_intra_alpha0")
    plt.clf()

    plt.scatter(hier_clusters, best_hier, s=.1)
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette Score")
    plt.savefig("/home/mwarr/hier_clus_min1_intra_alpha0")
    plt.clf()

    plt.hist(best_kmeans, bins=100)
    plt.xlabel("Silhouette Score")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.savefig("/home/mwarr/kmeans_sil_min1_intra_alpha0")
    plt.clf()

    plt.scatter(kmeans_clusters, best_kmeans, s=.1)
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette Score")
    plt.savefig("/home/mwarr/kmeans_clus_min1_intra_alpha0")
    plt.clf()

    plt.hist(best_db, bins=100)
    plt.xlabel("Silhouette Score")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.savefig("/home/mwarr/db_sil_min1_intra_alpha0")
    plt.clf()

    plt.scatter(db_clusters, best_db, s=.1)
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette Score")
    plt.savefig("/home/mwarr/db_clus_min1_intra_alpha0")
    plt.clf()

    plt.hist(best_spec, bins=100)
    plt.xlabel("Silhouette Score")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.savefig("/home/mwarr/spec_sil_min1_intra_alpha0")
    plt.clf()

    plt.scatter(spec_clusters, best_spec, s=.1)
    plt.xlabel("Number of clusters")
    plt.ylabel("Silhouette Score")
    plt.savefig("/home/mwarr/spec_clus_min1_intra_alpha0")
    plt.clf()

def get_stats(eval_file):
    no_score = 0
    score_sum = 0
    cluster_sum = 0
    count = 0
    with open(eval_file, "r") as file:
        file.readline() #read header
        for line in file:
            line_arr = line.split("\t")
            sil_score = float(line_arr[0])
            clusters = int(line_arr[1])
            if sil_score == -2:
                no_score += 1
            else:
                score_sum += sil_score
                cluster_sum += clusters
                count += 1
    print(f"Not clustered: {no_score}, Average sil_score: {round(score_sum / count, 2)}, Average clusters: {round(cluster_sum / float(count), 2)}")


distance_dir = Path("/home/mwarr/Data/Clustering/Distances_min1_intra_alpha0")
genome_dir = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes"
output_dir = "/home/mwarr/Data/Clustering/Evaluate"

lists = output_best_clusters(distance_dir, genome_dir, output_dir)
create_graphs(lists[0], lists[1], lists[2], lists[3])
