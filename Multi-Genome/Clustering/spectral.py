from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

def spec_cluster_get_sil(dist_matrix, k):
    #Subtract every element from 1 to make large numbers indicate more similar ACRs
    dist_matrix_inv = [[1 - x for x in list] for list in dist_matrix]
    clustering = SpectralClustering(n_clusters = k, random_state = 0, affinity="precomputed")
    labels = clustering.fit_predict(dist_matrix_inv)
    sil_score = silhouette_score(dist_matrix, labels, metric="precomputed")
    return (sil_score, labels)
