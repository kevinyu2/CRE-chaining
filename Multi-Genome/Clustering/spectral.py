from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score

def spec_cluster_get_sil(dist_matrix, k):
    clustering = SpectralClustering(n_clusters = k, n_init='auto', random_state = 0)
    labels = clustering.fit_predict(dist_matrix)
    sil_score = silhouette_score(dist_matrix, labels, metric="precomputed")
    return (sil_score, labels)