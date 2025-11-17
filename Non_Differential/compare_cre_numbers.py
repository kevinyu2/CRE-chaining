import sys
from matplotlib import pyplot as plt

'''
Helper function
<filename> should be formatted:
ACR: <name>
<motif>
...
ACR: <name>
<motif>
...

Returns a dictionary where the keys are ACR names and the values are number of motifs
'''
def motif_num_dict(filename):
    motif_num = {}
    with open(f"/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90_10/ref_set_motifs.txt") as in_file:
        line = in_file.readline()
        line_arr = line.split(" ")
        acr = line_arr[1].rstrip()
        for line in in_file:
            motifs = []
            while line != "" and line != "\n" and line[0] != "A":
                motifs.append(line.rstrip())
                line = in_file.readline()
                motif_num[acr] = len(motifs)
            
            if line != "" and line != "\n":
                line_arr = line.split(" ")
                acr = line_arr[1].rstrip()
    return motif_num

'''
<neg_type> should be 'rand' or 'up'

<motif_file> should be formatted:
ACR: <name>
<motif>
...
ACR: <name>
<motif>
...

Returns a dictionary with lengths as keys and number of ACRs with that length as values.
'''
def get_motifs_important_acrs(neg_type, motif_file):

    motif_num = motif_num_dict(motif_file)
    important_acrs = {}

    with open(f"{neg_type}_important_features.txt", "r") as in_file:
        for line in in_file:
            line_arr = line.split(" ")
            acr = (line_arr[0].rstrip())[:-4]
            num_motifs = motif_num_dict[acr]
            important_acrs[acr] = num_motifs
    
    count = {}
    for key, val in important_acrs.items():
        count[val] = count.get(val, 0) + 1
    return count

'''
Counts motifs in ACRs in <filename>.
<filename> should be a file formatted:
ACR: <name>
<motif>
...
ACR: <name>
<motif>
...

If rep_set=True, then only ACRs without _pos or _neg appended to them will be counted
Returns a dictionary with lengths as keys and number of ACRs with that length as values.
'''  
def get_motifs(filename, rep_set=False):
    acrs = {}
    with open(filename) as in_file:
        line = in_file.readline()
        line_arr = line.split(" ")
        acr = line_arr[1].rstrip()
        for line in in_file:
            motifs = []
            while line != "" and line != "\n" and line[0] != "A":
                motifs.append(line.rstrip())
                line = in_file.readline()
            if rep_set:
                if "pos" not in acr and "neg" not in acr:
                    acrs[acr] = motifs
            else:
                acrs[acr] = motifs
            
            if line != "" and line != "\n":
                line_arr = line.split(" ")
                acr = line_arr[1].rstrip()
    
    count = {}
    for key, val in acrs.items():
        count[len(val)] = count.get(len(val), 0) + 1
    return count

'''
Counts motifs in all the "large enough" clusters (i.e. whichever show up in the file)

<cluster_filename> should be formatted:
Cluster no. ####, Representative: <ACR>
<ACR>\t<ACR>\t...
...

<motif_file> should be formatted:
ACR: <name>
<motif>
...
ACR: <name>
<motif>
...

Returns a dictionary with lengths as keys and number of ACRs with that length as values.
'''
def get_motifs_large_cluster_acrs(cluster_filename, motif_file):
    motif_num = motif_num_dict(motif_file)

    count = {}
    with open(cluster_filename, "r") as cluster_file:
        for line in cluster_file:
            if line == "\n" or line[0:2] == "Cl":
                continue
            line_arr = line.rstrip().split("\t")
            for acr in line_arr:
                count[motif_num[acr]] = count.get(motif_num[acr], 0) + 1
    return count

'''
Counts motifs in representative ACRs.

<cluster_filename> should be formatted:
Cluster no. ####, Representative: <ACR>
<ACR>\t<ACR>\t...
...

<motif_file> should be formatted:
ACR: <name>
<motif>
...
ACR: <name>
<motif>
...

Returns a dictionary with lengths as keys and number of ACRs with that length as values.
'''
def get_motifs_representative(cluster_filename, motif_file):
    motif_num = motif_num_dict(motif_file)

    count = {}
    with open(cluster_filename) as cluster_file:
        for line in cluster_file:
            if line == "\n":
                continue
            if line[0:2] == "Cl":
                acr = line.split("Representative: ")[1].rstrip()
                count[motif_num[acr]] = count.get(motif_num[acr], 0) + 1
            else:
                continue
    return count

'''
<count> should be a dictionary with lengths as keys and number of ACRs with that length as values.
Returns a list of x values (number of CREs) and a list of corresponding y values (fraction of ACRs)
'''
def prepare_data_for_graph(count):
    x = []
    height = []
    total = 0
    for freq in count.values():
        total += freq
    for cre_count, freq in count.items():
        x.append(cre_count)
        height.append(round(freq / total, 2))
    return x, height

def plot(count1, count2, title1, title2, title, filename):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    x, height = prepare_data_for_graph(count1)
    ax1.bar(x, height)
    ax1.set_title(title1)

    x, height = prepare_data_for_graph(count2)
    ax2.bar(x, height)
    ax2.set_title(title2)

    # plt.xlim((0, 50))
    plt.ylim((0, .6))
    plt.suptitle(title)
    fig.supxlabel("Number of CREs")
    fig.supylabel("Fraction of ACRs")
    plt.savefig(f"experiment_1/plots/{filename}", dpi=300)

if __name__ == "__main__":
    count2 = get_motifs_large_cluster_acrs("/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/fix_clustering/clusters_0.55.tsv", "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90_10/all_motifs.txt")
    count1 = get_motifs_representative("/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/fix_clustering/clusters_0.55.tsv", "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/ml_exp_upstream_90_10/all_motifs.txt")
    plot(count1, count2, "Cluster Centers", "ACRs in Non-Trivial Clusters", "ACR CRE Counts (Threshold=.55)", "cluster_.55_acr_cre_count.png")
