# Count the number of exact CRE kmers, graph (just a quick test to see how viable kmers are)

import matplotlib.pyplot as plt

kmers = []

with open('/home/mwarr/Data/Chaining_one_local.tsv', 'r') as file:
    for line in file :
        line_arr = line.rstrip().split('\t')
        if float(line_arr[2]) / int(line_arr[3]) >= 4.999


plt.hist(kmers, s = 0.1)
plt.savefig('/home/kyu/local_one_scores.png')