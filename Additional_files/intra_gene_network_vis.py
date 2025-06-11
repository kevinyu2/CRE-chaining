import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

'''
Visualizes an intra-gene genome comparison as a graph
'''

ACR_FOLDER = '/home/mwarr/Chaining_min1_intra/'
ACR = 'Chr1_12966791to12968243.tsv'
mode = 'smallest'
min_weight = 5


###############


# Read 1: make a dict of max chain length
max_chain_lengths = {}

with open(ACR_FOLDER + ACR, 'r') as f:
    for line in f:
        line_arr = line.rstrip().split('\t')
        if line_arr[0] == line_arr[1] :
            max_chain_lengths[line_arr[0]] = int(line_arr[2])

# Read 2: Fill in graph
G = nx.Graph()
G.add_nodes_from(max_chain_lengths.keys())


with open(ACR_FILE, 'r') as f:
    for line in f:
        line_arr = line.rstrip().split('\t')
        if line_arr[0] != line_arr[1] :
            if mode == 'smallest' or mode == 's':
                if min(max_chain_lengths[line_arr[0]], max_chain_lengths[line_arr[1]]) - int(line_arr[2]) >= min_weight :
                    G.add_edge(line_arr[0], line_arr[1], weight = min(max_chain_lengths[line_arr[0]], max_chain_lengths[line_arr[1]]) - int(line_arr[2]))
            else :
                if max(max_chain_lengths[line_arr[0]], max_chain_lengths[line_arr[1]]) - int(line_arr[2]) >= min_weight :
                    G.add_edge(line_arr[0], line_arr[1], weight = max(max_chain_lengths[line_arr[0]], max_chain_lengths[line_arr[1]]) - int(line_arr[2]))
# Display graph
weights = [G[u][v]['weight'] for u, v in G.edges()]


cmap = plt.cm.Blues
norm = mpl.colors.Normalize(vmin=min(weights), vmax=max(weights))


# Draw graph
pos = nx.circular_layout(G) 
nx.draw(G, pos, with_labels=False, node_size=20, edge_color=weights, edge_cmap=cmap, width=0.5, node_color='lightgray')

# Draw custom labels from the dictionary
nx.draw_networkx_labels(G, pos, labels=max_chain_lengths, font_size=6)

if mode == 'smallest' or mode == 's': 
    plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm = norm), label='Smallest Max Possible - Chain Length')
else :
    plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm = norm), label='Largest Max Possible - Chain Length')

plt.axis('off')
plt.savefig(f'/home/kyu/{ACR}_network.png')