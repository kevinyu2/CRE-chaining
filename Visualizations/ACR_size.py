from pathlib import Path
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt


# Read motif size file
motif_df = pd.read_csv('/home/mwarr/motif_counts.tsv', sep = '\t')


motif_counts = []
acr_sizes = []
for i, row in motif_df.iterrows() :
    motif_counts.append(row['tomtom_motif_count'])
    start = int(row['genome'].split('to')[0].split('_')[-1])
    end = int(row['genome'].split('to')[1])
    acr_sizes.append(end - start)

# plt.figure()
# plt.hist(acr_sizes)
# plt.title("ACR size")
# plt.savefig("/home/mwarr/ACR_size.png")

plt.scatter(acr_sizes, motif_counts, s = 0.1)
plt.title("ACR Size vs Motif Count (Tomtom)")
plt.savefig('/home/kyu/motifs_v_size.png')