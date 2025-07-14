#Cluster a set of motifs
import numpy as np
import pandas as pd
import numpy as np
import re
from collections import defaultdict
# from biotite.sequence import NucleotideSequence
# from biotite.sequence.align import align_multiple

tomtom_out = '/home/kyu/tomtom_DAPv1.txt'
cluster_out = '/home/mwarr/Data/DAPv1_clusters.txt'
meme_file = '/home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme'

# tomtom_out = 'tomtom_DAPv1.txt'
# cluster_out = 'DAPv1_clusters.txt'
# meme_file = 'ArabidopsisDAPv1.meme'

CLUSTER_THRESH = 0.4
# meme_file = '/home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme'

#################################################################################
# CLUSTERING

# Load up the TOMTOM results
tomtom = pd.read_table(tomtom_out).rename(columns={'Query_ID': 'Query ID','Target_ID': 'Target ID', 'Optimal_offset' : 'Optimal offset',
                                                           'Query_consensus' : 'Query consensus', 'Target_consensus' : 'Target consensus'})
motifs = pd.unique(pd.concat([tomtom['Query ID'], tomtom['Target ID']]))

sim = tomtom.pivot_table(index='Query ID', columns='Target ID', values='E-value', fill_value=np.nan)
cols = sim.columns.intersection(motifs)
rows = sim.index.intersection(motifs)

sim = sim[cols].loc[rows]

x = sim.values

w = np.triu(x) +  np.triu(x, 1).T
v = np.tril(x) + np.tril(x, -1).T

sim.iloc[:,:] = np.nanmin(np.dstack([w, v]), axis=2)

sim.fillna(100, inplace=True)
sim = -np.log10(sim)
sim[np.isinf(sim)] = 10


from scipy.cluster.hierarchy import fcluster, linkage, dendrogram

Z = linkage(sim, method = 'complete', metric = 'correlation')

cl = fcluster(Z, CLUSTER_THRESH, criterion='distance')
o = dendrogram(Z, no_plot=True)['leaves']

print(f'Number of motif clusters: {max(cl)}')

# {Motif: cluster idx}
cluster_dict = {}
# {cluster idx : [motif, motif, ...]}
reverse_cluster_dict = defaultdict(list)

with open(cluster_out, 'w') as out :
    for idx, cnum in enumerate(cl) :
        out.write(f"{motifs[idx]}\tcluster_{cl[idx]}\n")
        cluster_dict[motifs[idx]] = f'cluster_{cl[idx]}'
        reverse_cluster_dict[f'cluster_{cl[idx]}'].append(motifs[idx])

###################################################################################
# PSSM ALIGNMENT

def parse_meme_pssms(filepath):
    with open(filepath) as f:
        lines = f.readlines()

    motifs = {}
    current_name = None
    current_pssm = []
    reading_matrix = False

    for line in lines:
        line = line.strip()

        # Detect new motif header
        if line.startswith("MOTIF"):
            parts = line.split()
            if len(parts) >= 2:
                current_name = parts[1]  # Use the first part as motif name
            current_pssm = []
            reading_matrix = False

        # Detect start of matrix
        elif line.startswith("letter-probability matrix"):
            reading_matrix = True
            current_pssm = []

        # Matrix lines
        elif reading_matrix:
            if re.match(r"^[0-9eE\.\-\s]+$", line):  # Match matrix rows
                row = list(map(float, line.split()))
                current_pssm.append(row)
            elif line.startswith("URL") or line == "":
                # End of current motif
                if current_pssm:
                    motifs[current_name] = np.array(current_pssm)
                    current_name = None
                    current_pssm = []
                    reading_matrix = False

    # Catch final motif (in case no URL after it)
    if current_name and current_pssm:
        motifs[current_name] = np.array(current_pssm)

    return motifs


pssms = (parse_meme_pssms(meme_file))

# Turn it into consensus sequence
def pssm_to_consensus(pssm, alphabet="ACGT"):
    idxs = np.argmax(pssm, axis=1)
    return ''.join(alphabet[i] for i in idxs)


from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import tempfile
import subprocess
import os

def run_msa(sequences):
    """
    Given a list of sequences (strings), return a BioPython alignment object.
    Requires MUSCLE installed and in PATH.
    """
    # 1. Write sequences to temp FASTA
    records = [SeqRecord(Seq(seq), id=f"seq{i}") for i, seq in enumerate(sequences)]
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".fasta") as fasta:
        SeqIO.write(records, fasta.name, "fasta")
        aligned_file = fasta.name + ".aln"

    muscle_path = os.path.expanduser("~/muscle")  # or absolute path
    cmd = f"{muscle_path} -align {fasta.name} -output {aligned_file}"

    # subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    result = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True  # ensures output is string, not bytes
    )

    # 3. Read alignment
    alignment = AlignIO.read(aligned_file, "fasta")

    # 4. Cleanup
    os.remove(fasta.name)
    os.remove(aligned_file)

    return alignment


def compute_consensus_pssm(msa, pssm_dict, id_to_motif_name):
    import numpy as np
    consensus_pssm = []
    pssm_positions = {}

    # Initialize pointer for each sequence
    for rec in msa:
        motif_name = id_to_motif_name[rec.id]
        pssm_positions[motif_name] = 0

    for col_idx in range(msa.get_alignment_length()):
        column_vector = []

        for rec in msa:
            base = rec.seq[col_idx]
            if base == "-":
                continue

            motif_name = id_to_motif_name[rec.id]
            pssm = pssm_dict[motif_name]
            pos = pssm_positions[motif_name]

            if pos < len(pssm):
                column_vector.append(pssm[pos])
                pssm_positions[motif_name] += 1

        if column_vector:
            avg_vector = np.mean(column_vector, axis=0).tolist()
        else:
            avg_vector = [0.0, 0.0, 0.0, 0.0]

        consensus_pssm.append(avg_vector)

    return consensus_pssm

pssm_out_dict = {}

for cluster_no, motif_list in reverse_cluster_dict.items() :
    curr_motif_consensus_sequences = []
    curr_motif_dict = {}

    for i, motif in enumerate(motif_list) :
        curr_motif_consensus_sequences.append(pssm_to_consensus(pssms[motif]))
        curr_motif_dict[f'seq{i}'] = motif

    if len(curr_motif_consensus_sequences) > 1 :
        msa = run_msa(curr_motif_consensus_sequences)
        pssm_out_dict[cluster_no] = compute_consensus_pssm(msa, pssms, curr_motif_dict)
    
    else :
        pssm_out_dict[cluster_no] = pssms[motif_list[0]]



def write_multiple_pssms_to_meme_file(pssm_dict, output_path, alphabet="ACGT"):
    with open(output_path, "w") as f:
        # Write header
        f.write("MEME version 4.4\n\n")
        f.write(f"ALPHABET= {alphabet}\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies (from /mnt/thumper-e1/home/shhuang/projects/dap/analysis.v4/gem07_rep_memechip03/ABI3VP1_tnt/AT5G18090_col_a/background):\n")
        f.write("A 0.30000 C 0.20000 G 0.20000 T 0.30000\n")

        # Write each PSSM
        for motif_name, pssm in pssm_dict.items():
            w = len(pssm)
            f.write(f"MOTIF {motif_name}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {w} nsites= {w} E= 0\n")
            for row in pssm:
                f.write("  " + "  ".join(f"{v:.6f}" for v in row) + "\n")
            f.write("\n")

write_multiple_pssms_to_meme_file(pssm_out_dict, '/home/mwarr/Data/DAPv1_clustered.meme')
