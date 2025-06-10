from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import Align
import sys

'''
Step 1 in the pipeline.

Creates a data file listing all pairs of motifs from the genomes which share a motif in the 
TOMTOM database. Also checks alignment and outputs a score.
'''

# Find pairs of motifs that each share a matched tomtom motif
# Dict is structured as: tomtom_motif : {acr_motif_1, acr_motif_2}
# Pairing is structured as a set of {(motif1, motif2)}, where each tuple is in sorted order
def create_motif_pairing() :
    print("######################", file=sys.stderr)
    print("Finding Motif Matches in Tomtom", file=sys.stderr)

    # Get all _xtreme folders
    search_dir = Path('../../projects/msu_nsf_pangenomics/pgrp/dACRxgenomes')
    xstreme_dirs = [p for p in search_dir.glob('*') if p.is_dir() and p.name.endswith('._xstreme')]

    motif_dict = defaultdict(set)
    motif_pairs = set()

    for x_dir in xstreme_dirs:
        tomtom_path = x_dir / 'meme_tomtom_out' / 'tomtom.tsv'
        if tomtom_path.exists():
            try : 
                tomtom_df = pd.read_csv(tomtom_path, sep = '\t')

                tomtom_df = tomtom_df.dropna(subset=['Target_consensus'])
                for i, row in tomtom_df.iterrows() :
                    motif_dict[row['Target_consensus']].add(row['Query_ID'])

                    # Add all motif pairs 
                    for motif in motif_dict[row['Target_consensus']] :
                        motif_list_sorted = sorted([motif, row['Query_ID']])
                        motif_pairs.add(tuple(motif_list_sorted))

            # Usually means the tsv file failed and is empty
            except :
                print(f"File failed: {tomtom_path}")

    # print(motif_pairs)

    print(f"Num. Motifs: {len(motif_dict.keys())}", file=sys.stderr)
    print(f"Num. Motif Pairs: {len(motif_pairs)}", file=sys.stderr)
    print("######################", file=sys.stderr)

    return motif_pairs

def align(motif_pairs):
    print("Aligning All Matches", file=sys.stderr)

    print("sequence_1\tsequence_2\tsequence_1_start\tsequence_2_start\tscore")
    aligner = Align.PairwiseAligner(scoring="blastn", mode="local")

    for pair in motif_pairs :
        alignments = aligner.align(pair[0], pair[1])
        best = alignments[0]

        seq1_aligned, seq2_aligned = best.aligned

        seq1_start = seq1_aligned[0][0]
        seq2_start = seq2_aligned[0][0]


        score = best.score
        score /= min(len(pair[0]), len(pair[1]))
        print(f"{pair[0]}\t{pair[1]}\t{seq1_start}\t{seq2_start}\t{score}")


# Driver
motif_pairs = create_motif_pairing()
align(motif_pairs)



