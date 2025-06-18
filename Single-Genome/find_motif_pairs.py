from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import Align
import sys
import itertools

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
    xstreme_dir = Path('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme')
    xstreme_dirs = [p for p in xstreme_dir.glob('*') if p.is_dir() and p.name.startswith('fimo_out')]

    motifs = []

    for x_dir in xstreme_dirs:
        fimo = x_dir / 'fimo.tsv'

        if fimo.exists():
            try : 
                fimo_df = pd.read_csv(fimo, sep = '\t')
                motif = fimo_df['motif_id'][1].split('-')[1]
                motifs.append(motif)

                


            # Usually means the tsv file failed and is empty
            except :
                print(f"File failed: {fimo}")




    return motifs

def align(motifs):
    print("Aligning All Motifs", file=sys.stderr)

    print("sequence_1\tsequence_2\tsequence_1_start\tsequence_2_start\tscore")
    aligner = Align.PairwiseAligner(scoring="blastn", mode="local")

    for a, b in itertools.combinations(motifs, 2):
        alignments = aligner.align(a, b)
        try :
            best = alignments[0]

            seq1_aligned, seq2_aligned = best.aligned

            seq1_start = seq1_aligned[0][0]
            seq2_start = seq2_aligned[0][0]


            score = best.score
            score /= min(len(a), len(b))
            if score >= 1:
                print(f"{a}\t{b}\t{seq1_start}\t{seq2_start}\t{score}")
        except: 
            print(f"Failed: {a}, {b}")

# Driver
motifs = create_motif_pairing()
align(motifs)



