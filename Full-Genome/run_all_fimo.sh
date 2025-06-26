#!/bin/bash

MAP_FILE="/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motifs.txt"
WANTED_LIST="/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/motifs.txt"
FASTA="/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa"
BGFILE="/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/background"
OUT_BASE="/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/fimo_outs/fimo_out"
BASE_PATH="/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme"

i=1
while IFS=$'\t' read -r motif_id xml_file alt_id; do
    full_xml_path="$BASE_PATH/$xml_file"
    out_dir=$(printf "%s_%03d" "$OUT_BASE" "$i")
    echo "Running FIMO for motif $alt_id (motif_id: $motif_id) from $full_xml_path → $out_dir"

    fimo --verbosity 1 \
         --oc "$out_dir" \
         --bgfile "$BGFILE" \
         --motif "$motif_id" \
         "$full_xml_path" \
         "$FASTA"

    ((i++))
done < "$MAP_FILE"
