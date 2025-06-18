from pathlib import Path
from matplotlib import pyplot as plt

def get_num_motifs(input_dir, output_dir):
    search_dir = Path(input_dir)
    ACR_dict = {}

    for folder in search_dir.glob("fimo*"):
        with open(f"{str(folder)}/fimo.tsv", 'r') as fimo_file:
            fimo_file.readline() #skip header
            for line in fimo_file:
                line_arr = line.split("\t")
                if len(line_arr) < 3: 
                    continue
                ACR = line_arr[2]
                ACR_dict[ACR] = ACR_dict.get(ACR, 0) + 1
    
    motif_num = list(ACR_dict.values())
    plt.xlim(0, 40)
    plt.hist(motif_num, bins=range(40))
    #plt.savefig("/home/mwarr/motif_hist.png")

    with open(f"{output_dir}/single_genome_motif_counts.tsv", 'w') as file:
        for ACR, num in ACR_dict.items():
            file.write(f"{ACR}\t{num}\n")
    

OUTPUT_DIR = "/home/mwarr/Data"
INPUT_DIR = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme"
get_num_motifs(INPUT_DIR, OUTPUT_DIR)