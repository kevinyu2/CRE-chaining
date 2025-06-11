from pathlib import Path
import time
from Bio import Align

def create_sequence_files(input_dir, output_dir):
    seq_dict = {}
    search_dir = Path(input_dir)
    start_time = time.time()
    for file in search_dir.glob("*dACRs.fa"):
        with open(file) as seq_file:
            line = seq_file.readline()
            while line and line[0] == ">":
                sequence_name = line[1:].strip()
                sequence = seq_file.readline().strip().upper()
                seq_dict[sequence_name] = sequence 
                line = seq_file.readline()
    print(f"Finished in {time.time() - start_time} seconds")

    count = 0
    start_time = time.time()
    for file in search_dir.glob("*group.txt"):
        with open(f"{output_dir}/{file.stem}.fa", "w") as output_file:
            with open(file) as group_file:
                for line in group_file:
                    sequence = seq_dict[line.strip()]
                    output_file.write(f">{line}{sequence}\n")
        count += 1
        print(f"Finished {count} files in {time.time() - start_time} seconds")

    

create_sequence_files("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes", "../ACR_sequences")