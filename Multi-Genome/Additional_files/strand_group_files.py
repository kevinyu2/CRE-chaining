from pathlib import Path
import time
from Bio import Align

def create_group_files(input_dir, output_dir):
    #create dictionary of ACRs and what strand they are on
    strand_dict = {}
    search_dir = Path(input_dir)
    start_time = time.time()
    for file in search_dir.glob("*dACRs.paf"):
        with open(file) as paf_file:
            for line in paf_file:
                line_arr = line.split("\t")
                if len(line_arr) < 2:
                    continue
                name = line_arr[0]
                start = line_arr[2]
                end = line_arr[3]
                strand = line_arr[4]
                strand_dict[f"{name}_{start}to{end}"] = strand


    start_time = time.time()
    for file in search_dir.glob("*.group.txt"):
        with open(f"{output_dir}/{file.stem}.txt", "w") as output_file:
            with open(file) as group_file:
                for line in group_file:
                    try:
                        strand = strand_dict[line.strip()]
                    except:
                        print(line)
                        print(file)
                        raise
                    output_file.write(f"{line.strip()}\t{strand}\n")

create_group_files("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes", "/home/mwarr/Strand_Group_Files")