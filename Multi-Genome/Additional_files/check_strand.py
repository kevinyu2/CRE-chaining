from pathlib import Path

def find_minus_strand(input_dir):
    search_dir = Path(input_dir)
    x_streme_folders = [folder for folder in search_dir.glob("*._xstreme") if folder.is_dir()]
    count = 0
    for folder in x_streme_folders:
        fimo_folders = [folder for folder in Path(f"{folder}").glob("fimo*") if folder.is_dir]
        for fimo in fimo_folders:
            with open(f"{folder}/{fimo.stem}/fimo.tsv") as file:
                for line in file:
                    line_arr = line.split("\t")
                    if len(line_arr) < 2:
                        continue
                    strand = line_arr[4]
                    if strand != "+":
                        print(f"Minus strand found in {fimo.stem} in {folder.stem}")
        count += 1
        print(f"Processed {count} x_streme folders")

def find_section_fna(file, start_ind, end_ind):
    START = 17054359 + start_ind - 1
    END = 17054359 + end_ind
    with open(file) as file:
        file.readline() # chromosome ID
        line = file.readline().strip()
        chromosome1 = line
        bp_count = len(line)
        while line[0] != ">":
            line = file.readline().strip()
            chromosome1 += line
            bp_count += len(line)
            if bp_count % 100000 == 0:
                print(bp_count)
            if bp_count > 17054559:
                break
    print(chromosome1[START:END])

def find_section_group(file, start_ind, end_ind, genome_ID):
    sequence = ""
    with open(file) as file:
        line = file.readline().strip()
        while line[0] == ">":
            if line[1:] == genome_ID:
                sequence = file.readline().strip()
                break
            file.readline()
            line = file.readline().strip()
    print(sequence[start_ind - 1 : end_ind])
        

#find_section_group("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/Chr1_15420017to15421784.group.fa", 74, 88, "IP-Sln-22_chr1_17054359to17054978")

#find_section_fna("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/GCA_946413935.1_IP-Sln-22.PacbioHiFiAssembly_genomic.fna", 74, 88)

#check_strand("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes")