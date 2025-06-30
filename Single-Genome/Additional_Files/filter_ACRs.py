from pathlib import Path
import os
from glob import glob
from matplotlib import pyplot as plot


def filter(INPUT_DIR, OUTPUT_DIR, ACR_length):
    short_ACRs = set()
    long_ACRs = set()

    for fimo_folder in INPUT_DIR.glob("fimo_out*"):
        new_folder_path = f"{str(OUTPUT_DIR)}/{fimo_folder.stem}"
        try:
            os.mkdir(new_folder_path)
        except Exception as e:
            print(f"Error: {e}")
        with open(f"{new_folder_path}/fimo.tsv", "w") as output_file:
            with open(f"{fimo_folder}/fimo.tsv", "r") as input_file:
                output_file.write(input_file.readline()) #copy the header
                for line in input_file:
                    line_arr = line.split("\t")
                    if len(line_arr) < 9:
                        output_file.write(line) #copy the line if it's not a data line
                        continue
                    #otherwise, check if the ACR is within the threshold
                    acr = line_arr[2]
                    if acr in short_ACRs:
                        output_file.write(line)
                    else:
                        ACR_start = int(acr[acr.index("_") + 1 : acr.index("to")])
                        ACR_end = int(acr[acr.index("to") + 2 : ])
                        if abs(ACR_end - ACR_start) < ACR_length + 100 and abs(ACR_end - ACR_start) > ACR_length - 200:
                            short_ACRs.add(acr)
                            output_file.write(line)
                        else:
                            long_ACRs.add(acr)

    print(len(short_ACRs))
    print("\n\n")
    print(len(long_ACRs))
    


def count_ACRs(INPUT_DIR):
    short_ACRs = set()
    long_ACRs = set()
    lengths = []
    pattern = os.path.join(INPUT_DIR, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob(pattern)

    for fimo_file in fimo_files:
        with open(fimo_file) as input_file:
            next(input_file) #skip header
            for line in input_file:
                line_arr = line.split("\t")
                if len(line_arr) < 9:
                    continue
                #check if the ACR is within the threshold
                acr = line_arr[2]
                ACR_start = int(acr[acr.index("_") + 1 : acr.index("to")])
                ACR_end = int(acr[acr.index("to") + 2 : ])
                length = abs(ACR_end - ACR_start)
                lengths.append(length)
                if length < 550 and length > 450:
                    short_ACRs.add(acr)
                else:
                    long_ACRs.add(acr)
    plot.hist(lengths, bins=1000)
    plot.title("ACR lengths for reference genome only")
    plot.xlabel("ACR lengths")
    plot.ylabel("Frequency")
    plot.savefig("/home/mwarr/ACR_lengths_ref.png")
        
    # print(len(short_ACRs))
    # print(short_ACRs)
    # print("\n\n")
    # print(len(long_ACRs))
    # print(long_ACRs)

INPUT_DIR = Path("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme/")
OUTPUT_DIR = Path("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/ACR-filter_rand/")

filter(INPUT_DIR, OUTPUT_DIR, 1000)
#count_ACRs("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/xstreme")