import sys
from pathlib import Path
import pandas as pd
import os
import time

'''
Step 1 in the pipeline.

Organizes the ATAC-seq and XStreme data. Creates a folder for each ACR. Each folder
contains a file for each genome which lists the motifs and corresponding locations
'''

search_dir = Path("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes")
root_output = "/home/mwarr/Preprocessed_Data"

x_streme_folders = [folder for folder in search_dir.glob("*._xstreme") if folder.is_dir()]

count = 0
start_time = time.time()

for x_treme_folder in x_streme_folders:
    #create folder corresponding with the ACR in the root output directory
    start_index = str(x_treme_folder).index("Chr")
    end_index = str(x_treme_folder).index("._xstreme")
    acr_name = str(x_treme_folder)[start_index:end_index]
    acr_folder_name = f"{root_output}/{acr_name}"
    try:
        os.mkdir(acr_folder_name)
    except FileExistsError:
        print(f"{acr_folder_name} already exists.", file = sys.stderr)

    #Create files and headers for each genome
    genome_list_file = f"{acr_name}.group.txt"
    with open(str(search_dir) + "/" + genome_list_file) as input_file:
        for line in input_file:
            with open(f"{acr_folder_name}/{line.strip()}.tsv", "w") as output_file:
                output_file.write("motif\tstart\tstop\tstrand\n")

     
    #Find all the motifs for each genome in this ACR and record
    for fimo_folder in x_treme_folder.glob("fimo*"):
        with open(f"{fimo_folder}/fimo.tsv", 'r') as read_file:
            #read header
            try:
                read_file.readline()
            except:
                print(f"{fimo_folder}/fimmo.tsv is empty")
            #iterate through each row in the table
            for line in read_file:
                line_arr = line.split("\t")
                if len(line_arr) < 5:
                    break
                #get only the sequence (the sequence is sometimes preceeded by a number)
                index = 0
                for char in line_arr[0]:
                    if char.isalpha():
                        break
                    else:
                        index += 1
                sequence = line_arr[0][index:]

                #check for data consistency between the group text file and the fimo files
                file = Path(f"{acr_folder_name}/{line_arr[2]}.tsv")
                if not file.is_file():
                    raise Exception(f"The file {file} does not exist. The genome name is not correctly"
                    +"listed in the group.txt file")
                with open(f"{acr_folder_name}/{line_arr[2]}.tsv", 'a') as write_file:
                    #write to file -- motif    start    stop
                    write_file.write(f"{sequence.strip()}\t{line_arr[3]}\t{line_arr[4]}\t{line_arr[5]}\n")
    count += 1
    print(f"Finished {count} folders in {time.time() - start_time}")

        

