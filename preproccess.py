import sys
from pathlib import Path

search_dir = Path("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes")
root_output = ""

x_streme_folders = [folder for folder in search_dir.glob("*._xstreme") if folder.is_dir()]

for folder in x_streme_folders:
    for item in folder:
        print(item)

