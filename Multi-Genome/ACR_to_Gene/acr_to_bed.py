from pathlib import Path
import sys


def acr_to_bed(output_file) :

    # Get all _xtreme folders
    search_dir = Path('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes')
    xstreme_dirs = [p for p in search_dir.glob('*') if p.is_dir() and p.name.endswith('._xstreme')]

    base_names = [p.name.replace('._xstreme', '') for p in xstreme_dirs]
    
    with open(output_file, 'w') as out:
        for name in base_names :
            stop = int(name.split('to')[1])
            start = int(name.split('to')[0].split('_')[1])
            
            # Swap if wrong order
            if stop < start :
                temp = stop
                stop = start
                start = temp
                print("swap")
            chromosome = name.split('_')[0]
            out.write(f"{chromosome}\t{start}\t{stop}\t{name}\n")

acr_to_bed("acr_regions.bed")