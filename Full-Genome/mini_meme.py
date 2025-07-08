with open('/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa', 'r') as input:
    chr_lines = {}
    current_chr = ""
    for line in input:
        if '>' in line :
            current_chr = line.rstrip().split('>')[1]
            chr_lines[current_chr] = []
        chr_lines[current_chr].append(line)

    for chr in {"Chr1", "Chr2", "Chr3", "Chr4", "Chr5"} :
        with open(f'/home/kyu/single_chrs/{chr}', 'w') as out:
            for line in chr_lines[chr] :
                out.write(line)

        
        
        # if curr_motif not in seen_motifs :

        