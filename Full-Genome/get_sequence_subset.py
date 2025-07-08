set_file = '/home/mwarr/Data/seta.txt'

set_to_use = set()

with open(set_file, 'r') as sf :
    for line in sf:
        set_to_use.add(line.rstrip())

with open('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/dACR.fa', 'r') as full :
    with open('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/seta.fa', 'w') as out :
        curr_acr = ""
        for line in full :
            if ">" in line :
                curr_acr = line.rstrip().split('>')[1] 
                if curr_acr in set_to_use :
                    out.write(line)
            else :
                if curr_acr in set_to_use :
                    out.write(line)
