import random
import time

'''
Finds all the genes which have an ACR within <threshold> distance upstream
of the gene and outputs the genes' positions and strands to a file
'''
def genes_near_ACRs(input_file_path, output_dir, threshold):
    far_gene_count = 0

    with open(f"{output_dir}/genes_threshold_{threshold}.tsv", "w") as gene_file:
        gene_file.write("Chromosome\tStart\tStop\tStrand")
        with open(f"{output_dir}/pseudo_threshold_{threshold}.tsv", "w") as pseudo_file:
            pseudo_file.write("Chromosome\tStart\tStop\tStrand")

            with open(input_file_path, "r") as input_file:
                for line in input_file:
                    line_arr = line.split("\t")
                    chromosome = line_arr[3]
                    gene_type = line_arr[5]
                    start = line_arr[6]
                    stop = line_arr[7]
                    strand = line_arr[9]
                    distance = abs(int(line_arr[-1]))
                    
                    if distance < threshold:
                        if gene_type == "gene":
                            gene_file.write(f"\n{chromosome}\t{start}\t{stop}\t{strand}")
                        elif gene_type == "pseudogene":
                            pseudo_file.write(f"\n{chromosome}\t{start}\t{stop}\t{strand}")
                        else:
                            print(f"Found {gene_type} gene type") #in case there are other gene types
                    else:
                        far_gene_count += 1
    print(f"There were {far_gene_count} far genes with threshold {threshold}")

'''
Creates a fasta file of the regions <region_length> base pairs upstream of the genes
in the input file. Genome file should be a fasta file of the entire genome.
Threshold is the distance from a gene which must contain an ACR in order for that gene 
to be included in the previous step.
'''
def upstream_region_fasta(input_file, genome_file, output_dir, region_length, threshold):
    genome_dict = {}
    start_time = time.time()
    with open(genome_file) as genome:
        chromosome = genome.readline().strip()
        while chromosome:
            chromosome = chromosome[1:]
            print(f"{chromosome}, time: {time.time() - start_time}")
            seq = ""
            line = genome.readline().strip()
            count = 0
            while line != "" and line[0] != ">":
                seq += line
                line = genome.readline().strip()
                count += len(line)
            genome_dict[chromosome] = seq
            chromosome = line
    print("Finished dictionary")
    
    with open(f"{output_dir}/upstream_regions_{region_length}_threshold{threshold}.fa", "w") as output:
        with open(input_file) as gene_pos:
            gene_pos.readline() #headers
            count = 0
            start_time = time.time()
            for line in gene_pos:
                count += 1
                line_arr = line.split("\t")
                chromosome = line_arr[0]
                if len(line_arr) < 4: #could be a newline at the end of the file
                    continue
                #Find upstream region
                if line_arr[3].strip() == "+": #plus strand
                    region_end = int(line_arr[1]) #beginning of gene
                    region_start = region_end - region_length + 1
                else: #minus strand
                    region_start = int(line_arr[2]) #beginning of gene (with respect to transcription direction)
                    region_end = region_start + region_length - 1
                
                #Output to file
                output.write(f">{chromosome}_{region_start}to{region_end}\n")
                sequence = genome_dict[chromosome][region_start - 1 : region_end]
                output.write(f"{sequence}\n")
                if count % 1000 == 0:
                    print(f"Finished {count} sequences in {time.time() - start_time} seconds")
            
            
def generate_test_files():
    #generate genome file
    bases = ['A', 'C', 'G', 'T']
    with open("test.fa", "w") as output:
        #generate sequence ID
        for i in range(1, 6):  
            output.write(f">Chr{i}\n")
            #generate random sequence
            seq = ""
            for i in range(50):
                seq += bases[random.randint(0, 3)]
            output.write(f"{seq}\n")
    #generate gene position file
    with open("gene_pos_test.tsv", "w") as output:
        output.write("header\n")
        for i in range(30):
            start_region = random.randint(10, 30)
            output.write(f"Chr{random.randint(1, 5)}\t{start_region}\t{start_region + 20}")
            if random.randint(0, 1) == 0:
                output.write("\t+\n")
            else:
                output.write("\t-\n")

#generate_test_files()
upstream_region_fasta("/home/mwarr/Data/nearby_genes/gene_locations/genes_threshold_2000.tsv", "/home/mwarr/Data/nearby_genes/tair10.fa", "./", 3000, 2000)          
