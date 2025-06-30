import random


def create_genome_dict(genome_file):
    #create genome dict {chr: seq, chr: seq, ...}
    genome_dict = {}
    with open(genome_file) as genome:
        chromosome = genome.readline().strip()
        while chromosome:
            chromosome = chromosome[1:]
            seq = ""
            line = genome.readline().strip()
            while line != "" and line[0] != ">":
                seq += line
                line = genome.readline().strip()
            genome_dict[chromosome] = seq
            chromosome = line


'''
Finds <num_regions> random regions of size <region_size> in a genome and outputs 
a fasta file.
Genome file should be a fasta file of the entire genome (each entry is a chromosome)
up_bound is the upper bound for the region size.
low_bound is the lower bound for the region size
'''
def random_region_fasta(genome_file, num_regions, low_bound, up_bound, output_file):
    genome_dict = create_genome_dict(genome_file)
#Write to a fasta file
    with open(output_file, "w") as output:
        for i in range(num_regions):
            #choose a random chromosome
            chromosome = random.randint(1, 5)
            seq = genome_dict[f"Chr{chromosome}"]
            region_size = random.randint(low_bound, up_bound)
            #get the max value for the range of random numbers to generate
            max_rand = len(seq) - region_size - 1
            #generate random region indices
            rand_start = random.randint(0, max_rand)
            rand_end = rand_start + region_size
            #find region
            region = seq[rand_start : rand_end]
            #write to fasta file
            output.write(f">Chr{chromosome}_{rand_start}to{rand_end}_rand\n")
            output.write(seq)
            #newline if this isn't the last line of the file
            if i != num_regions - 1:
                output.write("\n")




genome_file = "/home/mwarr/Data/random_regions/tair10.fa"
output_file = "/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/random_900-1100.fa"
random_region_fasta(genome_file, 2994, 900, 1100, output_file)
