import random

'''
Finds <num_regions> random regions of size <region_size> in a genome and outputs 
a fasta file.
Genome file should be a fasta file of the entire genome (each entry is a chromosome)
'''

def random_region_fasta(genome_file, num_regions, region_size, output_dir):
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

#Write to a fasta file
    with open(f"{output_dir}/random_length_{region_size}.fa", "w") as output:
        for i in range(num_regions):
            #choose a random chromosome
            chromosome = random.randint(1, 5)
            seq = genome_dict[f"Chr{chromosome}"]
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
output_dir = "/home/mwarr/Data/random_regions"
NUM_REGIONS =  35395 // 2 #number of genes near ACRS // 2
random_region_fasta(genome_file, NUM_REGIONS, 3000, output_dir)
