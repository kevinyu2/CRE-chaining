import random
from matplotlib import pyplot as plt
import time

'''
Creates a dictionary where the keys are the chromosome identifiers
for the genome and the values are the sequences for the given chromosome
'''
def create_genome_dict(genome_file):
    #create genome dict {chr: seq, chr: seq, ...}
    start_time = time.time()
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
            print(f"Finished {chromosome} after {time.time() - start_time} seconds", flush=True)
            chromosome = line
    return genome_dict


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
            output.write(region)
            #newline if this isn't the last line of the file
            if i != num_regions - 1:
                output.write("\n")

def random_region_size_file(genome_file, size_file, output_file):
    print("Began program", flush=True)
    start_time = time.time()
    genome_dict = create_genome_dict(genome_file)
    print(f"Genome dict created after {time.time() - start_time}", flush=True)
    sizes = []
    count = 0

    with open(size_file, "r") as input:
        count += 1
        if count % 500 == 0:
            print(f"Finished {count} sequences after {time.time() - start_time} seconds", flush=True)
        with open(output_file, "w") as output:
            for line in input:
                start = line[line.index("_") + 1 : line.index("to")]
                end = line[line.index("to") + 2 : ]
                region_size = abs(int(end) - int(start)) + 1
                sizes.append(region_size)

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
                output.write(f">Chr{chromosome}_{rand_start}to{rand_end - 1}_rand\n")
                output.write(region)
                output.write("\n")
    return sizes
    #plt.hist(sizes, bins=1000)
    #plt.savefig("/home/mwarr/ref_sizes.png")



if __name__ == "__main__":
    genome_file = "/home/mwarr/Data/random_regions/tair10.fa"
    output_file = "/home/mwarr/Data/random_regions/rand_match-size_b.fa"
    size_file = "/home/mwarr/Data/setb.txt"
    random_region_size_file(genome_file, size_file, output_file)
