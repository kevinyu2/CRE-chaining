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
Finds <num_regions> random regions of size between <low_bound> and <up_bound> (inclusive)
in a genome and outputs a fasta file.
Genome file should be a fasta file of the entire genome (each entry is a chromosome)
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

'''
Outputs a fasta file of random regions, where the sizes of the regions match 
the sizes of the regions in <reference_file>. <reference_file> should be a list of
regions in the form chr#_<start>to<end>.
'''
def random_region_match_sizes(genome_file, reference_file, output_file):
    print("Began program", flush=True)
    start_time = time.time()
    genome_dict = create_genome_dict(genome_file)
    print(f"Genome dict created after {time.time() - start_time}", flush=True)
    sizes = []
    count = 0

    with open(reference_file, "r") as input:
        count += 1
        if count % 500 == 0:
            print(f"Finished {count} sequences after {time.time() - start_time} seconds", flush=True)
        with open(output_file, "w") as output:
            for line in input:
                #get the size
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

'''
Takes in a fasta file and creates a text file where each line is the name of a 
region in the fasta file
'''
def make_random_list(input_fa, output_file):
    with open(input_fa) as input:
        with open(output_file, "w") as output:
            for line in input:
                if line[0] == ">":
                    output.write(f"{line[1:].strip()}\n")


if __name__ == "__main__":
    make_random_list("/home/mwarr/Data/One_Genome/random_regions/setb_10_match_sizes.fa", "/home/mwarr/Data/One_Genome/random_regions/setb_rand.txt")

