from random_regions import random_region_size_file
import random

def generate_genome_test():
    bases = ['A', 'C', 'T', 'G']
    with open("genome_test.fa", "w") as genome:
        for i in range(1, 6):
            genome.write(f">Chr{i}\n")
            for j in range(10000):
                genome.write(bases[random.randint(0, 3)])
            if i != 5:
                genome.write("\n")

def generate_ACR_sizes():
    with open("acr_sizes_test.txt", "w") as file:
        for i in range(1000):
            end = random.randint(10, 200)
            file.write(f"Chr1_1to{end}\n")

if __name__ == "__main__":
    #generate_genome_test()
    #generate_ACR_sizes()
    sizes = random_region_size_file("genome_test.fa", "acr_sizes_test.txt", "test.out")
    #recreate genome_dict
    genome_dict = {}
    with open("genome_test.fa") as genome:
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
    with open("test.out") as file:
        count = 0
        for line in file:
            if line[0] == ">":
                chr = line[1 : 5]
                start = int(line[line.index("_") + 1 : line.index("to")])
                end = int(line[line.index("to") + 2 : line.index("_rand")])
                length = abs(int(end) - int(start)) + 1
            elif line[0] != "\n":
                assert(len(line.strip()) == sizes[count])
                assert(sizes[count] == length)
                assert(genome_dict[chr][start : end + 1] == line.strip())
                count += 1
    print("Finished test ")

        
