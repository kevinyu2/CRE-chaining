from pathlib import Path
import time
from matplotlib import pyplot as plt
import random

'''
Step 4 in the pipeline.

Takes in the files with the chain length data. Creates files with the frequencies of 
chain length, number of anchors, and ratio of chain length to number of anchors. Creates
histograms with this data.
'''

def output_frequency_files(input_dir, output_dir):
    chain_length = {}
    anchor_length = {}
    ratios = {}
    input_dir = Path(input_dir)
    count = 0
    start_time = time.time()
    for file_name in input_dir.glob("*"):
        count += 1
        if count % 100 == 0:
            print(f"{count} files parsed in {time.time() - start_time}", flush=True)
        with open(file_name) as acr_file:
            for line in acr_file:
                line_arr = line.split("\t")
                #count frequencies of each stat
                if len(line_arr) < 2:
                    continue
                chain_len = round(float(line_arr[2]))
                chain_length[chain_len] = chain_length.get(chain_len, 0) + 1
                anchor_length[int(line_arr[3])] = anchor_length.get(int(line_arr[3]), 0) + 1
                ratio = float(line_arr[2]) / float(line_arr[3])
                ratios[round(ratio, 2)] = ratios.get(round(ratio, 2), 0) + 1
    
    #write to files to save information
    with open(f"{output_dir}/chain_length_frequencies.txt", "w") as file:
        for key, value in chain_length.items():
            file.write(f"{key}\t{value}\n")
    
    with open(f"{output_dir}/anchor_length_frequencies.txt", "w") as file:
        for key, value in anchor_length.items():
            file.write(f"{key}\t{value}\n")
    
    with open(f"{output_dir}/ratio_frequencies.txt", "w") as file:
        for key, value in ratios.items():
            file.write(f"{key}\t{value}\n")

    
    return (chain_length, anchor_length, ratios)

def input_frequency_file(filename):
    frequency_dict = {}
    with open(filename) as file:
        for line in file:
            line_arr = line.split("\t")
            if (len(line_arr) > 1):
                frequency_dict[float(line_arr[0])] = float(line_arr[1])
    return frequency_dict


def make_histogram(frequency_dict, filename, width, xaxis, title):
        print("Generating graph")
        items_sorted = sorted(frequency_dict.items())
        num, freq = zip(*items_sorted)
        plt.bar(num, freq, width=width)
        plt.yscale("log")
        plt.title(title)
        plt.xlabel(xaxis)
        plt.ylabel("Frequency")
        plt.savefig(filename)
        plt.clf()


def generate_test_data(filename, num_lines):
    with open(filename, "w") as file:
        for i in range(num_lines):
            file.write(f"hi\thi\t{random.randint(1, 10)}")
            file.write(f"\t{random.randint(10, 15)}\n")



dictionaries = output_frequency_files("../Chaining_min1_intra", "../Frequencies_min1_intra")
chain_length = dictionaries[0]
anchor_length = dictionaries[1]
ratios = dictionaries[2]
#chain_length = input_frequency_file("../Frequencies_min1/chain_length_frequencies.txt")
#anchor_length = input_frequency_file("../Frequencies_min1/anchor_length_frequencies.txt")
#ratios = input_frequency_file("../Frequencies_min1/ratio_frequencies.txt")

make_histogram(chain_length, 'chain_min1_intra.png', 1, "Chain Length", "Chain Length Frequencies (Intra-ACR)")
make_histogram(anchor_length, 'anchor_min1_intra.png', 1, "Number of Anchors", "Frequencies of the Number of Anchors (Intra-ACR)")
make_histogram(ratios, 'ratios_min1_intra.png', .01, "Ratio of Chain Length to Number of Anchors", "Ratio Frequencies (Intra-ACR)")
