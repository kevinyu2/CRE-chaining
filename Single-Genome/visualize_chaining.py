from pathlib import Path
from matplotlib import pyplot as plt

def output_frequency_files(file_name, output_dir):
    chain_length = {}
    anchor_length = {}
    ratios = {}
    with open(file_name) as acr_file:
            for line in acr_file:
                line_arr = line.split("\t")
                #count frequencies of each stat
                if len(line_arr) < 2:
                    continue
                chain_len = round(float(line_arr[2]))
                chain_length[chain_len] = chain_length.get(chain_len, 0) + 1
                ratio = float(line_arr[2]) / float(line_arr[3])
                ratios[round(ratio, 2)] = ratios.get(round(ratio, 2), 0) + 1
    
    
    #write to files to save information
    with open(f"{output_dir}/chain_length_frequencies.txt", "w") as file:
        for key, value in chain_length.items():
            file.write(f"{key}\t{value}\n")
    
    with open(f"{output_dir}/ratio_frequencies.txt", "w") as file:
        for key, value in ratios.items():
            file.write(f"{key}\t{value}\n")
    
    return (chain_length, ratios)

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

dictionaries = output_frequency_files("/home/mwarr/Data/One_Genome/Chaining_one_local.tsv", "/home/mwarr/Data/One_Genome")
chain_length = dictionaries[0]
ratios = dictionaries[1]


make_histogram(chain_length, '../Visualizations/Visuals/Single_Genome/chain_one_local.png', 1, "Chain Length", "One Genome Chain Length Frequencies, Local Chaining")
make_histogram(ratios, '../Visualizations/Visuals/Single_Genome/ratios_one_local.png', .01, "Ratio of Chain Length to Number of Anchors", "One Genome Ratio Frequencies, Local Chaining")
