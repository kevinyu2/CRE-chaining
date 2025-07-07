import random
from collections import defaultdict
from matplotlib import pyplot as plt #type:ignore
from pathlib import Path
import os
import time
import numpy as np #type:ignore

'''
Helper function. Returns the length of an ACR in the format Chr#_#####to#####
'''
def get_length(item):
    start = int(item[item.index("_") + 1 : item.index("to")])
    if "rand" in item:
        end = int(item[item.index("to") + 2 : item.index("_rand")])
    else:
        end = int(item[item.index("to") + 2 : ])
    length = abs(end - start) + 1
    return length

'''
Helper function. Takes in a file with a list of ACRs and creates a set
'''
def create_ref_set(input_path):
    ref_set = set()
    with open(input_path, "r") as file:
        for line in file:
            line = line.strip()
            if line != "":
                ref_set.add(line)
    return ref_set

'''
Takes in a file with a list of ACRs and creates a set. Only includes ACRs with lengths between
<low> and <high>
'''
def create_ref_set_filter(input_path, low, high):
    ref_set = set()
    with open(input_path, "r") as file:
        for line in file:
            line = line.strip()
            if line != "":
                length = get_length(line)
                if length >= low and length <= high:
                    ref_set.add(line)
    return ref_set

'''
Takes in a chaining file with all random regions and ACRs pairwise chained. Also takes in the
set of "reference" ACRs, i.e the ACRs which should be chained against the random regions and other
ACRs.
Non_ref_num is the number of non-reference and random regions there should be.

To restrict the dictionaries to only ACRs of a certain length, pass in parameters for low (lower bound
on the length) and high (upper bound on the length)
'''
def create_dicts(chaining_file, reference_set, non_ref_num, low=None, high=None):
    #For every non-reference ACR, this dictionary contains a list of scores of that ACR
    #chained with every reference ACR
    ACR_dict = defaultdict(list)
    #For every random sequence, this dictionary contains a list of scores of that random sequence
    #chained with every reference ACR
    rand_dict = defaultdict(list)

    with open(chaining_file, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 3:
                continue
            item1 = line_arr[0]
            item2 = line_arr[1]
            score = float(line_arr[2])

            if low != None: #skip if outside the bounds for the length
                len1 = get_length(item1)
                len2 = get_length(item2)
                if len1 < low or len1 > high:
                    continue
                if len2 < low or len2 > high:
                    continue

            if item1[-4:] == "rand":
                if item2 in reference_set:
                    rand_dict[item1].append(score)
            elif item1 in reference_set:
                if item2[-4:] == "rand":
                    rand_dict[item2].append(score)
                elif item2 not in reference_set:
                    ACR_dict[item2].append(score)
            else:
                if item2 in reference_set:
                    ACR_dict[item1].append(score)
        
    #Account for any reference ACRs which have chaining scores of 0 with some ACR / random region
    for lst in ACR_dict.values():
        missing_values = len(reference_set) - len(lst)
        lst.extend([0 for _ in range(missing_values)])
        assert(len(lst) == len(reference_set))
    for lst in rand_dict.values():
        missing_values = len(reference_set) - len(lst)
        lst.extend([0 for _ in range(missing_values)])
        assert(len(lst) == len(reference_set))
    
    #Account for any non-reference and random regions with chaining scores of 0 with every ref ACR
    for i in range(len(ACR_dict), non_ref_num):
        ACR_dict[f"unknown_{i}"].extend([0 for i in range(len(reference_set))])
    assert(len(ACR_dict) == non_ref_num)
    
    for i in range(len(rand_dict), non_ref_num):
        rand_dict[f"unknown_rand_{i}"].extend([0 for i in range(len(reference_set))])
    assert(len(rand_dict) == non_ref_num)
    
    return(ACR_dict, rand_dict)

'''
Outputs score frequencies to a file
'''
def output_all_score_freq(ACR_dict, rand_dict, output_dir):
    ACR = {}
    rand = {}
    for score_list in ACR_dict.values():
        for item in score_list:
            ACR[item] = ACR.get(item, 0) + 1
    for score_list in rand_dict.values():
        for item in score_list:
            rand[item] = rand.get(item, 0) + 1
    
    with open(f"{output_dir}/ACR_vs_ACR_all_freq.tsv", "w") as file:
        for score in ACR:
            file.write(f"{score}\t{ACR[score]}\n")
    
    with open(f"{output_dir}/rand_vs_ACR_all_freq.tsv", "w") as file:
        for score in rand:
            file.write(f"{score}\t{rand[score]}\n")

'''
Does some operation on each non-reference ACR / random region to get some number (e.g. max, average,
min, etc.). Outputs the frequencies of this number to a file.  
'''
def output_score_freq(ACR_dict, rand_dict, output_dir, list_op, op_name):
    ACR = {}
    rand = {}
    for score_list in ACR_dict.values():
        value = list_op(score_list)
        ACR[value] = ACR.get(value, 0) + 1
    for score_list in rand_dict.values():
        value = list_op(score_list)
        rand[value] = rand.get(value, 0) + 1
    
    with open(f"{output_dir}/ACR_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in ACR:
            file.write(f"{score}\t{ACR[score]}\n")
    
    with open(f"{output_dir}/rand_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in rand:
            file.write(f"{score}\t{rand[score]}\n")

'''
Reads frequencies from files and creates two side-by-side plots.
File1 should contain ACR chaining with ACR data.
File2 should contain ACR chaining with random data.
'''
def create_plots(file1, file2, title, x_label, max_index, min_index=0):
    data1 = [0 for i in range((max_index - min_index) + 1)]
    acr_outliers = 0
    rand_outliers = 0

    with open(file1, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data1[score-min_index] += freq
            except IndexError:
                acr_outliers += freq

    data2 = [0 for i in range((max_index - min_index) + 1)]
    with open(file2, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data2[score-min_index] += freq
            except IndexError:
                rand_outliers += freq

    try:
        assert(sum(data1) + acr_outliers == sum(data2) + rand_outliers)
    except:
        print(sum(data1) + acr_outliers)
        print(sum(data2) + rand_outliers)
        raise

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)

    axes[0].bar(range(min_index, max_index+1), data1, color='blue')
    axes[0].set_title("ACRs aligned with ACRs")
    axes[0].set_ylabel("Frequency")
    axes[0].set_xlabel(x_label)

    axes[1].bar(range(min_index, max_index+1), data2, color='green')
    axes[1].set_title("ACRs aligned with random regions")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlabel(x_label)

    #plt.axis('equal')
    plt.suptitle(title)
    plt.tight_layout()
    file_name = title.replace(" ", "_")
    plt.savefig(f"/home/mwarr/{file_name}.png")
    print(f"ACR outliers {acr_outliers}")
    print(f"Random outliers {rand_outliers}")

def get_average(file_path):
    sum = 0
    num_regions = 0
    with open(file_path) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                break
            sum += float(line_arr[0]) * int(line_arr[1])
            num_regions += int(line_arr[1])
    return sum / num_regions

def create_average_graph(input_dir_base):
    options = ["all", "max"]
    chain_type = ["local", "global"]
    labels = ["All, Local", "Max, Local", "All, Global", "Max, Global"]
    ACR_data = []
    rand_data = []
    for type in chain_type:
        for opt in options:
            ACR_data.append(get_average(f"{input_dir_base}/{type}/ACR_vs_ACR_{opt}_freq.tsv"))
            rand_data.append(get_average(f"{input_dir_base}/{type}/rand_vs_ACR_{opt}_freq.tsv"))
    print(ACR_data)
    print(rand_data)
    width = .4
    x = np.arange(len(labels))
    plt.bar(x - width / 2, ACR_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.title("Comparison of ACR Chaining Score Averages")
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.legend()
    plt.savefig("/home/mwarr/averages_all-sizes.png")

def create_average_subgraphs(input_dir_base):
    options = ["all", "max"]
    chain_type = ["local", "global"]
    labels = ["All, Local", "Max, Local", "All, Global", "Max, Global"]
    ACR_data = []
    rand_data = []
    for type in chain_type:
        for opt in options:
            ACR_data.append(get_average(f"{input_dir_base}/{type}/ACR_vs_ACR_{opt}_freq.tsv"))
            rand_data.append(get_average(f"{input_dir_base}/{type}/rand_vs_ACR_{opt}_freq.tsv"))
    width = .4
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    x = np.arange(2)

    
    axes[0].bar(x - width / 2, ACR_data[2:], label="ACR", width=width)
    axes[0].bar(x + width / 2, rand_data[2:], label="Rand", width=width)
    axes[0].set_title("Global Alignment Scores")
    plt.sca(axes[0])
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.xticks(x, labels[0:2])
    plt.ylim(11000, 12500)
    plt.legend()

    axes[1].bar(x - width / 2, ACR_data[0:2], label="ACR", width=width)
    axes[1].bar(x + width / 2, rand_data[0:2], label="Rand", width=width)
    axes[1].set_title("Local Alignment Scores")
    plt.sca(axes[1])
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.xticks(x, labels[2:])

    
    plt.legend()
    plt.savefig("/home/mwarr/averages_all-sizes.png")


def driver_all():
    # ref_set = create_ref_set("/home/mwarr/Data/seta.txt")
    # dicts = create_dicts("/home/mwarr/Data/alignment/global/alignment_scores_glob.tsv", ref_set, 15648)
    # print(f"Non reference ACRs: {len(dicts[0])}")
    # print(f"Random regions: {len(dicts[1])}")
    # output_all_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global")
    # output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", max, "max")
    # output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", lambda lst: sum(lst) / float(len(lst)), "avg")
    create_plots("/home/mwarr/Data/alignment/global/ACR_vs_ACR_all_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_all_freq.tsv", "All Global Alignment Scores", "Score", 11000, 10000)
    create_plots("/home/mwarr/Data/alignment/global/ACR_vs_ACR_avg_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_avg_freq.tsv", "Average Global Alignment Score", "Average Score", 11000, 10000)
    create_plots("/home/mwarr/Data/alignment/global/ACR_vs_ACR_max_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_max_freq.tsv", "Max Global Alignment Score", "Max Score", 11000, 10000)

def driver_filter():
    for i in range(650, 2050, 100):
        start_time = time.time()
        dir = f"/home/mwarr/Data/One_Genome/local_freq_{i}-{i+100}"
        os.mkdir(dir)
        ref_set = create_ref_set_filter("/home/mwarr/Data/seta.txt", i, i+100)
        non_ref_set = create_ref_set_filter("/home/mwarr/Data/setb.txt", i, i+100)
        dicts = create_dicts("/home/mwarr/Data/One_Genome/Chaining_one_acr_rand_loc.tsv", ref_set, len(non_ref_set), i, i+100)
        output_all_score_freq(dicts[0], dicts[1], dir)
        output_score_freq(dicts[0], dicts[1], dir, max, "max")
        output_score_freq(dicts[0], dicts[1], dir, lambda lst: sum(lst) / float(len(lst)), "avg")
        create_plots(f"{dir}/ACR_vs_ACR_all_freq.tsv", f"{dir}/rand_vs_ACR_all_freq.tsv", f"All Local Chaining Scores Lengths {i}-{i+100}", "Score")
        create_plots(f"{dir}/ACR_vs_ACR_avg_freq.tsv", f"{dir}/rand_vs_ACR_avg_freq.tsv", f"Average Local Chaining Score Lengths {i}-{i+100}", "Average Score")
        create_plots(f"{dir}/ACR_vs_ACR_max_freq.tsv", f"{dir}/rand_vs_ACR_max_freq.tsv", f"Highest Local Chaining Score Lengths {i}-{i+100}", "Max Score")
        print(f"Finished {i}-{i+100} in {time.time() - start_time} seconds", flush=True)


if __name__ == "__main__":
    create_average_subgraphs("/home/mwarr/Data/alignment")