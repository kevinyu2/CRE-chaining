import random
from collections import defaultdict
from matplotlib import pyplot as plt #type:ignore
from pathlib import Path
import os
import time
import numpy as np #type:ignore
from useful_functions import *


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

            if low != None: 
                len1 = get_length(item1)
                len2 = get_length(item2)
                #skip if outside the bounds for the length
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
        if type(value) == list:
            for val in value:
                ACR[val] = ACR.get(val, 0) + 1
        else:
            ACR[value] = ACR.get(value, 0) + 1
    for score_list in rand_dict.values():
        value = list_op(score_list)
        if type(value) == list:
            for val in value:
                rand[val] = rand.get(val, 0) + 1
        else:
            rand[value] = rand.get(value, 0) + 1
    
    with open(f"{output_dir}/ACR_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in ACR:
            file.write(f"{score}\t{ACR[score]}\n")
    
    with open(f"{output_dir}/rand_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in rand:
            file.write(f"{score}\t{rand[score]}\n")

'''
Gets the average chaining score from a file with the frequencies of chaining scores
'''
def get_average(file_path, exclude_zero=False):
    sum = 0
    num_regions = 0
    with open(file_path) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                break
            if exclude_zero and float(line_arr[0]) == 0:
                continue
            sum += float(line_arr[0]) * int(line_arr[1])
            num_regions += int(line_arr[1])
    return sum / num_regions

'''
Creates a bar graph showing the average ACR v ACR score and ACR v random score,
as well as the average max ACR v ACR score and ACR v random score. Shows both
local and global chaining score averages.
'''
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

'''
Creates a graph similar to the function above (create_average_graph). However, it will
plot local and global on different subplots so that the axes and scale can be different 
for each.
'''
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

    
    axes[0].bar(x - width / 2, ACR_data[0:2], label="ACR", width=width)
    axes[0].bar(x + width / 2, rand_data[0:2], label="Rand", width=width)
    axes[0].set_title("Local Alignment Scores")
    plt.sca(axes[0])
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.xticks(x, labels[0:2])
    plt.legend()

    axes[1].bar(x - width / 2, ACR_data[2:4], label="ACR", width=width)
    axes[1].bar(x + width / 2, rand_data[2:4], label="Rand", width=width)
    axes[1].set_title("Global Alignment Scores")
    plt.sca(axes[1])
    plt.xlabel("Data Averaged and Chain Type")
    plt.xticks(x, labels[2:4])
    plt.ylim(11000, 12500)
    plt.legend()
    
    fig.suptitle("Comparison of Alignment Score Averages")
    plt.savefig("/home/mwarr/averages_alignment.png")

def average_graph_top(input_dir_base, top_num):
    acr_data = []
    rand_data = []
    labels = ["Max", "2nd-Highest", "3rd-Highest"]
    for i in range(1, top_num + 1):
        if i > 3:
            labels.append(f"{i}th-Highest")
        acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_{i}-highest_freq.tsv"))
        rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_{i}-highest_freq.tsv"))
    # labels.append(f"Top-{top_num}-Highest")
    # acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_top-5-avg_freq.tsv"))
    # rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_top-5-avg_freq.tsv"))
    print(acr_data)
    print(rand_data)
    width = .4
    x = np.arange(len(labels))
    plt.figure(figsize=(15, 6))
    plt.bar(x - width / 2, acr_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.title("Average Top Alignment Scores, Local")
    plt.xlabel("Data Averaged")
    plt.ylabel("Average Score")
    plt.legend()
    plt.savefig("/home/mwarr/align_exp2_loc_original.png")

'''
Function to be passed into output_score_freq as <list_op>
'''
def top_avg(lst):
    lst.sort(reverse=True)
    sum = 0
    for i in range(10):
        sum += lst[i]
    return sum / 10

def driver_all_summary():
    # ref_set = create_ref_set("/home/mwarr/Data/seta.txt")
    # dicts = create_dicts("/home/mwarr/Data/alignment/global/alignment_scores_glob.tsv", ref_set, 15648)
    # print(f"Non reference ACRs: {len(dicts[0])}")
    # print(f"Random regions: {len(dicts[1])}")
    # output_all_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global")
    # output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", max, "max")
    # output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", lambda lst: sum(lst) / float(len(lst)), "avg")
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_all_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_all_freq.tsv", "All Global Alignment Scores", "Score", 11000, 10000)
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_avg_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_avg_freq.tsv", "Average Global Alignment Score", "Average Score", 11000, 10000)
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_max_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_max_freq.tsv", "Max Global Alignment Score", "Max Score", 11000, 10000)

def driver_filter_summary():
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
        create_histograms(f"{dir}/ACR_vs_ACR_all_freq.tsv", f"{dir}/rand_vs_ACR_all_freq.tsv", f"All Local Chaining Scores Lengths {i}-{i+100}", "Score")
        create_histograms(f"{dir}/ACR_vs_ACR_avg_freq.tsv", f"{dir}/rand_vs_ACR_avg_freq.tsv", f"Average Local Chaining Score Lengths {i}-{i+100}", "Average Score")
        create_histograms(f"{dir}/ACR_vs_ACR_max_freq.tsv", f"{dir}/rand_vs_ACR_max_freq.tsv", f"Highest Local Chaining Score Lengths {i}-{i+100}", "Max Score")
        print(f"Finished {i}-{i+100} in {time.time() - start_time} seconds", flush=True)

def driver_all_top():
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    dicts = create_dicts("/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/local/alignment_top5_local.tsv", ref_set, 3130)
    print(f"Non reference ACRs: {len(dicts[0])}")
    print(f"Random regions: {len(dicts[1])}")
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/local/local_freq"
    for i in range(1, 6):
        output_score_freq(dicts[0], dicts[1], base_dir, lambda lst: sorted(lst)[-(i)], f"{i}-highest")
        
    # for i in range(1, 6):
    #     create_histograms(f"{base_dir}/ACR_vs_ACR_{i}-highest_freq.tsv", f"{base_dir}/rand_vs_ACR_{i}-highest_freq.tsv", f"{i}-highest Local", 
    #                  f"{i}-highest Chaining Score", 100)
        
    # output_score_freq(dicts[0], dicts[1], base_dir, top_avg, f"top-5-avg")
    average_graph_top("/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/local/local_freq", 5)
    # create_histograms(f"{base_dir}/ACR_vs_ACR_top-10-avg_freq.tsv", f"{base_dir}/rand_vs_ACR_top-5-avg_freq.tsv", f"Top 10 Averaged Chaining Score Frequencies Global", 
    #                  "Average Over Top 10 Chaining Scores", 100)


if __name__ == "__main__":
    driver_all_top()