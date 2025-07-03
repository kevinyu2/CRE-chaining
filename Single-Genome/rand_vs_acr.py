import random
from collections import defaultdict
from matplotlib import pyplot as plt
from pathlib import Path
import os
import time

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
            if len(line_arr) < 4:
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
def create_plots(file1, file2, title, x_label):
    data1 = [0 for i in range(75)]
    ACR_count = 0
    rand_count = 0
    with open(file1, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            index = round(float(line_arr[0]))
            try:
                data1[index] += int(line_arr[1])
            except IndexError:
                ACR_count += 1
            
    data2 = [0 for i in range(75)]
    with open(file2, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            index = round(float(line_arr[0]))
            try:
                data2[index] += int(line_arr[1])
            except IndexError:
                rand_count += 1
    assert(sum(data1) == sum(data2))
    
    #loop through data arrays backwards to find max index which is not 0
    for i in range(len(data1) - 1, -1, -1):
        if data1[i] != 0 or data2[i] != 0:
            max_index = i + 1
            data1 = data1[: max_index]
            data2 = data2[: max_index]
            break

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)

    axes[0].bar(range(max_index), data1, color='blue')
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    axes[0].set_xlabel(x_label)

    axes[1].bar(range(max_index), data2, color='green')
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlabel(x_label)

    #plt.axis('equal')
    plt.suptitle(title)
    plt.tight_layout()
    file_name = title.replace(" ", "_")
    plt.savefig(f"/home/mwarr/{file_name}.png")
    print(f"ACR outliers: {ACR_count}")
    print(f"random outliers: {rand_count}")

def get_averages(input_dir):
    options = ["all", "avg", "max"]
    region_type = ["ACR", "rand"]
    for opt in options:
        print(f"-------{opt}--------")
        for type in region_type:
            for file_path in Path(input_dir).glob(f"{type}*{opt}_freq.tsv"):
                sum = 0
                num_regions = 0
                with open(file_path) as file:
                    for line in file:
                        line_arr = line.split("\t")
                        if len(line_arr) < 2:
                            break
                        sum += float(line_arr[0]) * int(line_arr[1])
                        num_regions += int(line_arr[1])
                
                print(f"Average for {type}: {sum / num_regions}")

def driver_all():
    ref_set = create_ref_set("/home/mwarr/Data/seta.txt")
    dicts = create_dicts("/home/mwarr/Data/One_Genome/Chaining_one_acr_rand_loc.tsv", ref_set, 15648)
    print(f"Non reference ACRs: {len(dicts[0])}")
    print(f"Random regions: {len(dicts[1])}")
    output_all_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome/local_freq_all")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome/local_freq_all", max, "max")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome/local_freq_all", lambda lst: sum(lst) / float(len(lst)), "avg")
    create_plots("/home/mwarr/Data/One_Genome/local_freq_all/ACR_vs_ACR_all_freq.tsv", "/home/mwarr/Data/One_Genome/local_freq_all/rand_vs_ACR_all_freq.tsv", "All Chaining Scores", "Score")
    create_plots("/home/mwarr/Data/One_Genome/local_freq_all/ACR_vs_ACR_avg_freq.tsv", "/home/mwarr/Data/One_Genome/local_freq_all/rand_vs_ACR_avg_freq.tsv", "Average Chaining Score", "Average Score")
    create_plots("/home/mwarr/Data/One_Genome/local_freq_all/ACR_vs_ACR_max_freq.tsv", "/home/mwarr/Data/One_Genome/local_freq_all/rand_vs_ACR_max_freq.tsv", "Highest Chaining Score", "Max Score")

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
    driver_filter()