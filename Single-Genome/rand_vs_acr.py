import random
from collections import defaultdict
from matplotlib import pyplot as plt
import time

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
Driver function.

Takes in two files:
1. pairwise upstream genes chaining 
2. pairwise upstream + random chaining

Randomly chooses <reference_size> of the ACR regions as 'references'. Returns a dictionary with
every non-reference ACR and a list of the chaining scores with all reference ACRs.
Returns a second dictionary with every random region and the list of chaining scores
with all reference ACRs.
'''
def create_dicts(all_ACR_file, ACR_rand_file, reference_size):
    start_time = time.time()
    #For every non-reference ACR, this dictionary contains a list of scores of that ACR
    #chained with every reference ACR
    ACR_dict = defaultdict(list)
    #For every random sequence, this dictionary contains a list of scores of that random sequence
    #chained with every reference ACR
    rand_dict = defaultdict(list)

    ACR_ref = set() #set of random ACR references

    #populate ACR dict (and choose random ACRs to be references)
    with open(all_ACR_file, "r") as ACR_file:
            print("Before reading lines")
            lines = ACR_file.read().split("\n")
            print(f"Read in ACR lines after {time.time() - start_time} seconds")
            random.shuffle(lines)
            print(f"Shuffled lines after {time.time() - start_time} seconds")
            for line in lines:
                line_arr = line.split("\t")
                if len(line_arr) < 4:
                    continue
                item1 = line_arr[0]
                item2 = line_arr[1]
                score = float(line_arr[2])
                if item1 in ACR_ref:
                    if item2 in ACR_ref: #both are references
                        continue
                    else:
                        ACR_dict[item2].append(score) #item 1 is reference; item 2 is not
                elif item2 in ACR_ref: #item 2 is reference; item1 is not
                    ACR_dict[item1].append(score)
                elif len(ACR_ref) < reference_size: #neither are references yet
                    if item1 != item2: #Don't count chaining with itself
                        ACR_dict[item2].append(score)
                    ACR_ref.add(item1)
    #populate rand dict
    with open(ACR_rand_file) as rand_file:
        print(f"Began rand dict after {time.time() - start_time} seconds")
        for line in rand_file:
            line_arr = line.split("\t")
            if len(line_arr) < 4:
                continue
            item1 = line_arr[0]
            item2 = line_arr[1]
            score = float(line_arr[2])
            if item1[-4:] == "rand": #item 1 is random
                if item2 in ACR_ref: #item 2 is a reference
                    rand_dict[item1].append(score)
            elif item1 in ACR_ref: #item 1 is a reference
                if item2[-4:] == "rand": #item 2 is random
                    rand_dict[item2].append(score)
    print(f"Reference size: {len(ACR_ref)}")
    return (ACR_dict, rand_dict)

'''
Driver function
Takes in a chaining file with all random regions and ACRs pairwise chained. Also takes in the
set of "reference" ACRs, i.e the ACRs which should be chained against the random regions and other
ACRs.
Non_ref_num is the number of non-reference and random regions there should be
'''
def create_dicts_ref_set(chaining_file, reference_set, non_ref_num):
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
    data1 = [0 for i in range(200)]
    ACR_count = 0
    rand_count = 0
    with open(file1, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            index = round(float(line_arr[0]))
            try:
                data1[index] = int(line_arr[1])
            except IndexError:
                ACR_count += 1
            
    data2 = [0 for i in range(200)]
    with open(file2, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            index = round(float(line_arr[0]))
            try:
                data2[index] = int(line_arr[1])
            except IndexError:
                rand_count += 1
    
    max_axis = max(max(data1), max(data2))
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].bar(range(200), data1, color='blue')
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    axes[0].set_xlabel(x_label)

    axes[1].bar(range(200), data2, color='green')
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlabel(x_label)

    plt.suptitle(title)
    plt.tight_layout()
    file_name = title.replace(" ", "_")
    plt.savefig(f"./{file_name}.png")
    print(f"ACR outliers: {ACR_count}")
    print(f"random outliers: {rand_count}")

if __name__ == "__main__":
    ref_set = create_ref_set()
    dicts = create_dicts_ref_set("/home/mwarr/Data/One_Genome/Chaining_one_rand_loc.tsv", ref_set, 15648)
    print(f"Non reference ACRs: {len(dicts[0])}")
    print(f"Random regions: {len(dicts[1])}")
    output_all_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome", max, "max")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/One_Genome", lambda lst: sum(lst) / float(len(lst)), "avg")
    create_plots("/home/mwarr/Data/One_Genome/ACR_vs_ACR_all_freq.tsv", "/home/mwarr/Data/One_Genome/rand_vs_ACR_all_freq.tsv", "All Chaining Scores", "Score")
    create_plots("/home/mwarr/Data/One_Genome/ACR_vs_ACR_avg_freq.tsv", "/home/mwarr/Data/One_Genome/rand_vs_ACR_avg_freq.tsv", "Average Chaining Score", "Average Score")
    create_plots("/home/mwarr/Data/One_Genome/ACR_vs_ACR_max_freq.tsv", "/home/mwarr/Data/One_Genome/rand_vs_ACR_max_freq.tsv", "Highest Chaining Score", "Max Score")