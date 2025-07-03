import random
from rand_vs_acr import *
import os
import sys

'''
Testing rand_vs_acr.py
'''

def generate_ref_file(size):
    with open("ref_set_test.txt", "w") as file:
        for i in range(size):
            start = random.randint(1000, 3000)
            size = random.randint(500, 1500)
            file.write(f"Chr{random.randint(1, 6)}_{start}to{start+size}\n")

def generate_chaining_file(length, ref_set):
    ref_set = list(ref_set)
    try:
        os.remove("ACR_v_rand_test.tsv")
    except:
        print("Probably file isn't there")
    with open("ACR_v_rand_test.tsv", "w") as file: 
        for i in range(length):
            try:
                acr_options = [random_ACR_name(), random_ACR_name(), random_ACR_name() + "_rand", random_ACR_name() + "_rand", 
                            ref_set[i], ref_set[i + 1]]
            except:
                acr_options = [random_ACR_name(), random_ACR_name(), random_ACR_name() + "_rand", random_ACR_name() + "_rand", 
                            ref_set[i], ref_set[i - 1]]

            option_1 = acr_options[random.randint(0, 5)]
            option_2 = option_1
            while option_1 == option_2:
                option_2 = acr_options[random.randint(0, 5)]
            
            file.write(f"{option_1}\t{option_2}\t{random.randint(5, 25)}\t{random.randint(5, 25)}\n")

def random_ACR_name():
    start = random.randint(1000, 3000)
    size = random.randint(1000, 1500)
    name = f"Chr{random.randint(1, 5)}_{start}to{start+size}"
    return name

def rand_v_acr_dict_all():
    generate_ref_file(1000)
    ref_set = create_ref_set("ref_set_test.txt")
    generate_chaining_file(500, ref_set)
    dicts = create_dicts("ACR_v_rand_test.tsv", ref_set, 300)

    acr_count = 0
    rand_count = 0
    for key in dicts[0].keys():
        if "unknown" not in key:
            acr_count += 1

    for key in dicts[1].keys():
        if "unknown" not in key:
            rand_count += 1

    acr_act_count = 0
    rand_act_count = 0
    with open("ACR_v_rand_test.tsv") as file:
        for line in file:
            line_arr = line.split("\t")
            if "rand" in line_arr[0]:
                if "Chr" in line_arr[1]:
                    rand_act_count += 1
            elif "Chr" in line_arr[0]:
                if "rand" in line_arr[1]:
                    rand_act_count += 1
                elif "Chr" not in line_arr[1]:
                    acr_act_count += 1
            else:
                if "Chr" in line_arr[1]:
                    acr_act_count +=1

    assert(acr_count == acr_act_count)
    assert(rand_count == rand_act_count)

def rand_v_acr_dict_filter():
    generate_ref_file(1000)
    ref_set = create_ref_set_filter("ref_set_test.txt", 1000, 1500)
    generate_chaining_file(100, ref_set)
    dicts = create_dicts("ACR_v_rand_test.tsv", ref_set, 300, 1000, 1500)
    acr_count = 0
    rand_count = 0
    for key in dicts[0].keys():
        if "unknown" not in key:
            acr_count += 1

    for key in dicts[1].keys():
        if "unknown" not in key:
            rand_count += 1

    acr_act_count = 0
    rand_act_count = 0
    with open("ACR_v_rand_test.tsv") as file:
        for line in file:
            line_arr = line.split("\t")
            line_arr[0] = line_arr[0].strip()
            line_arr[1] = line_arr[0].strip()
            if "rand" in line_arr[0]:
                if line_arr[1] in ref_set:
                    rand_act_count += 1
            elif line_arr[0] in ref_set:
                if "rand" in line_arr[1]:
                    rand_act_count += 1
                elif line_arr[1] not in ref_set:
                    acr_act_count += 1
            else:
                if line_arr[1] in ref_set:
                    acr_act_count +=1
    print(acr_count)
    print(acr_act_count)
    print(rand_count)
    print(rand_act_count)
    assert(acr_count == acr_act_count)
    assert(rand_count == rand_act_count)
    
def get_total_freq(input_file):
    count = 0
    with open(input_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                break
            count += int(line_arr[1])
    return count

def test_frequencies():
    ref_set = create_ref_set("ref_set_test.txt")
    generate_chaining_file(500, ref_set)
    dicts = create_dicts("ACR_v_rand_test.tsv", ref_set, 300)
    output_all_score_freq(dicts[0], dicts[1], ".")
    output_score_freq(dicts[0], dicts[1], ".", max, "max")
    output_score_freq(dicts[0], dicts[1], ".", lambda lst: sum(lst) / float(len(lst)), "avg")
    # create_plots("ACR_vs_ACR_all_freq.tsv", "rand_vs_ACR_all_freq.tsv", "All Chaining Scores", "Score")
    # create_plots("ACR_vs_ACR_avg_freq.tsv", "rand_vs_ACR_avg_freq.tsv", "Average Chaining Score", "Average Score")
    # create_plots("ACR_vs_ACR_max_freq.tsv", "rand_vs_ACR_max_freq.tsv", "Highest Chaining Score", "Max Score")

def ref_set_filter_test():
    size = 1000
    generate_ref_file(size)
    ref_set = create_ref_set_filter("ref_set_test.txt", 501, 501)
    print(ref_set)
    print(len(ref_set))

rand_v_acr_dict_filter()











