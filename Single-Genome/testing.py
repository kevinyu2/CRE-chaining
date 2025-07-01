import random
from rand_vs_acr import *
import os
import sys

'''
Testing rand_vs_acr.py
'''

def generate_ref_file():
    with open("ref_set_test.txt", "w") as file:
        for i in range(1000):
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
            acr_options = [random_ACR_name(), random_ACR_name(), random_ACR_name() + "_rand", random_ACR_name() + "_rand", 
                        ref_set[i], ref_set[i + 1]]
        
            option_1 = acr_options[random.randint(0, 5)]
            option_2 = option_1
            while option_1 == option_2:
                option_2 = acr_options[random.randint(0, 5)]
            
            score = random.randint(5, 25)
            if "rand" in option_1:
                if "Chr" in option_2:
                    score = 30
            elif "Chr" in option_1:
                if "rand" in option_2:
                    score = 30
                else:
                    score = 50
            elif "Chr" in option_2:
                score = 50
            
            file.write(f"{option_1}\t{option_2}\t{score}\t{random.randint(5, 25)}\n")


def random_ACR_name():
    name = ""
    for i in range(10):
        name += chr(random.randint(97, 122))
    return name

def main_rand_v_acr_test():
    ref_set = create_ref_set("ref_set_test.txt")
    generate_chaining_file(500, ref_set)
    dicts = create_dicts_ref_set("ACR_v_rand_test.tsv", ref_set, 300)

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

    print(f"ACR count: {acr_count}, ACR_act count: {acr_act_count}")
    print(f"rand count: {rand_count}, rand_act count: {rand_act_count}")

def test_frequencies():
    ref_set = create_ref_set("ref_set_test.txt")
    generate_chaining_file(500, ref_set)
    dicts = create_dicts_ref_set("ACR_v_rand_test.tsv", ref_set, 300)
    output_all_score_freq(dicts[0], dicts[1], ".")
    output_score_freq(dicts[0], dicts[1], ".", max, "max")
    output_score_freq(dicts[0], dicts[1], ".", lambda lst: sum(lst) / float(len(lst)), "avg")
    # create_plots("ACR_vs_ACR_all_freq.tsv", "rand_vs_ACR_all_freq.tsv", "All Chaining Scores", "Score")
    # create_plots("ACR_vs_ACR_avg_freq.tsv", "rand_vs_ACR_avg_freq.tsv", "Average Chaining Score", "Average Score")
    # create_plots("ACR_vs_ACR_max_freq.tsv", "rand_vs_ACR_max_freq.tsv", "Highest Chaining Score", "Max Score")

test_frequencies()








