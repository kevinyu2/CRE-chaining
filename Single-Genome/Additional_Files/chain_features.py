from matplotlib import pyplot as plt
import sys
sys.path.append("../")
from useful_functions import *
from rand_vs_acr import *
from collections import defaultdict

'''
Takes in a chaining file with all random regions and ACRs pairwise chained. Also takes in the
set of "reference" ACRs, i.e the ACRs which should be chained against the random regions and other
ACRs.
Non_ref_num is the number of non-reference and random regions there should be.

Returns two dictionaries with the number of anchors and reference ACR lengths
(for exact contents and structure, see in-line comment)

Note: Unlike the similar function in rand_vs_acr.py, this function does not account for chain
lengths of 0
'''
def create_dicts_mult_feat(chaining_file, reference_set):
    #For every non-reference ACR, this dictionary contains a list of anchor numbers and lengths for that ACR
    #chained with every reference ACR
    #{non-ref-acr: [(score, anchor_num, ref-acr-len), (.., ..)...]}
    ACR_dict = defaultdict(list)
    #For every random sequence, this dictionary contains a list of anchor numbers and lengths of that random sequence
    #chained with every reference ACR
    #{random_reg: [(score, anchor_num, ref-acr-len), (.., ..)...]}
    rand_dict = defaultdict(list)

    with open(chaining_file, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 4:
                continue
            item1 = line_arr[0]
            item2 = line_arr[1]
            score = float(line_arr[2])
            anchor_num = int(line_arr[3])
            item1_len = get_length(item1)
            item2_len = get_length(item2)

            if item1[-4:] == "rand":
                if item2 in reference_set:
                    rand_dict[item1].append((score, anchor_num, item2_len))
            elif item1 in reference_set:
                if item2[-4:] == "rand":
                    rand_dict[item2].append((score, anchor_num, item1_len))
                elif item2 not in reference_set:
                    ACR_dict[item2].append((score, anchor_num, item1_len))
            else:
                if item2 in reference_set:
                    ACR_dict[item1].append((score, anchor_num, item2_len))
    
    return(ACR_dict, rand_dict)

'''
Wrapper function to create function for output_score_freq.
'''
def list_op_anchor(chain_score):
    '''
    Function to be passed into output_score_freq. Returns a list of all the anchor numbers
    of the regions with a score of <chain_score> if and only if the 5th highest score is <chain_score>. 
    Otherwise returns 0.
    '''
    def anchor_num(lst):
        score_lst = [item[0] for item in lst]
        if sorted(score_lst)[-5] == chain_score: #5th highest score
            anchor_lst = []
            for item in lst:
                if item[0] == chain_score:
                    anchor_lst.append(item[1])
            return anchor_lst
        return 0
    return anchor_num

'''
Function to be passed into output_score_freq. Returns a list of all the reference region lengths
of the regions with a score of 10 if and only if the 5th highest score is 10. Otherwise returns 0.
'''
def list_op_len_ref(chain_score):
    def len_ref(lst):
        score_lst = [item[0] for item in lst]
        if sorted(score_lst)[-5] == chain_score: #5th highest score is 10
            len_lst = []
            for item in lst:
                if item[0] == chain_score:
                    len_lst.append(item[2])
            return len_lst
        return 0
    return len_ref

def non_ref_len_freq(ACR_dict, rand_dict, output_dir, op_name):
    ACR = {}
    rand = {}
    for key in ACR_dict.keys():
        length = get_length(key)
        ACR[length] = ACR.get(length, 0) + 1
    for key in rand_dict.keys():
        length = get_length(key)
        rand[length] = rand.get(length, 0) + 1
    output_freq_to_file(f"{output_dir}/ACR_vs_ACR_{op_name}", ACR)
    output_freq_to_file(f"{output_dir}/rand_vs_ACR_{op_name}", rand)

def average_graph_constant_chain(input_dir_base, op_name, ylabel, title):
    acr_data = []
    rand_data = []
    labels = [i for i in range(2, 20)]
    for i in range(2, 20):
        acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_{op_name}_{i}_freq.tsv", True))
        rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_{op_name}_{i}_freq.tsv", True))
    width = .4
    x = np.arange(len(labels))
    plt.figure(figsize=(15, 6))
    plt.bar(x - width / 2, acr_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.ylim(870, 890)
    plt.title(title)
    plt.xlabel("5th-Highest Chain Length")
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(f"/home/mwarr/chainging_exp2_{op_name}_glob_original.png")

def driver():
    base = "/home/mwarr/Data/One_Genome/other_features"
    # ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    # dicts = create_dicts_mult_feat("/home/mwarr/Data/One_Genome/experiment2_10-90/Chaining_one_acr_rand_10-90_glob.tsv", ref_set)
    # for i in range(2, 20):
    #     output_score_freq(dicts[0], dicts[1], base, list_op_anchor(i), f"anchor_num_{i}")
    #     output_score_freq(dicts[0], dicts[1], base, list_op_len_ref(i), f"ref_len_{i}") 
    #     non_ref_len_freq(dicts[0], dicts[1], base, f"non-ref-len_{i}")   
    # average_graph_constant_chain(f"{base}/anchor_num", "anchor_num", "Average Number of Anchors", "Average Number of Anchors for a Given 5th-Highest Chain Length, Global")
    # average_graph_constant_chain(f"{base}/ref_len", "ref_len", "Average Reference ACR Length", "Average Reference ACR Length for a Given 5th-Highest Chain Length, Global")
    average_graph_constant_chain(f"{base}/non-ref_len", "non-ref-len", "Average Test Region Length", "Average Test Region Length for a Given 5th-Highest Chain Length, Global")

if __name__ == "__main__":
    driver()