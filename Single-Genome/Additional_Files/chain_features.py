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
def create_dicts(chaining_file, reference_set):
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

def list_op_anchor(lst):
    score_lst = [item[0] for item in lst]
    if max(score_lst) == 10:
        return lst[score_lst.index(10)][1]

def list_op_len_ref(lst):
    score_lst = [item[0] for item in lst]
    if max(score_lst) == 10:
        return lst[score_lst.index(10)][2]


if __name__ == "__main__":
    base = "/home/mwarr/Data/One_Genome/other_features"
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/setb_10.txt")
    dicts = create_dicts("/home/mwarr/Data/One_Genome/experiment2_10-90/Chaining_one_acr_rand_10-90_glob.tsv", ref_set)
    output_score_freq(dicts[0], dicts[1], base, list_op_anchor, "anchor_num")
    output_score_freq(dicts[0], dicts[1], base, list_op_len_ref, "ref_len")
    create_histograms(f"{base}/ACR_vs_ACR_anchor_num_freq.tsv", f"{base}/rand_vs_ACR_anchor_num_freq.tsv", "Anchor number", "", 500)
    create_histograms(f"{base}/ACR_vs_ACR_ref_len_freq.tsv", f"{base}/rand_vs_ACR_ref_len_freq.tsv", "Region Length", "", 500)
