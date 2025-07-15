from collections import defaultdict
from useful_functions import *
from rand_vs_acr import *

'''
Populates <chain_dict> with pairs of regions and a score for that pair. Returns
the min and max chain scores.
'''
def create_score_dict(dict, filepath):
    max_chain = None
    min_chain = None
    #Put the chaining file into a dictionary: {(reg1, reg2): chain_score, ...}
    with open(filepath) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            chain_score = float(line_arr[2])
            reg1 = min(line_arr[0], line_arr[1])
            reg2 = max(line_arr[0], line_arr[1])
            dict[(reg1, reg2)] = chain_score
            if max_chain == None or chain_score > max_chain:
                max_chain = chain_score
            if min_chain == None or chain_score < min_chain:
                min_chain = chain_score
    return (min_chain, max_chain)

'''
Creates a new score file in the same format as the chaining and alignment files.
Normalizes each score with respect to the max chain/alignment score.
Flips the normalized alignment score so that better alignment = lower score.
Weights the alignment score at <align_frac> and the chain score at 1 - <align_frac>.
Outputs the two regions and score to a line of the file.
'''
def align_chain_scores(align_file, chain_file, output_file, align_frac):
    chain_dict = {}
    min_chain, max_chain = create_score_dict(chain_dict, chain_file)
    
    #Put the alignment score file into a list of tuples (region1, region2, score) and
    #find the max alignment score
    align_lst = []
    max_align = 0
    min_align = 10000
    with open(align_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            align_score = float(line_arr[2])
            reg1 = min(line_arr[0], line_arr[1])
            reg2 = max(line_arr[0], line_arr[1])
            align_lst.append((reg1, reg2, align_score))
            if  align_score > max_align: 
                max_align = align_score
            if align_score < min_align:
                min_align = align_score
    
    chain_frac = 1 - align_frac
    with open(output_file, "w") as file:
        for reg1, reg2, align_score in align_lst:
            chain_score = chain_dict.get((reg1, reg2), 0)
            try:
                align_score_norm = align_frac * (1 - ((align_score - min_align) / (max_align - min_align)))
                chain_score_norm = chain_frac * ((chain_score - min_chain)/ (max_chain - min_chain))
            except ZeroDivisionError:
                print("Scores are all the same value")
                raise
            file.write(f"{reg1}\t{reg2}\t{align_score_norm + chain_score_norm}\n")

'''
Creates a dictionary for ACR test regions and for random regions in the form:
{region: [chain_score, chain_score, ...], ...}
Excludes any chain scores where the normalized alignment score is above <thresh_frac>
'''
def exclude_high_align(align_file, chain_file, ref_set, thresh_frac):
    align_dict = {}
    min_align, max_align = create_score_dict(align_dict, align_file)
    thresh = (thresh_frac * (max_align - min_align)) + min_align # Turn normalized value back to regular value
    return create_dicts(chain_file, ref_set, 3130, align_dict=align_dict, thresh=thresh)

def align_chain_driver(type, type_short, align_frac):
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90"
    align_file = f"{base_dir}/alignment/{type}/alignment_90-10_{type_short}.tsv"
    chain_file = f"{base_dir}/Chaining_one_acr_rand_10-90_{type_short}.tsv"
    output_file = f"{base_dir}/chain_and_align/align-{align_frac}_chain_{type_short}_scores.tsv"

    align_chain_scores(align_file, chain_file, output_file, align_frac)

def exclude_high_align_driver(type, type_short, thresh):
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90"
    dicts = exclude_high_align(f"{base_dir}/alignment/{type}/alignment_90-10_{type_short}.tsv", f"{base_dir}/Chaining_one_acr_rand_10-90_{type_short}.tsv", ref_set, thresh)
    for i in range(1, 6):
        lst_op = lambda lst: sorted(lst)[-i]
        output_score_freq(dicts[0], dicts[1], f"{base_dir}/chain_and_align/exclude_high_align", lst_op, f"exclude_{thresh}")

if __name__ == "__main__":
    exclude_high_align_driver("global", "glob", .05)
    exclude_high_align_driver("local", "loc", .05)
    # for frac in [.25, .5, .75]:
    #     align_chain_driver("global", "glob", frac)
    #     align_chain_driver("local", "loc", frac)



            
