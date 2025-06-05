# Should be tsv of motif1, motif2, and score. Motif1 and motif2 are sorted in alphabetical order

def read_motif_pairs(min_score):
    SCORES_FILE_NAME = "scores_norm.tsv"

    pairwise_score_dict = {}

    with open(SCORES_FILE_NAME) as pair_scores:
        for line in pair_scores :
            line_arr = line.rstrip().split('\t')
            if float(line_arr[2]) > min_score:
                pairwise_score_dict[(line_arr[0], line_arr[1])] = float(line_arr[2])
    return pairwise_score_dict

# Genome 1 is 
def find_anchors(pairwise_score_dict, genome_1, genome_2, acr_1, acr_2):
    
        
