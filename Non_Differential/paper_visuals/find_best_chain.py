
ref_set_path = "/home/mwarr/Data/paper_visualizing/ref_set.txt"
acr_pred_path = "/home/mwarr/Data/paper_visualizing/random_true_pred_histgbc.txt"
chain_file = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_50-50/all_chain.tsv"
out_path = "/home/mwarr/Data/paper_visualizing/random_weighted.txt"

# Put ACR predictions and ref ACRs into sets

ref_set = set()

with open(ref_set_path) as in_file:
    for line in in_file:
        ref_set.add(line.rstrip())

acr_pred = set()

with open(acr_pred_path) as in_file:
    for line in in_file:
        acr_pred.add(line.rstrip())

print("Finished loading predictions and ref set", flush=True)

#find max chain score

chain_scores = []
with open(chain_file) as in_file:
    for line in in_file:
        line_arr = line.split("\t")
        if len(line_arr) < 3:
            continue

        curr_chain = 0
        curr_ref = ""
        curr_pred = ""
        if line_arr[0].rstrip() in ref_set and line_arr[1].rstrip() in acr_pred:
            curr_ref = line_arr[0].rstrip()
            curr_pred = line_arr[1].rstrip()
            curr_chain = float(line_arr[2].rstrip())
        elif line_arr[1].rstrip() in ref_set and line_arr[0].rstrip() in acr_pred:
            curr_ref = line_arr[1].rstrip()
            curr_pred = line_arr[0].rstrip()
            curr_chain = float(line_arr[2].rstrip())

        if curr_chain > 0:
            chain_scores.append((curr_chain, curr_ref, curr_pred))

chain_scores = sorted(chain_scores, reverse=True)

# output
with open(out_path, "w") as out:
    for i in range(20):
        out.write(f"{chain_scores[i][0]}\treference: {chain_scores[i][1]}\tpredicted: {chain_scores[i][2]}\n")
