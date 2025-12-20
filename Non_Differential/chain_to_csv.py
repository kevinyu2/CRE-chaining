import csv
import time

def load_region_list(path):
    """Load a list of region names, one per line."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]

def load_scores(path):
    """
    Load the main score table.
    File format (tab-separated):
        region_1   region_2   score   other_number
    Return a dictionary:
        scores[(region1, region2)] = score (float)
        scores[(region2, region1)] = score (float)   # symmetric lookup
    """
    scores = {}
    with open(path) as f:
        line_count = 0
        start_time = time.time()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            r1, r2, score_str, _ = parts
            score = float(score_str)

            # Record both directions for convenience
            scores[(r1, r2)] = score
            scores[(r2, r1)] = score

            line_count += 1
    
    return scores

def write_matrix(output_path, row_regions, col_regions, scores):
    """
    Make a matrix (row_regions x col_regions) of scores and save as CSV.
    Missing entries become 0.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow(["non_rep_region"] + col_regions + ["acr_label"])
        
        # Rows
        for r in row_regions:
            row = [r]
            for c in col_regions:
                row.append(scores.get((r, c), 0))
            if r[-3:] == "neg":
                row.append(0)
            elif r[-3:] == "pos":
                row.append(1)
            else:
                raise Exception("No ACR label")
            writer.writerow(row)

def make_datasets(score_file, train_file, test_file, validate_file, reference_file, out_dir):
    # Load files
    print("Beginning program", flush=True)
    scores = load_scores(score_file)
    print("Scores loaded", flush=True)
    train_regions = load_region_list(train_file)
    test_regions = load_region_list(test_file)
    validate_regions = load_region_list(validate_file)
    reference_regions = load_region_list(reference_file)
    print("Regions loaded", flush=True)

    # Produce outputs
    print("outputing train", flush=True)
    write_matrix(f"{out_dir}/train.csv", train_regions, reference_regions, scores)
    print("outputing test", flush=True)
    write_matrix(f"{out_dir}/test.csv", test_regions, reference_regions, scores)
    print("outputing val", flush=True)
    write_matrix(f"{out_dir}/validate.csv", validate_regions, reference_regions, scores)
    print("done", flush=True)

input_dir = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50"
make_datasets(f"{input_dir}/all_chain.tsv", f"{input_dir}/train_regions.txt", f"{input_dir}/test_regions.txt",
              f"{input_dir}/validate_regions.txt", f"{input_dir}/ref_set.txt", input_dir)