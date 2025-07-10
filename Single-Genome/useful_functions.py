from matplotlib import pyplot as plt #type:ignore

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
Creates a dictionary where the keys are the chromosome identifiers
for the genome and the values are the sequences for the given chromosome
'''
def create_genome_dict(genome_file):
    #create genome dict {chr: seq, chr: seq, ...}
    genome_dict = {}
    with open(genome_file) as genome:
        chromosome = genome.readline().strip()
        while chromosome:
            chromosome = chromosome[1:]
            seq = ""
            line = genome.readline().strip()
            while line != "" and line[0] != ">":
                seq += line
                line = genome.readline().strip()
            genome_dict[chromosome] = seq
            chromosome = line
    return genome_dict


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

def output_freq_to_file(filename, freq_dict):
    with open(filename, "w") as file:
        for val, freq in freq_dict.items():
            file.write(f"{val}\t{freq}\n")

'''
Reads frequencies from files and creates two side-by-side plots.
<file1> should contain ACR chaining with ACR frequency data.
<file2> should contain ACR chaining with random frequency data.
Plots will only include values between <min_index> and <max_index>.

'''
def create_histograms(file1, file2, title, x_label, max_index, equal_freq=False, min_index=0):
    data1 = [0 for i in range((max_index - min_index) + 1)]
    acr_outliers = 0
    rand_outliers = 0

    with open(file1, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data1[score-min_index] += freq
            except IndexError:
                acr_outliers += freq

    data2 = [0 for i in range((max_index - min_index) + 1)]
    with open(file2, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data2[score-min_index] += freq
            except IndexError:
                rand_outliers += freq

    #loop through data arrays backwards to find max index which is not 0
    for i in range(len(data1) - 1, -1, -1):
        if data1[i] != 0 or data2[i] != 0:
            max_index = i
            data1 = data1[: max_index + 1]
            data2 = data2[: max_index + 1]
            break

    if equal_freq:
        try:
            assert(sum(data1) + acr_outliers == sum(data2) + rand_outliers)
        except:
            print(sum(data1) + acr_outliers)
            print(sum(data2) + rand_outliers)
            raise

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)

    axes[0].bar(range(min_index, max_index+1), data1, color='blue')
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    axes[0].set_xlabel(x_label)

    axes[1].bar(range(min_index, max_index+1), data2, color='green')
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlabel(x_label)

    #plt.axis('equal')
    plt.suptitle(title)
    plt.tight_layout()
    file_name = title.replace(" ", "_")
    print("Saving figure")
    plt.savefig(f"/home/mwarr/{file_name}.png")
    print(f"ACR outliers {acr_outliers}")
    print(f"Random outliers {rand_outliers}")


if __name__ == "__main__":
    base_dir = "/home/mwarr/Data/One_Genome/other_features/counts"
    create_histograms(f"{base_dir}/ACR_vs_ACR_count_5-highest_freq.tsv", f"{base_dir}/rand_vs_ACR_count_5-highest_freq.tsv",
                      "Frequency of Count of 5th Highest Chain Scores", "Count of 5th Highest Chain Scores", 200, True)