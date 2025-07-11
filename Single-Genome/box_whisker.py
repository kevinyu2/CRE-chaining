from matplotlib import pyplot as plt
import random

'''
Takes in a frequency file and creates a list of all the values
'''
def format_data(input_file, ignore_zeros=True, max_value=None):
    data = []
    with open(input_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            value = float(line_arr[0])
            freq = int(line_arr[1])
            if ignore_zeros and value == 0:
                continue
            if max_value != None and value > max_value:
                continue
            extension = [value for _ in range(freq)]
            data.extend(extension)
    return data

'''
Generates the box plot and saves it to /home/mwarr/{filename}
'''
def create_box_plot(data, xlabels, title, xtitle, ytitle, filename, colors=None):
    plt.figure(figsize=(15, 8))
    bplot = plt.boxplot(data, labels=xlabels, sym="o", patch_artist=True, 
                        whiskerprops=dict(linewidth=2), medianprops=dict(linewidth=2, color="red"),
                        boxprops=dict(linewidth=2), capprops=dict(linewidth=2))
    color_rand = (random.random(), random.random(), random.random())
    for ind, box in enumerate(bplot["boxes"]):
        if colors == None:
            box.set_facecolor(color_rand)
            if ind % 2 == 1: #regenerate every 2 boxes
                color_rand = (random.random(), random.random(), random.random())
        else:
            box.set_facecolor(colors[ind // 2])
    plt.title(title, fontweight='bold', pad=15)
    plt.xlabel(xtitle, fontweight='bold', labelpad=15)
    plt.ylabel(ytitle, fontweight='bold', labelpad=15)
    plt.savefig(f"/home/mwarr/{filename}")

'''
Reads in the frequency files from <base_dir> which contain frequencies of some feature of
top 5 scores and generates a box plot.
'''
def top_5_all(base_dir, op, title, xtitle, ytitle, filename, max_value=None):
    data = []
    labels = []
    for i in range(1, 6):
        acr_data = format_data(f"{base_dir}/ACR_vs_ACR_{op}_{i}-highest_freq.tsv", max_value)
        rand_data = format_data(f"{base_dir}/rand_vs_ACR_{op}_{i}-highest_freq.tsv", max_value)
        data.append(acr_data)
        labels.append(f"{i}-highest, \nACR")
        data.append(rand_data)
        labels.append(f"{i}-highest, \nrandom")
    colors = ["yellow", "blue", "green", "pink", "purple"]
    create_box_plot(data, labels, title, xtitle, ytitle, filename, colors)

'''
Reads in the frequency files from <base_dir> which contain the frequencies of some
feature given a chain length.
'''
def fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename, max_value=None):
    data = []
    labels = []
    for i in range(2, 20):
        acr_data = format_data(f"{base_dir}/ACR_vs_ACR_{op}_{i}_freq.tsv", max_value)
        rand_data = format_data(f"{base_dir}/rand_vs_ACR_{op}_{i}_freq.tsv", max_value)
        data.append(acr_data)
        labels.append(f"{i}, \nACR")
        data.append(rand_data)
        labels.append(f"{i}, \nrand")
    create_box_plot(data, labels, title, xtitle, ytitle, filename)

def driver_counts():
    base_dir = "/home/mwarr/Data/One_Genome/other_features_glob/counts"
    op = "count"
    title = "Distribution of the Number of Top Scores, Global"
    xtitle = "Score Category and Region Type"
    ytitle = "Number of Top Scores"
    filename = "box_counts_glob.png"
    top_5_all(base_dir, op, title, xtitle, ytitle, filename)

def driver_scores():
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90/loc_freq"
    op = "score"
    title = "Distribution of Top Scores, Local"
    xtitle = "Score Category and Region Type"
    ytitle = "Score"
    filename = "box_scores_loc.png"
    top_5_all(base_dir, op, title, xtitle, ytitle, filename)

def driver_anchor_num():
    base_dir = "/home/mwarr/Data/One_Genome/other_features_glob/anchor_num"
    op = "anchor_num"
    title = "Distribution of Number of Anchors for a 5th-Highest Chain Length"
    xtitle = "5th-Highest Chain Length"
    ytitle = "Number of Anchors"
    filename = "box_anchor_num_glob.png"
    fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename)

def driver_ref_len():
    base_dir = "/home/mwarr/Data/One_Genome/other_features_glob/ref_len"
    op = "ref_len"
    title = "Distribution of the Reference Region Length of the 5th-Highest Chain Score"
    xtitle = "5th-Highest Chain Score"
    ytitle = "Reference Region Length"
    filename = "box_ref_len_glob.png"
    fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename)

def driver_combine_score():
    data = []
    base_dir = "/home/mwarr/Data/One_Genome/other_features_glob/combine_score"
    data.append(format_data(f"{base_dir}/ACR_vs_ACR_combine_freq.tsv"))
    data.append(format_data(f"{base_dir}/rand_vs_ACR_combine_freq.tsv"))
    labels = ["ACR", "Random"]
    #create_box_plot(data, labels, "Distribution of Combined Score, Global", "Type of Region", "Score", "box_combined_glob.png")
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)
    axes[0].hist(data[0], bins=100)
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    
    axes[1].hist(data[1], bins=100)
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")

    plt.suptitle("Combined Score Frequencies, Global")
    plt.tight_layout()
    plt.savefig("/home/mwarr/combined_freq_glob.png")

if __name__ == "__main__":
    driver_scores()
