from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier #type: ignore
from sklearn.metrics import accuracy_score, f1_score #type: ignore
from sklearn.neighbors import KNeighborsClassifier #type:ignore
from sklearn.svm import LinearSVC #type:ignore
from sklearn.gaussian_process import GaussianProcessClassifier #type:ignore
from sklearn.tree import DecisionTreeClassifier #type:ignore
from sklearn.neural_network import MLPClassifier #type:ignore
from sklearn.naive_bayes import GaussianNB #type:ignore
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis #type:ignore
from matplotlib import pyplot as plt #type:ignore


import random
from collections import defaultdict

'''
<rand_acr_list> should be the list of all the genomic regions, with each name on its own line.
Random regions should have _rand at the end of the name. It will be assumed that all other regions
are ACRs.

Randomly splits the regions into TEST and TRAIN sets.

Outputs a tuple. The first item is a set with the names in the training set. The second item is a set
with the names in the test set.
'''
def get_region_sets(rand_acr_list):
    region_list = []
    with open(rand_acr_list, "r") as file:
        for line in file:
            region_list.append(line.strip())
    random.shuffle(region_list)
    train_set = set(region_list[ : int(len(region_list) / 2)])
    test_set = set(region_list[int(len(region_list) / 2) : ])
    return train_set, test_set

'''
Takes in a file with all the names of the reference ACRs on their own line.
Outputs a list of the names of the reference ACRs and a dictionary mapping the names
to their index in the list.
'''
def rep_list(rep_list_file):
    rep_list = []
    index_dict = {}
    with open(rep_list_file) as file:
        line_count = 0
        for line in file:
            rep_list.append(line.strip())
            index_dict[line.strip()] = line_count
            line_count += 1
    return rep_list, index_dict

'''
<score_file> should be a file with the train and set regions chained with the representative
regions. A pairwise chaining file should not be passed in here.

<train_set> and <test_set> should be sets containing the names of the respective regions

<rep_list> should be a list of the names of the representative regions
<reg_dict> should be a dictinoary mapping the names of the representative regions to their index
in rep_list
Outputs a tuple. First item is a dictionary mapping every training region to a list of all chaining scores
with representative regions. Second item is a dictionary mapping every test region to a list of all chaining scores
with representative regions.
'''
def get_chain_dicts(score_file, train_set, test_set, rep_list, rep_dict):
    test_dict = defaultdict(lambda: [0 for _ in range(len(rep_list))])
    train_dict = defaultdict(lambda:[0 for _ in range(len(rep_list))])
    with open(score_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 3:
                continue
            reg1 = line_arr[0].strip()
            reg2 = line_arr[1].strip()
            score = float(line_arr[2].strip())

            if reg1 in train_set:
                index = rep_dict[reg2]
                train_dict[reg1][index] = score
            elif reg1 in test_set:
                index = rep_dict[reg2]
                test_dict[reg1][index] = score
            else: #reg1 is a representative region
                index = rep_dict[reg1]
                if reg2 in train_set:
                    train_dict[reg2][index] = score
                if reg2 in test_set:
                    test_dict[reg2][index] = score
    return train_dict, test_dict

'''
<train_dict> should be a dictionary with all of the training regions mapped to a list of scores
from chaining with the representative regions
<test_dict> should be a dictionary with all of the test regions mapped to a list of scores
from chaining with the representative regions
<settings> should be a dictionary with arguments for RandomForestClassifier
Outputs a tuple. The first value is a list with the predicted classes and the second value is a list
of the actual classes. 0 indicates a random region and 1 indicates an ACR.

A classifier may be passed in. This will default to RandomForestClassifier. If <classifier> is set
to 'control', random predictions will be made.
'''          
def train_predict(train_dict, test_dict, fraction=.5, classifier=RandomForestClassifier()):
    feature_train = []
    labels_train = []
    for region, feature_list in train_dict.items():
        feature_train.append(feature_list)
        if "rand" in region:
            labels_train.append(0)
        else:
            labels_train.append(1)
    
    feature_test = []
    labels_test = []
    threshold = fraction * len(test_dict)
    acr_count = 0
    to_delete = []
    
    for region, feature_list in test_dict.items():
        if "rand" in region:
            labels_test.append(0)
            feature_test.append(feature_list)
        else:
            if acr_count > threshold: #we have included enough ACRs
                to_delete.append(region) #save for deletion
            else:
                labels_test.append(1)
                feature_test.append(feature_list)
                acr_count += 1
    if classifier == "control":
        return [random.randint(0, 1) for _ in range(len(labels_test))], labels_test
    else:
        classifier.fit(feature_train, labels_train)
        return classifier.predict(feature_test), labels_test

'''
Calculates and returns a dictionary of true and false positives,
as well as true and false negatives.
'''
def get_stats_for_eval(predict_list, actual_list):
    true_pos = 0
    false_pos = 0
    true_neg = 0
    false_neg = 0
    for i in range(len(predict_list)):
        if predict_list[i] == 1:
            if actual_list[i] == 1:
                true_pos += 1
            else:
                false_pos += 1
        else:
            if actual_list[i] == 1:
                false_neg += 1
            else:
                true_neg += 1
    return {"true_pos": true_pos, "false_pos": false_pos, "true_neg": true_neg, "false_neg": false_neg}

'''
<score_file> should be a file with the train and set regions chained with the representative
regions. A pairwise chaining file should not be passed in here.
<rep_list_file> should be a file with the names of all the reference ACRs, each on their own line.

This function randomly creates a test and train set, trains a model based on the train set feature
vectors consisting of chaining scores to some central ACRs, predicts the region type of the regions in
test_set, and then returns a dictionary of performance statistics, as well as a dictionary of control statistics
(random predictions)
'''
def process(score_file, rep_list_file, fraction=.5, classifier=RandomForestClassifier()):
    rand_acr_file = "/home/mwarr/Data/One_Genome/exp3_ml/setb_and_rand_list.txt"
    train_set, test_set = get_region_sets(rand_acr_file)
    rep_lst, index_dict = rep_list(rep_list_file)
    train_dict, test_dict = get_chain_dicts(score_file, train_set, test_set, rep_lst, index_dict)
    predict, actual = train_predict(train_dict, test_dict, fraction=fraction, classifier=classifier)
    return predict, actual

'''
Makes a Matplotlib heatmap with <table_data> and saves it to <output_file>
'''
def create_heatmap(table_data, xlabels, ylabels, output_file, title):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.imshow(table_data)
    ax.set_xticks(range(len(xlabels)), labels=xlabels, rotation=45)
    ax.set_yticks(range(len(ylabels)), labels=ylabels)
    ax.tick_params(axis='both', labelsize=15)
    for i in range(len(xlabels)):
        for j in range(len(ylabels)):
            ax.text(i, j, table_data[j][i], fontsize=15, ha='center')

    plt.title(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(output_file)

'''
Driver for experiment where <fraction> of TEST regions are ACRs and <1-fraction> are not.
Outputs results to <output_file> in the format:
<type>
actual <accuracy> <f1>
'''
def driver_data_file(output_file, fraction=.5):
    base = "/home/mwarr/Data/One_Genome/exp3_ml"
    score_files = ["Chaining_local_rep_weight.tsv", "Chaining_global_rep_weight.tsv", "Chaining_global_consensus.tsv", "Chaining_local_consensus.tsv"]
    classifiers = [KNeighborsClassifier(), LinearSVC(), GaussianProcessClassifier(), RandomForestClassifier(), DecisionTreeClassifier(),
                   AdaBoostClassifier(), MLPClassifier(), GaussianProcessClassifier(), GaussianNB(), QuadraticDiscriminantAnalysis(), "control"]
    
    with open(output_file, "w") as file:
        for classifier in classifiers:
            file.write(f"{str(classifier)}\n")
            for score_file in score_files:
                if "rep" in score_file:
                    predict, actual = process(f"{base}/{score_file}", f"{base}/seta_rep_clustered_loc_list.txt", fraction=fraction, classifier=classifier)
                else:
                    predict, actual = process(f"{base}/{score_file}", f"{base}/seta_consensus_DAPv1_clustered_loc_list.txt", fraction=fraction, classifier=classifier)
                
                #get evaluation metrics
                accuracy = accuracy_score(actual, predict)
                f1 = f1_score(actual, predict)
                
                file.write(f"{score_file[:-4]}\n")
                file.write(f"accuracy\t{round(accuracy, 2)}\tf1\t{round(f1, 2)}\n")
                #file.write(f"{str(get_stats_for_eval(predict, actual))}\n")
            file.write("\n\n")

'''
Driver for experiment where <fraction> of TEST regions are ACRs and <1-fraction> are not.
Outputs 2 tables to <output_base> as .png files. One file has the accuracy scores, and the
other has f1 scores.
'''
def driver_heatmap(input_file, output_base, fraction=.5):
    table_f1 = []
    table_acc = []

    type = ["local rep", "global rep", "global consensus", "local consensus"]
    algo = []
    algo_count = 0
    with open(input_file, "r") as file:
        for line in file:
            algo_line = line.strip()
            if len(algo_line) < 2:
                continue
            if "control" in algo_line:
                algo.append(algo_line)
            else:
                algo.append(algo_line[: -2])
            row_f1 = []
            row_acc = []
            #for each type
            for _ in range(4):
                file.readline() #type line
                scores = file.readline().split("\t")
                acc = float(scores[1].strip())
                f1 = float(scores[3].strip())
                row_f1.append(f1)
                row_acc.append(acc)
            table_f1.append(row_f1)
            table_acc.append(row_acc)
            algo_count += 1

    create_heatmap(table_f1, type, algo, f"{output_base}/f1_heat_{fraction}.png", f"F1 scores, {round(fraction, 2) *100}% ACRs")
    create_heatmap(table_acc, type, algo, f"{output_base}/acc_heat_{fraction}.png", f"Accuracies, {round(fraction, 2) *100}% ACRs")
    
      
if __name__ == "__main__":
    #driver_data_file("/home/mwarr/Data/One_Genome/exp3_ml/results/90-10_algo_compare.txt", fraction=.1)
    driver_heatmap("/home/mwarr/Data/One_Genome/exp3_ml/results/50-50_algo_compare.txt", "/home/mwarr/Data/One_Genome/exp3_ml/results", fraction=.5)