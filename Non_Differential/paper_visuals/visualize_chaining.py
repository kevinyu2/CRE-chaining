#Chat GPT
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib
import numpy as np

def build_motif_color_map(test_motifs, rep_motifs):
    unique_motifs = sorted(set(test_motifs) | set(rep_motifs))
    n = len(unique_motifs)

    cmap = matplotlib.colormaps.get_cmap("tab20")

    motif_to_color = {
        motif: cmap(i / max(n - 1, 1))
        for i, motif in enumerate(unique_motifs)
    }
    return motif_to_color

def parse_motif_file(filepath):
    test_motifs = []
    rep_motifs = []

    current_region = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line == "test_region":
                current_region = "test"
            elif line == "rep_region":
                current_region = "rep"
            elif line.startswith("motif_"):
                motif_id = line.split("_")[1]
                if current_region == "test":
                    test_motifs.append(motif_id)
                elif current_region == "rep":
                    rep_motifs.append(motif_id)

    return test_motifs, rep_motifs

def get_anchors(test_motifs, rep_motifs):
    anchors = []
    for i, motif1 in enumerate(test_motifs):
        for j, motif2 in enumerate(rep_motifs):
            if motif1 == motif2:
                anchors.append((i, j))
    return anchors


def plot_motif_connections(test_motifs, rep_motifs, chain_pairs, cut=0):
    fig, ax = plt.subplots(figsize=(7.5, 2))

    # y positions
    y_test = 0
    y_rep = 1

    # x positions (left to right)
    x_test = ([i*2 for i in range(len(test_motifs))])[::-1]
    x_rep = list(range(len(rep_motifs)))[::-1]

    # Draw matching lines
    for i, motif in enumerate(test_motifs):
        for j, rep_motif in enumerate(rep_motifs):
            if j < cut:
                continue
            if motif == rep_motif and (i, j) in chain_pairs:
                ax.plot(
                    [x_test[i], x_rep[j]],
                    [y_test, y_rep],
                    linewidth=1.5,
                    alpha=1,
                    color="black",
                    zorder=2
                )
            elif motif == rep_motif:
                ax.plot(
                    [x_test[i], x_rep[j]],
                    [y_test, y_rep],
                    linewidth=1,
                    alpha=0.5,
                    color="lightgrey",
                    zorder=1
                )
    
    motif_to_color = build_motif_color_map(test_motifs, rep_motifs)

    for i, motif in enumerate(test_motifs):
        ax.scatter(
            x_test[i],
            y_test,
            s=50,
            color=motif_to_color[motif],
            zorder=3
        )

    for i, motif in enumerate(rep_motifs):
        if i < cut:
            continue
        ax.scatter(
            x_rep[i],
            y_rep,
            s=50,
            color=motif_to_color[motif],
            zorder=3
        )

    # Labels and formatting
    ax.set_yticks([y_test, y_rep])
    ax.set_yticklabels(["Predicted ACR", "Known ACR"])
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(-1, max(len(x_test), len(x_rep) - cut))
    ax.set_xticks([])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    plt.tight_layout()
    plt.savefig("chaining_top_weighted_cut.png", dpi=300)


#Code from Geeks For Geeks -- geeksforgeeks.org/segment-tree-set-2-range-maximum-query-node-update/
from math import ceil, log

# A utility function to get the
# middle index of given range.
def getMid(s, e):
    return s + (e - s) // 2

# /* A recursive function to get the sum of
    # values in given range of the array.
    # The following are parameters for this
    # function.
    #
    # st     -> Pointer to segment tree
    # node     -> Index of current node in
    #             the segment tree .
    # ss & se -> Starting and ending indexes
    #             of the segment represented
    #             by current node, i.e., st[node]
    # l & r -> Starting and ending indexes
    #             of range query */
def MaxUtil(st, ss, se, l, r, node):

    # If segment of this node is completely
    # part of given range, then return
    # the max of segment
    if (l <= ss and r >= se):
        return st[node]

    # If segment of this node does not
    # belong to given range
    if (se < l or ss > r):
        return -1

    # If segment of this node is partially
    # the part of given range
    mid = getMid(ss, se)

    return max(MaxUtil(st, ss, mid, l, r,
                       2 * node + 1),
               MaxUtil(st, mid + 1, se, l,
                       r, 2 * node + 2))

#
# /* A recursive function to update the nodes which
# have the given index in their range. The following
# are parameters st, ss and se are same as defined
# above index -> index of the element to be updated.*/
def updateValue(arr, st, ss, se, index, value, node):
    if (index < ss or index > se):
        print("Invalid Input")
        return

    if (ss == se):

        # update value in array and in segment tree
        arr[index] = value
        st[node] = value
    else:
        mid = getMid(ss, se)

        if (index >= ss and index <= mid):
            updateValue(arr, st, ss, mid, index,
                        value, 2 * node + 1)
        else:
            updateValue(arr, st, mid + 1, se,
                        index, value, 2 * node + 2)

        st[node] = max(st[2 * node + 1],
                       st[2 * node + 2])
    return

# Return max of elements in range from
# index l (query start) to r (query end).
def getMax(st, n, l, r):

    # Check for erroneous input values
    if (l < 0 or r > n - 1 or l > r):
        print("Invalid Input")
        return -1

    return MaxUtil(st, 0, n - 1, l, r, 0)

# A recursive function that constructs Segment
# Tree for array[ss..se]. si is index of
# current node in segment tree st
def constructSTUtil(arr, ss, se, st, si):

    # If there is one element in array, store
    # it in current node of segment tree and return
    if (ss == se):
        st[si] = arr[ss]
        return arr[ss]

    # If there are more than one elements, then
    # recur for left and right subtrees and
    # store the max of values in this node
    mid = getMid(ss, se)

    st[si] = max(constructSTUtil(arr, ss, mid, st,
                                 si * 2 + 1),
                 constructSTUtil(arr, mid + 1, se,
                                 st, si * 2 + 2))

    return st[si]
#
# /* Function to construct segment tree from given array.
# This function allocates memory for segment tree.*/
def constructST(arr, n):

    # Height of segment tree
    x = ceil(log(n, 2))

    # Maximum size of segment tree
    max_size = 2 * pow(2, x) - 1

    # Allocate memory
    st = [0]*max_size

    # Fill the allocated memory st
    constructSTUtil(arr, 0, n - 1, st, 0)

    # Return the constructed segment tree
    return st

# This code is contributed by mohit kumar 29


# Do argsort but break ties by taking the later one in the array first
# This is done by just adding a small amount less than 1 to each in the array,
# Giving a greater value to the earlier ones (so argsort prioritizes later ones)
def argsort_reverse_ties(arr) :
    len_arr = len(arr)

    # Increment each so argsort does it in reverse tiebreaker
    for i in range(len_arr) :
        arr[i] += float((len_arr - i - 1) / len_arr)

    return np.argsort(arr)

#MEMs should be a list of tuples. Each entry in the list is an MEM. For each MEM, 
#there is a tuple where the first entry is the index of the match in the first sequence
#and the second entry is the index of the match in the second sequence (a, c)
def chain(mems):
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    mems_len = len(mems)

    #get list of first values (a)
    mems_a = [element[0] for element in mems]

    # Order in terms of a
    order = argsort_reverse_ties(mems_a)

    # Array for the segment tree
    arr = [0]*mems_len
    maxChainST = constructST(arr, mems_len)

    for index in order :
        # Get the max length chain up to the current c
        maxPrev = getMax(maxChainST, mems_len, 0, index)
        #update arr[index] to maxPrev + 1
        updateValue(arr, maxChainST, 0, mems_len - 1, index, maxPrev + 1, 0)

    curr_len = max(arr)
    a_ind = []
    for i in range(len(arr) - 1, -1, -1):
        if arr[i] == curr_len:
            a_ind.append(i)
            curr_len -= 1

    #iterate through arr backwards. Start at 7 and find the next smallest number
    #save those indices (a's which are part of chain; test regions)
    #for each index, find that number in 'order' and save the index that number is at (c's which are part of chain; rep regions)
    return a_ind

if __name__ == "__main__":
    filepath = "/home/mwarr/Data/paper_visualizing/results/top_weighted_motifs.txt"
    test_motifs, rep_motifs = parse_motif_file(filepath)
    anchors = get_anchors(test_motifs, rep_motifs)
    anchor_ind = chain(anchors)
    chain_pairs = set()
    for ind in anchor_ind:
        chain_pairs.add(anchors[ind])
    plot_motif_connections(test_motifs, rep_motifs, chain_pairs, cut=50)

