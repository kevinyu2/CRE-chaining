from collections import defaultdict
import os
from glob import glob
import heapq

class Node:
    def __init__(self, start, stop, label):
        self.start = start
        self.stop = stop
        self.label = label
        self.left = None
        self.right = None

    def __repr__(self):
        return f"Node({self.start}, {self.stop}, '{self.label}')"

class BST:
    def __init__(self):
        self.root = None

    def insert(self, start, stop, label):
        new_node = Node(start, stop, label)
        if self.root is None:
            self.root = new_node
        else:
            self._insert(self.root, new_node)

    def _insert(self, current, new_node):
        if new_node.start < current.start:
            if current.left is None:
                current.left = new_node
            else:
                self._insert(current.left, new_node)
        else:
            if current.right is None:
                current.right = new_node
            else:
                self._insert(current.right, new_node)

    def inorder(self):
        return self._inorder(self.root)

    def _inorder(self, node):
        if node is None:
            return []
        return self._inorder(node.left) + [node] + self._inorder(node.right)

    def search_le(self, target):
        return self._search_le(self.root, target, None)

    def _search_le(self, node, target, best):
        if node is None:
            return best

        if node.start == target:
            return node
        elif node.start < target:
            # Node is a candidate; go right to look for closer one
            return self._search_le(node.right, target, node)
        else:
            # Go left to look for a smaller one
            return self._search_le(node.left, target, best)

def get_starts_and_stops(acr_set_file) :
    acr_trees = defaultdict(BST)

    with open(acr_set_file, 'r') as asf:
        for line in asf :
            start = int(line.rstrip().split('_')[1].split('to')[0])
            stop = int(line.rstrip().split('_')[1].split('to')[1])
            chr = line.split('_')[0]
            acr_trees[chr].insert(start, stop, line.rstrip())

    return acr_trees


def get_motif_loc_dict(data_dir, acr_bsts) :
    # Holds where each motif is located {acr: {motif: [loc, ..., loc], motif:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))

    print("Filling Motif Dict")

    pattern = os.path.join(data_dir, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob(pattern)

    for fimo_file in fimo_files:
        
        # Fill in the dict
        with open(fimo_file, "r") as f:
            # Ignore first line
            next(f)
            for line in f: 
                line_arr = line.rstrip().split('\t')
                # The file ends tsv info early
                if len(line_arr) < 5: 
                    break
                node = acr_bsts[line_arr[2]].search_le(int(line_arr[3]))
                if node != None :
                    if node.stop >= int(line_arr[4]) :
                        acr = node.label
                        if int(line_arr[3]) not in motif_loc_dict[acr][line_arr[0]] :
                            motif_loc_dict[acr][line_arr[0]].append(int(line_arr[3]))
    
    # Remove duplicates from repeated sequences (basically remove overlaps)
    for acr, single_acr_dict in motif_loc_dict.items() :
        for motif, loc_list in single_acr_dict.items() :
            motif_len = len(motif)
            loc_list.sort()
            to_remove = set()
            for i in range(len(loc_list) - 1) :
                if loc_list[i + 1] - loc_list[i] <= motif_len :
                    to_remove.add(loc_list[i])
            single_acr_dict[motif] = [x for x in loc_list if x not in to_remove]

    return motif_loc_dict




def merge_dict_lists(list_dict):
    heap = []
    result = []

    # Initialize the heap with the first element of each list
    for key, lst in list_dict.items():
        if lst:  # Skip empty lists
            heapq.heappush(heap, (lst[0], key, 0))  # (value, source_key, index_in_list)

    while heap:
        value, key, idx = heapq.heappop(heap)
        result.append((value, key))

        # If there are more items in this list, push the next one
        if idx + 1 < len(list_dict[key]):
            next_value = list_dict[key][idx + 1]
            heapq.heappush(heap, (next_value, key, idx + 1))

    return result

def print_sequences(motif_loc_dict, outfile) :
    with open(outfile, 'w') as out :
        for acr, single_acr_dict in motif_loc_dict.items() :
            out.write(f"ACR: {acr}\n")
            result_list = merge_dict_lists(single_acr_dict)
            for motif_loc_pair in result_list :
                out.write(f'{motif_loc_pair[1]}\n')

acr_bsts = get_starts_and_stops('/home/mwarr/Data/seta.txt')

mld = get_motif_loc_dict("/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_ArabidopsisDAPv1/", acr_bsts)

print_sequences(mld, '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/acr_motif_sequences_ArabidopsisDAPv1.txt')