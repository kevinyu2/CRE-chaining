# from segment_tree import *
import bisect

import numpy as np

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
    mems_len = len(mems)
    if mems_len == 0 :
        return 0
    
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    mems_a = [element[0] for element in mems]

    # Order in terms of a
    order = argsort_reverse_ties(mems_a)

    return lengthOfLIS(order)

# From a leetcode submission
def lengthOfLIS(nums) :

    sub = []
    for x in nums:
        if len(sub) == 0 or sub[-1] < x:
            sub.append(x)
        else:
            idx = bisect.bisect_left(sub, x)  # Find the index of the first element >= x
            sub[idx] = x  # Replace that number with x

    return len(sub)