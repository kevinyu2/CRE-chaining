from segment_tree import *
import bisect

import numpy as np

'''
Performs co-linear chainging with two different implementations. One implementation allows
for weights to be specified on the anchors. This implementation uses a binary indexed tree. The
other implementation is much faster and reduces the problem to longest increasing subsequence.
'''

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


# mems are [(position1, position2, weight)]
def chain_weighted(mems):
    mems_len = len(mems)
    if mems_len <= 1 :
        return mems_len
    
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    mems_a = [element[0] for element in mems]

    # Order in terms of a
    order = argsort_reverse_ties(mems_a)

    bit = PrefixMaxBIT(mems_len)
    for i in order:
        maxPrev = bit.query(i + 1)
        bit.update(i + 1, maxPrev + mems[i][2])
    
    return bit.query(mems_len)

###############################
# BIT

class PrefixMaxBIT:
    def __init__(self, size):
        self.n = size
        self.tree = [float(0)] * (self.n + 1)  # 1-based index

    def update(self, i, val):
        while i <= self.n:
            self.tree[i] = max(self.tree[i], val)
            i += i & -i

    def query(self, i):
        res = float('-inf')
        while i > 0:
            res = max(res, self.tree[i])
            i -= i & -i
        return res

# mems = [(1, 1, 4), (2, 1, 5), (3, 3, 1.5), (4, 10, 9), (5, 5, 6), (6, 6, 5)]
# print(chain_weighted(mems) == 17.5)

# mems = [(3, 5, 10), (4, 7, 20), (5, 6, 1), (4, 5, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 1, 4)]
# print(chain_weighted(mems) == 35)

# mems = [(4, 5, 1), (2, 3, 1), (3, 4, 1), (5, 1, 7)]
# print(chain_weighted(mems) == 7)