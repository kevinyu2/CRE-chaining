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
def argsort_reverse_ties(arr) :
    len_arr = len(arr)

    # Increment each so argsort does it in reverse tiebreaker
    for i in range(len_arr) :
        arr[i] += float((len_arr - i - 1) / len_arr)

    return np.argsort(arr)

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


# mems are [(position1, position2, weight)]
def chain_weighted(mems):
    mems_len = len(mems)
    if mems_len <= 1 :
        return mems_len
    
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    arr = np.array(mems) 
    mems_a = arr[:, 0] 

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



# mems are numbered based on how many motifs apart (i.e not counting base pair location, but each motif is 1, 2, etc)
# match, mismatch, and gap scores like alignment
# Returns best local chain in the form of chain_score, chain_len (where the len is the best among all with the best chain score)
def chain_local(mems, match, mismatch, gap) :
    mems_len = len(mems)
    if mems_len == 0 :
        return 0
    
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    mems_a = [element[0] for element in mems]

    # Order in terms of a
    order = argsort_reverse_ties(mems_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((mems_len, 2))

    for mem_idx in order :
        max_score = match
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = match + score_list[j][0] + gap * abs(mems[mem_idx][0] - mems[j][0] - mems[mem_idx][1] + mems[j][1]) + mismatch * (min(mems[mem_idx][0] - mems[j][0], mems[mem_idx][1] - mems[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length
    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])


# mems are numbered based on how many motifs apart (i.e not counting base pair location, but each motif is 1, 2, etc)
# mismatch and gap scores like alignment. No match score is a multiplier to the weight.
# Returns best local chain in the form of chain_score, chain_len (where the len is the best among all with the best chain score)
def chain_local_weighted(mems, match, mismatch, gap) :
    mems_len = len(mems)
    if mems_len == 0 :
        return 0
    
    #sort mems by second value (but tiebreak by the first value, descending)
    mems.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    mems_a = [element[0] for element in mems]

    # Order in terms of a
    order = argsort_reverse_ties(mems_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((mems_len, 2))

    for mem_idx in order :
        max_score = match * mems[mem_idx][2]
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = match * mems[mem_idx][2] + score_list[j][0] + gap * abs(mems[mem_idx][0] - mems[j][0] - mems[mem_idx][1] + mems[j][1]) + mismatch * (min(mems[mem_idx][0] - mems[j][0], mems[mem_idx][1] - mems[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length

    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])

########################
# CALL THESE FUNCTIONS

def chain_driver(mems, is_weighted):
    if is_weighted:
        negative_mems = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in mems]
        return max(chain_weighted(mems), chain_weighted(negative_mems))

    else :
        negative_mems = [(anchor[0], -1 * anchor[1]) for anchor in mems]
        return (max(chain(mems), chain(negative_mems)))


def chain_local_driver(mems, match, mismatch, gap, is_weighted):
    if is_weighted:
        negative_mems = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in mems]
        return max(chain_local_weighted(mems, match, mismatch, gap)[0], chain_weighted(negative_mems, match, mismatch, gap)[0])

    else :
        negative_mems = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in mems]
        return (max(chain_local(mems, match, mismatch, gap)[0], chain_local(negative_mems, match, mismatch, gap)[0]))


# mems = [(0,0), (4,5), (2,8), (5,8), (1, 6), (6,12), (4,9)]
# print(chain_local(mems, 4, -2, -1))

# mems = [(0,0,2), (4,5,1), (2,8,3), (5,8,2), (1, 6,1.5), (6,12,3.2), (4,9,1.2)]
# print(chain_local_weighted(mems, 2, -2, -1))


# mems = [(1, 1, 4), (2, 1, 5), (3, 3, 1.5), (4, 10, 9), (5, 5, 6), (6, 6, 5)]
# print(chain_driver(mems, True))
# # == 17.5

# mems = [(3, 5, 10), (4, 7, 20), (5, 6, 1), (4, 5, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 1, 4)]
# print(chain_driver(mems, True))
# # == 35

# mems = [(4, 5, 1), (2, 3, 1), (3, 4, 1), (5, 1, 7)]
# print(chain_driver(mems, True))
# # == 8

# mems = [(1, 7), (2, 6), (3, 5), (4, 4)]
# print(chain_driver(mems, False))
# # == 4

# mems = [(1,3), (4,4), (5,3), (6,2)]
# print(chain_driver(mems, False))
# # == 3