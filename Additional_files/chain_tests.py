from chaining import chain
import random 
from collections import defaultdict
import time

'''
Tests our chaining file by comparing the output with longest common subsequence (should be equivalent)
'''

def generate_str(length):
    rand_str = ""
    for i in range(length):
        rand_str += chr(random.randint(1, 127))
    return rand_str

def lcs(X, Y): 
    # find the length of the strings 
    m = len(X) 
    n = len(Y) 

    # declaring the array for storing the dp values 
    L = [[None]*(n + 1) for i in range(m + 1)] 

    """Following steps build L[m + 1][n + 1] in bottom up fashion 
    Note: L[i][j] contains length of LCS of X[0..i-1] 
    and Y[0..j-1]"""
    for i in range(m + 1): 
        for j in range(n + 1): 
            if i == 0 or j == 0 : 
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: 
                L[i][j] = L[i-1][j-1]+1
            else: 
                L[i][j] = max(L[i-1][j], L[i][j-1]) 

    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return L[m][n] 

def lcs_test() :
    rand_str_1 = generate_str(random.randint(2000, 3000))
    rand_str_2 = generate_str(random.randint(2000, 3000))

    start_time = time.time()
    str1_dict = defaultdict(list)
    for i in range(len(rand_str_1)) :
        str1_dict[rand_str_1[i]].append(i)

    anchors = []
    for i in range(len(rand_str_2)) :
        for str1_index in str1_dict[rand_str_2[i]] :
            anchors.append((str1_index, i))
    print(f"Anchor time: {time.time() - start_time}")
    start_time = time.time()
    dp_arr = chain(anchors)
    chain_max = max(dp_arr)
    print(f"chaining: {time.time() - start_time} (num anchors: {len(anchors)})")
    start_time = time.time()
    lcs_max = lcs(rand_str_1, rand_str_2)
    print(f"lcs: {time.time() - start_time}")
    assert(chain_max == lcs_max)


lcs_test()