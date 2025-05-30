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