"""
Binary Indexed Tree / Fenwick Tree
https://www.hackerearth.com/practice/notes/binary-indexed-tree-made-easy-2/
https://www.topcoder.com/community/data-science/data-science-tutorials/binary-indexed-trees/
https://www.youtube.com/watch?v=v_wj_mOAlig
https://www.youtube.com/watch?v=kPaJfAUwViY
"""


def bit_update(index, value, array, bi_tree):
	"""
	Updates the binary indexed tree with the given value
	:param index: index at which the update is to be made
	:param value: the new element at the index
	:param array: the input array
	:param bi_tree: the array representation of the binary indexed tree
	:return: void
	"""

	while index < len(array):
		bi_tree[index] += value
		index += index & -index


def bit_get_max(index, bi_tree):
	"""
	Calculates the max of the elements from the beginning to the index
	:param index: index till which the sum is to be calculated
	:param bi_tree: the array representation of the binary indexed tree
	:return: max of the elements from beginning till index
	"""
	ans = 0

	while index > 0:
		ans = max(ans, bi_tree[index])
		index -= index & -index

	return ans


