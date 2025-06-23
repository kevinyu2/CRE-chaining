# Quick check to verify the results of parallelziation did not change anything
# import pickle

# parallel = pickle.load('test_unordered.pkl')
# non_parallel = pickle.load('test_unpar.pkl')

# # Check if keys are the same
# keys1 = set(parallel.keys())
# keys2 = set(non_parallel.keys())

# print("Pairs Not in Parallel:")
# for i in keys2 - keys1 :
#     print(i)
# print("Pairs Not in Non Parallel:")
# for i in keys1 - keys2 :
#     print(i)




# for pair, anchors in non_parallel :
#     anchor_set = set(anchors)
#     if pair in keys1 :
#         par_set = set(parallel[pair])
#         if par_set != anchor_set :
#             print(f"Mismatch: {pair}")
#             print(anchor_set)
#             print(par_set)

soln_dict = {}

with open("/home/mwarr/Data/Chaining_one.tsv", 'r') as f:
    for line in f:
        line_arr = line.rstrip().split('\t')
        soln_dict[(line_arr[0], line_arr[1])] = (line_arr[2], line_arr[3])

print("Checking")
with open("./Chaining_one_par_f.txt", 'r') as f:
    for line in f:
        line_arr = line.rstrip().split('\t')
        if soln_dict[(line_arr[0], line_arr[1])] != (line_arr[2], line_arr[3]) :
            print(f"{line_arr}, {soln_dict[(line_arr[0], line_arr[1])]}")
        # soln_dict[(line_arr[0], line_arr[1])] = (line_arr[2], line_arr[3])