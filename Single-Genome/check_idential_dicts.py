# Quick check to verify the results of parallelziation did not change anything
import pickle

parallel = pickle.load('test_unordered.pkl')
non_parallel = pickle.load('test_unpar.pkl')

# Check if keys are the same
keys1 = set(parallel.keys())
keys2 = set(non_parallel.keys())

print("Pairs Not in Parallel:")
for i in keys2 - keys1 :
    print(i)
print("Pairs Not in Non Parallel:")
for i in keys1 - keys2 :
    print(i)




for pair, anchors in non_parallel :
    anchor_set = set(anchors)
    if pair in keys1 :
        par_set = set(parallel[pair])
        if par_set != anchor_set :
            print(f"Mismatch: {pair}")
            print(anchor_set)
            print(par_set)

    