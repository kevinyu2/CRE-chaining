import random
import matplotlib.pyplot as plt
import numpy as np
from chaining import *

# Makes num_anchors number of anchors, each a random number from 0 to len_sequence
def make_random_anchors(num_anchors, len_sequence) :
    anchors = set()
    while len(anchors) < num_anchors :
        start1 = random.randint(0, len_sequence)
        start2 = random.randint(0, len_sequence)

        anchors.add((start1, start2))
    return list(anchors)

def make_permutation(num) :
    perm = random.sample(range(num), num)
    return perm

# average_chain_lens = []
average_lis_len = []
# stdev_chain_lens = []
stdev_lis_lens = []
x = range(2, 200)
for i in x :

    chain_len = []
    lis_len = []

    for j in range(10000) :
        # chain_len.append(chain(make_random_anchors(i, 1000)))
        lis_len.append(lengthOfLIS(make_permutation(i)))
    
    # average_chain_lens.append(np.mean(chain_len))
    # stdev_chain_lens.append(np.std(chain_len))
    average_lis_len.append(np.mean(lis_len))
    stdev_lis_lens.append(np.std(lis_len))

from scipy.optimize import curve_fit

def model(x, a, b):
    return a * np.power(x, b)

params, _ = curve_fit(model, x, stdev_lis_lens)
a_fit, b_fit = params

# Fit the curve
c_fit = params[0]
print(f"Fitted parameter constant = {a_fit}")
print(f"Fitted power = {b_fit}")

# Generate fitted curve points for smooth plotting
y_fit = []
for xval in x :
    y_fit.append(xval ** b_fit * a_fit)

with open('/home/kyu/LIS_len.tsv', 'w') as f:
    f.write("Size\tLIS_len_avg\tLIS_len_stdev\n")
    for i in range(len(x)) :
        f.write(f"{x[i]}\t{average_lis_len[i]}\t{stdev_lis_lens[i]}\n")

plt.figure()
plt.plot(x, stdev_lis_lens, label = 'Stdev LIS Length')
plt.plot(x, y_fit, label = 'Model')
plt.legend()
plt.xlabel("Permutation Size")
plt.title(f"LIS Length Stdev Vs Model (cx^(1/6), c = {c_fit})")
plt.savefig("/home/kyu/model_LIS_stdev.png")

# The theory suggested model for expectation
y_fit = []

for xval in x :
    y_fit.append(2 * np.sqrt(xval) - 1.771088 * (xval ** (1/6)))
plt.figure()
plt.plot(x, average_lis_len, label = 'Avg LIS Length')
plt.plot(x, y_fit, label = 'Model')
plt.legend()
plt.xlabel("Permutation Size")
plt.title(f"LIS Length Stdev Vs Model")
plt.savefig("/home/kyu/model_LIS_exp.png")


# plt.plot(x, [i / (num + 2) for num, i in enumerate(average_chain_lens)], label = 'Average Chain Ratio')
'''
plt.figure()
x = range(2, 200)
plt.plot(x, average_chain_lens, label = 'Average Chain Length')
plt.plot(x, average_lis_len, label = 'Average LIS Length')
# plt.plot(x, [i / (num + 2) for num, i in enumerate(average_chain_lens)], label = 'Average Chain Ratio')

plt.xlabel("Number of Anchors")
plt.legend()
plt.savefig('/home/kyu/Expected_Chain.png')

plt.figure()
x = range(2, 200)
plt.plot(x, stdev_chain_lens, label = 'Stdev Chain Length')
plt.plot(x, stdev_lis_lens, label = 'Stdev LIS Length')
# plt.plot(x, [i / (num + 2) for num, i in enumerate(average_chain_lens)], label = 'Average Chain Ratio')

plt.xlabel("Number of Anchors")
plt.legend()
plt.savefig('/home/kyu/Stdev_Chain.png')


'''