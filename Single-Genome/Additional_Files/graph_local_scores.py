from glob import glob
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt

scores = []
nums = []

with open('/home/mwarr/Data/Chaining_one_local.tsv', 'r') as file:
    for line in file :
        line_arr = line.rstrip().split('\t')
        scores.append(float(line_arr[2]))
        nums.append(int(line_arr[3]))
        


plt.scatter(nums, scores, s = 0.1)
plt.savefig('/home/kyu/local_one_scores.png')