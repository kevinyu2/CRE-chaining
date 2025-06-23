from glob import glob
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt

avg_scores = []
avg_nums = []
for chain_file in tqdm(Path('/home/mwarr/Data/Chaining_min1_local').rglob("*")):
    with open(chain_file, 'r') as f :
        total_score = 0
        total_num = 0
        total_count = 0
        for line in f :
            genome1, genome2, score, num = line.rstrip().split('\t')
            total_count += 1
            total_num += float(num)
            total_score += float(score)

        avg_scores.append(total_score / total_count)
        avg_nums.append(total_num / total_count)

plt.scatter(avg_nums, avg_scores, s = 0.1)
plt.savefig('/home/kyu/local_multi_scores.png')