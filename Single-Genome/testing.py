import random
from rand_vs_upstream import *
import os

#testing rand_vs_upstream.py
def generate_chaining_file(length):
    try:
        os.remove("ACR_v_ACR_test.tsv")
        os.remove("ACR_v_rand_test.tsv")
    except:
        print("Probably file isn't there")
        
    references = set()
    
    for i in range(length):
        name = random_ACR_name()
        name2 = random_ACR_name()
        random_name = random_ACR_name()
        references.add(name)
        with open("ACR_v_ACR_test.tsv", "a") as file:
            file.write(f"{name}\t{name2}\t{random.randint(2, 10)}\t{random.randint(2, 10)}\n")
        with open("ACR_v_rand_test.tsv", "a") as file:
            file.write(f"{name}\t{random_name + "_rand"}\t{random.randint(2, 10)}\t{random.randint(2, 10)}\n")
    return references


def random_ACR_name():
    name = ""
    for i in range(10):
        name += chr(random.randint(97, 122))
    return name

# references = generate_chaining_file(20)
dicts = create_dicts("ACR_v_ACR_test.tsv", "ACR_v_rand_test.tsv", 20)
output_score_freq(dicts[0], dicts[1], ".", lambda lst: sum(lst) / len(lst), "avg")