import matplotlib.pyplot as plt
from pathlib import Path


intra = True

if not intra:
    distance_file = '/home/mwarr/Data/Clustering/Distances_min1_alpha0.tsv'


    distances = []
    with open(distance_file, 'r') as f:
        for line in f:
            line_arr = line.split('\t')
            distances.append(float(line_arr[2]))

    plt.hist(distances)
    plt.yscale('log')
    plt.savefig('/home/mwarr/Distances_min1_alpha0.png')

else :
    distance_folder = Path('/home/mwarr/Data/Clustering/Distances_min1_intra_alpha50')

    num_files = 0


    distances = []
    for intra_file in distance_folder.glob('*'):

        num_files += 1

        # if num_files == 100 :
        #     break

        with open(intra_file, 'r') as f:
            skip_file = False
            for line in f:
                line_arr = line.split('\t')
                if len(line_arr) < 3:
                    skip_file = True
                    break
                distances.append(float(line_arr[2]))
            if skip_file:
                continue

    plt.hist(distances)
    plt.savefig('/home/mwarr/Distances_alpha50_intra')
        



