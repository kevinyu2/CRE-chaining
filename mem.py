from collections import defaultdict

# This algorithm finds matching k-mer start locations then discounts the ones that get exended by tracking with a 2D array
def find_mems(k, seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    # Add all kmers of seq 1 into a dict
    k_dict = defaultdict(list)
    for i in range(0, n - k + 1):
        kmer = seq1[i:i+k]
        k_dict[kmer].append(i)

    # Matrix to track k-mers
    precense_matrix = [[0]*m]*n

    # List of mems
    mem_starts = []

    for j in range(0, m - k + 1):
        kmer = seq2[j:j+k]
        matches = k_dict[kmer]
        for i in matches:
            precense_matrix[i][j] = 1
            if j != 0 and i != 0 :
                if precense_matrix[i-1][j-1] == 0:
                    mem_starts.append((i, j))
            else :
                mem_starts.append((i, j))
    return mem_starts
    
print(find_mems(3, "ACGAT", "ACGAACGATTACG"))

print(find_mems(3, "CGGCGG", "CGGCGG"))
