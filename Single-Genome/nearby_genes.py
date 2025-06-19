
def genes_near_ACRs(input_file_path, output_dir, threshold):
    far_gene_count = 0

    with open(f"{output_dir}/genes_threshold_{threshold}.tsv", "w") as gene_file:
        gene_file.write("Chromosome\tStart\tStop\tStrand")
        with open(f"{output_dir}/pseudo_threshold_{threshold}.tsv", "w") as pseudo_file:
            pseudo_file.write("Chromosome\tStart\tStop\tStrand")

            with open(input_file_path, "r") as input_file:
                for line in input_file:
                    line_arr = line.split("\t")
                    chromosome = line_arr[3]
                    gene_type = line_arr[5]
                    start = line_arr[6]
                    stop = line_arr[7]
                    strand = line_arr[9]
                    distance = abs(int(line_arr[-1]))
                    
                    if distance < threshold:
                        if gene_type == "gene":
                            gene_file.write(f"\n{chromosome}\t{start}\t{stop}\t{strand}")
                        elif gene_type == "pseudogene":
                            pseudo_file.write(f"\n{chromosome}\t{start}\t{stop}\t{strand}")
                        else:
                            print(f"Found {gene_type} gene type")
                    else:
                        far_gene_count += 1
    print(f"There were {far_gene_count} far genes with threshold {threshold}")

genes_near_ACRs("/home/mwarr/Data/nearby_genes/all_dACRs_nearby_genes_iu.bed", "/home/mwarr/Data/nearby_genes", 2000)
