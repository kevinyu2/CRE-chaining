def get_parent_from_attrs(attrs_str):
    attrs = attrs_str.split(";")
    for attr in attrs:
        if attr.startswith("Parent="):
            return attr.split("=", 1)[1].split('.')[0]
    return None

def get_id_from_attrs(attrs_str):
    attrs = attrs_str.split(";")
    for attr in attrs:
        if attr.startswith("ID="):
            return attr.split("=", 1)[1]
    return None

def gff_to_bed(input_gff, output_bed, tss_thresh, use_fiveputr) :
    with open(input_gff) as gff, open(output_bed, 'w') as out:
        # {gene_id : [promoter start, promoter stop, chrom, gene id, strand]}
        gene_start_stop_dict = {}

        for line in gff:
            # Housekeeping
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs_str = parts
            if feature == 'gene':
                gene_id = get_id_from_attrs(attrs_str)
                # print(attrs_str)

                # The strand determins where is upstream
                if strand == "+" :
                    gene_start_stop_dict[gene_id] = [max(0, int(start) - tss_thresh - 1), int(start) - 1, chrom, gene_id, strand]
                else :
                    gene_start_stop_dict[gene_id] = [int(end), int(end) + tss_thresh, chrom, gene_id, strand]
            
            if feature == 'five_prime_UTR' and use_fiveputr:
                gene_id = get_parent_from_attrs(attrs_str)
                if strand == "+" :
                    gene_start_stop_dict[gene_id][1] = end
                else :
                    gene_start_stop_dict[gene_id][0] = start
        
        for key in gene_start_stop_dict.keys() :
            # chrom, start, stop, id, 0, strand
            gene_tuple = gene_start_stop_dict[key]
            out.write(f"{gene_tuple[2]}\t{gene_tuple[0]}\t{gene_tuple[1]}\t{gene_tuple[3]}\t0\t{gene_tuple[4]}\n")


gff_to_bed("TAIR10_GFF3_genes.gff", "tair10promoter.bed", 2000, True)