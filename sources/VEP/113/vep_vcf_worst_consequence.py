#!/usr/bin/env python3

import sys
import itertools

cons_key = "VEP_Consequence"
impact_key = "VEP_IMPACT"
gene_key = "VEP_SYMBOL"
csn_key = "VEP_CSN"
sift_key = "VEP_SIFT"
polyphen_key = "VEP_PolyPhen"
biotype_key = "VEP_BIOTYPE"

new_cons_key = "VEP_Worst_Consequence"
new_impact_key = "VEP_Worst_Impact"
new_gene_key = "VEP_Worst_Gene"
new_csn_key = "VEP_Worst_CSN"
new_sift_key = "VEP_Worst_SIFT"
new_polyphen_key = "VEP_Worst_PolyPhen"
new_biotype_key = "VEP_Worst_BIOTYPE"

impact_ranked = {'HIGH': 1, 'MODERATE':2, 'LOW':3, 'MODIFIER':4}

# see: https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
consequences = '''transcript_ablation
splice_acceptor_variant
splice_donor_variant
stop_gained
frameshift_variant
stop_lost
start_lost
transcript_amplification
feature_elongation
feature_truncation
inframe_insertion
inframe_deletion
missense_variant
protein_altering_variant
splice_donor_5th_base_variant
splice_region_variant
splice_donor_region_variant
splice_polypyrimidine_tract_variant
incomplete_terminal_codon_variant
start_retained_variant
stop_retained_variant
synonymous_variant
coding_sequence_variant
mature_miRNA_variant
5_prime_UTR_variant
3_prime_UTR_variant
non_coding_transcript_exon_variant
intron_variant
NMD_transcript_variant
non_coding_transcript_variant
coding_transcript_variant
upstream_gene_variant
downstream_gene_variant
TFBS_ablation
TFBS_amplification
TF_binding_site_variant
regulatory_region_ablation
regulatory_region_amplification
regulatory_region_variant
intergenic_variant
sequence_variant'''

#cons_ranked = [(int(x[0]),x[1]) for x in [y.split(' ') for y in consequences.split('\n')]]
cons_ranked = {}
for rank,name in enumerate(consequences.split('\n')):
    #rank, name = y.split(' ')
    cons_ranked[name] = int(rank)

def parse_pos(line):
    cols = line.strip().split('\t')
    outcols = cols[:7]
    
    info = cols[7]
    
    info_vals = info.split(';')
    info_out = []

    genes_vals = []
    csns_vals = []
    poly_vals = []
    sift_vals = []
    impact_vals = []
    cons_vals = []
    biotype_vals = []

    gene_idx = -1

    # The values are a modified CSV list. Meaning, if there is only one unique value, then only one value will be returned (and there won't be a comma).
    # This is the same for consequence and impact too. But, if there is only one value for those, and there are multiple transcripts, we will return the 
    # first CSN/Gene/etc... as the first transcript returned, should be the the longest/higher priority

    for ival in info_vals:
        if '=' not in ival:
            info_out.append(ival)
        else:
            info_out.append(ival)
            left, right = ival.split('=')
            if left == gene_key:
                genes_vals = right.split(',')
            elif left == csn_key:
                csns_vals = right.split(',')
            elif left == polyphen_key:
                poly_vals = right.split(',')
            elif left == sift_key:
                sift_vals = right.split(',')
            elif left == impact_key:
                impact_vals = right.split(',')
            elif left == cons_key:
                cons_vals = right.split(',')
            elif left == biotype_key:
                biotype_vals = right.split(',')
            else:
                # skip all others for priority?
                pass

    if impact_vals and cons_vals:
        merged = []
        for idx, (impact, cons) in enumerate(zip(impact_vals, cons_vals)):
            for cons2 in cons.split("&"):
                merged.append((impact_ranked[impact], cons_ranked[cons2], idx, impact, cons2))
            
        merged = sorted(merged)
        worst_idx = merged[0][2]
        
        info_out.append('%s=%s' % (new_impact_key, merged[0][3]))
        info_out.append('%s=%s' % (new_cons_key, merged[0][4]))
        
        if genes_vals:
            # If there is only one gene, it won't be comma delimited
            newgene = ''
            if len(genes_vals) == 1:
                newgene = genes_vals[0]
            elif len(genes_vals) > worst_idx:
                newgene = genes_vals[worst_idx]

            if newgene:
                info_out.append('%s=%s' % (new_gene_key, newgene))

        if biotype_vals:
            # If there is only one biotype, it won't be comma delimited
            newbiotype = ''
            if len(biotype_vals) == 1:
                newbiotype = biotype_vals[0]
            elif len(biotype_vals) > worst_idx:
                newbiotype = biotype_vals[worst_idx]

            if newbiotype:
                info_out.append('%s=%s' % (new_biotype_key, newbiotype)) 

        if csns_vals and len(csns_vals) > worst_idx and csns_vals[worst_idx]:
            info_out.append('%s=%s' % (new_csn_key, csns_vals[worst_idx]))
        if poly_vals and len(poly_vals) > worst_idx and poly_vals[worst_idx]:
            info_out.append('%s=%s' % (new_polyphen_key, poly_vals[worst_idx]))
        if sift_vals and len(sift_vals) > worst_idx and sift_vals[worst_idx]:
            info_out.append('%s=%s' % (new_sift_key, sift_vals[worst_idx]))


    # if genes and gene_idx > -1:
    #     if len(genes) > 1:
    #         info_out.append('%s=%s' % (new_gene_key, genes[gene_idx]))
    #     else:
    #         info_out.append('%s=%s' % (new_gene_key, genes[0]))

    # if csns and gene_idx > -1:
    #     if len(csns) > 1:
    #         info_out.append('%s=%s' % (new_csn_key, csns[gene_idx] if csns[gene_idx] else '.'))
    #     else:
    #         info_out.append('%s=%s' % (new_csn_key, csns[0] if csns[0] else '.'))

    outcols.append(';'.join(info_out))
    outcols.extend(cols[8:])

    return outcols



def parse_main(fin=sys.stdin, fout=sys.stdout):

    printed_header = False

    for line in fin:
        if line[0] == '#' and line[1] == '#':
            fout.write(line)
        elif line[0] == '#':
            if not printed_header:
                printed_header = True
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst Consequence">\n' % (new_cons_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst Impact">\n' % (new_impact_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst Gene">\n' % (new_gene_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst SIFT">\n' % (new_sift_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst PolyPhen">\n' % (new_polyphen_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst CSN">\n' % (new_csn_key))
                fout.write('##INFO=<ID=%s,Number=1,Type=String,Description="VEP Worst BIOTYPE">\n' % (new_biotype_key))
            fout.write(line)
        else:
            fout.write('%s\n' % '\t'.join(parse_pos(line)))
    

if __name__ == '__main__':
    parse_main()
