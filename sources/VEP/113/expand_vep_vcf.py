#!/usr/bin/env python3

import sys

prefix="VEP_"
vep_info_name = "CSQ"


def parse_info_format(line):
###INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|...">

    vals = []
    buf = ''
    idx = line.find('Format: ') + len('Format: ')

    while line[idx] != '"':
        if line[idx] == '|':
            vals.append(buf)
            buf = ''
        else:
            buf += line[idx]
        idx += 1

    if buf:
        vals.append(buf)
    return vals


def parse_pos(vep_cols, line):
    cols = line.strip().split('\t')
    outcols = cols[:7]
    
    info = cols[7]
    
    info_vals = info.split(';')
    info_out = []
    for ival in info_vals:
        if '=' not in ival:
            info_out.append(ival)
        else:
            left, right = ival.split('=')
            if left != vep_info_name:
                info_out.append(ival)
            else:
                vep_list = right.split(',')
                vep_out = {}
                for vlval in vep_list:
                    vep_vals = vlval.split('|')
                    for i,vval in enumerate(vep_vals):
                        if not vep_cols[i] in vep_out:
                            vep_out[vep_cols[i]] = [vval,]
                        else:
                            vep_out[vep_cols[i]].append(vval)

                for k in vep_out:
                    diff = False
                    empty = True
                    for v in vep_out[k]:
                        if v:
                            empty = False
                        if v != vep_out[k][0]:
                            diff = True
                        
                    if not empty:
                        if not diff:
                            info_out.append('%s%s=%s' % (prefix, k, vep_out[k][0]))
                        else:
                            info_out.append('%s%s=%s' % (prefix, k, ','.join(vep_out[k])))




    outcols.append(';'.join(info_out))
    outcols.extend(cols[8:])

    return outcols



def parse_main(fin=sys.stdin, fout=sys.stdout):

    vep_cols = []
    for line in fin:
        if line[0] == '#':
            if line.startswith('##INFO=<ID=%s,' % vep_info_name):
                vep_cols = parse_info_format(line)
                for col in vep_cols:
                    if col.endswith('AF'):
                        fout.write('##INFO=<ID=%s%s,Number=.,Type=Float,Description="VEP Annotation %s">\n' % (prefix, col, col))
                    else:
                        fout.write('##INFO=<ID=%s%s,Number=.,Type=String,Description="VEP Annotation %s">\n' % (prefix, col, col))
            else:
                fout.write(line)

        else:
            fout.write('%s\n' % '\t'.join(parse_pos(vep_cols, line)))
    

if __name__ == '__main__':
    parse_main()
