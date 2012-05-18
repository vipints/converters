#!/usr/bin/python
"""Program to convert condensed state alignment file to 
COL-0 coordinate state.

Usage: condense_convert.py in.scafold_map in.mummer_map 
"""
import re, sys 

def grep_scaf(fname):
    """read condense mapping file and create a mapping db 
       between condensed contigs and original scafolds.
    """
    cont_scaf = dict()
    scafh=open(fname)
    for line in scafh:
        line=line.strip('\n\r').split('\t')
        if re.match(r"^#", line[0]) or line[1]=="NSPACER":
            continue
        if line[0] in cont_scaf:
            cont_scaf[line[0]][(int(line[2]), int(line[3]))]=line[1]
        else:
            cont_scaf[line[0]]={(int(line[2]), int(line[3])):line[1]}
    scafh.close()
    return cont_scaf

def grep_col(fname):
    """read mummer whole genome alignment file and fetch 
       the scafold to main genome contig mapping relation.
    """
    scaf_chr=dict()
    fh=open(fname)
    for line in fh:
        line=line.strip('\n\r').split('\t')
        if re.match(r"^#", line[0]):
            continue
        if line[2] in scaf_chr:
            scaf_chr[line[2]][(int(line[4]), int(line[4])+int(line[6])-1)]=(line[1], (int(line[3]), int(line[3])+int(line[5])-1), line[7])
        else:
            scaf_chr[line[2]]={(int(line[4]), int(line[4])+int(line[6])-1):(line[1], (int(line[3]), int(line[3])+int(line[5])-1), line[7])}
    fh.close()
    return scaf_chr

def __main__():
    try:
        scafname = sys.argv[1]
        mumfname = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    maps_scaf=grep_scaf(scafname)
    maps_chr=grep_col(mumfname)
    # strand ? at final stage 

if __name__=="__main__":__main__()
