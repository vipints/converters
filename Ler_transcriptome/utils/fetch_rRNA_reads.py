#!/usr/bin/env python
"""
Program to fetch reads falling on ribosomal RNA region.
Resulting file in BAM format.

Usage: samtools view -h in.bam | grep -e "^@" | fetch_rRNA_reads.py in_rRNA.bed in.bam -f| samtools view -bS - > out.bam 
"""
import re, sys 
import time
import pysam 

def get_Feature(fname):
    """Extract genome annotation information from a provided GFF file.
    """
    ribo_featdb=dict()
    fh=open(fname, 'rU')
    for line in fh:
        line=line.strip('\n\r').split(' ')
        if line[0] in ribo_featdb:
            ribo_featdb[line[0]][(int(line[1]), int(line[2]))]=line[3]
        else:
            ribo_featdb[line[0]]={(int(line[1]), int(line[2])):line[3]}
    fh.close()
    return ribo_featdb

def AlignGenerator(fbam, rib_db):
    """Parsing alignment file, extracting the information.
    """
    samfile = pysam.Samfile(fbam, 'rb') 
    outsam = pysam.Samfile('-', 'w', template = samfile)
    for rec in samfile.fetch():
        for details, strand_info in rib_db[samfile.getrname(rec.rname)].items():
            if details[0]-42 <= int(rec.pos) and int(rec.pos) <= details[1]+5:
                outsam.write(rec)
                break
    samfile.close()
    outsam.close()

if __name__=="__main__":
    try:
        anno_file=sys.argv[1]
        bamf=sys.argv[2]
        header=sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
    ribo_featdb=get_Feature(anno_file) 
    if header=='-f':
        for line in sys.stdin:
            line=line.strip()
            print line 
    AlignGenerator(bamf, ribo_featdb)
