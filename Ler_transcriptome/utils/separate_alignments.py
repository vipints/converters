#!/usr/bin/python 
"""This program used to separate the read alignments according to the feature type.
Usage: python separate_alignments.py in.bam STDOUT SAM file 
"""
import re, sys 
import pysam 

def BAMParse(fname):
    """filter BAM file 
    """
    bfh = pysam.Samfile(fname, 'rb') 
    outsam = pysam.Samfile('-', 'w', template = bfh)
    for rec in bfh.fetch():
        if len(rec.seq)==40:
            if rec.qlen==40:
                outsam.write(rec)
        elif len(rec.seq)==42:
            if rec.qlen==42:
                outsam.write(rec)
    bfh.close()
    outsam.close()

if __name__=="__main__":
    try:
        fbam = sys.argv[1]
    except:
        print 'Incorrect number of arguments provided:'
        print __doc__
        sys.exit(-1)
    BAMParse(fbam)
