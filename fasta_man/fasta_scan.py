#!/usr/bin/python
"""Gives an overall statistics on contigs in a FASTA file
USAGE: feature_scan.py in.fasta
"""
from Bio import SeqIO
import re, sys 
from operator import itemgetter

def __main__():

    try:
        fa_name = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    seq_info = dict()
    fah = open(fa_name)
    for rec in SeqIO.parse(fah, "fasta"):
        seq_info[rec.id] = len(rec.seq)
    fah.close()
    print '# of FASTA entries : ', len(seq_info)
    for long_one in sorted(seq_info.items(), key=itemgetter(1), reverse=True):
        print 'Long feature (bp): ', long_one[1]
        break
    for short_one in sorted(seq_info.items(), key=itemgetter(1)):
        print 'Short feature (bp): ', short_one[1]
        break
    flength = 0 
    for ele in sorted(seq_info.items(), key=itemgetter(1)):
        flength += ele[1]
    print 'Average length of FASTA entry (bp): ', (flength/len(seq_info))

if __name__=="__main__":__main__()
