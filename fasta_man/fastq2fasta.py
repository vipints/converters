#!/usr/bin/env python
"""
Create a FASTA file from FASTQ
Usage:fastq2fasta.py in.fastq out.fasta
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def __main__():
    try:
        fq_in=open(sys.argv[1], "rU")
        fa_out=open(sys.argv[2], 'w+')
    except:
        print __doc__
        sys.exit(-1)
    for rec in SeqIO.parse(fq_in, 'fastq'):
        SeqIO.write([rec], fa_out, "fasta") 
    fq_in.close()
    fa_out.close()

if __name__=="__main__": __main__()
