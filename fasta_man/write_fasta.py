#!/usr/bin/env python
"""Write a FASTA sequence output for a specified chromosome.
Usage: write_fasta.py in.fasta chr_id 
"""
import re, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def __main__():
    try:
        inh = open(sys.argv[1], "rU")
        chr_id = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    seq_rec = dict()
    for rec in SeqIO.parse(inh, 'fasta'):
        if rec.id==chr_id: 
            seq_rec[rec.id]=rec
    fasta_out = open(chr + '_seq.fa', 'w')
    SeqIO.write([seq_rec[chr_id]], fasta_out, "fasta")
    fasta_out.close()

if __name__=="__main__": __main__()
