#!/usr/bin/env python 
"""
Convert a GFF and associated FASTA file into GenBank format.

Usage: gff_to_genbank.py in.gff in.fasta out.gb
"""
import sys, os 
from Core import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def __main__():
    try:
        gff_fname = sys.argv[1]
        fasta_fname = sys.argv[2]
        gb_fname = sys.argv[3]
    except: 
        print __doc__
        sys.exit(-1)
    fasta_rec = SeqIO.to_dict(SeqIO.parse(fasta_fname, "fasta", generic_dna))
    gff_rec = GFF.parse(gff_fname, fasta_rec)
    SeqIO.write(gff_rec, gb_fname, "genbank")

if __name__=="__main__":__main__()
