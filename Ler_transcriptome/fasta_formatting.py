#!/usr/bin/python
# Formatting FASTA file, it takes a unformatted FASTA file and results a Nice formatted FASTA file.

import re, sys, os 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

def __main__():
    
    try:
        fasta_raw = sys.argv[1]
    except: 
        sys.stderr.write('Access denied for the input FASTA file ! terminating\n')
        sys.exit(-1)
    fah = open(fasta_raw, 'rU')
    for rec in SeqIO.parse(fah, 'fasta'):
        rna_genome = rec.seq.tomutable()
        i = 0 
        for nt in rna_genome:
            if nt == 'U':
                rna_genome[i] = 'T'
            i += 1
        rna_genome = 100*'N' + rna_genome + 100*'N' ## adding a non-nucleotide clips at both ends of the sequence.

        gseq = SeqRecord(rna_genome, id=rec.id, description='A_thaliana_ribosomal_RNA')
        print gseq.format("fasta")
        #break
    fah.close()

if __name__=="__main__":__main__()
