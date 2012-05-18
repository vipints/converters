#!/usr/bin/python 
"""FASTQ file pre-processing, before alignment to check for any invalid 
read entries which were not filtered out.

"""
import re, sys
from Bio import SeqIO

outfq=open(sys.argv[2], 'w')
for rec in SeqIO.parse(sys.argv[1], "fastq"):
    if len(rec.seq)<90:
        continue
    SeqIO.write(rec, outfq, 'fastq')
outfq.close()
