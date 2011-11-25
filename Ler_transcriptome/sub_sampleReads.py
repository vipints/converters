#!/usr/bin/python 
"""
Take a sub sample of reads from FASTQ file. 
"""

import re, sys, random
from Bio import SeqIO

def __main__():

    try:
        fastqFile = sys.argv[1]
        subCnt = int(sys.argv[2]) # sub sample count
    except:
        sys.stderr.write('Access denied for a FASTQ format file !\n')
        sys.exit(-1)
    
    if fastqFile == '-':
        ffh = sys.stdin
    else:
        ffh = open(fastqFile)
    read_cnt = 0 
    for rec in SeqIO.parse(ffh, 'fastq'):
        read_cnt += 1
    ffh.close()
    print 'Number of reads in FASTQ: ', read_cnt
    assert subCnt <= read_cnt, str(subCnt) + ' (sub-sample count) should be less than total read count ' + str(read_cnt)

    try:
        accept_prob = (1.0*subCnt)/read_cnt
    except:
        accept_prob = 1
    #print accept_prob

    cnt, sub_cnt = 0, 0
    subFile = open(str(subCnt)+ '_N_' + fastqFile, 'w+')
    ffh = open(fastqFile)
    for rec in SeqIO.parse(ffh, 'fastq'):
        rnb = random.random()
        cnt += 1
        if rnb <= accept_prob:
            sub_cnt += 1 
            subFile.write(rec.format("fastq"))
        if subCnt == sub_cnt:break
    ffh.close()
    subFile.close()

    print 'Number of reads scanned: ', cnt
    print 'Number of reads in: ', sub_cnt

if __name__ == "__main__":__main__()
