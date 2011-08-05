#!/usr/bin/python
# Count the number of occurrences of a motif in drosophila introns near the splice sites (distance less than 200nt)
# The motif is YGCY where at least one Y has to be a T. [TGCC, TGCT, CGCT]

import sys, re 
from Bio import SeqIO

def intDB(intron_file):

    intdb = dict()
    fh = open(intron_file, 'rU')
    for line in fh:
        line = line.strip('\n\r').split('\t')
        if re.match(r'#', line[0]) or re.match(r'>', line[0]):continue
        if len(line) == 1:continue
        if '' in line:continue
        assert len(line) == 9, '\t'.join(line)
        if line[2] == 'intron': 
            if line[0] in intdb:
                intdb[line[0]].append((int(line[3]), int(line[4]), line[6]))
            else:
                intdb[line[0]] = [(int(line[3]), int(line[4]), line[6])]
    fh.close()
    return intdb

def getSeq(gseq, cid, features):

    mcnt = 0 
    fah = open(gseq, 'rU')
    for rec in SeqIO.parse(fah, 'fasta'):
        if rec.id == cid:
            for fp in features:
                motif_seq = rec.seq[fp[0]-1:fp[1]]
                motif_seq = motif_seq.tomutable() 
                motif_seq = str(motif_seq)
                if len(motif_seq) < 200:
                    imseq = motif_seq.split('TGCC')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1 
                    imseq = motif_seq.split('TGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1 
                    imseq = motif_seq.split('CGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1 
                else:
                    dsp = motif_seq[0:200]
                    asp = motif_seq[len(motif_seq)-200:]
                    imseq = dsp.split('TGCC')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    imseq = dsp.split('TGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    imseq = dsp.split('CGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    imseq = asp.split('TGCC')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    imseq = asp.split('TGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    imseq = asp.split('CGCT')
                    if len(imseq) != 1:
                        mcnt += len(imseq)-1
                    
    fah.close()
    return mcnt

def __main__():
    
    try:
        intron_info = sys.argv[1]
        fasta_seq = sys.argv[2]
    except:
        sys.stderr.write("Provide Intron information and FASTA file with genome sequence \n")
        sys.exit(-1)
    
    introndb = intDB(intron_info)
    for contig, feature in introndb.items():
        motif_cnter = getSeq(fasta_seq, contig, feature)
        print contig, len(feature), motif_cnter

if __name__ == '__main__':__main__()
