#!/usr/bin/python
"""
Create a Table in which decide the read is aligned to Sense or Antisense or undecide orientation.

ex: Table representation
Library name    feature:# of reads  Sense   Antisense   Doubt
Ler-1_BR-1  591117  574892  5454    10771

Usage:
samtools view rRNA_alignment.bam | FeatureSummary.py -f ~/tmp/Ath/rRNA_orientation.txt
"""

def infoFeature(feature):

    featdb = dict()
    fh = open(feature, 'rU') 
    for line in fh:
        line = line.strip('\n\r').split('\t')
        featdb[line[0]] = line[1]
    fh.close()
    return featdb

def AlignGenerator(fh):

    infodb = dict()
    for line in fh:
        line = line.strip('\n\r').split('\t')
        if len(line) < 12 :continue 
        NM = re.search(r'NM:i:(.+)', line[11]).group(1)
        orient = None
        for field in line:
            if re.search(r'XS:A:*', field):
                orient = re.search(r'XS:A:(.+)', field).group(1)
                break
        if line[0] in infodb:
            infodb[line[0]].append((orient, int(NM), line[2]))
        else:
            infodb[line[0]] = [(orient, int(NM), line[2])]
    fh.close()
    return infodb
    
def getSmallMMIndex(MM):
    
    previous, first_iter = 0, 0
    for item in MM:
        if first_iter == 0:
            previous = item
            first_iter = 1
            continue
        if item < previous:
            previous = item
    return MM.index(previous) 

def StrandDecision(alignDB, featdb):
    
    doubt, ant_sense_st, sense_st = 0, 0, 0
    for rid, rinfo in alignDB.items():
        strand, contig, MM = [], [], []
        for ele in rinfo:
            strand.append(ele[0])
            MM.append(ele[1])
            contig.append(ele[2])
        Cstrand = list(set(strand))
        if len(Cstrand) == 1: # SENSE/ANTI-SENSE
            sns_fl, antsns_fl = 0, 0 
            for cid in contig:
                if strand[0] == featdb[cid]:
                    sns_fl = 1
                else:
                    antsns_fl = 1 
            if sns_fl == 1 and antsns_fl == 1:
                if len(set(MM))==1:
                    doubt += 1   
                else: # get the minimum mis-match orientation 
                    sidx = getSmallMMIndex(MM)
                    if strand[sidx] == featdb[contig[sidx]]:
                        sense_st += 1
                    elif strand[sidx] != featdb[contig[sidx]]:
                        ant_sense_st += 1
                    else:
                        doubt +=1 
            elif sns_fl == 1:
                sense_st += 1
            elif antsns_fl == 1:
                ant_sense_st += 1 
        else: # Undecide 
            if len(set(MM))==1:
                doubt += 1   
            else:
                sidx = getSmallMMIndex(MM)
                if strand[sidx] == featdb[contig[sidx]]:
                    sense_st += 1
                elif strand[sidx] != featdb[contig[sidx]]:
                    ant_sense_st += 1
                else:
                    doubt += 1   
    return doubt, ant_sense_st, sense_st

import sys, re

if __name__ == "__main__":

    try:
        if sys.argv[1] == '-f':
            bfh = sys.stdin
        ref_genome = sys.argv[2]
    except:
        print 'Incorrect number of arguments:'
        print __doc__
        sys.exit(-1)

    # get the feature orientations
    featdb = infoFeature(ref_genome)
    # get alignment data 
    alignDB = AlignGenerator(bfh)
    # figure out strandness of read and feature 
    doubt, ant_sense_st, sense_st = StrandDecision(alignDB, featdb)
    xls_line = [str(len(alignDB)),
            str(sense_st),
            str(ant_sense_st),
            str(doubt)
            ]
    print '\t'.join(xls_line)

