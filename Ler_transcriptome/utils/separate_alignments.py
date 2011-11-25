#!/usr/bin/python 
"""This program used to separate the read alignments according to the feature type.
Usage:
samtools view BAM | python separate_alignments.py -f <feature annotation in GFF>
"""
import re, sys 

def ParseAnno(feat_file):
    
    candidates = dict()
    afh = open(feat_file)
    for line in afh:
        line = line.strip('\n\r').split('\t')
        if re.search(r'^#', line[0]):
            continue
        fid = None
        if line[2] == 'pseudogene':
            for ele in line[-1].split(';'):
                if re.search(r'ID=', ele):
                    fid = re.search(r'ID=(.+)', ele).group(1)
                    break
            candidates[fid]=1
    return candidates
    afh.close()

def BAMParse(bfh, featdb):
    
    PSDfh = open('Psedo.sam', 'w')
    STRfh = open('Oth.sam', 'w')
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        if line[2] in featdb:
            PSDfh.write('\t'.join(line)+'\n')
        else:
            STRfh.write('\t'.join(line)+'\n')
    PSDfh.close()
    STRfh.close()
    bfh.close()

if __name__=="__main__":
    try:
        if sys.argv[1] == '-f':
            bfh = sys.stdin
        anno_file = sys.argv[2]
    except:
        print 'Incorrect number of arguments provided:'
        print __doc__
        sys.exit(-1)
    
    featdb = ParseAnno(anno_file)
    for fid in featdb:
        print fid
        break
    BAMParse(bfh, featdb)
    #print len(featdb)
