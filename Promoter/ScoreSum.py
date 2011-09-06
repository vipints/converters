#!/usr/bin/python
"""
Normalize the promoter prediction score to [0-1] interval on each chromosome. The resulting BED file will be used in pppBenchmark program to evaluate the promoter prediction status.

Usage:
ScoreSum.py <BED file>
"""

import re, sys 
from operator import itemgetter

def BEDreader(fname):
    
    bed_score = []
    bfh = open(fname)
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        assert len(line) == 5, '\t'.join(line)
        bed_score.append(float(line[3]))
    bfh.close()
    return bed_score

def __main__():
    
    try:
        plus_fname = sys.argv[1]
    except:
        print 'Incorrect Argument supplied'
        print __doc__
        sys.exit(-1)

    score = BEDreader(plus_fname)
    min_score, max_score = None, None
    score.sort()
    for fs in score:
        min_score=fs
        if min_score != -42.0:
            break
    max_score = score[-1]
    bfh = open(plus_fname)
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        assert len(line) == 5, '\t'.join(line)
        if float(line[3]) == -42.0:
            line[3] = min_score
        norm_score=(float(line[3])-min_score)/(max_score-min_score)
        bline = [line[0],
            line[1],
            line[2],
            str(round(norm_score, 4)),
            line[-1]
            ]
        print '\t'.join(bline)
    bfh.close()

if __name__=="__main__":__main__()
