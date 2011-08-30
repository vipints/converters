#!/usr/bin/python
"""
Normalize the promoter prediction score to [0-1] interval on each chromosome.
"""

import re, sys 
from operator import itemgetter

def BEDreader(fname):
    
    bed_score = dict()
    bfh = open(fname)
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        assert len(line) == 5, '\t'.join(line)
        if (line[0], line[4]) in bed_score:
            bed_score[(line[0], line[4])][(int(line[1]), int(line[2]))] = float(line[3])
        else:
            bed_score[(line[0], line[4])] = {(int(line[1]), int(line[2])): float(line[3])}
    bfh.close()
    return bed_score

def __main__():
    
    try:
        plus_fname = sys.argv[1]
    except:
        sys.stderr.write('Access denied for BED file! \n')
        sys.exit(-1)

    score = BEDreader(plus_fname)
    for fid, fls in score.items():
        min_score, max_score = None, None
        for fp, fs in sorted(fls.items(), key=itemgetter(1), reverse=True):
            max_score=fs 
            break
        for fp, fs in sorted(fls.items(), key=itemgetter(1)):
            min_score=fs 
            if min_score != -42.0: 
                break
        for fp, fs in sorted(fls.items()):
            if fs == -42.0:
                fs=min_score
            norm_score=(fs-min_score)/(max_score-min_score)
            #print norm_score
            bline = [fid[0],
                    str(fp[0]),
                    str(fp[1]),
                    str(round(norm_score, 4)),
                    fid[1]
                    ]
            print '\t'.join(bline)

if __name__=="__main__":__main__()
