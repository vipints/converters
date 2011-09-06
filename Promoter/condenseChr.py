#!/usr/bin/python
"""
Condense Contig locations in a BED file to a bin of 500 nts location on a single chromosome. The score of that bin will be the maximum score obtained in the 500 nts bin.

Usage:
condenseChr.py -f <STDIN of a BED file> 
"""
import sys 

if __name__=="__main__":

    try:
        if sys.argv[1] == '-f':
            wigfh = sys.stdin
    except:
        print 'Incorrect argument supplied'
        print __doc__
        sys.exit()

    win_score = []
    wind_cnt, condenseStart, LastPos = 0, 0, 0
    for line in wigfh:
        line = line.strip('\n\r').split('\t')
        line[3] = round(float(line[3]), 4)
        win_score.append(line[3])
        if wind_cnt == 499: # change the nucleotide count according to the bin size 
            win_score.sort()
            dense_line = [line[0],
                    str(condenseStart),
                    str(condenseStart+1),
                    str(win_score[-1]),
                    line[4]
                    ]
            print '\t'.join(dense_line)
            condenseStart +=500 # change the begnining line position according to the bin size.
            wind_cnt = 0
            win_score = []
            continue
        wind_cnt += 1
        LastPos = int(line[1])
    wigfh.close()
   
    # to capture the final bin 
    if wind_cnt != 0:
        win_score.sort()
        dense_line = [line[0],
                    str(LastPos),
                    str(LastPos+1),
                    str(win_score[-1]),
                    line[4]
                    ]
        print '\t'.join(dense_line)
