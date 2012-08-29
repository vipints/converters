#!/usr/bin/env python 
"""
Objective:
Best rRNA fragments filtered from the BLAT alignment of genbank entries and silvia rRNA database. 
Genome sequence of A. thaliana undergone a repeat masking step with RepeatMasker program and the types of 
repeats are identified along the genome. 
Based on the presence of rRNA fragment, TE family members are classified as rRNA+ and rRNA- group. 
Look at the right and left region of each group of members and see whether we can clearly distinguish any 
difference in repeat pattern which we can connect to the slippage or transposition process.

Usage:  TE family members with location on the genome, 
        Repeat type information from RepeatMasker result file
"""

import sys, re 

def __main__():
    try:
        tefname = sys.argv[1]
        repeatfname = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)

if __name__ == "__main__":
    __main__()
