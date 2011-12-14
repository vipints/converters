#!/usr/bin/python
"""Program to read a binary file.
Usage: read_binary.py in.bin 
"""
import sys, struct

def read_int(fn):
    """
    Read integer from a binary file.
    """
    region = []
    fhd = open(fn, 'rb') 
    while True:
        rec = fhd.read(4)
        if len(rec) != 4:
            break
        (pos,) = struct.unpack('i', rec)
        region.append(pos) 
    fhd.close() 
    return region

if __name__=="__main__":
    try:
        fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    locations = read_int(fname)
    print locations[0], locations[-1]
