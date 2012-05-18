#!/usr/bin/python 
"""
Write a BAM file content to a plain text format. 
Usage: python bam_to_txt.py in.bam > out.txt  
"""
import re, sys 
import pysam 

def AnnoJ_BAM(bam_read):
    for rec in sam_reader.fetch():
        orient=None
        for atb in rec.tags:
            if atb[0]=='XS':
                orient=atb[1]
        print rec.qname, orient, rec.pos, rec.qlen, sam_reader.getrname(rec.rname), rec.query

def sort_BAM_TAGS(bam_read):
    """
    """
    for line in sys.stdin:
        print line.strip('\n\r')
    outsam=pysam.Samfile("-", "w", template = bam_read)
    for rec in bam_read.fetch():
        new_cont=rec.tags
        new_cont.insert(0, rec.tags[-1])
        new_cont.pop()
        rec.tags=new_cont
        outsam.write(rec)
    outsam.close()

def __main__():

    try:
        bam_file=sys.argv[1]
    except:
        print __doc__
        sys.exit(1)
    sam_reader = pysam.Samfile(bam_file, "rb")
    sort_BAM_TAGS(sam_reader)
    sam_reader.close()

if __name__=="__main__":__main__()
