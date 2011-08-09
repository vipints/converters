#!/usr/bin/python
"""
Find the reads mapped to rRNA, Coding region, TE and other structural RNA regions, both individually and in a shared way to plot as Venn diagram.
"""
import re, sys, os

def AlignReader(fname):
    
    fbh = os.popen('/fml/ag-raetsch/share/software/samtools//samtools view ' + fname)
    for line in fbh:
        line = line.strip('\r\n').split('\t')
        yield line[0]
    fbh.close()

def __main__():

    try:
        ribo_bam = sys.argv[1]
        cg_bam = sys.argv[2]
        te_bam = sys.argv[3]
    except:
        sys.stderr.write('Access denied for BAM files \n')
        sys.exit(-1)
    
    Generator = AlignReader(ribo_bam)
    ribo_read = dict()
    for rid in Generator:
        ribo_read[rid] = 1
    cg_read = dict()
    Generator = AlignReader(cg_bam)
    for rid in Generator:
        cg_read[rid] = 1
    te_read = dict() 
    Generator = AlignReader(te_bam)
    for rid in Generator:
        te_read[rid] = 1

    te_reads = [fid for fid in te_read]
    te_reads = set(tuple(te_reads))
    cg_reads = [fid for fid in cg_read]
    cg_reads = set(tuple(cg_reads))
    ribo_reads = [fid for fid in ribo_read]
    ribo_reads = set(tuple(ribo_reads))

    te_cg_cnt = len(te_reads & cg_reads)
    te_ribo_cnt = len(te_reads & ribo_reads)
    cg_ribo_cnt = len(cg_reads & ribo_reads)
    te_ribo_cg_cnt =  len(cg_reads & ribo_reads & te_reads)

    print 'Individual counts: rRNA, TE, Coding genes\t', len(ribo_reads), len(te_reads), len(cg_reads)
    print 'TE and Coding genes\t', te_cg_cnt
    print 'TE and rRNA\t', te_ribo_cnt
    print 'Coding gene and rRNA\t', cg_ribo_cnt
    print 'TE, Coding gene, rRNA\t', te_ribo_cg_cnt

    print 'TE bin: ', len(te_reads) - (te_ribo_cnt + te_cg_cnt + te_ribo_cg_cnt)
    print 'rRNA bin: ', len(ribo_reads) - (te_ribo_cnt + cg_ribo_cnt + te_ribo_cg_cnt)
    print 'Coding gene: ', len(cg_reads) - (te_cg_cnt + cg_ribo_cnt + te_ribo_cg_cnt)

if __name__=="__main__":__main__()
