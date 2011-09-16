#!/usr/bin/python
"""Program to find the reads mapped in Yeast and Human genome. 
Usage: 
shared_reads.py <YEAST BAM> <HUMAN BAM>
"""
import sys, re, os

def BAMreader(bam_file, reads):

    bfh = os.popen('/fml/ag-raetsch/share/software/samtools//samtools view '+ bam_file)
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        reads[line[0]]=1
    bfh.close()
    return reads

def BAMWriter(bam_file, get_out, org):
    
    bfh = os.popen('/fml/ag-raetsch/share/software/samtools//samtools view '+ bam_file)
    fname = re.search(r'(.*)\.bam$', bam_file).group(1)
    samfh = open(fname+'_uniq.sam', 'w') 
    for line in bfh:
        line = line.strip('\n\r').split('\t')
        if not line[0] in get_out:
            samfh.write('\t'.join(line)+'\n')
    samfh.close()
    bfh.close()
    if org == 'yeast':
        os.system("/fml/ag-raetsch/share/software/samtools//samtools import /fml/ag-raetsch/nobackup2/projects/sequencing_runs/yeast/reads/SR_2010/genome/yeast_genome.fa.fai " + fname + "_uniq.sam " + fname + "_uniq.bam")
    elif org == 'human':
        os.system("/fml/ag-raetsch/share/software/samtools//samtools import /fml/ag-raetsch/nobackup/projects/rgasp.2/genomes/human/hg19/hg19.fa.fai " + fname + "_uniq.sam " + fname + "_uniq.bam")
    os.system("rm -fr " + fname + "_uniq.sam")

if __name__=="__main__":
    try:
        yeast_bam = sys.argv[1]
        human_bam = sys.argv[2]
    except:
        print 'Incomplete arguments:'
        print __doc__
        sys.exit(-1)
    common_reads = dict()
    common_reads_yeast = BAMreader(yeast_bam, common_reads)
    common_reads = dict()
    common_reads_human = BAMreader(human_bam, common_reads)
    both_cnt = 0
    com = dict()
    for rid in common_reads_yeast:
        if rid in common_reads_human:
            both_cnt += 1
            com[rid] = 1
    print yeast_bam, human_bam
    print both_cnt
    # get out the common reads in the genome
    BAMWriter(yeast_bam, com, 'yeast')
    BAMWriter(human_bam, com, 'human')
