#!/usr/bin/python
"""
Program to Count the Number of alignments of reads aligned to each annotated features in the genome and Create a SAM file with reads mapped to inter-genic region.

Usage:
TableFeature.py in.bam in.gff out.sam
"""

def get_Feature(fname):
    """Extract genome annotation information from a provided GFF file.
    """
    from Core import GFF
    te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb = dict(), dict(), dict(), dict(), dict()
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrM', 'ChrC']: ## change according to the GFF file TODO According to the file type add chromosome number automatically.
        print cid
        fh = open(fname, 'rU')
        limit_info = dict(gff_id=[cid], gff_source=['TAIR10']) ## change the source TODO According to the file type add source flag automatically 
        for rec in GFF.parse(fh, limit_info=limit_info):
            for each_rec in rec.features:
                if each_rec.type=='gene':
                    for child in each_rec.sub_features:
                        if child.type=='mRNA':
                            if cid in cg_featdb:
                                cg_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=each_rec.strand
                            else:
                                cg_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)):
                                    each_rec.strand}
                        elif each_rec.sub_features[0].type in ['miRNA', 'ncRNA', 'snoRNA', 'snRNA', 'tRNA']:
                            if cid in oth_featdb:
                                oth_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=each_rec.strand
                            else:
                                oth_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)):
                                    each_rec.strand}
                        elif each_rec.sub_features[0].type=='rRNA':
                            if cid in ribo_featdb:
                                ribo_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=each_rec.strand
                            else:
                                ribo_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)): 
                                    each_rec.strand}
                elif each_rec.type=='pseudogene':
                    if cid in psd_featdb:
                        psd_featdb[cid][(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=each_rec.strand
                    else:
                        psd_featdb[cid]={(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                each_rec.strand}
                elif each_rec.type in ['transposable_element', 'transposable_element_gene']:
                    if cid in te_featdb:
                        te_featdb[cid][(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=each_rec.strand
                    else:
                        te_featdb[cid]={(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                each_rec.strand}
        fh.close()
    return te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb

def AlignGenerator(fbam):
    """Parsing alignment file, extracting the information.
    """
    samdb = dict()
    smap={'+' : 1, '-' : -1}
    samfile = pysam.Samfile(fbam, 'rb') 
    for rec in samfile.fetch():
        for cont in rec.tags:
            if cont[0]=='NM':NM=cont[1]
            if cont[0]=='XS':
                orient=cont[1]
                break
        if rec.qname in samdb:
            samdb[rec.qname].append((samfile.getrname(rec.rname), rec.pos, smap[orient], NM))
        else:
            samdb[rec.qname]=[(samfile.getrname(rec.rname), rec.pos, smap[orient], NM)]
    samfile.close()
    return samdb
    
def readScaner(alignDB, cg_featdb, oth_featdb, ribo_featdb, te_featdb, psd_featdb):
    """Find the location of reads in the genome especially in the annotated region.
    """
    read_dist, rd_cnt, no_amnt = [], 0, 0
    cnt_share, cnt_cg, cnt_te, cnt_ribo, cnt_pd, cnt_oth, cnt_int = 0, 0, 0, 0, 0, 0, 0
    shared_reads, inter_gene_reads, CG_reads = dict(), dict(), dict()
    #Oth_reads, Ribo_reads, TE_reads, Psd_reads = dict(), dict(), dict(), dict()
    for rid, rinfo in alignDB.items():
        #read_strand=[]
        ribXrd, cgXrd, othXrd, teXrd, pdXrd = 0, 0, 0, 0, 0
        ribX, cgX, othX, teX, pdX = 0, 0, 0, 0, 0
        for rdet in rinfo:
            #read_strand.append(rdet[2]) # storing the strand information of each alignment 
            if rdet[0] in cg_featdb:
                for details, strand_info in cg_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        cgX+=1
                        if rdet[2]==strand_info:
                            cgXrd=1
                        else:
                            cgXrd=-1
                        break
            if rdet[0] in oth_featdb:
                for details, strand_info in oth_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        othX+=1
                        if rdet[2]==strand_info:
                            othXrd=1
                        else:
                            othXrd=-1
                        break
            if rdet[0] in ribo_featdb:
                for details, strand_info in ribo_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        ribX+=1
                        if rdet[2]==strand_info:
                            ribXrd=1
                        else:
                            ribXrd=-1
                        break
            if rdet[0] in te_featdb:
                for details, strand_info in te_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        teX+=1
                        if rdet[2]==strand_info:
                            teXrd=1
                        else:
                            teXrd=-1
                        break
            if rdet[0] in psd_featdb:
                for details, strand_info in psd_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        pdX+=1
                        if rdet[2]==strand_info:
                            pdXrd=1
                        else:
                            pdXrd=-1
                        break
        #print ribXrd, cgXrd, othXrd, teXrd, pdXrd
        no_amnt += len(rinfo)
        if (no_amnt%1000)==0:print no_amnt
        # Counting starts from here ..
        if len(set([ribXrd, cgXrd, othXrd, teXrd, pdXrd]))==1: # intergenic read 
            inter_gene_reads[rid]=1
            cnt_int += len(rinfo)
        elif [ribXrd, cgXrd, othXrd, teXrd, pdXrd].count(0)<=3: # multiple alignments spanning across different features
            #shared_reads[rid]=1
            cnt_share +=len(rinfo)
        elif cgXrd !=0:
            #CG_reads[rid]=1
            cnt_cg+=cgX
            if cgX != len(rinfo):
                cnt_int +=(len(rinfo)-cgX)
        elif othXrd !=0:
            #Oth_reads[rid]=1
            cnt_oth+=othX
            if othX != len(rinfo):
                cnt_int +=(len(rinfo)-othX)
        elif teXrd !=0:
            #TE_reads[rid]=1
            cnt_te+=teX
            if teX != len(rinfo):
                cnt_int +=(len(rinfo)-teX)
        elif pdXrd !=0:
            #Psd_reads[rid]=1
            cnt_pd+=pdX
            if pdX != len(rinfo):
                cnt_int +=(len(rinfo)-pdX)
        elif ribXrd !=0:
            #Ribo_reads[rid]=1
            cnt_ribo+=ribX
            if ribX != len(rinfo):
                cnt_int+=(len(rinfo)-ribX)

    pline = ['Nr. of Alignments',
            str(len(alignDB)), 
            str(no_amnt), 
            str(cnt_cg), 
            str(cnt_pd), 
            str(cnt_oth), 
            str(cnt_te), 
            str(cnt_ribo), 
            str(cnt_int), 
            str(cnt_share)
            ]
    print '\t'.join(pline)
    return inter_gene_reads

import sys, re
from numpy import *
import scipy.io as sio
import pysam 

if __name__ == "__main__":

    import time
    print time.asctime( time.localtime(time.time()) )
    try:
        bamf = sys.argv[1]
        anno_file = sys.argv[2]
        samf = sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
    # get the feature orientations
    print 'Annotation file parsing started...'
    te_db, psd_db, cg_db, oth_db, ribo_db=get_Feature(anno_file)
    print '...Parsing completed successfully'
    # get read alignment data 
    print 'Getting read alignment information from BAM file...'
    read_db = AlignGenerator(bamf)
    print '...Alignments are stored in internal db'
    # figure out each reads location on the genome 
    print 'Figure out the read location on genome...'
    read_distribution = readScaner(read_db, cg_db, oth_db, ribo_db, te_db, psd_db)
    print '...Read coverage identified'
    print 'Writing SAM file with Intergenic reads...'
    samfile = pysam.Samfile(bamf, 'rb')
    outsam = pysam.Samfile(samf, 'w', template = samfile)
    for rec in samfile.fetch():
        if rec.qname in read_distribution: 
            outsam.write(rec)
    outsam.close()
    samfile.close()
    print '...SAM file ready', samf
    print time.asctime( time.localtime(time.time()) )
