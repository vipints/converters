#!/usr/bin/python
"""
Create a matrix representation of reads aligned to the features in the genome.

Columns are:
[ID] [No. of alignments] [read orientation] [TE] [rRNA] [CG] [OtherRNA] etc.
1   10  -1  0   0   1

Strand Orientation are noted as 1/-1/0 decoded as +/-/cannot decide
Which Orientation the read mapped to the feature as  0/1/-1/2 decoded as no match/sense/anti-sense/sense+anti-sense

Usage:
matrix_gen.py in.gff in.bam out.mat
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
    for rid, rinfo in alignDB.items():
        read_strand=[]
        ribXrd, cgXrd, othXrd, teXrd, pdXrd = [], [], [], [], []
        for rdet in rinfo:
            read_strand.append(rdet[2]) # storing the strand information of each alignment 
            if rdet[0] in cg_featdb:
                for details, strand_info in cg_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        if rdet[2]==strand_info:
                            cgXrd.append(1)
                        else:
                            cgXrd.append(-1)
                        break
            if rdet[0] in oth_featdb:
                for details, strand_info in oth_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        if rdet[2]==strand_info:
                            othXrd.append(1)
                        else:
                            othXrd.append(-1)
                        break
            if rdet[0] in ribo_featdb:
                for details, strand_info in ribo_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        if rdet[2]==strand_info:
                            ribXrd.append(1)
                        else:
                            ribXrd.append(-1)
                        break
            if rdet[0] in te_featdb:
                for details, strand_info in te_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        if rdet[2]==strand_info:
                            teXrd.append(1)
                        else:
                            teXrd.append(-1)
                        break
            if rdet[0] in psd_featdb:
                for details, strand_info in psd_featdb[rdet[0]].items():
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                        if rdet[2]==strand_info:
                            pdXrd.append(1)
                        else:
                            pdXrd.append(-1)
                        break
        #print ribXrd, cgXrd, othXrd, teXrd, pdXrd
        no_amnt += len(rinfo)
        rd_cnt += 1
        ind_rd = [rd_cnt, len(rinfo)] # read id mapped to integer starting from 1 to ..
        read_strand = list(set(read_strand))
        if len(read_strand)==2:
            ind_rd.append(0)
        else:
            ind_rd.append(read_strand[0])
        for xq in [cgXrd, ribXrd, teXrd, othXrd, pdXrd]: # strand information are processed for matrix creation
            xq = list(set(xq))
            if len(xq)==0:
                ind_rd.append(0) # No match in the feature 
            elif len(xq)==2:
                ind_rd.append(2) # Sense + anti-sense
            else:
                ind_rd.append(xq[0]) # Sense / anti-sense
        ind_rd = tuple(ind_rd)
        read_dist.append(ind_rd)
        if (no_amnt%1000)==0:print no_amnt
    return read_dist

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
        #samf = sys.argv[3]
        mat_file = sys.argv[3]
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
    """print 'Writing SAM file with Intergenic reads...'
    samfile = pysam.Samfile(bamf, 'rb')
    outsam = pysam.Samfile(samf, 'w', template = samfile)
    for rec in samfile.fetch():
        if rec.qname in read_distribution: 
            outsam.write(rec)
    outsam.close()
    samfile.close()
    print '...SAM file ready', samf"""
    # re-arrange and write into a matrix form 
    cmp_type=array(read_distribution)
    sio.savemat(mat_file, {'rdist':cmp_type})
    #TODO Saving intergenic reads, Creating mat file 
    #print 'Alignments were saved in matrix form, finished successfully !'
    print time.asctime( time.localtime(time.time()) )
