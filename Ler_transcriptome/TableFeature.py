#!/usr/bin/python
"""
Create a matrix representation of reads aligned to the features in the genome.

Columns are:
[ID] [No. of alignments] [read orientation] [TE] [rRNA] [CG] [OtherRNA] etc.
1   10  -1  0   0   1

Strand Orientation are noted as 1/-1/0 decoded as +/-/cannot decide
Which Orientation the read mapped to the feature as  0/1/-1 decoded as no match/sense/anti-sense

Usage:
samtools view BAM.bam | TableFeature.py -f <anno db> <mat file>
"""

def get_Feature(fname):
    """Extract genome annotation information from provided GFF file.
    """

    te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb = dict(), dict(), dict(), dict(), dict()
    from Core import GFF
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrM', 'ChrC']: ## change according to the GFF file TODO According to the file type add chromosome number automatically.
        print cid
        fh = open(fname, 'rU')
        limit_info = dict(gff_id=[cid], gff_source=['TAIR10']) ## change the source TODO According to the file type add source flag automatically 
        for rec in GFF.parse(fh, limit_info=limit_info):
            for each_rec in rec.features:
                if each_rec.type=='gene':
                    if each_rec.sub_features[0].type=='mRNA':
                        if cid in cg_featdb:
                            cg_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            cg_featdb[cid]=[(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='miRNA':
                        if cid in oth_featdb:
                            oth_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            oth_featdb[cid]=[(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='ncRNA':
                        if cid in oth_featdb:
                            oth_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            oth_featdb[cid]=[(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='snoRNA':
                        if cid in oth_featdb:
                            oth_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            oth_featdb[cid]=[(int(each_rec.location._start.position),
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='snRNA':
                        if cid in oth_featdb:
                            oth_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            oth_featdb[cid]=[(int(each_rec.location._start.position),
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='tRNA':
                        if cid in oth_featdb:
                            oth_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position),
                                    each_rec.strand))
                        else:
                            oth_featdb[cid]=[(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                    if each_rec.sub_features[0].type=='rRNA':
                        if cid in ribo_featdb:
                            ribo_featdb[cid].append((int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand))
                        else:
                            ribo_featdb[cid]=[(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position), 
                                    each_rec.strand)]
                elif each_rec.type=='pseudogene':
                    if cid in psd_featdb:
                        psd_featdb[cid].append((int(each_rec.location._start.position), 
                                int(each_rec.location._end.position), 
                                each_rec.strand))
                    else:
                        psd_featdb[cid]=[(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position), 
                                each_rec.strand)]
                elif each_rec.type=='transposable_element' or each_rec.type=='transposable_element_gene':
                    if cid in te_featdb:
                        te_featdb[cid].append((int(each_rec.location._start.position), 
                                int(each_rec.location._end.position), 
                                each_rec.strand))
                    else:
                        te_featdb[cid]=[(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position), 
                                each_rec.strand)]
        fh.close()
    return te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb

def AlignGenerator(fh):

    infodb = dict()
    strand_map={'+' : 1, '-' : -1}
    for line in fh:
        line = line.strip('\n\r').split('\t')
        if len(line) < 12 :continue 
        NM = re.search(r'NM:i:(.+)', line[11]).group(1)
        orient = None
        for field in line:
            if re.search(r'XS:A:*', field):
                orient = re.search(r'XS:A:(.+)', field).group(1)
                break
        if line[0] in infodb:
            infodb[line[0]].append((line[2], int(line[3]), strand_map[orient], int(NM)))
        else:
            infodb[line[0]] = [(line[2], int(line[3]), strand_map[orient], int(NM))]
    fh.close()
    return infodb
    

def readScaner(alignDB, cg_featdb, oth_featdb, ribo_featdb, te_featdb, psd_featdb):
    
    read_dist, rd_cnt, no_amnt = [], 1, 0
    cnt_cg, cnt_te, cnt_ribo, cnt_pd, cnt_oth, cnt_int = 0, 0, 0, 0, 0, 0 
    for rid, rinfo in alignDB.items():
        ind_rd = []
        ind_rd.append(rd_cnt) # read id mapped to intger starting from 1 to ..
        rd_cnt += 1
        ind_rd.append(len(rinfo)) # number of alignments 
        read_strand=[]
        ribXrd, cgXrd, othXrd, teXrd, pdXrd = 0, 0, 0, 0, 0
        ribX, cgX, othX, teX, pdX = 0, 0, 0, 0, 0
        no_amnt += len(rinfo)
        for rdet in rinfo:
            read_strand.append(rdet[2])
            if rdet[0] in cg_featdb:
                for details in cg_featdb[rdet[0]]:
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+95:
                        cgX += 1
                        if rdet[2]==details[2]:
                            cgXrd=1
                        else:
                            cgXrd=-1
                        break
            if rdet[0] in oth_featdb:
                for details in oth_featdb[rdet[0]]:
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+95:
                        othX += 1 
                        if rdet[2]==details[2]:
                            othXrd=1
                        else:
                            othXrd=-1
                        break
            if rdet[0] in ribo_featdb:
                for details in ribo_featdb[rdet[0]]:
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+95:
                        ribX +=1 
                        if rdet[2]==details[2]:
                            ribXrd=1
                        else:
                            ribXrd=-1
                        break
            if rdet[0] in te_featdb:
                for details in te_featdb[rdet[0]]:
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+95:
                        teX += 1
                        if rdet[2]==details[2]:
                            teXrd=1
                        else:
                            teXrd=-1
                        break
            if rdet[0] in psd_featdb:
                for details in psd_featdb[rdet[0]]:
                    if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+95:
                        pdX += 1
                        if rdet[2]==details[2]:
                            pdXrd=1
                        else:
                            pdXrd=-1
                        break
        read_strand = list(set(read_strand))
        if len(read_strand)==2:
            ind_rd.append(0)
        else:
            ind_rd.append(read_strand[0])
        ind_rd.append(cgXrd)
        ind_rd.append(ribXrd)
        ind_rd.append(teXrd)
        ind_rd.append(othXrd)
        ind_rd.append(pdXrd)
        if cgXrd != 0:  ## counting based on number of alignments of a read. 
            cnt_cg += cgX
        if ribXrd != 0:
            cnt_ribo += ribX
        if teXrd != 0:
            cnt_te += teX 
        if othXrd != 0:
            cnt_oth += othX
        if pdXrd !=0:
            cnt_pd += pdX
        if ribXrd ==0 and cgXrd ==0 and teXrd ==0 and othXrd == 0 and pdXrd ==0:
            cnt_int += len(rinfo)
            print rid, rinfo
        ind_rd = tuple(ind_rd)
        read_dist.append(ind_rd)

    pline = [str(rd_cnt-1),
            str(no_amnt),
            str(cnt_cg),
            str(cnt_ribo),
            str(cnt_te),
            str(cnt_oth),
            str(cnt_pd),
            str(cnt_int)]
    print '\t'.join(pline)
    return read_dist

import sys, re
from numpy import *
import scipy.io as sio

if __name__ == "__main__":

    try:
        if sys.argv[1]=='-f': 
            bfh = sys.stdin
        anno_file = sys.argv[2]
        mat_file = sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
    # get the feature orientations
    te_db, psd_db, cg_db, oth_db, ribo_db=get_Feature(anno_file)
    print 'Annotation file parsed successfully...'
    # get read alignment data 
    alignDB = AlignGenerator(bfh)
    print 'Read alignments are stored in internal db...'
    # figure out each reads location on the genome 
    read_distribution = readScaner(alignDB, cg_db, oth_db, ribo_db, te_db, psd_db)
    print 'checked the read distribution'
    print 'started writing to mat file'
    # re-arrange and write into a matrix form 
    cmp_type=array(read_distribution)
    sio.savemat(mat_file, {'rdist':cmp_type})
