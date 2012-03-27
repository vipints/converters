#!/usr/bin/env python 
"""
Program to find the number of reads aligned to each TE element region in 
A. thaliana genome. 
Usage: fetch_BAM.py in.bam in.gff3
"""
import sys, re, pysam 

def get_alignment(bfname):
    """Read BAM file and count the number of alignments for each read.
    """
    aln_cnt=dict()
    bam_read = pysam.Samfile(bfname, 'rb') 
    for rec in bam_read.fetch():
        if rec.qname in aln_cnt:
            aln_cnt[rec.qname]+=1
        else:
            aln_cnt[rec.qname]=1
    bam_read.close()
    return aln_cnt

def get_annotation(fname):
    """Parse genome annotation from GFF3 file.
    """
    from Core import GFF
    te_featdb, teg_featdb, ribo_featdb, oth_featdb=dict(), dict(), dict(), dict()
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrM', 'ChrC']: ## change according to the GFF file TODO According to the file type add chromosome number automatically.
        print cid
        fh = open(fname, 'rU')
        limit_info = dict(gff_id=[cid], gff_source=['TAIR10']) ## change the source TODO According to the file type add source flag automatically 
        for rec in GFF.parse(fh, limit_info=limit_info):
            for each_rec in rec.features:
                if each_rec.type=='gene':
                    for child in each_rec.sub_features:
                        if child.type=='rRNA':
                            if cid in ribo_featdb:
                                ribo_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=each_rec.strand
                            else:
                                ribo_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)): 
                                    each_rec.strand}
                        elif child.type in ['mRNA', 'miRNA', 'ncRNA', 'snoRNA', 'snRNA', 'tRNA']:
                            if cid in oth_featdb:
                                oth_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=each_rec.strand
                            else:
                                oth_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)):
                                    each_rec.strand}
                elif each_rec.type=='transposable_element_gene':
                    if cid in teg_featdb:
                        teg_featdb[cid][(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=each_rec.strand
                    else:
                        teg_featdb[cid]={(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                each_rec.strand}
                elif each_rec.type=='transposable_element':
                    if cid in te_featdb:
                        te_featdb[cid][(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=each_rec.strand
                    else:
                        te_featdb[cid]={(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                each_rec.strand}
                elif each_rec.type=='pseudogene':
                    if cid in oth_featdb:
                        oth_featdb[cid][(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=each_rec.strand
                    else:
                        oth_featdb[cid]={(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                each_rec.strand}
        fh.close()
    """ribo_featdb['Chr3'][(14199917, 14203578)]=1 ## Unannotated rRNA from A thaliana genome. 
    ribo_featdb['Chr2'][(5784, 9683)]=1 ## Unannotated rRNA from A thaliana genome. 
    ribo_featdb['Chr2'][(2821, 3704)]=1 #
    ribo_featdb['Chr3'][(14196614, 14197675)]=1 #
    ribo_featdb['Chr3'][(14194052, 14194611)]=1 #
    ribo_featdb['Chr3'][(14199498, 14199751)]=1 #
    ribo_featdb['Chr3'][(14195564, 14195739)]=1 #
    ribo_featdb['ChrM'][(11426, 11883)]=-1 #
    ribo_featdb['ChrM'][(364594, 365124)]=1 #"""
    return te_featdb, teg_featdb, ribo_featdb, oth_featdb

def get_region_alignment(bam_fh, ctg_id, start, stop):
    """Get the read alignment from each TE region based on the request.
    """
    samdb=dict()
    smap={'+' : 1, '-' : -1}
    for each_align in bam_fh.fetch(ctg_id, start, stop):
        for cont in each_align.tags:
            if cont[0]=='NM':
                NM=cont[1]
            elif cont[0]=='XS':
                orient=cont[1]
                break
        if each_align.qname in samdb:
            samdb[each_align.qname].append((bam_fh.getrname(each_align.rname), each_align.pos, smap[orient], NM, each_align.qlen))
        else:
            samdb[each_align.qname]=[(bam_fh.getrname(each_align.rname), each_align.pos, smap[orient], NM, each_align.qlen)]
    bamdb = [(fid, finfo) for fid, finfo in samdb.items()]
    return bamdb

def get_ribo_reads(bam_fh, rib_fdb):
    """Take the reads from ribosomal RNA region
    """
    ribo_db=dict()
    smap={'+' : 1, '-' : -1}
    for contig_nb, rib_info in rib_fdb.items():
        for cod, fstd in rib_info.items(): 
            for each_align in bam_fh.fetch(contig_nb, cod[0], cod[1]):
                for cont in each_align.tags:
                    if cont[0]=="XS":
                        orient=cont[1]
                        break
                if each_align.qname in ribo_db:
                    ribo_db[each_align.qname].append((fstd, smap[orient]))
                else:
                    ribo_db[each_align.qname]=[(fstd, smap[orient])]
    return ribo_db

if __name__=="__main__":
    try:
        bamfname=sys.argv[1] 
        gffname=sys.argv[2]
        contig_id=sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
    read_db=get_alignment(bamfname) # get alignment count for each read
    te_fdb, teg_fdb, rib_fdb, oth_fdb=get_annotation(gffname) # fetch genome annotation 
    bam_reader=pysam.Samfile(bamfname, "rb") # BAM file handle for querying different feature interval
    ribo_read_db=get_ribo_reads(bam_reader, rib_fdb) # get rRNA region reads 
    # fetch read alignments for each TE 
    # TODO Write a sam file with unique reads set. 
    if contig_id in te_fdb: 
        for features in te_fdb[contig_id].items(): # just one TE is processing at a time.
            unq_s, unq_a, all_s, all_a, te_rib_s, te_rib_a=0, 0, 0, 0, 0, 0
            align_db=get_region_alignment(bam_reader, contig_id, features[0][1], features[0][2])# fetch the read alignment corresponding to this feature
            if len(align_db)>0: # feature has a coverage 
                for each_item in align_db: # aligned reads 
                    rid, rinfo=each_item
                    if read_db[rid]>1: # multiple alignments
                        std=[aln[2] for aln in rinfo]# alignment orientation 
                        std=list(set(std))
                        if len(std)>1: # multiple alignment fall in same region different orientation. 
                            print rid
                            all_s+=0.5
                            all_a+=0.5
                        else:
                            if features[1]==std[0]: # sense/anti-sense orientation
                                all_s+=1
                            else:
                                all_a+=1
                        if rid in ribo_read_db: # TODO different combination check whether the reads are sharing from rRNA 
                            if features[1]==std[0] and ribo_read_db[rid][0]==ribo_read_db[rid][1]:
                                te_rib_s+=1
                            if features[1]!=std[0] and ribo_read_db[rid][0]!=ribo_read_db[rid][1]:
                                te_rib_a+=1
                    else: # check whether this unique alignment is in a shared annotated feature region from GFF.
                        orientation=None
                        if rinfo[0][0] in rib_fdb:
                            for region, strand in rib_fdb[rinfo[0][0]].items():
                                if region[0]<=rinfo[0][1]+rinfo[0][4] and rinfo[0][1]+rinfo[0][4]<=region[1]: # + aligned read length 
                                    orientation=strand
                                    break
                        if orientation:
                            if features[1]==orientation: # same orientation two different feature annotation hard to decide 
                                continue 
                            else:
                                if features[1]==rinfo[0][2]: # only considering the sense direction, anti-sense we cannot make decision when there is something present in the opposite direction
                                    unq_s+=1
                            continue 
                        if rinfo[0][0] in oth_fdb:
                            for region, strand in oth_fdb[rinfo[0][0]].items():
                                if region[0]<=rinfo[0][1]+rinfo[0][4] and rinfo[0][1]+rinfo[0][4]<=region[1]: # + aligned read length 
                                    orientation=strand
                                    break
                        if orientation:
                            if features[1]==orientation: 
                                continue 
                            else:
                                if features[1]==rinfo[0][2]:
                                    unq_s+=1
                            continue 
                        if features[1]==rinfo[0][2]:
                            unq_s+=1
                        else:
                            unq_a+=1
                #break # one covered feature 
            #print features[0][0], unq_s, unq_a, all_s, all_a, te_rib_s, te_rib_a 
    bam_reader.close()
