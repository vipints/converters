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

def process_elements(contig_id, read_db, te_fdb, ribo_read_db, bam_reader, rib_fdb, oth_fdb, teg_fdb, te_marker):
    """fetch read alignments for each TE 
    """
    if contig_id in te_fdb:
        cnt=1
        unq_s_rid, unq_a_rid=dict(), dict()
        for features in te_fdb[contig_id].items(): # just one TE is processing at a time.
            unq_s, unq_a, all_s, all_a, te_rib_s, te_rib_a=0, 0, 0, 0, 0, 0
            ste_a_rib, ate_s_rib=0, 0
            align_db=get_region_alignment(bam_reader, contig_id, features[0][1], features[0][2])# fetch the read alignment corresponding to this feature
            sys.stderr.write(str(cnt)+" "+features[0][0]+"\n")
            cnt+=1
            if len(align_db)>0: # feature has a coverage 
                for each_item in align_db: # aligned reads 
                    rid, rinfo=each_item
                    if read_db[rid]>1: # multiple alignments
                        std=[aln[2] for aln in rinfo] # alignment orientation 
                        std=list(set(std))
                        if len(std)>1:continue # different orientation alignment for a read in single location 
                        if features[1]==std[0]: # sense/anti-sense orientation
                            all_s+=1
                        else:
                            all_a+=1
                        if rid in ribo_read_db: # check whether the reads are sharing with rRNA 
                            if features[1]==std[0] and ribo_read_db[rid][0][0]==ribo_read_db[rid][0][1]: 
                                te_rib_s+=1 # sense reads to TE & rRNA
                            elif features[1]!=std[0] and ribo_read_db[rid][0][0]!=ribo_read_db[rid][0][1]:
                                te_rib_a+=1 # Antisense read to TE & rRNA
                            elif features[1]==std[0] and ribo_read_db[rid][0][0]!=ribo_read_db[rid][0][1]:
                                ste_a_rib+=1 # sense TE antisense rRNA
                            elif features[1]!=std[0] and ribo_read_db[rid][0][0]==ribo_read_db[rid][0][1]:
                                ate_s_rib+=1 # antisense TE sense rRNA
                    else: # check whether this unique alignment is in a shared annotated feature region from GFF.
                        if features[0][0] in te_marker:
                             if te_marker[features[0][0]]:
                                if features[1]==te_marker[features[0][0]]:
                                    continue
                                else:
                                    if features[1]==rinfo[0][2]:
                                        unq_s+=1
                                        all_s+=1
                                        unq_s_rid[rid]=1
                             else:
                                if features[1]==rinfo[0][2]:
                                    unq_s+=1
                                    all_s+=1
                                    unq_s_rid[rid]=1
                                else:
                                    unq_a+=1
                                    all_a+=1
                                    unq_a_rid[rid]=1
                #break # one covered feature 
            print features[0][0], unq_s, unq_a, all_s, all_a, te_rib_s, te_rib_a, ste_a_rib, ate_s_rib 
        return unq_s_rid, unq_a_rid  

def writeSAM(bam_reader, sense_reads, asense_reads, contig_id):

    outsam_sense=pysam.Samfile(contig_id+'_sense.sam', 'w', template=bam_reader)
    outsam_asense=pysam.Samfile(contig_id+'_antisense.sam', 'w', template=bam_reader)
    for rec in bam_reader.fetch():
        if rec.qname in sense_reads:
            outsam_sense.write(rec)
        elif rec.qname in asense_reads:
            outsam_asense.write(rec)
    outsam_sense.close()
    outsam_asense.close()

def overlap_feature(te_fdb, teg_fdb, rib_fdb, oth_fdb, chr_id):
    """Indexing the TE based on whether they are sharing in genome with annotated regions.
    """
    te_index=dict()
    for element in te_fdb[chr_id].items(): 
        orientation=None
        for region, strand in oth_fdb[chr_id].items():
            if (element[0][1]<=region[0] and element[0][2]>=region[1]) or (element[0][1]>=region[0] and element[0][1]<= region[1]) or (element[0][2]>=region[0] and element[0][2]<=region[1]):
                orientation=strand
                break
        if orientation:
            te_index[element[0][0]]=orientation
        else:
            if chr_id in rib_fdb:
                for region, strand in rib_fdb[chr_id].items():
                    if (element[0][1]<=region[0] and element[0][2]>=region[1]) or (element[0][1]>=region[0] and element[0][1]<= region[1]) or (element[0][2]>=region[0] and element[0][2]<=region[1]):
                        orientation=strand
                        break
            te_index[element[0][0]]=orientation
    return te_index

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
    te_marker=overlap_feature(te_fdb, teg_fdb, rib_fdb, oth_fdb, contig_id)
    bam_reader=pysam.Samfile(bamfname, "rb") # BAM file handle for querying different feature interval
    ribo_read_db=get_ribo_reads(bam_reader, rib_fdb) # get rRNA region reads 
    sense_reads, asense_reads=process_elements(contig_id, read_db, te_fdb, ribo_read_db, bam_reader, rib_fdb, oth_fdb, teg_fdb, te_marker)
    writeSAM(bam_reader, sense_reads, asense_reads, contig_id)
    bam_reader.close()
