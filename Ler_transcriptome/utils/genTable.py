#!/usr/bin/env python 
"""
Generate a table with annotated TE's from Col-0 genome and calculate the 
read distribution in sense/anti-sense/uniq locations. 

Usage: genTable.py in.bam ribo.bam in.gff3
"""
import os 
import sys, re, pysam 
import multiprocessing
import collections
import itertools
import time 

class MapReduce(object):
    
    def __init__(self, map_func, reduce_func, num_workers=None):
        self.map_func = map_func
        self.reduce_func = reduce_func
        self.pool = multiprocessing.Pool(num_workers)
    
    def partition(self, mapped_values):
        partitioned_data = collections.defaultdict(list)
        for key, value in mapped_values:
            partitioned_data[key].append(value)
        return partitioned_data.items()
    
    def __call__(self, inputs, chunksize=1):
        map_responses = self.pool.map(self.map_func, inputs, chunksize=chunksize)
        partitioned_data = self.partition(itertools.chain(*map_responses))
        reduced_values = self.pool.map(self.reduce_func, partitioned_data)
        return reduced_values

def get_alignment(bfname):
    """Read BAM file and count the number of alignments for each read.
    """
    if not os.path.exists(bfname+ ".bai"):
        pysam.index(bfname)
    aln_cnt=dict()
    bam_read=pysam.Samfile(bfname, 'rb') 
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
    smap={1:'+',  -1:'-'}
    from Core import GFF
    te_featdb, teg_featdb, ribo_featdb, oth_featdb=dict(), dict(), dict(), dict()
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrM', 'ChrC']: 
        sys.stderr.write(cid+" .....\n")
        fh=open(fname, 'rU')
        limit_info=dict(gff_id=[cid], gff_source=['TAIR10']) 
        for rec in GFF.parse(fh, limit_info=limit_info):
            for each_rec in rec.features:
                if each_rec.type=='gene':
                    for child in each_rec.sub_features:
                        if child.type in ['mRNA', 'miRNA', 'ncRNA', 'snoRNA', 'snRNA', 'tRNA']:
                            if cid in oth_featdb:
                                oth_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=smap[each_rec.strand]
                            else:
                                oth_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)):
                                    smap[each_rec.strand]}
                        elif child.type in ['rRNA']:
                            if cid in ribo_featdb:
                                ribo_featdb[cid][(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position))]=smap[each_rec.strand]
                            else:
                                ribo_featdb[cid]={(int(each_rec.location._start.position), 
                                    int(each_rec.location._end.position)):
                                    smap[each_rec.strand]}
                elif each_rec.type=='transposable_element_gene':
                    if cid in teg_featdb:
                        teg_featdb[cid][(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=smap[each_rec.strand]
                    else:
                        teg_featdb[cid]={(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                smap[each_rec.strand]}
                elif each_rec.type=='transposable_element':
                    if cid in te_featdb:
                        te_featdb[cid][(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=smap[each_rec.strand]
                    else:
                        te_featdb[cid]={(each_rec.id, int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                smap[each_rec.strand]}
                elif each_rec.type=='pseudogene':
                    if cid in oth_featdb:
                        oth_featdb[cid][(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position))]=smap[each_rec.strand]
                    else:
                        oth_featdb[cid]={(int(each_rec.location._start.position), 
                                int(each_rec.location._end.position)): 
                                smap[each_rec.strand]}
        fh.close()
    return te_featdb, teg_featdb, ribo_featdb, oth_featdb

def get_region_alignment(bam_fh, ctg_id, start, stop):
    """Get the read alignment from each TE region based on the request.
    """
    samdb=dict()
    for each_align in bam_fh.fetch(ctg_id, start, stop):
        for cont in each_align.tags:
            if cont[0]=='NM':
                NM=cont[1]
            elif cont[0]=='XS':
                orient=cont[1]
                break
        if each_align.qname in samdb:
            samdb[each_align.qname].append((bam_fh.getrname(each_align.rname), each_align.pos, orient, NM, each_align.qlen))
        else:
            samdb[each_align.qname]=[(bam_fh.getrname(each_align.rname), each_align.pos, orient, NM, each_align.qlen)]
    bamdb = [(fid, finfo) for fid, finfo in samdb.items()]
    return bamdb

def get_ribo_reads(rib_fname):
    """Take the reads from ribosomal RNA region
    """
    ribo_db=dict()
    bam_read = pysam.Samfile(rib_fname, 'rb') 
    for each_align in bam_read.fetch():
        for cont in each_align.tags:
            if cont[0]=="XS":
                orient=cont[1]
                break
        ribo_db[each_align.qname]=('+', orient)
    bam_read.close()
    return ribo_db

def ReduceFN(sub_items):
    """reduce the result and make a 
    """
    fid, values=sub_items
    return [(fid, values)]  

def MapFN(items):
    """calculate the type of read distribution on each TE 
    Unique sense |Unique antisense |All sense |All antisense |TE rRNA sense |TE rRNA antisense |Sense TE antisense rRNA |Antisense TE sense rRNA
    """
    features, strand=items
    unq_s, unq_a, te_rib_s, te_rib_a=0, 0, 0, 0
    all_s, all_a, ste_a_rib, ate_s_rib=0, 0, 0, 0, 
    contig_id='Chr'+re.search(r'AT(.+)TE*', features[0]).group(1)
    align_db=get_region_alignment(bam_reader, contig_id, features[1], features[2])
    sys.stderr.write(features[0]+"\n")
    if len(align_db)>0: # feature has a coverage 
        for each_item in align_db: # aligned reads 
            rid, rinfo=each_item
            if read_db[rid]>1: # multiple alignments
                std=[aln[2] for aln in rinfo] # alignment orientation 
                std=list(set(std))
                if len(std)>1:continue # different orientation alignment for a read in single location 
                if strand==std[0]: # sense/anti-sense orientation
                    all_s+=1
                else:
                    all_a+=1
                if rid in ribo_read_db: # check whether the reads are sharing with rRNA 
                    if strand==std[0] and ribo_read_db[rid][0]==ribo_read_db[rid][1]: 
                        te_rib_s+=1 # sense reads to TE & rRNA
                    elif strand!=std[0] and ribo_read_db[rid][0]!=ribo_read_db[rid][1]:
                        te_rib_a+=1 # Antisense read to TE & rRNA
                    elif strand==std[0] and ribo_read_db[rid][0]!=ribo_read_db[rid][1]:
                        ste_a_rib+=1 # sense TE antisense rRNA
                    elif strand!=std[0] and ribo_read_db[rid][0]==ribo_read_db[rid][1]:
                        ate_s_rib+=1 # antisense TE sense rRNA
            else: # check whether this unique alignment is in a shared annotated feature region from GFF.
                if features[0] in te_marker:
                    if te_marker[features[0]]:
                        if strand==te_marker[features[0]]:
                            continue
                        else:
                            if strand==rinfo[0][2]:
                                unq_s+=1
                                all_s+=1
                    else:
                        if strand==rinfo[0][2]:
                            unq_s+=1
                            all_s+=1
                        else:
                            unq_a+=1
                            all_a+=1
    return [(features[0], (unq_s, unq_a, all_s, all_a, te_rib_s, te_rib_a, ste_a_rib, ate_s_rib))] 

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
                            if features[1]==std[0] and ribo_read_db[rid][0]==ribo_read_db[rid][1]: 
                                te_rib_s+=1 # sense reads to TE & rRNA
                            elif features[1]!=std[0] and ribo_read_db[rid][0]!=ribo_read_db[rid][1]:
                                te_rib_a+=1 # Antisense read to TE & rRNA
                            elif features[1]==std[0] and ribo_read_db[rid][0]!=ribo_read_db[rid][1]:
                                ste_a_rib+=1 # sense TE antisense rRNA
                            elif features[1]!=std[0] and ribo_read_db[rid][0]==ribo_read_db[rid][1]:
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

    outsam_sense=pysam.Samfile(contig_id+'_S.sam', 'w', template=bam_reader)
    outsam_asense=pysam.Samfile(contig_id+'_AS.sam', 'w', template=bam_reader)
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

def std_ribo_reads(bamfname, ribo_read_db, rib_fdb):
    """get the reads mapping to the ribosomal region defined in the feature annotations.
    """
    bam_fh = pysam.Samfile(bamfname, 'rb') 
    for cid, loq in rib_fdb.items():
        for locations, strand in loq.items(): 
            for each_align in bam_fh.fetch(cid, locations[0], locations[1]):
                for cont in each_align.tags:
                    if cont[0]=="XS":
                        orient=cont[1]
                        break
                ribo_read_db[each_align.qname]=(strand, orient)
    bam_fh.close()     
    return ribo_read_db

def __main__():
    try:
        bamfname=sys.argv[1] 
        rib_bamf=sys.argv[2]
        gffname=sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
    ## read frequency table 
    global read_db
    read_db=get_alignment(bamfname)
    sys.stderr.write(str(len(read_db))+' mapped reads stored in db\n')
    ## fetch genome annotation 
    te_fdb, teg_fdb, rib_fdb, oth_fdb=get_annotation(gffname) 
    sys.stderr.write('genome annotation stored in db\n')
    ## reads from ribosomal RNA region (extra rRNA genome)
    global ribo_read_db
    ribo_read_db=get_ribo_reads(rib_bamf)
    sys.stderr.write(str(len(ribo_read_db))+' reads from rRNA\n')
    ## reads to annotated ribosomal RNA from feature description  
    ribo_read_db=std_ribo_reads(bamfname, ribo_read_db, rib_fdb)
    sys.stderr.write(str(len(ribo_read_db))+' reads from rRNA\n')
    ## execution in 5 different process 
    global bam_reader
    bam_reader=pysam.Samfile(bamfname, "rb") # BAM file handle for querying different feature interval
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']: # check for TE's in Chr1..5
        global te_marker
        te_marker=overlap_feature(te_fdb, teg_fdb, rib_fdb, oth_fdb, cid)
        sys.stderr.write(cid+' feature indexing done\n')
        JobMap=MapReduce(MapFN, ReduceFN, 3)
        ResData=JobMap(te_fdb[cid].items())
        for element in ResData:
            xq=[str(ele) for ele in element[0][1][0]]
            xq.insert(0, element[0][0])
            print ' '.join(xq)
    #    sense_reads, asense_reads=process_elements(cid, read_db, te_fdb, ribo_read_db, bam_reader, rib_fdb, oth_fdb, teg_fdb, te_marker)
    #    writeSAM(bam_reader, sense_reads, asense_reads, cid)
    bam_reader.close()

if __name__=="__main__":
    __main__()
