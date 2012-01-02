#!/usr/bin/python
"""
Program to calculate the statistics of reads aligned to different annotated 
features in a genome. This is a time-cost effective method of implementation 
using the module called multiprocessing, along with MapReduce functionality. 

Usage: distributed_feature_stats.py in.bam in.gff3
"""
import re, sys 
import time
import pysam 
import multiprocessing
import collections
import itertools

class MapReduce(object):
    
    def __init__(self, map_func, reduce_func, num_workers=None):
        """
        map_func
          Function to map inputs to intermediate data. Takes as
          argument one input value and returns a tuple with the key
          and a value to be reduced.
        
        reduce_func
          Function to reduce partitioned version of intermediate data
          to final output. Takes as argument a key as produced by
          map_func and a sequence of the values associated with that
          key.
         
        num_workers
          The number of workers to create in the pool. Defaults to the
          number of CPUs available on the current host.
        """
        self.map_func = map_func
        self.reduce_func = reduce_func
        self.pool = multiprocessing.Pool(num_workers)
    
    def partition(self, mapped_values):
        """Organize the mapped values by their key.
        Returns an unsorted sequence of tuples with a key and a sequence of values.
        """
        partitioned_data = collections.defaultdict(list)
        for key, value in mapped_values:
            partitioned_data[key].append(value)
        return partitioned_data.items()
    
    def __call__(self, inputs, chunksize=1):
        """Process the inputs through the map and reduce functions given.
        
        inputs
          An iterable containing the input data to be processed.
        
        chunksize=1
          The portion of the input data to hand to each worker.  This
          can be used to tune performance during the mapping phase.
        """
        map_responses = self.pool.map(self.map_func, inputs, chunksize=chunksize)
        partitioned_data = self.partition(itertools.chain(*map_responses))
        reduced_values = self.pool.map(self.reduce_func, partitioned_data)
        return reduced_values

def map_fn(items):
    """Map function checking each reads location at the genome. 
    """
    rid, rinfo = items
    sys.stdout.write(rid);sys.stdout.flush()
    sys.stdout.write('\r' + ' '*len(rid))
    sys.stdout.flush()
    read_strand=[]
    ribX, cgX, othX, teX, pdX = 0, 0, 0, 0, 0
    intXrd, ribXrd, cgXrd, othXrd, teXrd, pdXrd = 0, 0, 0, 0, 0, 0
    for rdet in rinfo:
        read_strand.append(rdet[2])
        [ribXrdx, cgXrdx, othXrdx, teXrdx, pdXrdx] = [0, 0, 0, 0, 0]
        if rdet[0] in cg_featdb:
            for details, strand_info in cg_featdb[rdet[0]].items():
                if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+10:
                    cgX+=1
                    cgXrdx=1
                    if rdet[2]==strand_info:
                        cgXrd=1
                    else:
                        cgXrd=-1
                    break
        if rdet[0] in oth_featdb:
            for details, strand_info in oth_featdb[rdet[0]].items():
                if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                    othX+=1
                    othXrdx=1
                    if rdet[2]==strand_info:
                        othXrd=1
                    else:
                        othXrd=-1
                    break
        if rdet[0] in ribo_featdb:
            for details, strand_info in ribo_featdb[rdet[0]].items():
                if details[0]-95 <= rdet[1] and rdet[1] <= details[1]+10:
                    ribX+=1
                    ribXrdx =1
                    if rdet[2]==strand_info:
                        ribXrd=1
                    else:
                        ribXrd=-1
                    break
        if rdet[0] in te_featdb:
            for details, strand_info in te_featdb[rdet[0]].items():
                if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                    teX+=1
                    teXrdx=1
                    if rdet[2]==strand_info:
                        teXrd=1
                    else:
                        teXrd=-1
                    break
        if rdet[0] in psd_featdb:
            for details, strand_info in psd_featdb[rdet[0]].items():
                if details[0]-95 <= rdet[1] and rdet[1] <= details[1]:
                    pdX+=1
                    pdXrdx=1
                    if rdet[2]==strand_info:
                        pdXrd=1
                    else:
                        pdXrd=-1
                    break
        if [ribXrdx, cgXrdx, othXrdx, teXrdx, pdXrdx].count(0)==5:
            intXrd = 1
    xq = []
    if [ribXrd, cgXrd, othXrd, teXrd, pdXrd].count(0)==5: # intergenic read 
        xq = [('int_NR', 1), ('int_NA', len(rinfo))] 
    elif [ribXrd, cgXrd, othXrd, teXrd, pdXrd, intXrd].count(0)<=4: # multiple alignments spanning across different features
        xq = [('sh_NR', 1), ('sh_NA', len(rinfo))]
    elif cgXrd !=0:
        xq = [('cg_NR', 1), ('cg_NA', cgX)]
    elif othXrd !=0:
        xq = [('oth_NR', 1), ('oth_NA', othX)]
    elif teXrd !=0:
        xq = [('te_NR', 1), ('te_NA', teX)]
    elif pdXrd !=0:
        xq = [('pd_NR', 1), ('pd_NA', pdX)]
    elif ribXrd !=0:
        xq = [('rib_NR', 1), ('rib_NA', ribX)]
    return xq

def reduce_fn(sub_items):
    """Reduce function collecting each feature count
    """
    ftype, count_ftype = sub_items 
    return [(ftype, sum(count_ftype))]

def get_Feature(fname):
    """Extract genome annotation information from a provided GFF file.
    """
    from Core import GFF
    global te_featdb; global psd_featdb; global cg_featdb; global oth_featdb; global ribo_featdb
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
    ribo_featdb['Chr3'][(14199917, 14203578)]=1 ## Unannotated rRNA from A thaliana genome. 
    ribo_featdb['Chr2'][(5784, 9683)]=1 ## Unannotated rRNA from A thaliana genome. 
    ribo_featdb['Chr2'][(2821, 3704)]=1 #
    ribo_featdb['Chr3'][(14196614, 14197675)]=1 #
    ribo_featdb['Chr3'][(14194052, 14194611)]=1 #
    ribo_featdb['Chr3'][(14199498, 14199751)]=1 #
    ribo_featdb['Chr3'][(14195564, 14195739)]=1 #
    ribo_featdb['ChrM'][(11426, 11883)]=-1 #
    ribo_featdb['ChrM'][(364594, 365124)]=1 #"""
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
    bamdb = [(fid, finfo) for fid, finfo in samdb.items()]
    return bamdb
    
if __name__ == "__main__":
    try:
        bamf = sys.argv[1]
        anno_file = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    print time.asctime( time.localtime(time.time()) )
    te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb = dict(), dict(), dict(), dict(), dict() # declaring variable 
    te_featdb, psd_featdb, cg_featdb, oth_featdb, ribo_featdb = get_Feature(anno_file) # parse function for getting annotation 
    read_db = AlignGenerator(bamf) # parse function for getting read alignment information 
    master_job = MapReduce(map_fn, reduce_fn, 31) # create an object such that the jobs are distributing in 31 CPU's 
    results = master_job(read_db) # start the core the job
    print 
    for element in sorted(results):
        print element
    print time.asctime( time.localtime(time.time()) )
