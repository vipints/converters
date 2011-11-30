#!/usr/bin/python 
"""From BAM file find out the genomic location with a big cluster of reads 
Usage: cluster_region.py in.bam
"""

import re, sys, os 
from bx.intervals.cluster import ClusterTree
import collections
import pysam 
import time
from operator import itemgetter

def process_bam(bam_file):
    """Extract contents from BAM file
    """
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    for rec in sam_reader.fetch():
        yield rec.qname, sam_reader.getrname(rec.rname), rec.pos, rec.pos+rec.qlen
    sam_reader.close()

if __name__ == '__main__':

    try:
        bam_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1) 
    print time.asctime( time.localtime(time.time()) )
    cluster_distance = 20 # here considering reads falling in 20 nts distance share the same cluster   
    cluster_trees = collections.defaultdict(lambda:ClusterTree(cluster_distance, 2))
    rid_map, cnt = dict(), 0
    align_db=process_bam(bam_file)
    for rid, match_id, start, end in align_db:
        if not rid in rid_map: # get rid of long read id with small numbers 
            cnt +=1
            rid_map[rid] = cnt 
        cluster_trees[match_id].insert(start, end, rid_map[rid])
    for contig, ctree in sorted(cluster_trees.items()):
        freq_db = dict()
        for start, end, read_ids in ctree.getregions():
            freq_db[(start, end)] = len(read_ids)
        cand = 0 
        for region, rcnt in sorted(freq_db.items(), key=itemgetter(1), reverse=True):
            print contig + ':' + str(region[0]) + '..' + str(region[1]) + '\t' + str(rcnt)
            cand += 1
            if cand == 5:break ## Top 5 candidates from intergenic regions
        #break
    print time.asctime( time.localtime(time.time()) )
