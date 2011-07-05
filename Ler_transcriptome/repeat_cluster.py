#!/usr/bin/python 

## Description: Unmapped reads from Ler genome mapped against COL-0 genome, Now analyse whether the reads are aligned to 50bp or 78bp repeats of COL-0 genome

import re, sys
from bx.intervals.cluster import ClusterTree
import collections

def repeat_parse(score_file):
    
    infh = open(score_file, 'rU')
    for line in infh:
        line = line.strip().split(' ')
        if re.search(r'variableStep', line[0]):
            contig = re.search(r'chrom=(.+)', line[1]).group(1)
            if contig== 'Chr2':break
            continue
        yield contig, int(line[0]), int(line[0]), int(line[1])
    infh.close()
   
def sam_reader(sam_file):

    reg_cov = dict()
    infh = open(sam_file, 'rU')
    for line in infh:
        line = line.strip().split('\t')
        try:
            if line[2] in reg_cov:
                reg_cov[line[2]].append(int(line[3]))
            else:
                reg_cov[line[2]] = [int(line[3])]
        except:
            continue
    infh.close()
    return reg_cov

if __name__ == '__main__':
    
    try:
        repeat_file_1 = sys.argv[1]
        alignment_file = sys.argv[2]
    except:
        print 'Provide repeat region file in WIG format, Alignment file in SAM format'
        sys.exit(-1)

    cluster_distance = 1 
    repeat_regions_50 = collections.defaultdict(lambda:ClusterTree(cluster_distance, 2))

    repeat_generator = repeat_parse(repeat_file_1)
    for match_id, start, end, score in repeat_generator:
        repeat_regions_50[match_id].insert(start, end, score)
    #print 'Number of clusters: ' + str(len(repeat_regions_50))
    location_db = sam_reader(alignment_file)
    repeat_cnt = 0
    for chrom, cluster_tree in repeat_regions_50.items():
        for start, end, scores in cluster_tree.getregions():
            for rloc in location_db[chrom]:
                if (rloc >= start and rloc <= end) or (rloc+80 >= start and rloc+80 <= end):repeat_cnt += 1
    print repeat_cnt 
