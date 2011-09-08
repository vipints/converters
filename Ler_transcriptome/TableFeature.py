#!/usr/bin/python
"""
Create a Table of reads which aligned to the different genomic features.

ex: Table representation

Read ID [+/-] [No. of alignment]    [TE]    [rRNA]  [CG]    [OtherRNA]
4546697989  [+] [2] [0] [1]

Usage:
samtools view BAM_MM_0-1-2_CG/TE/rRNA/OtherRNA | TableFeature.py -f <MAT file name>

This will give a complete picture of reads aligned to different genome. 
Also make a association map of reads between different genome.
"""

class Table(object):
    def __init__(self):
        self.rows = []
        self.nameMaps = dict()

    def addRow(self, row):
        self.rows.append(row)
        self.nameMap[row['name']] = row

    def getRow(self, name):
        return self.nameMap[name]


def infoFeature(feature):

    featdb = dict()
    fh = open('/fml/ag-raetsch/home/vipin/tmp/rRNA_Ath_seq/' + feature, 'rU') ## care about the path in this line 
    for line in fh:
        line = line.strip('\n\r').split('\t')
        featdb[line[0]] = line[1]
    fh.close()
    return featdb

def AlignGenerator(fh):

    infodb = dict()
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
            infodb[line[0]].append((orient, int(NM), line[2]))
        else:
            infodb[line[0]] = [(orient, int(NM), line[2])]
    fh.close()
    return infodb
    
def init_table():
    table = dict(id = '',
        count = None,
        orient = None,
        CG = [],
        CG_MM = [],
        TE = [],
        TE_MM = [], 
        Oth = [], 
        Oth_MM = [], 
        ribo = [],
        ribo_MM = []
        )
    return table 

def readScaner(alignDB, cg_featdb, oth_featdb, ribo_featdb, te_featdb):
    
    read_dist = [] 
    for rid, rinfo in alignDB.items():
        data_table = init_table()
        data_table['id'] = rid
        data_table['count'] = len(rinfo)
        read_strand = [] 
        for rdet in rinfo:
            read_strand.append(rdet[0])
            if rdet[2] in cg_featdb:
                if rdet[0] == cg_featdb[rdet[2]]:
                    data_table['CG'].append(0)
                    data_table['CG_MM'].append(rdet[1])
                else:
                    data_table['CG'].append(1)
                    data_table['CG_MM'].append(rdet[1])
            elif rdet[2] in oth_featdb:
                if rdet[0] == oth_featdb[rdet[2]]:
                    data_table['Oth'].append(0)
                    data_table['Oth_MM'].append(rdet[1])
                else:
                    data_table['Oth'].append(1)
                    data_table['Oth_MM'].append(rdet[1])
            elif rdet[2] in ribo_featdb:
                if rdet[0] == ribo_featdb[rdet[2]]:
                    data_table['ribo'].append(0)
                    data_table['ribo_MM'].append(rdet[1])
                else:
                    data_table['ribo'].append(1)
                    data_table['ribo_MM'].append(rdet[1])
            elif rdet[2] in te_featdb:
                if rdet[0] == te_featdb[rdet[2]]:
                    data_table['TE'].append(0)
                    data_table['TE_MM'].append(rdet[1])
                else:
                    data_table['TE'].append(1)
                    data_table['TE_MM'].append(rdet[1])
        read_strand = list(set(read_strand))

        strand_types = np.zeros((len(read_strand),), dtype=np.object)
        for i in range(len(read_strand)):
            strand_types[i]= np.array(read_strand[i])
        data_table['orient'] = strand_types

        read_dist.append(data_table)
    return read_dist

import sys, re
import numpy as np 
import scipy.io as sio

if __name__ == "__main__":

    try:
        if sys.argv[1] == '-f':
            bfh = sys.stdin
        mat_file = sys.argv[2]
    except:
        print 'Incorrect number of alignments supplied'
        print __doc__
        sys.exit(-1)

    fnames = ['CG_orientations.txt', 'Other_RNA_orientations.txt', 'rRNA_orientation.txt', 'TE_orientations.txt']

    # get the feature orientations
    cg_featdb = infoFeature(fnames[0])
    oth_featdb = infoFeature(fnames[1])
    ribo_featdb = infoFeature(fnames[2])
    te_featdb = infoFeature(fnames[3])

    # get alignment data 
    alignDB = AlignGenerator(bfh)

    # figure out each reads locations on different genome
    read_distribution = readScaner(alignDB, cg_featdb, oth_featdb, ribo_featdb, te_featdb)

    # write into a matlab struture.
    sio.savemat(mat_file, {'read_dist':read_distribution}, oned_as='row')

