#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Vipin T. Sreedharan
# Copyright (C) 2009-2011 Friedrich Miescher Laboratory of the Max Planck Society
#
# Description: This program is used to get genome annotation from a valid GFF3 formated file.
# 

import re, sys, time
import numpy as np
import scipy.io
from Core import GFF

class ordered_dict(dict):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self._order = self.keys()

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key in self._order:
            self._order.remove(key)
        self._order.append(key)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._order.remove(key)

    def order(self):
        return self._order[:]

    def ordered_items(self):
        return [(key,self[key]) for key in self._order]


def OrganizePacket(singlegene):
    """
    Writing into matlab needs a compressed way of representation of feature descriptions.
    """
    comp_exon = np.zeros((len(singlegene['exons']),), dtype=np.object)
    for i in range(len(singlegene['exons'])):
        comp_exon[i]= np.array(singlegene['exons'][i])
    singlegene['exons'] = comp_exon
    comp_transcripts = np.zeros((len(singlegene['transcripts']),), dtype=np.object)
    for i in range(len(singlegene['transcripts'])):
        comp_transcripts[i]= np.array(singlegene['transcripts'][i])
    singlegene['transcripts'] = comp_transcripts
    return singlegene

def init_gene():
    """Defining the gene information"""
    gene_info = dict(
          id = '',
          name = '',
          source = '',
          strand = '',
          chr = '',
          transcripts = [],
          exons = [],
          is_alt_spliced = 0
          )
    return gene_info

def ParserCommonType(fname, source):

    gff_type=['transcript', 'exon']
    source_type = zip([source] * len(gff_type), gff_type)
    filter_type = dict(gff_source_type = source_type, gff_id = ['I', 'II', 'III', 'IV', 'V', 'X'])
    gene_cnt, genes_data = 1, []
    gid = ordered_dict()
    gfh = open(fname)
    for rec in GFF.parse(gfh, limit_info=filter_type):
        for feature in rec.features:
            if feature.qualifiers['gene'][0] in gid:
                gene = gid[feature.qualifiers['gene'][0]]
                gene['is_alt_spliced'] = 1
                gene['transcripts'].append(feature.id)
                exon_pos = [] 
                for xlevel in feature.sub_features:
                    exon_pos.append([xlevel.location._start.position + 1, xlevel.location._end.position])
                if orient == '-':
                    if exon_pos != [] and len(exon_pos) != 1:
                        if exon_pos[0][0] > exon_pos[-1][0]: 
                            exon_pos.reverse()
                gene["exons"].append(exon_pos)
                gid[feature.qualifiers['gene'][0]]=gene
            else:
                gene = init_gene()
                gene['id'] = gene_cnt
                gene['name'] = feature.qualifiers['gene'][0]
                gene['chr'] = rec.id
                gene['source'] = feature.qualifiers['source'][0]
                orient = None
                if feature.strand == 1:
                    orient = '+'
                else:
                    orient = '-'
                gene['strand'] = orient
                gene['is_alt_spliced'] = 0
                gene['transcripts'].append(feature.id)
                exon_pos = [] 
                for xlevel in feature.sub_features:
                    exon_pos.append([xlevel.location._start.position + 1, xlevel.location._end.position])
                if orient == '-':
                    if exon_pos != [] and len(exon_pos) != 1:
                        if exon_pos[0][0] > exon_pos[-1][0]: 
                            exon_pos.reverse()
                gene["exons"].append(exon_pos)
                gene_cnt += 1
                gid[feature.qualifiers['gene'][0]]=gene
    gfh.close()
    for ent in gid.ordered_items():
        cand = OrganizePacket(ent[1])
        genes_data.append(cand)
    return genes_data

if __name__=="__main__":
    
    try:
        gff_fname = sys.argv[1]
        sname = sys.argv[2]
        mat_fname = sys.argv[3]
    except:
        print 'Access denied for GFF file'
        sys.exit(-1)
    #for source in ["miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "protein_coding", "pseudogene", "rRNA", "snRNA", "snoRNA"]:
    #for source in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']:
    genes_list = ParserCommonType(gff_fname, sname)
    scipy.io.savemat(mat_fname, mdict={'gene_models':genes_list}, format='5', oned_as='row')
