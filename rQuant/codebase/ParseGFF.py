#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010 Vipin T Sreedharan
# Copyright (C) 2010 Max Planck Society
#

"""
Description: Parse genome annotation from a GFF3 (a tab delimited format for storing sequence features and annotations:
http://www.sequenceontology.org/gff3.shtml) 
file and create gene struct which can be used for rQuant downstream processing.
"""

import re, sys, os
import scipy.io

def CreateExon(strand_p, five_p_utr, cds_cod, three_p_utr):
    """Create exon cordinates from UTR's and CDS region"""
     
    exon_pos = []
    if strand_p == '+':        
        utr5_start, utr5_end = 0, 0
        if five_p_utr != []:
            utr5_start = five_p_utr[-1][0] 
            utr5_end = five_p_utr[-1][1]
        cds_5start = cds_cod[0][0]
        cds_5end = cds_cod[0][1]  
        jun_exon = []
        if cds_5start-utr5_end == 0 or cds_5start-utr5_end == 1:
            jun_exon = [utr5_start, cds_5end]    
        if len(cds_cod) == 1:
            five_prime_flag = 0
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                five_prime_flag = 1
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []: 
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            if utr3_start-cds_5end == 0 or utr3_start-cds_5end == 1:
                jun_exon = [cds_5start, utr3_end]
            three_prime_flag = 0
            if jun_exon != []: 
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
                three_prime_flag = 1
            if five_prime_flag == 1 and three_prime_flag == 1:
                exon_pos.append([utr5_start, utr3_end])
            if five_prime_flag == 1 and three_prime_flag == 0:
                exon_pos.append([utr5_start, cds_5end])
                cds_cod = cds_cod[:-1]
            if five_prime_flag == 0 and three_prime_flag == 1:
                exon_pos.append([cds_5start, utr3_end])
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
        else:    
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []:
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            cds_3start = cds_cod[-1][0]
            cds_3end = cds_cod[-1][1]
            if utr3_start-cds_3end == 0 or utr3_start-cds_3end == 1:       
                jun_exon = [cds_3start, utr3_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
    elif strand_p == '-':
        utr3_start, utr3_end = 0, 0        
        if three_p_utr != []:
            utr3_start = three_p_utr[-1][0]
            utr3_end = three_p_utr[-1][1]
        cds_3start = cds_cod[0][0]
        cds_3end = cds_cod[0][1]
        jun_exon = []
        if cds_3start-utr3_end == 0 or cds_3start-utr3_end == 1:
            jun_exon = [utr3_start, cds_3end]  
        if len(cds_cod) == 1:    
            three_prime_flag = 0
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                three_prime_flag = 1
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]
            if utr5_start-cds_3end == 0 or utr5_start-cds_3end == 1:
                jun_exon = [cds_3start, utr5_end]
            five_prime_flag = 0
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
                five_prime_flag = 1
            if three_prime_flag == 1 and five_prime_flag == 1:
                exon_pos.append([utr3_start, utr5_end])
            if three_prime_flag == 1 and five_prime_flag == 0:
                exon_pos.append([utr3_start, cds_3end])
                cds_cod = cds_cod[:-1]
            if three_prime_flag == 0 and five_prime_flag == 1:
                exon_pos.append([cds_3start, utr5_end])        
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
        else:
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr3 in three_p_utr:
                exon_pos.append(utr3)   
            if jun_exon != []:
                exon_pos.append(jun_exon)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]    
            cds_5start = cds_cod[-1][0]
            cds_5end = cds_cod[-1][1]
            if utr5_start-cds_5end == 0 or utr5_start-cds_5end == 1:
                jun_exon = [cds_5start, utr5_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            if jun_exon != []:
                exon_pos.append(jun_exon)    
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
    return exon_pos

def init_gene():
    """Initializing the gene structure"""
    
    gene_details = dict(chr = '', exons = [], gene_info = {}, id = '', is_alt_spliced = 0, name = '', source = '', start = '', stop = '', strand = '', transcripts = [])
    return gene_details

def FeatureValueFormat(singlegene):
    """Make feature value compactable to write in a .mat format"""

    ## based on the feature set including for rQuant process each genes selected feature values. 
    import numpy as np
    comp_exon = np.zeros((len(singlegene['exons']),), dtype=np.object)
    for i in range(len(singlegene['exons'])):
        comp_exon[i]= np.array(singlegene['exons'][i])
    singlegene['exons'] = comp_exon
    comp_transcripts = np.zeros((len(singlegene['transcripts']),), dtype=np.object)
    for i in range(len(singlegene['transcripts'])):
        comp_transcripts[i] = np.array(singlegene['transcripts'][i])
    singlegene['transcripts'] = comp_transcripts
    return singlegene 

def CreateGeneModels(genes_cmpt, transcripts_cmpt, exons_cmpt, utr3_cmpt, utr5_cmpt, cds_cmpt):
    """Creating Coding/Non-coding gene models from parsed GFF objects."""

    gene_counter, gene_models = 1, []
    for gene_entry in genes_cmpt: ## Figure out the genes and transcripts associated feature 
        if gene_entry in transcripts_cmpt:
            gene = init_gene() ## gene section related tags
            gene['id'] = gene_counter
            gene['name'] = gene_entry[1]
            gene['chr'] = genes_cmpt[gene_entry]['chr']
            gene['source'] = genes_cmpt[gene_entry]['source']
            gene['start'] = genes_cmpt[gene_entry]['start']
            gene['stop'] = genes_cmpt[gene_entry]['stop']
            gene['strand'] = genes_cmpt[gene_entry]['strand']
            if gene['strand'] != '+' and gene['strand'] != '-': gene['strand'] = '.' # Strand info not known replaced with a dot symbol instead of None, ?, . etc.
            general_info = dict()
            ## TODO add more gene related information from attribute column of GFF file based on the reserved key words
            if 'Name' in genes_cmpt[gene_entry]:general_info['Name'] = genes_cmpt[gene_entry]['Name']
            if 'Note' in genes_cmpt[gene_entry]:general_info['Note'] = genes_cmpt[gene_entry]['Note']
            if 'Alias' in genes_cmpt[gene_entry]:general_info['Alias'] = genes_cmpt[gene_entry]['Alias']
            if general_info == {}:general_info['ID'] = gene_entry[1]
            gene['gene_info'] = general_info
            if len(transcripts_cmpt[gene_entry]) > 1:gene['is_alt_spliced'] = 1
            for tids in transcripts_cmpt[gene_entry]: ## transcript section related tags 
                gene['transcripts'].append(tids['ID'])
                exon_cod = []
                if len(exons_cmpt) != 0: ## rQuant requires only exon coordinates of the transcripts 
                    if (gene['chr'], tids['ID']) in exons_cmpt:
                        for feat_exon in exons_cmpt[(gene['chr'], tids['ID'])]:exon_cod.append([feat_exon['start'], feat_exon['stop']])
                else: ## build exon coordinates from UTR3, UTR5 and CDS
                    utr5_pos, cds_pos, utr3_pos = [], [], []
                    if (gene['chr'], tids['ID']) in utr5_cmpt:
                        for feat_utr5 in utr5_cmpt[(gene['chr'], tids['ID'])]:utr5_pos.append([feat_utr5['start'], feat_utr5['stop']])
                    if (gene['chr'], tids['ID']) in cds_cmpt:
                        for feat_cds in cds_cmpt[(gene['chr'], tids['ID'])]:cds_pos.append([feat_cds['start'], feat_cds['stop']])
                    if (gene['chr'], tids['ID']) in utr3_cmpt:
                        for feat_utr3 in utr3_cmpt[(gene['chr'], tids['ID'])]:utr3_pos.append([feat_utr3['start'], feat_utr3['stop']])
                    exon_cod = CreateExon(gene['strand'], utr5_pos, cds_pos, utr3_pos) 
                ## generalize the coordinate system for exons, GFF file may contain ascending or descending order.
                if gene['strand'] == '-':
                    if exon_cod != [] and len(exon_cod) != 1:
                        if exon_cod[0][0] > exon_cod[-1][0]: exon_cod.reverse()
                if exon_cod: gene['exons'].append(exon_cod)
            ## make a compact form of features in each gene struct to write into .mat format.
            gene = FeatureValueFormat(gene)
            gene_counter += 1
            gene_models.append(gene)
    return gene_models    

def GFFParse(gff_file):
    """Parsing GFF file based on feature"""

    genes, transcripts, exons, utr3, utr5, cds = {}, {}, {}, {}, {}, {}
    gff_handle = open(gff_file, "rU")
    for gff_line in gff_handle:
        gff_line = gff_line.strip('\n\r').split('\t')
        if not gff_line:continue
        if re.match(r'#', gff_line[0]) or re.match(r'>', gff_line[0]):continue
        if len(gff_line) == 1:continue ## GFF files with genome sequence in FASTA at the end 
        assert (len(gff_line)==9), '\t'.join(gff_line)
        if gff_line[3] == '' or gff_line[4] == '' or gff_line[-1] == '':sys.stdout.write('Warning: invalid GFF line\t' + '\t'.join(gff_line) + '\n');continue
        if gff_line[2] == 'gene' or gff_line[2] == 'pseudogene':
            gid, gene_info = None, dict()
            gene_info['start'] = int(gff_line[3])
            gene_info['stop'] = int(gff_line[4])
            gene_info['chr'] = gff_line[0]
            gene_info['source'] = gff_line[1]
            gene_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'ID':gid=attr[1];continue 
                gene_info[attr[0]] = attr[1]
            genes[(gff_line[0], gid)] = gene_info
        elif gff_line[2] == 'mRNA' or gff_line[2] == 'transcript' or gff_line[2] == 'ncRNA' or gff_line[2] == 'miRNA' or gff_line[2] == 'pseudogenic_transcript' or gff_line[2] == 'rRNA' or gff_line[2] == 'snoRNA' or gff_line[2] == 'snRNA' or gff_line[2] == 'tRNA' or gff_line[2] == 'scRNA': # TODO Include non coding transcripts 
            gid, mrna_info = None, dict() 
            mrna_info['start'] = int(gff_line[3])
            mrna_info['stop'] = int(gff_line[4])
            mrna_info['chr'] =  gff_line[0]
            mrna_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':gid=attr[1];continue
                mrna_info[attr[0]] = attr[1]
            if (gff_line[0], gid) in transcripts:
                transcripts[(gff_line[0], gid)].append(mrna_info)
            else:
                transcripts[(gff_line[0], gid)] = [mrna_info]
        elif gff_line[2] == 'exon':
            tids, exon_info = None, dict()
            exon_info['start'] = int(gff_line[3])
            exon_info['stop'] = int(gff_line[4])
            exon_info['chr'] =  gff_line[0]
            exon_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tids=attr[1];continue
                exon_info[attr[0]] = attr[1]
            for tid in tids.split(','):
                if (gff_line[0], tid) in exons:
                    exons[(gff_line[0], tid)].append(exon_info)
                else:
                    exons[(gff_line[0], tid)] = [exon_info]
        elif gff_line[2] == 'five_prime_UTR':
            utr5_info, tids = dict(), None
            utr5_info['start'] = int(gff_line[3])
            utr5_info['stop'] = int(gff_line[4])
            utr5_info['chr'] =  gff_line[0]
            utr5_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tids=attr[1];continue
                utr5_info[attr[0]] = attr[1]
            for tid in tids.split(','):
                if (gff_line[0], tid) in utr5:
                    utr5[(gff_line[0], tid)].append(utr5_info)
                else:
                    utr5[(gff_line[0], tid)] = [utr5_info]
        elif gff_line[2] == 'CDS':
            cds_info, tids = dict(), None
            cds_info['start'] = int(gff_line[3])
            cds_info['stop'] = int(gff_line[4])
            cds_info['chr'] =  gff_line[0]
            cds_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tids=attr[1];continue
                cds_info[attr[0]] = attr[1]
            for tid in tids.split(','):
                if (gff_line[0], tid) in cds:
                    cds[(gff_line[0], tid)].append(cds_info)
                else:
                    cds[(gff_line[0], tid)] = [cds_info]
        elif gff_line[2] == 'three_prime_UTR':
            utr3_info, tids = dict(), None
            utr3_info['start'] = int(gff_line[3])
            utr3_info['stop'] = int(gff_line[4])
            utr3_info['chr'] =  gff_line[0]
            utr3_info['strand'] = gff_line[6]
            for attr in gff_line[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tids=attr[1];continue
                utr3_info[attr[0]] = attr[1]
            for tid in tids.split(','):
                if (gff_line[0], tid) in utr3:
                    utr3[(gff_line[0], tid)].append(utr3_info)
                else:
                    utr3[(gff_line[0], tid)] = [utr3_info]
    gff_handle.close()
    return genes, transcripts, exons, utr3, utr5, cds

def __main__():
    """This function provides a best way to extract genome feature information from a GFF3 file for the rQuant downstream processing"""

    try:
        gff_file = sys.argv[1]
        gff_source = sys.argv[2]
        gff_contigs = sys.argv[3]
        mat_file = sys.argv[4]
    except:
        sys.stderr.write('Check Input Requirements:\n\t1. Genome annotation in GFF3\n\t2. Gene model source (all- all source gene models/specify source -help:2nd column in GFF file) ex: Coding_transcript\n\t3. Contig(s) (all- all contigs/comma seperated list of contigs) ex: II,IV\n\t4. Result file name ex: gene_models.mat\n')
        sys.exit(-1)
    
    ## 1. Check for user supplied contig present in GFF file 
    contigs = []
    first_col = os.popen('cut -f 1 ' + gff_file + ' | sort | uniq')
    for entries in first_col:
        entries = entries.strip( '\r\n' )
        if not entries or entries[0] == '#':continue
        contigs.append(entries)
    first_col.close() 
    contig_flag = 0
    for gff_chr in gff_contigs.split(','):
        if not gff_chr in contigs:
            if not gff_chr == 'all':sys.stdout.write('Warning: Contig -' + gff_chr + '- NOT present in GFF file\n');contig_flag +=1
    if len(gff_contigs.split(',')) == contig_flag:sys.stderr.write('Error: Absence of provided CONTIG(s) in GFF file, Cannot continue, Program terminating.\n');sys.exit(-1)        
    ## 2. Check for the Source in GFF file 
    source, feature = [], []
    twoThree_col = os.popen('cut -f 2 ' + gff_file + ' | sort | uniq')
    for entries in twoThree_col:
        entries = entries.strip( '\r\n' ).split( '\t' )
        if not entries or re.match(r'#', entries[0]):continue
        source.append(entries[0])
    twoThree_col.close()
    source_flag = 0
    for src in gff_source.split(','):
        if not src in source:
            if not src == 'all':sys.stdout.write('Warning: Source -' + src + '- NOT present in GFF file.\n');source_flag +=1
    if len(gff_source.split(','))==source_flag:sys.stderr.write('Error: Absence of provided SOURCE(s) in GFF file, Cannot continue, Program terminating.\n');sys.exit(-1)        
    ## 3. Creating gene models based on the following types. 
    type_col = os.popen('cut -f 3 ' + gff_file + ' | sort | uniq')
    for types in type_col:
        types = types.strip('\r\n')
        if not types or re.match(r'#', types):continue
        feature.append(types)
    type_col.close()
    if len(feature)<3:sys.stderr.write('Error: Creating Gene models need at least 3 features like gene, mRNA, exon. Cannot continue, Program terminating.\n');sys.exit(-1)    
    ## Parse content according to Parent Child relationship   
    genes, transcripts, exons, utr3, utr5, cds = GFFParse(gff_file)
    ## Creating gene struct from parsed GFF content 
    gene_models = CreateGeneModels(genes, transcripts, exons, utr3, utr5, cds)
    # TODO Write to matlab struct instead of Cell arrays.
    # Saving genome annotation features to a .mat file 
    scipy.io.savemat(mat_file, mdict={'genes':gene_models}, format='5', oned_as='row')

if __name__=='__main__':__main__()
