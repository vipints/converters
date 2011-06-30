#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010 Vipin T Sreedharan, Friedrich Miescher Laboratory of the Max Planck Society
# Copyright (C) 2010 Max Planck Society
#
# Description : Convert a genome annotation in GFF3 format to UCSC 12 column Wiggle BED format. BED format typically represents the transcript models. 

import re, sys

def WriteBED(tinfo, einfo):
    
    for contig_id, features in tinfo.items():
        for tid, tloc in features.items():
            if tid in einfo: # get corresponding exon info 
                exon_cnt, exon_len, exon_cod, fex, rstart = 0, '', '', 0, None
                if tloc[-1] == '-':
                    if einfo[tid][0][1] > einfo[tid][-1][1]:einfo[tid].sort()
                for ex_ele in einfo[tid]:
                    if ex_ele[0] != contig_id:continue
                    exon_cnt += 1
                    exon_len += str(int(ex_ele[2])-int(ex_ele[1])+1) + ','
                    if fex == 0: # calculate the relative exon start 
                        exon_cod += '0,'
                        fex = 1
                        rstart = int(ex_ele[1])
                    else:
                        exon_cod += str(int(ex_ele[1])-rstart) + ','
                if exon_len: # display bed line
                    print contig_id + '\t' + tloc[0] + '\t' + tloc[1] + '\t' + tid + '\t' + tloc[2] + '\t' + tloc[-1] + '\t' + tloc[0] + '\t' + tloc[1] + '\t0\t' + str(exon_cnt) + '\t' + exon_len + '\t' + exon_cod 

def ParseAnno(gff_fh):
   
    tinfo, einfo = dict(), dict()
    for gff_line in gff_fh:
        gff_line = gff_line.strip('\n\r').split('\t')
        if not gff_line:continue ## not considering an empty line 
        if re.match(r'#', gff_line[0]) or re.match(r'>', gff_line[0]):continue ## not considering commented and FASTA header lines from GFF
        if len(gff_line) == 1: continue ## not considering if FASTA sequence along with GFF
        assert len(gff_line) == 9, '\t'.join(gff_line) ## a valid GFF line contains only 9 tab-delimited fields
        if '' in gff_line:continue ## empty fields in any line ? 
        if gff_line[2] == 'transcript' or gff_line[2] == 'scRNA' or gff_line[2] == "mRNA" or gff_line[2] == 'ncRNA' or gff_line[2] == 'miRNA' or gff_line[2] == 'rRNA' or gff_line[2] == 'snoRNA' or gff_line[2] == 'snRNA' or gff_line[2] == 'tRNA' or gff_line[2] == 'pseudogenic_transcript':
            col9 = gff_line[-1].split(';')
            tid = None
            for ele in col9:
                if re.search(r'ID=', ele):tid = re.search(r'ID=(.+)', ele).group(1);break
            if gff_line[0] in tinfo:
                tinfo[gff_line[0]][tid] = (gff_line[3], gff_line[4], gff_line[5], gff_line[6])
            else:
                tinfo[gff_line[0]] = {tid:(gff_line[3], gff_line[4], gff_line[5], gff_line[6])}
        if gff_line[2] == 'exon':
            col9 = gff_line[-1].split(';')
            pid = None
            for ele in col9:
                if re.search(r'Parent=', ele):pid = re.search(r'Parent=(.+)', ele).group(1);break
            if pid in einfo:
                einfo[pid].append((gff_line[0], int(gff_line[3]), int(gff_line[4])))
            else:
                einfo[pid] = [(gff_line[0], int(gff_line[3]), int(gff_line[4]))]
    gff_fh.close()
    return tinfo, einfo

if __name__ == "__main__": 

    try:
        gff_fh = open(sys.argv[1], 'rU')
    except:
        sys.stderr.write('\nGFF format file fail to open, Cannot continue...\nUSAGE: gff3_to_bed_converter.py <gff format file> > <bed format file>\n\n')
        sys.exit(-1)
    
    ## get transcript annotation
    tinfo, einfo = ParseAnno(gff_fh)
    ## write into bed format 
    WriteBED(tinfo, einfo)
