#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010 Vipin T Sreedharan
# Copyright (C) 2010 Friedrich Miescher Laboratory of the Max Planck Society
#
# Description : Convert a GTF format file to GFF3 format

import re, sys
from operator import itemgetter

def addCDSphase(strand, cds):
    """Add CDS phase to the CDS exons"""
    
    cds_region, cds_flag = [], 0 
    if strand == '+':
        for cdspos in cds:
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:
                xy = 0
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1: 
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1 
    elif strand == '-':
        cds.reverse()
        for cdspos in cds: 
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:  
                xy = 0 
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1:
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1
        cds_region.reverse()
    return cds_region

def buildUTR(cc, ec, strand):
    """Build UTR regions from CDS and exon coordiantes"""
    
    utr5, utr3 = [], []
    if strand == '+':
        cds_s = cc[0][0]
        for ex in ec:
            if ex[0] <= cds_s and cds_s <= ex[1]:
                if ex[0] != cds_s:utr5.append((ex[0], cds_s-1))
                break
            else:
                utr5.append(ex)
        cds_e = cc[-1][1]
        for i in range(len(ec)):
            i += 1
            if ec[-i][0] <= cds_e and cds_e <= ec[-i][1]:
                if ec[-i][1] != cds_e:utr3.append((cds_e +1, ec[-i][1]))
                break
            else:
                utr3.append(ec[-i]) 
        utr3.reverse()
    elif strand == '-':
        cds_s = cc[-1][1]
        for i in range(len(ec)):
            i += 1
            if ec[-i][0] <= cds_s and cds_s <= ec[-i][1]:
                if ec[-i][1] != cds_s:utr5.append((cds_s+1, ec[-i][1]))
                break
            else:
                utr5.append(ec[-i])
        utr5.reverse()
        cds_e = cc[0][0] 
        for ex in ec:
            if ex[0] <= cds_e and cds_e <= ex[1]:
                if ex[0] != cds_e:utr3.append((ex[0], cds_e-1))
                break
            else:
                utr3.append(ex)
    return utr5, utr3

def GFFWriter(gtf_file_content):
    """Write into GFF3 format"""

    for contig, contig_info in sorted(gtf_file_content.items()): # first level, chromosome
        for feature, details in contig_info.items(): # second level, gene
            source, gene_start, gene_stop, transcript_details, orient, gname = dict(), [], [], dict(), None, None
            for ftid, tinfo in details.items(): # third level, transcripts 
                tinfo['exon'].sort() # generalize the coordinate system in ascending order. 
                tinfo['CDS'].sort()
                if tinfo['exon']:gene_start.append(tinfo['exon'][0][0])
                if tinfo['exon']:gene_stop.append(tinfo['exon'][-1][1])
                if not gene_start:continue
                orient = tinfo['info'][1]
                if tinfo['info'][0] in source:
                    source[tinfo['info'][0]] += 1
                else:
                    source[tinfo['info'][0]] = 1
                gname = tinfo['info'][3]   
                if len(tinfo['CDS']) == 0: # non coding transcript 
                    transcript_details[ftid] = dict(info = tinfo['info'], exon = tinfo['exon'], tpe = 'transcript')
                else:
                    if tinfo['stop_codon']: # in GTF stop codon are seperated from CDS, here we are adding the stop codon to CDS region based on the strand.
                        if orient == '+':
                            if tinfo['stop_codon'][0][0]-tinfo['CDS'][-1][1] == 1:
                                tinfo['CDS'][-1] = (tinfo['CDS'][-1][0], tinfo['stop_codon'][0][1])
                            else:
                                tinfo['CDS'].append(tinfo['stop_codon'][0])
                        if orient == '-':
                            if tinfo['CDS'][0][0]-tinfo['stop_codon'][0][1] == 1:
                                tinfo['CDS'][0] = (tinfo['stop_codon'][0][0], tinfo['CDS'][0][1])
                            else:
                                tinfo['CDS'].insert(0, tinfo['stop_codon'][0])
                    if tinfo['exon']:
                        utr5, utr3 = buildUTR(tinfo['CDS'], tinfo['exon'], orient) # getting UTR info from CDS and exon.
                        transcript_details[ftid] = dict(info = tinfo['info'], exon = tinfo['exon'], utr5 = utr5, utr3 = utr3, cds = tinfo['CDS'], tpe = 'mRNA')
            
            if gene_start and gene_stop:# displying Gene, transcript and subfeatures
                gene_start.sort()
                gene_stop.sort()
                if gname == None:gname = feature # assign gene name as gene id, if not defined 
                if len(source) == 1: # to get apt source for gene feature 
                    for src in source:break
                else:
                    for src, freq in sorted(d.items(), key=itemgetter(1), reverse=True):break
                print contig + '\t' + src + '\tgene\t' + str(gene_start[0]) + '\t' + str(gene_stop[0]) + '\t.\t' + orient + '\t.\tID=' + feature + ';Name=' + gname
                for dtid, dinfo in transcript_details.items():
                    if dinfo['info'][4]:
                        print contig + '\t' + dinfo['info'][0] + '\t' + dinfo['tpe'] + '\t' + str(dinfo['exon'][0][0]) + '\t' + str(dinfo['exon'][-1][1]) + '\t' + dinfo['info'][2] + '\t' + orient + '\t.\tID=' + dtid + ';Parent=' + feature + ';Name=' + dinfo['info'][4]
                    else:
                        print contig + '\t' + dinfo['info'][0] + '\t' + dinfo['tpe'] + '\t' + str(dinfo['exon'][0][0]) + '\t' + str(dinfo['exon'][-1][1]) + '\t' + dinfo['info'][2] + '\t' + orient + '\t.\tID=' + dtid + ';Parent=' + feature
                    if 'utr5' in dinfo:
                        for ele in dinfo['utr5']:print contig + '\t' + dinfo['info'][0] + '\tfive_prime_UTR\t' + str(ele[0]) + '\t' + str(ele[1]) + '\t.\t' + orient + '\t.\tParent=' + dtid 
                    if 'cds' in dinfo:
                        cds_w_phase = addCDSphase(orient, dinfo['cds'])
                        for ele in cds_w_phase:print contig + '\t' + dinfo['info'][0] + '\tCDS\t' + str(ele[0]) + '\t' + str(ele[1]) + '\t.\t' + orient + '\t' + str(ele[-1]) +'\tParent=' + dtid
                    if 'utr3' in dinfo:
                        for ele in dinfo['utr3']:print contig + '\t' + dinfo['info'][0] + '\tthree_prime_UTR\t' + str(ele[0]) + '\t' + str(ele[1]) + '\t.\t' + orient + '\t.\tParent=' + dtid
                    if 'exon' in dinfo:
                        for ele in dinfo['exon']:print contig + '\t' + dinfo['info'][0] + '\texon\t' + str(ele[0]) + '\t' + str(ele[1]) + '\t.\t' + orient + '\t.\tParent=' + dtid


def get_GTF_info(gtf_file):
    """Parse informations from GTF file"""

    gtf_content = dict()
    gtf_fh = open(gtf_file, 'rU')
    recall = None
    for gtf_line in gtf_fh:
        gtf_line = gtf_line.strip('\n\r').split('\t')
        if re.match(r'#', gtf_line[0]):continue
        if re.match(r'>', gtf_line[0]):continue
        if len(gtf_line) == 1:continue
        assert len(gtf_line) == 9, '\t'.join(gtf_line)
        if gtf_line[2] == 'start_codon':continue

        attributes = gtf_line[-1].split(';')
        gid, tid, pid, gname, tname, ex_cnt = None, None, None, None, None, None
        for ele in attributes:
            if re.search(r'^\s?$', ele):continue
            if re.match(r'^\s', ele):ele = ele.lstrip()
            try:
                tag, name = re.search(r'^(\w+?)\s\"?(.+)\"?', ele).group(1, 2)
            except:
                tag, name = ele.split(' ')
            name = name.strip('"')
            if re.search(r'^(gene_id|geneid|name)$', tag, re.IGNORECASE):gid = name;continue
            if re.search(r'^(transcript_id|transcriptId)$', tag, re.IGNORECASE):tid = name;continue
            if re.search(r'^(protein_id|proteinid)$', tag, re.IGNORECASE):pid = name;continue
            if re.search(r'^(gene_name|genename)$', tag, re.IGNORECASE): gname= name;continue
            if re.search(r'^(transcript_name|transcriptname)$', tag, re.IGNORECASE):tname= name;continue
            if re.search(r'^(exonNumber|exon_number)$', tag, re.IGNORECASE):ex_cnt = name;continue

        if tid == None:
            if gtf_line[2] == 'CDS':tid = recall # JGI Joint Genome Institute GTF files dont have transcript ID for CDS line
        if tid == pid == None:continue # stop_codon icluded in CDS coordinates of JGI GTF files, moreover stop_codon lines dont have any transcript identifications.
        if gid == tid:gid = 'Gene:' + str(gid);tid = 'Transcript:' + str(tid) # UCSC gene and transcript ID are similar, differentaiting with Gene: and Transcript: tag in the begning 
        
        if gtf_line[0] in gtf_content: # existing chromosome
            if gid in gtf_content[gtf_line[0]].keys(): # existing gene 
                if tid in gtf_content[gtf_line[0]][gid].keys(): # existing transcript
                    if gtf_line[2] == 'exon':gtf_content[gtf_line[0]][gid][tid]['exon'].append((int(gtf_line[3]), int(gtf_line[4])))
                    elif gtf_line[2] == 'CDS':gtf_content[gtf_line[0]][gid][tid]['CDS'].append((int(gtf_line[3]), int(gtf_line[4])))
                    elif gtf_line[2] == 'stop_codon':gtf_content[gtf_line[0]][gid][tid]['stop_codon'].append((int(gtf_line[3]), int(gtf_line[4])))
                else: # new transcript 
                    if gtf_line[2] == 'exon':gtf_content[gtf_line[0]][gid][tid] = dict(exon = [(int(gtf_line[3]), int(gtf_line[4]))], CDS = [], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])
                    elif gtf_line[2] == 'CDS':gtf_content[gtf_line[0]][gid][tid] = dict(exon = [], CDS = [(int(gtf_line[3]), int(gtf_line[4]))], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])
                    elif gtf_line[2] == 'stop_codon':gtf_content[gtf_line[0]][gid][tid] = dict(exon = [], CDS = [], stop_codon = [(int(gtf_line[3]), int(gtf_line[4]))], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])
            else: # new gene  
                if gtf_line[2] == 'exon':gtf_content[gtf_line[0]][gid] = {tid : dict(exon = [(int(gtf_line[3]), int(gtf_line[4]))], CDS = [], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}
                elif gtf_line[2] == 'CDS':gtf_content[gtf_line[0]][gid] = {tid : dict(exon = [], CDS = [(int(gtf_line[3]), int(gtf_line[4]))], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}
                elif gtf_line[2] == 'stop_codon':gtf_content[gtf_line[0]][gid] = {tid : dict(exon = [], CDS = [], stop_codon = [(int(gtf_line[3]), int(gtf_line[4]))], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}
        else: # new chromosome
            if gtf_line[2] == 'exon':gtf_content[gtf_line[0]] = {gid : {tid : dict(exon = [(int(gtf_line[3]), int(gtf_line[4]))], CDS = [], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}}
            elif gtf_line[2] == 'CDS':gtf_content[gtf_line[0]] = {gid : {tid : dict(exon = [], CDS = [(int(gtf_line[3]), int(gtf_line[4]))], stop_codon = [], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}}
            elif gtf_line[2] == 'stop_codon':gtf_content[gtf_line[0]] = {gid : {tid : dict(exon = [], CDS = [], stop_codon = [(int(gtf_line[3]), int(gtf_line[4]))], info = [gtf_line[1], gtf_line[6], gtf_line[5], gname, tname])}}
        recall = tid

    gtf_fh.close()
    return gtf_content

def __main__():

    try:
        gtf_file = sys.argv[1]
    except:
        sys.stderr.write('\nGTF format file fail to open, Cannot continue...\n\tUSAGE: gtf_to_gff3_converter.py <file in gtf> > <file in GFF3>\n\n')
        sys.exit(-1)
    
    gtf_file_content = get_GTF_info(gtf_file)
    print '##gff-version 3'
    GFFWriter(gtf_file_content)

if __name__ == "__main__": __main__()
