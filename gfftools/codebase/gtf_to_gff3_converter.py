#!/usr/bin/env python
"""
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Written (W) 2009-2011 Vipin T Sreedharan
Copyright (C) 2009-2011 Friedrich Miescher Laboratory of the Max Planck Society

Description : Convert data in Gene Transfer Format (GTF) to Generic Feature Format Version 3 (GFF3) file.
"""

import re, sys
from optparse import OptionParser

def addCDSphase(strand, cds):
    """Calculate CDS phase and add to the CDS exons
    """
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
    """Build UTR regions from a given set of CDS and exon coordiantes of a gene
    """
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

def GFFWriter(gff_fh, gtf_file_content):
    """Write into GFF3 
    """
    gff_fh.write('##gff-version 3\n')
    for contig, contig_info in sorted(gtf_file_content.items()): # first level, chromosome
        for feature, details in contig_info.items(): # second level, gene
            gene_start, gene_stop, transcript_details, orient, gname = [], [], dict(), None, None
            for ftid, tinfo in details.items(): # third level, transcripts 
                tinfo['exon'].sort() # generalize the coordinate system in ascending order. 
                tinfo['CDS'].sort()
                if tinfo['exon']:
                    gene_start.append(tinfo['exon'][0][0])
                    gene_stop.append(tinfo['exon'][-1][1])
                if not gene_start:
                    continue
                orient = tinfo['info'][1]
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
                if gname == None:
                    gname = feature # assign gene name as gene id, if not defined 
                pline = [str(contig[0]),
                        contig[1],
                        'gene',
                        str(gene_start[0]),
                        str(gene_stop[0]),
                        '.',
                        orient,
                        '.',
                        'ID=' + feature + ';Name=' + gname]
                gff_fh.write('\t'.join(pline) + '\n')
                for dtid, dinfo in transcript_details.items():
                    if dinfo['info'][4]:
                        pline = [str(contig[0]),
                                dinfo['info'][0],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][2],
                                orient,
                                '.',
                                'ID=' + dtid + ';Parent=' + feature + ';Name=' + dinfo['info'][4]]
                        gff_fh.write('\t'.join(pline) + '\n')
                    else:
                        pline = [str(contig[0]),
                                dinfo['info'][0],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][2],
                                orient,
                                '.',
                                'ID=' + dtid + ';Parent=' + feature]
                        gff_fh.write('\t'.join(pline) + '\n')
                    if 'utr5' in dinfo:
                        for ele in dinfo['utr5']:
                            pline = [str(contig[0]),
                                    dinfo['info'][0], 
                                    'five_prime_UTR',
                                    str(ele[0]), 
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n')
                    if 'cds' in dinfo:
                        cds_w_phase = addCDSphase(orient, dinfo['cds'])
                        for ele in cds_w_phase:
                            pline = [str(contig[0]),
                                    dinfo['info'][0],
                                    'CDS',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    str(ele[-1]),
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n')
                    if 'utr3' in dinfo:
                        for ele in dinfo['utr3']:
                            pline = [str(contig[0]),
                                    dinfo['info'][0],
                                    'three_prime_UTR',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n')
                    if 'exon' in dinfo:
                        for ele in dinfo['exon']:
                            pline = [str(contig[0]),
                                    dinfo['info'][0],
                                    'exon',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n')
    gff_fh.close()                            

def getGTFcontent(gtf_file):
    """Parse informations from a GTF file
    """
    gtf_content, recall = dict(), None
    for gtf_line in gtf_file:
        gtf_line = gtf_line.strip('\n\r').split('\t')
        if re.match(r'#|>', gtf_line[0]) or len(gtf_line)==1: # pass all commented and FASTA lines included in the GTF 
            continue 
        assert len(gtf_line) == 9, '\t'.join(gtf_line)
        gid, tid, attb_tag = None, None, dict()
        for attb in gtf_line[-1].split(';'):
            if re.search(r'^\s?$', attb):
                continue
            attb = re.sub('"', '', attb).strip() # trimming out quotes over id and white space in between the fields 
            attb = attb.split()
            if re.search(r'^(gene_id|geneid|name)$', attb[0], re.IGNORECASE): # TODO growing list of standard key words describing the identifiers.  
                gid = attb[1]
            elif re.search(r'^(transcript_id|transcriptId)$', attb[0], re.IGNORECASE):
                tid = attb[1]
            else:
                attb_tag[attb[0]] = attb[1]
        if gid == tid: # UCSC genome browser GTF files were having similar gene and transcript identifier 
            gid = 'Gene:'+str(gid) 
            tid = 'Transcript:'+str(tid)
        if tid == None: # JGI Joint Genome Institute GTF file dont have transcript ID for CDS line
            tid = recall 
        exon, cds, sp_cod, st_cod = [], [], [], []
        if re.search(r'^exon$', gtf_line[2], re.IGNORECASE): # TODO include GFF2 features name such that the line should be included for GFF3 contruction. lines like intron/utr's 
            exon = [(int(gtf_line[3]), int(gtf_line[4]))]
        elif re.search(r'^CDS$', gtf_line[2], re.IGNORECASE):
            cds = [(int(gtf_line[3]), int(gtf_line[4]))]
        elif re.search(r'^(stop_codon|stop-codon|stopcodon)$', gtf_line[2], re.IGNORECASE):
            sp_cod = [(int(gtf_line[3]), int(gtf_line[4]))]
        elif re.search(r'^(start_codon|start-codon|startcodon)$', gtf_line[2], re.IGNORECASE):
            st_cod = [(int(gtf_line[3]), int(gtf_line[4]))]
        else: # other lines are not required to contruct 
            continue
        if gtf_line[0] in gtf_content: # adding to existing chromosome
            if (gid, gtf_line[1]) in gtf_content[gtf_line[0]].keys(): # adding to existing gene 
                if tid in gtf_content[gtf_line[0]][(gid, gtf_line[1])].keys(): # adding to existing transcript
                    if exon:
                        gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid]['exon'].append(exon[0])
                    elif cds:
                        gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid]['CDS'].append(cds[0])
                    elif sp_cod:    
                        gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid]['sp_cod'].append(sp_cod[0])
                    elif st_cod:
                        gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid]['st_cod'].append(st_cod[0])
                else: # inserting new transcript
                    gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid] = dict(exon = exon, 
                                                            CDS = cds, 
                                                            sp_cod = sp_cod, 
                                                            st_cod = st_cod, 
                                                            info = [gtf_line[6], gtf_line[5], attb_tag])
            else: # inserting new gene 
                gtf_content[gtf_line[0]][(gid, gtf_line[1])] = {tid : dict(exon = exon, 
                                                    CDS = cds,
                                                    sp_cod = sp_cod, 
                                                    st_cod = st_cod,
                                                    info = [gtf_line[6], gtf_line[5], attb_tag])}
        else: # inserting new chromosome identifier 
            gtf_content[gtf_line[0]] = {(gid, gtf_line[1]) : {tid : dict(exon = exon, 
                                            CDS = cds,
                                            sp_cod = sp_cod, 
                                            st_cod = st_cod,
                                            info = [gtf_line[6], gtf_line[5], attb_tag])}}
        recall = tid
    gtf_file.close()
    return gtf_content

def stop_err(fmsg):
    """Error message
    """
    sys.stderr.write('%s\n' % fmsg)
    sys.exit(-1)

def __main__():
    cmd_arg = OptionParser()
    cmd_arg.add_option('', '-q', dest='query_file', help='Query file in Gene transfer format (GTF)')
    cmd_arg.add_option('', '-o', dest='result_file', help='Output file name with GFF3 extension')
    options, args = cmd_arg.parse_args()
    if len(sys.argv) < 2:
        cmd_arg.print_help()
        sys.exit(-1)
    try:
        gtf_fh = open(options.query_file, 'rU') 
    except Exception, erm:
        stop_err('Error reading query file ' + str(erm))
    try:
        gff_fh = open(options.result_file, 'w')
    except Exception, erm:
        stop_err('Error writing result file ' + str(erm))
    gtf_file_content = getGTFcontent(gtf_fh)
    #GFFWriter(gff_fh, gtf_file_content)

if __name__ == "__main__": __main__()
