#!/usr/bin/env python
"""
Convert data in Gene Transfer Format (GTF) to Generic Feature Format Version 3 (GFF3) file.

Written (W) 2009-2012 Vipin T Sreedharan
Copyright (C) 2009-2012 Friedrich Miescher Laboratory of the Max Planck Society

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

def GFFWriter(gff_fh, gtf_file_cont):
    """Write into GFF3 
    """
    gff_fh.write('##gff-version 3\n')
    for contig, contig_info in sorted(gtf_file_cont.items()): # chromosome 
        for feature, details in contig_info.items(): # gene with source 
            gene_start, gene_stop, tnames, gnames, transcript_details = [], [], dict(), None, dict()
            for ftid, tinfo in details.items(): # transcripts 
                tinfo['exon'].sort() # coordinate system in ascending order. 
                tinfo['CDS'].sort()
                if tinfo['exon']:
                    gene_start.append(tinfo['exon'][0][0])
                    gene_stop.append(tinfo['exon'][-1][1])
                if not gene_start:
                    continue
                orient = tinfo['info'][0]
                tnames[ftid]=tinfo['info'][-1]
                gnames=tinfo['info'][-2]
                if len(tinfo['CDS']) == 0: # non coding transcript 
                    transcript_details[ftid] = dict(info = tinfo['info'], 
                                                exon = tinfo['exon'], 
                                                tpe = 'transcript')
                else:
                    if tinfo['sp_cod']: # stop codon are seperated from CDS, add the coordinates based on strand 
                        if orient == '+':
                            if tinfo['sp_cod'][0][0]-tinfo['CDS'][-1][1] == 1:
                                tinfo['CDS'][-1] = (tinfo['CDS'][-1][0], tinfo['sp_cod'][0][1])
                            else:
                                tinfo['CDS'].append(tinfo['sp_cod'][0])
                        if orient == '-':
                            if tinfo['CDS'][0][0]-tinfo['sp_cod'][0][1] == 1:
                                tinfo['CDS'][0] = (tinfo['sp_cod'][0][0], tinfo['CDS'][0][1])
                            else:
                                tinfo['CDS'].insert(0, tinfo['sp_cod'][0])
                    if tinfo['exon']:
                        utr5, utr3 = buildUTR(tinfo['CDS'], tinfo['exon'], orient) # getting UTR info from CDS and exon.
                        transcript_details[ftid] = dict(info = tinfo['info'], 
                                                    exon = tinfo['exon'], 
                                                    utr5 = utr5, 
                                                    utr3 = utr3, 
                                                    cds = tinfo['CDS'], 
                                                    tpe = 'mRNA')
            if gene_start and gene_stop: # displying Gene, transcript and subfeatures
                gene_start.sort()
                gene_stop.sort()
                if gnames == None:
                    gnames = feature[0] # assign gene name as gene id, if not defined 
                pline = [str(contig),
                        feature[1],
                        'gene',
                        str(gene_start[0]),
                        str(gene_stop[-1]),
                        '.',
                        orient,
                        '.',
                        'ID=' + feature[0] + ';Name=' + gnames]
                gff_fh.write('\t'.join(pline) + '\n') # writing gene line 
                for dtid, dinfo in transcript_details.items():
                    if dinfo['info'][3]:
                        pline = [str(contig),
                                feature[1],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][1],
                                orient,
                                '.',
                                'ID=' + str(dtid) + ';Parent=' + feature[0] + ';Name=' + str(dinfo['info'][3]) ]
                                #'ID=' + str(dtid) + ';Parent=' + feature[0] + ';Name=' + str(dinfo['info'][3]) +';Type='+dinfo['info'][4]]
                    else:
                        pline = [str(contig),
                                feature[1],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][1],
                                orient,
                                '.',
                                'ID=' + dtid + ';Parent=' + feature[0]]
                    gff_fh.write('\t'.join(pline) + '\n') # writing transcript line 
                    if 'utr5' in dinfo:
                        for ele in dinfo['utr5']:
                            pline = [str(contig),
                                    feature[1], 
                                    'five_prime_UTR',
                                    str(ele[0]), 
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n') # writing 5 prime UTR line 
                    if 'cds' in dinfo:
                        cds_w_phase = addCDSphase(orient, dinfo['cds'])
                        for ele in cds_w_phase:
                            pline = [str(contig),
                                    feature[1],
                                    'CDS',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    str(ele[-1]),
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n') # writing CDS line 
                    if 'utr3' in dinfo:
                        for ele in dinfo['utr3']:
                            pline = [str(contig),
                                    feature[1],
                                    'three_prime_UTR',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n') # writing 3 prime UTR line 
                    if 'exon' in dinfo:
                        intron_start = 0
                        for xq, ele in enumerate(dinfo['exon']):
                            pline = [str(contig),
                                    feature[1],
                                    'exon',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n') # writing exon line 
                            if xq ==0:
                                intron_start = ele[1]+1 
                                continue 
                            pline = [str(contig),
                                    feature[1],
                                    'intron',
                                    str(intron_start),
                                    str(ele[0]-1),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            gff_fh.write('\t'.join(pline) + '\n') # writing intron line 
                            intron_start = ele[1]+1 
                            #if len(dinfo['exon'])-1 == xq:
                            #    continue
                             
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
        if '' in gtf_line:
            continue
        if re.search(r'^(start_codon|start-codon|startcodon)$', gtf_line[2], re.IGNORECASE):
            continue
        gid, tid, gname, tname, ttype = None, None, None, None, None
        for attb in gtf_line[-1].split(';'):
            if re.search(r'^\s?$', attb):
                continue
            attb = re.sub('"', '', attb).strip() # trimming out quotes over id and white space in between the fields 
            attb = attb.split()
            if re.search(r'^(gene_id|geneid|name)$', attb[0], re.IGNORECASE): # TODO growing list of standard key words describing the identifiers.  
                gid = attb[1]
            elif re.search(r'^(transcript_id|transcriptId)$', attb[0], re.IGNORECASE):
                tid = attb[1]
            elif re.search(r'^(gene_name|genename)$', attb[0], re.IGNORECASE):
                gname = attb[1]
            elif re.search(r'^(transcript_name|transcriptname)$', attb[0], re.IGNORECASE):
                tname = attb[1]
            elif re.search(r'^(transcript_type)$', attb[0], re.IGNORECASE):
                ttype = attb[1]
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
        else: # other lines are not required to contruct GFF3 lines 
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
                else: # inserting new transcript
                    gtf_content[gtf_line[0]][(gid, gtf_line[1])][tid] = dict(exon = exon, 
                                                            CDS = cds, 
                                                            sp_cod = sp_cod, 
                                                            info = [gtf_line[6], gtf_line[5], gname, tname, ttype])
            else: # inserting new gene 
                gtf_content[gtf_line[0]][(gid, gtf_line[1])] = {tid : dict(exon = exon, 
                                                    CDS = cds,
                                                    sp_cod = sp_cod, 
                                                    info = [gtf_line[6], gtf_line[5], gname, tname, ttype])}
        else: # inserting new chromosome identifier 
            gtf_content[gtf_line[0]] = {(gid, gtf_line[1]) : {tid : dict(exon = exon, 
                                            CDS = cds,
                                            sp_cod = sp_cod, 
                                            info = [gtf_line[6], gtf_line[5], gname, tname, ttype])}}
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
    cmd_arg.add_option('', '-o', dest='result_file', help='Output file in Generic feature format version 3 (GFF3)')
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
    GFFWriter(gff_fh, gtf_file_content)

if __name__ == "__main__": __main__()
