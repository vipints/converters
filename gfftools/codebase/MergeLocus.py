#!/usr/bin/env python
"""
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Written (W) 2010 Vipin T Sreedharan
Copyright (C) 2010 - 2012 Friedrich Miescher Laboratory of the Max Planck Society
Description : to merge same transcripts in single loci and define as an alternative spliced form for the gene.
Usage: MergeLocus.py in.gff3 out.gff3
"""

def display_content(final_dict):
    """displaying the summary from GFF file
    """
    print "\tUnique combination of Source(s), Feature type(s) and corresponding count:"
    for sftype, cnt in sorted(final_dict['gff_source_type'].items()):
        if sftype[1] == 'gene':
            print '\t' + str(cnt) + '\t' + str(sftype[0]) + ', '+ str(sftype[1])

def available_limits(gff_file):
    """Figure out the available feature types from the given GFF file
    """
    gff_handle = open(gff_file, 'rU')
    filter_info = dict(gff_id = [0], gff_source_type = [1, 2],
                gff_source = [1], gff_type = [2])
    cur_limits = dict()
    for filter_key in filter_info.keys():
        cur_limits[filter_key] = collections.defaultdict(int)
    for line in gff_handle:
        if line.strip('\n\r')[0] != "#":
            parts = [p.strip() for p in line.split('\t')]
            if len(parts) == 1 and re.search(r'\w+', parts[0]):continue # GFF files with FASTA sequence together 
            assert len(parts) == 9, line
            for filter_key, cur_indexes in filter_info.items():
                cur_id = tuple([parts[i] for i in cur_indexes])
                cur_limits[filter_key][cur_id] += 1
    gff_handle.close()
    final_dict = dict()
    for key, value_dict in cur_limits.items(): # get rid of the default dicts
        if len(key) == 1:key = key[0]
        final_dict[key] = dict(value_dict)
    return final_dict

def GFFWriter(merged_info, genes, transcripts, exons, utr5, cds, utr3, out_file):
    """Write GFF3 file with merged feature description
    """
    out_fh = open(out_file, 'w')
    for ginfo, regions in merged_info.items():
        gene_cnt = 1
        for interval, features in sorted(regions.items()):# master gene feature
            pline = [ginfo[0], 
                    ginfo[1],
                    'gene',
                    str(interval[0]),
                    str(interval[1]),
                    '.',
                    ginfo[2],
                    '.',
                    'ID=Gene_' + ginfo[0] + ginfo[2] + '_' + str(gene_cnt).zfill(5) + ';Name=Gene_' + ginfo[0] + ginfo[2] +'_' + str(gene_cnt).zfill(5)]
            out_fh.write('\t'.join(pline) + '\n')
            for geneid in features:# corresponding transcript info 
                if geneid in transcripts:
                    for tinfo in transcripts[geneid]:# transcript feature line 
                        pline = [ginfo[0],
                                ginfo[1],
                                tinfo['type'],
                                str(tinfo['start']),
                                str(tinfo['stop']),
                                '.',
                                ginfo[2],
                                '.',
                                'ID=' + tinfo['ID']+ ';Parent=Gene_' + ginfo[0] + ginfo[2] + '_' + str(gene_cnt).zfill(5)]
                        out_fh.write('\t'.join(pline) + '\n')
                        if tinfo['ID'] in utr5:# check for 5 prime UTR 
                            for u5info in utr5[tinfo['ID']]:
                                pline = [ginfo[0],
                                        ginfo[1],
                                        'five_prime_UTR',
                                        str(u5info['start']),
                                        str(u5info['stop']),
                                        '.',
                                        ginfo[2],
                                        '.',
                                        'Parent=' + tinfo['ID']]
                                out_fh.write('\t'.join(pline) + '\n')
                        if tinfo['ID'] in cds:# check for CDS 
                            for cdsinfo in cds[tinfo['ID']]:
                                pline = [ginfo[0],
                                        ginfo[1],
                                        'CDS',
                                        str(cdsinfo['start']),
                                        str(cdsinfo['stop']),
                                        '.',
                                        ginfo[2],
                                        '.',
                                        'Parent=' + tinfo['ID']]
                                out_fh.write('\t'.join(pline) + '\n')
                        if tinfo['ID'] in utr3:# check for 3 prime UTR
                            for u3info in utr3[tinfo['ID']]:
                                pline = [ginfo[0],
                                        ginfo[1],
                                        'three_prime_UTR',
                                        str(u3info['start']),
                                        str(u3info['stop']),
                                        '.',
                                        ginfo[2],
                                        '.',
                                        'Parent=' + tinfo['ID']]
                                out_fh.write('\t'.join(pline) + '\n')
                        if tinfo['ID'] in exons:# check for exons
                            for exinfo in exons[tinfo['ID']]:
                                pline = [ginfo[0],
                                        ginfo[1],
                                        'exon',
                                        str(exinfo['start']),
                                        str(exinfo['stop']),
                                        '.',
                                        ginfo[2],
                                        '.',
                                        'Parent=' + tinfo['ID']]
                                out_fh.write('\t'.join(pline) + '\n')
            gene_cnt += 1
    out_fh.close()

def UniqLoci(genes, transcripts, exons):
    """determine unique location where features annotated multiple times
    """
    uniq_loci = dict()
    for gid, parts in genes.items():
        gene_info = (parts['chr'], parts['source'], parts['strand'])
        if gene_info in uniq_loci: # same contig, orientation, source: look for merging transcripts based on the nearby location
            if (int(parts['start']), int(parts['stop'])) in uniq_loci[gene_info].keys(): ## similar transcripts will catch here (start and stop are same may be exon, CDS or intron content may vary) 
                uniq_loci[gene_info][(int(parts['start']), int(parts['stop']))].append(gid)
            else: # heuristic approach to include closely related region on a single master loci.
                got_a_range = 0
                for floc in uniq_loci[gene_info].keys():# look whether it lies closely to any intervel which is already defined  
                    if (floc[1]-parts['start']) < 150 or (parts['stop']-floc[0]) < 150:continue ## TODO boundary spanning length in same orientation for genes of each species will be great.
                    if floc[0] <= parts['start'] and parts['start'] < floc[1]: # the start of the new candidate is inside of any of the already defined interval ?
                        non_coding = 0
                        try: # check for small transcript whether they belong to a existing one or a new non-coding candidate.  
                            if len(transcripts[gid]) == 1: 
                                if len(exons[transcripts[gid][0]['ID']]) == 1:non_coding = 1
                            if non_coding == 0:
                                if parts['stop'] > floc[1]:# making global gene coordinate from individual transcript model
                                    entries = uniq_loci[gene_info][floc]
                                    del uniq_loci[gene_info][floc] # remove the existing interval, here we got a longer downstream position from the candidate
                                    entries.append(gid)
                                    uniq_loci[gene_info][(floc[0], parts['stop'])] = entries
                                else:
                                    uniq_loci[gene_info][floc].append(gid)
                            else:# create a new interval for non-coding type entry 
                                uniq_loci[gene_info][(parts['start'], parts['stop'])] = [gid] 
                            got_a_range = 1
                            break
                        except: # dont have any transcripts or exons defined.
                            break
                    elif floc[0] < parts['stop'] and parts['stop'] <= floc[1]: # the stop of the new candidate is inside of any of the pre-defined interval ? the candidate seems to be from more upstream
                        non_coding = 0
                        try:
                            if len(transcripts[gid]) == 1: 
                                if len(exons[transcripts[gid][0]['ID']]) == 1:non_coding = 1
                            if non_coding == 0:
                                entries = uniq_loci[gene_info][floc]
                                del uniq_loci[gene_info][floc] # remove the existing interval, here we got a upstream position from which the candidate transcribing 
                                entries.append(gid)
                                uniq_loci[gene_info][(int(parts['start']), floc[1])] = entries
                            else: # create a new interval for non-coding type entry 
                                uniq_loci[gene_info][(parts['start'], parts['stop'])] = [gid] 
                            got_a_range = 1
                            break
                        except:
                            break
                    elif floc[0] > parts['start'] and floc[1] < parts['stop']: # whether the whole feature floc region (--) resides in the candidate location (----------) ? 
                        non_coding = 0 # here the candidate seems to be longer than the pre-defined interval, check all entries from the pre-defined interval whether it is a small region, any chance as non-coding.
                        try:
                            for features in uniq_loci[gene_info][floc]:
                                if len(transcripts[features]) == 1:
                                    if len(exons[transcripts[features][0]['ID']]) == 1:non_coding = 1
                            if non_coding == 1: # create a new interval for non coding 
                                uniq_loci[gene_info][(parts['start'], parts['stop'])] = [gid]
                            else: # append the existing transcript cluster, here change the interval position based on the candidate location 
                                entries = uniq_loci[gene_info][floc]
                                del uniq_loci[gene_info][floc] # remove the existing interval, here we got a longer upstream and downstream region.
                                entries.append(gid)
                                uniq_loci[gene_info][(parts['start'], parts['stop'])] = entries
                            got_a_range = 1
                            break
                        except:
                            break
                ## or create a new interval ?? 
                if got_a_range == 0:
                    uniq_loci[gene_info][(parts['start'], parts['stop'])] = [gid]
        else:
            uniq_loci[gene_info] = {(int(parts['start']), int(parts['stop'])): [gid]}
    return uniq_loci

def ParseGFF(gff_file):
    """feature extraction from provided GFF file
    """
    gff_handle = open(gff_file, 'rU')
    genes, transcripts, exons, utr5, cds, utr3 = dict(), dict(), dict(), dict(), dict(), dict()
    for gff_line in gff_handle:
        parts = gff_line.strip('\n\r').split('\t')
        if gff_line[0] == '#' or gff_line[0] == '>':continue
        if len(parts) == 1:continue ## Few data centers create GFF files with FASTA sequence together 
        assert len(parts) == 9, '\t'.join(line)
        if parts[2] == 'gene':
            gene_info = dict()
            gene_info['start'] = int(parts[3])
            gene_info['stop'] = int(parts[4])
            gene_info['chr'] = parts[0]
            gene_info['source'] = parts[1]
            gene_info['strand'] = parts[6]
            gid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue ## GFF line may end with a ';' symbol
                attr = attr.split('=')
                if attr[0] == 'ID':gid=attr[1];continue 
                gene_info[attr[0]] = attr[1]
            if gid != '': genes[gid] = gene_info
        if parts[2] in ['mRNA', 'transcript', 'ncRNA', 'tRNA', 'snRNA', 'scRNA', 'snoRNA', 'snlRNA', 'rRNA', 'miRNA']:
            mrna_info = dict() 
            mrna_info['start'] = int(parts[3])
            mrna_info['stop'] = int(parts[4])
            mrna_info['chr'] =  parts[0]
            mrna_info['strand'] = parts[6]
            mrna_info['type'] = parts[2]
            gid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue ## GFF line may end with a ';' symbol
                attr = attr.split('=')
                if attr[0] == 'Parent':gid=attr[1];continue
                mrna_info[attr[0]] = attr[1]
            if gid in transcripts:
                transcripts[gid].append(mrna_info)
            else:
                transcripts[gid] = [mrna_info]
        if parts[2] == 'exon':
            exon_info = dict()
            exon_info['start'] = int(parts[3])
            exon_info['stop'] = int(parts[4])
            exon_info['chr'] =  parts[0]
            exon_info['strand'] = parts[6]
            tid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue ## GFF line may end with a ';' symbol
                attr = attr.split('=')
                if attr[0] == 'Parent':tid=attr[1];continue
                exon_info[attr[0]] = attr[1]
            if tid in exons:
                exons[tid].append(exon_info)
            else:
                exons[tid] = [exon_info]
        if parts[2] == 'five_prime_UTR':
            utr5_info = dict()
            utr5_info['start'] = int(parts[3])
            utr5_info['stop'] = int(parts[4])
            utr5_info['chr'] =  parts[0]
            utr5_info['strand'] = parts[6]
            tid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue ## GFF line may end with a ';' symbol
                attr = attr.split('=')
                if attr[0] == 'Parent':tid=attr[1];continue
                utr5_info[attr[0]] = attr[1]
            if tid in utr5:
                utr5[tid].append(utr5_info)
            else:
                utr5[tid] = [utr5_info]
        if parts[2] == 'CDS':
            cds_info = dict()
            cds_info['start'] = int(parts[3])
            cds_info['stop'] = int(parts[4])
            cds_info['chr'] =  parts[0]
            cds_info['strand'] = parts[6]
            tid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tid=attr[1];continue
                cds_info[attr[0]] = attr[1]
            if tid in cds:
                cds[tid].append(cds_info)
            else:
                cds[tid] = [cds_info]
        if parts[2] == 'three_prime_UTR':
            utr3_info = dict()
            utr3_info['start'] = int(parts[3])
            utr3_info['stop'] = int(parts[4])
            utr3_info['chr'] =  parts[0]
            utr3_info['strand'] = parts[6]
            tid = ''
            for attr in parts[-1].split(';'):
                if attr == '':continue
                attr = attr.split('=')
                if attr[0] == 'Parent':tid=attr[1];continue
                utr3_info[attr[0]] = attr[1]
            if tid in utr3:
                utr3[tid].append(utr3_info)
            else:
                utr3[tid] = [utr3_info]
    gff_handle.close()
    return genes, transcripts, exons, utr5, cds, utr3

import re, sys, time 
import collections

if __name__=='__main__':
    try:
        gff_file = sys.argv[1]
        out_file = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
    stime = time.asctime( time.localtime(time.time()) )
    print '-------------------------------------------------------'
    print 'MergeLoci started on ' + stime
    print '-------------------------------------------------------'
    print '--------'
    print 'Level: 1- ' + 'Reading GFF file: ' + gff_file
    print '--------'
    print '--------'
    print 'Level: 2- ' + 'BEFORE processing, Merging feature distribution in GFF file'
    print '--------'
    # initial feature distribution in file
    final_dict = available_limits(gff_file)
    display_content(final_dict)
    # determine the whole content from GFF file 
    genes, transcripts, exons, utr5, cds, utr3 = ParseGFF(gff_file)
    print '--------'
    print 'Level: 3- ' + 'Start merging feature(s) from similar locations...'
    print '--------'
    # determine the same gene loci on specific chromosome based on the same source
    merged_regions = UniqLoci(genes, transcripts, exons)
    print '\tDone.'
    print '--------'
    print 'Level: 4- ' + 'Writing merged feature annotation to GFF format...' 
    print '--------'
    # write new GFF file with merged loci information for gene feature
    GFFWriter(merged_regions, genes, transcripts, exons, utr5, cds, utr3, out_file)
    print '\tDone.'
    # after processing display the feature distribution in the result file 
    print '--------'
    print 'Level: 5- ' + 'Merged feature(s) summary from GFF file' 
    print '--------'
    final_dict = available_limits(out_file)
    display_content(final_dict)
    etime = time.asctime( time.localtime(time.time()) )
    print '-------------------------------------------------------'
    print 'MergeLoci finished at ' + etime
    print '-------------------------------------------------------'
