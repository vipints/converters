#!/usr/bin/env python
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Vipin T Sreedharan, Friedrich Miescher Laboratory
# Copyright (C) 2009-2011 Max Planck Society
#
# Description : Provide a mapping of parent to child relationships in a Generic Feature Format Version 3 file. 

import re, sys 
import collections
import urllib
import time 

def _gff_line_map(line):
    """Parses a line of GFF into a dictionary.
    Given an input line from a GFF file, this:
        - breaks it into component elements
        - determines the type of attribute (flat, parent, child or annotation)
        - generates a dictionary of GFF info 
    """
    gff3_kw_pat = re.compile("\w+=")
    def _split_keyvals(keyval_str):
        """Split key-value pairs in a GFF2, GTF and GFF3 compatible way.

        GFF3 has key value pairs like:
          count=9;gene=amx-2;sequence=SAGE:aacggagccg
        GFF2 and GTF have:           
          Sequence "Y74C9A" ; Note "Clone Y74C9A; Genbank AC024206"
          name "fgenesh1_pg.C_chr_1000003"; transcriptId 869
        """
        quals = collections.defaultdict(list)
        if keyval_str is None:
            return quals
        # ensembl GTF has a stray semi-colon at the end
        if keyval_str[-1] == ';':
            keyval_str = keyval_str[:-1]
        # GFF2/GTF has a semi-colon with at least one space after it.
        # It can have spaces on both sides; wormbase does this.
        # GFF3 works with no spaces.
        # Split at the first one we can recognize as working
        parts = keyval_str.split(" ; ")
        if len(parts) == 1:
            parts = keyval_str.split("; ")
            if len(parts) == 1:
                parts = keyval_str.split(";")
        # check if we have GFF3 style key-vals (with =)
        is_gff2 = True
        if gff3_kw_pat.match(parts[0]):
            is_gff2 = False
            key_vals = [p.split('=') for p in parts]
        # otherwise, we are separated by a space with a key as the first item
        else:
            pieces = []
            for p in parts:
                # fix misplaced semi-colons in keys in some GFF2 files
                if p and p[0] == ';':
                    p = p[1:]
                pieces.append(p.strip().split(" "))
            key_vals = [(p[0], " ".join(p[1:])) for p in pieces]
        for key, val in key_vals:
            # remove quotes in GFF2 files
            if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
                val = val[1:-1] 
            if val:
                quals[key].extend(val.split(','))
            # if we don't have a value, make this a key=True/False style
            # attribute
            else:
                quals[key].append('true')
        for key, vals in quals.items():
            quals[key] = [urllib.unquote(v) for v in vals]
        return quals, is_gff2

    def _nest_gff2_features(gff_parts):
        """Provide nesting of GFF2 transcript parts with transcript IDs.

        exons and coding sequences are mapped to a parent with a transcript_id
        in GFF2. This is implemented differently at different genome centers
        and this function attempts to resolve that and map things to the GFF3
        way of doing them.
        """
        # map protein or transcript ids to a parent
        for transcript_id in ["transcript_id", "transcriptId", "proteinId"]:
            try:
                gff_parts["quals"]["Parent"] = \
                        gff_parts["quals"][transcript_id]
                break
            except KeyError:
                pass
        # case for WormBase GFF -- everything labelled as Transcript or CDS
        for flat_name in ["Transcript", "CDS"]:
            if gff_parts["quals"].has_key(flat_name):
                # parent types
                if gff_parts["type"] in [flat_name]:
                    if not gff_parts["id"]:
                        gff_parts["id"] = gff_parts["quals"][flat_name][0]
                        gff_parts["quals"]["ID"] = [gff_parts["id"]]
                # children types
                elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                        "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                        "start_codon"]:
                    gff_parts["quals"]["Parent"] = gff_parts["quals"][flat_name]
                break

        return gff_parts

    line = line.strip()
    if line == '':return [('directive', line)] # sometimes the blank lines will be there 
    if line[0] == '>':return [('directive', '')] # sometimes it will be a FATSA header
    if line[0] == "#":
        return [('directive', line[2:])]
    elif line:
        parts = line.split('\t')
        if len(parts) == 1 and re.search(r'\w+', parts[0]):return [('directive', '')] ## GFF files with FASTA sequence together 
        assert len(parts) == 9, line
        gff_parts = [(None if p == '.' else p) for p in parts]
        gff_info = dict()
            
        # collect all of the base qualifiers for this item
        quals, is_gff2 = _split_keyvals(gff_parts[8])

        gff_info["is_gff2"] = is_gff2

        if gff_parts[1]:quals["source"].append(gff_parts[1])
        gff_info['quals'] = dict(quals)

        # if we are describing a location, then we are a feature
        if gff_parts[3] and gff_parts[4]:
            gff_info['type'] = gff_parts[2]
            gff_info['id'] = quals.get('ID', [''])[0]
            
            if is_gff2:gff_info = _nest_gff2_features(gff_info)
            # features that have parents need to link so we can pick up
            # the relationship
            if gff_info['quals'].has_key('Parent'):
                final_key = 'child'
            elif gff_info['id']:
                final_key = 'parent'
            # Handle flat features
            else:
                final_key = 'feature'
        # otherwise, associate these annotations with the full record
        else:
            final_key = 'annotation'
        return [(final_key, gff_info)]
    
def parent_child_id_map(gff_handle):
    """Provide a mapping of parent to child relationships in the file.

    Gives a dictionary of parent child relationships:

    keys -- tuple of (source, type) for each parent
    values -- tuple of (source, type) as children of that parent"""

    # collect all of the parent and child types mapped to IDs
    parent_sts = dict()
    child_sts = collections.defaultdict(list)

    for line in gff_handle:
        line_type, line_info = _gff_line_map(line)[0]
        if (line_type == 'parent' or (line_type == 'child' and line_info['id'])):
            parent_sts[line_info['id']] = (line_info['quals']['source'][0], line_info['type'])
        if line_type == 'child':
            for parent_id in line_info['quals']['Parent']:
                child_sts[parent_id].append((line_info['quals']['source'][0], line_info['type']))
    gff_handle.close()
     
    # generate a dictionary of the unique final type relationships
    pc_map = collections.defaultdict(list)
    for parent_id, parent_type in parent_sts.items():
        for child_type in child_sts[parent_id]:
            pc_map[parent_type].append(child_type)
    pc_final_map = dict()
    for ptype, ctypes in pc_map.items():
        unique_ctypes = list(set(ctypes))
        unique_ctypes.sort()
        pc_final_map[ptype] = unique_ctypes

    # Check for Parent Child relations
    level1, level2, level3, sec_level_mis = {}, {}, {}, {}
    for etype, fchild in pc_final_map.items():
        level2_flag = 0
        for kp, vp in pc_final_map.items():
            if etype in vp:level2_flag = 1; level2[etype] = 1 # check for second level features
        if level2_flag == 0: # first level features
            level1[etype] =1 
            for eachfch in fchild: # perform a check for all level1 objects values were defined as level2 keys.  
                if not eachfch in pc_final_map.keys(): # figure out the missing level2 objects  
                    if etype in sec_level_mis:
                        sec_level_mis[etype].append(eachfch)
                    else:
                        sec_level_mis[etype]=[eachfch]
        if level2_flag == 1:level3[str(fchild)] =1  # taking third level features 
    # disply the result 
    print
    print 'Found following level features from the given GFF file'
    print
    if level1==level2==level3=={} and sec_level_mis == {}:
        print 'FIRST level feature(s):'
        source_type = dict()
        gff_handle = open(gff_handle.name, 'rU')
        for line in gff_handle:
            line = line.strip('\n\r')
            if line[0] == '#': continue
            parts = line.split('\t')
            if parts[-1] == '':parts.pop()
            assert len(parts) == 9, line
            source_type[(parts[1], parts[2])] = 1
        gff_handle.close()
        for ele in source_type:print '\t' + str(ele)
        print
    else:
        print 'FIRST level feature(s):'
        for ele in level1: print '\t' + str(ele)
        print 
        print 'SECOND level feature(s):'
        for ele in level2: print '\t' + str(ele)
        print 
        print 'THIRD level feature(s):'
        for ele in level3:print '\t' + str(ele[1:-1])
        print
        # wrong way mapped feature mapping description 
        for wf, wfv in sec_level_mis.items():
            if wf[1]=='gene':
                print '**FYI**'
                print 'GFF Parsing modules from publicly available packages like Bio-python, Bio-perl etc. are heavily dependent on feature identifier mapping.' 
                print 'Here few features seems to be wrongly mapped to its child, which inturn cause problems while extracting the annotation based on feature identifier.'
                print 
                for ehv in wfv:
                    if ehv[1]=='exon' or ehv[1]=='intron' or ehv[1]=='CDS' or ehv[1]=='three_prime_UTR' or ehv[1]=='five_prime_UTR':
                        print 'Error in ID mapping: Level1 feature ' + str(wf) + ' maps to Level3 feature ' + str(ehv)
                print 

if __name__=='__main__':

    try:
        gff_handle = open(sys.argv[1], 'rU')
    except:
        sys.stderr.write("Access denied for a GFF file, Please check the query file\nUSAGE: gff_examiner.py <gff file>\n")
        sys.exit(-1)

    parent_child_id_map(gff_handle)
