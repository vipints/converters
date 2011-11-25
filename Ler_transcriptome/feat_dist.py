#!/usr/bin/python
# What is the status of Left and Right reads aligned to the genome, whether one belongs to rRNA and the other belongs to TE region. Find a statistics on aligned Paired-end reads.

import re, sys, os 

def AnnoDB(gff_file):

    tinfo, ribo = dict(), dict()
    gff_fh = open(gff_file, 'rU')
    for gff_line in gff_fh:
        gff_line = gff_line.strip('\n\r').split('\t')
        if re.match(r'#', gff_line[0]) or re.match(r'>', gff_line[0]) or len(gff_line) == 1:continue ## not considering commented and FASTA  lines from GFF
        if '' in gff_line:continue ## empty fields in any line ? 
        assert len(gff_line) == 9, '\t'.join(gff_line) ## a valid GFF line contains only 9 tab-delimited fields
        if gff_line[2] == 'transposable_element' or gff_line[2] == 'transposable_element_gene':
            if gff_line[0] in tinfo:
                tinfo[gff_line[0]].append((int(gff_line[3]), int(gff_line[4])))
            else:
                tinfo[gff_line[0]] = [(int(gff_line[3]), int(gff_line[4]))]
        if gff_line[2] == 'rRNA':
            if gff_line[0] in ribo:
                ribo[gff_line[0]].append((int(gff_line[3]), int(gff_line[4])))
            else:
                ribo[gff_line[0]] = [(int(gff_line[3]), int(gff_line[4]))]
    gff_fh.close()
    return tinfo, ribo

def AlignGenerator(bam_file):

    read_loc = dict()
    bfh = os.popen('/fml/ag-raetsch/share/software/samtools/samtools view ' + bam_file)
    for bline in bfh:
        bline = bline.strip('\n\r').split('\t')
        if bline[2] in read_loc:
            if bline[0] in read_loc[bline[2]].keys():
                read_loc[bline[2]][bline[0]].append(int(bline[3]))
            else:
                read_loc[bline[2]][bline[0]] = [int(bline[3])]
        else:
            read_loc[bline[2]] = {bline[0]:[int(bline[3])]}
    bfh.close()
    return read_loc

def DecideAlignment(pos, fAnno):
    
    dflag = 0 
    for featloc in fAnno:
        if featloc[0]-100 <= pos and pos <= featloc[1]:dflag = 1
    return dflag

def __main__():

    import time 
    print time.asctime( time.localtime(time.time()) )
    try:
        left_reads = sys.argv[1]
        right_reads = sys.argv[2]
        anno_file = sys.argv[3]
    except:
        sys.stderr.write('Please provide, Left reads and Right reads in BAM file + annotation in GFF\n')
        sys.exit(-1)

    teAnno, riboAnno = AnnoDB(anno_file) 
    #print len(teAnno['Chr3'])
    #print len(riboAnno['Chr3'])
    Leftq = AlignGenerator(left_reads)
    #print len(Leftq['Chr3'])
    Rightq = AlignGenerator(right_reads)
    #print len(Rightq['Chr3'])
    # check for the aligned region on the genome
    cnt_r1t2, cnt_t1r2 = 0, 0 
    for contig, q1 in Leftq.items():
        for fid, floc in q1.items():
            #print fid
            cid = re.sub(r'/1', r'/2', fid)
            #print cid
            try:
                for fp in floc:
                    #print fp 
                    r1flag = DecideAlignment(fp, riboAnno[contig])
                    t1flag = DecideAlignment(fp, teAnno[contig])
                    #print rflag 
                for cloc in Rightq[contig][cid]:
                    #print cloc
                    t2flag = DecideAlignment(fp, teAnno[contig])
                    r2flag = DecideAlignment(fp, riboAnno[contig])
                    #print tflag
                if r1flag == 1 and t2flag == 1: cnt_r1t2 += 1  
                if t1flag == 1 and r2flag == 1: cnt_t1r2 += 1  
            except:
                continue
            #break
    print 'L: rRNA, R: TE\t', cnt_r1t2
    print 'L: TE, R: rRNA\t', cnt_t1r2

    print time.asctime( time.localtime(time.time()) )

if __name__ == "__main__":__main__()
