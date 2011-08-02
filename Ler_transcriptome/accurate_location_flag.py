#!/usr/bin/python
# DESC : Find the reads aligned to rRNA, Mitochondrial, Chloroplast chromosome regions of the Columbia-0 genome. This is with data set from Lisa Sequencing runs OLD.
import os, sys, re

def AnnoDB(gff_file):

    tinfo = dict()
    gff_fh = open(gff_file, 'rU')
    for gff_line in gff_fh:
        gff_line = gff_line.strip('\n\r').split('\t')
        if re.match(r'#', gff_line[0]) or re.match(r'>', gff_line[0]) or len(gff_line) == 1:continue ## not considering commented and FASTA  lines from GFF
        if '' in gff_line:continue ## empty fields in any line ? 
        assert len(gff_line) == 9, '\t'.join(gff_line) ## a valid GFF line contains only 9 tab-delimited fields

        if gff_line[2] == 'rRNA':
            tid = None
            for ele in gff_line[-1].split(';'):
                if re.search(r'ID=', ele):
                    tid = re.search(r'ID=(.+)', ele).group(1)
                    break
            if gff_line[0] in tinfo:
                tinfo[gff_line[0]].append((int(gff_line[3]), int(gff_line[4])))
            else:
                tinfo[gff_line[0]] = [(int(gff_line[3]), int(gff_line[4]))]
    gff_fh.close()
    return tinfo

def AlignGenerator(bam_file):
    
    bfh = os.popen('/fml/ag-raetsch/share/software/samtools/samtools view ' + bam_file)
    for bline in bfh:
        bline = bline.strip('\n\r').split('\t')
        yield bline
        #break
    bfh.close()

def __main__():
    import time 
    print time.asctime( time.localtime(time.time()) )
    try:
        bam_file = sys.argv[1]
        gff_file = sys.argv[2]
    except:
        print 'Missing BAM and GFF file'
        sys.exit(-1)
    anno_db = AnnoDB(gff_file)
    #print len(anno_db)
    contig = ['Chr1', 'Chr4', 'Chr5']
    rrna_fh = open('rRNA_COL.sam', 'w')
    chrc_fh = open('ChrC_COL.sam', 'w')
    chrm_fh = open('ChrM_COL.sam', 'w')
    generator = AlignGenerator(bam_file)
    for line in generator:
        if line[2] in contig:
            #print '\t'.join(line)
            continue
        if line[2] == 'Chr2' and 'NM:i:0' in line:
            for floc in anno_db[line[2]]:
                if int(line[3]) >= (floc[0]-40) and  int(line[3]) <= floc[1]:
                    line.append('XF:i:0')
                    rrna_fh.write('\t'.join(line) + '\n')
                    break
        elif line[2] == 'Chr3' and 'NM:i:0' in line:
            for floc in anno_db[line[2]]:
                if int(line[3]) >= (floc[0]-40) and  int(line[3]) <= floc[1]:
                    line.append('XF:i:0')
                    rrna_fh.write('\t'.join(line) + '\n')
                    break
        elif line[2] == 'ChrC' and 'NM:i:0' in line:
            mflag = 0 
            for floc in anno_db[line[2]]:
                if int(line[3]) >= (floc[0]-40) and  int(line[3]) <= floc[1]:
                    mflag = 1
                    break
            if mflag == 1: 
                line.append('XF:i:60')
            else:
                line.append('XF:i:6')
            chrc_fh.write('\t'.join(line) + '\n')
        elif line[2] == 'ChrM' and 'NM:i:0' in line:
            mflag = 0 
            for floc in anno_db[line[2]]:
                if int(line[3]) >= (floc[0]-40) and  int(line[3]) <= floc[1]:
                    mflag = 1
                    break
            if mflag == 1: 
                line.append('XF:i:70')
            else:
                line.append('XF:i:7')
            chrm_fh.write('\t'.join(line) + '\n')
    rrna_fh.close()
    chrc_fh.close()
    chrm_fh.close()
    print time.asctime( time.localtime(time.time()) )
if __name__=='__main__':__main__()
