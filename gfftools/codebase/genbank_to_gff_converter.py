#!/usr/bin/env python
"""Convert genbank to gff.
USAGE:genbank_to_gff_converter.py in.gbk > file.gff
"""
import sys
from Bio import SeqIO

def gbk2gff(fname):
    """
    Converting data from genbank to GFF
    """
    for rec in SeqIO.parse(fname, "genbank"):
        strand, chr_id, mol_type, gname = None, None, None, None
        gene_start, gene_stop = 0, 0
        db_xref, cds, exon = [], [], []
        for ind_rec in rec.features:
            if ind_rec.type=='source':
                chr_id = ind_rec.qualifiers['chromosome'][0]
                mol_type = ind_rec.qualifiers['mol_type'][0]
            if ind_rec.type=='gene':
                if ind_rec.strand>0: strand = '+'
                else: strand = '-'
                gene_start = ind_rec.location._start.position+1
                gene_stop = ind_rec.location._end.position
                db_xref = ind_rec.qualifiers['db_xref']
                gname = ind_rec.qualifiers['gene'][0]
            if ind_rec.type=='exon':
                exon.append((ind_rec.location._start.position+1, ind_rec.location._end.position))
            if ind_rec.type=='CDS':
                cds.append((ind_rec.location._start.position+1, ind_rec.location._end.position))
        line = [str(chr_id), 
                'source',
                'gene',
                str(gene_start),
                str(gene_stop),
                '.',
                strand,
                '.',
                'ID='+str(gname)+';Note='+','.join(db_xref)]
        print '\t'.join(line) # gene line 
        line = [str(chr_id),
                'source',
                mol_type,
                str(gene_start),
                str(gene_stop),
                '.',
                strand,
                '.',
                'ID='+str(gname)+'.a;Parent='+str(gname)] # TODO introducing multiple transcript id's, get a better understanding of genbank structure.
        print '\t'.join(line) # transcript line 
        for ex in sorted(cds):
            line = [str(chr_id),
                    'source',
                    'CDS',
                    str(ex[0]),
                    str(ex[1]),
                    '.',
                    strand,
                    '.',
                    'Parent='+str(gname)+'.a']
            print '\t'.join(line) # CDS line 
        for ex in sorted(exon):
            line = [str(chr_id),
                    'source',
                    'exon',
                    str(ex[0]),
                    str(ex[1]),
                    '.',
                    strand,
                    '.',
                    'Parent='+str(gname)+'.a']
            print '\t'.join(line) # exon line 

if __name__=='__main__': 
    try:
        gb_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    gbk2gff(gb_fname)
