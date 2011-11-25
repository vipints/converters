#!/usr/bin/python
# Formatting FASTA file, it takes a unformatted FASTA file and results a Nice formatted FASTA file.
# Take a GFF file and create FASTA file according to the features present in it.

import re, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

def GFFParse(fname):

    featuredb, tegenedb = dict(), dict()
    fh = open(fname, 'rU')
    for line in fh:
        line = line.strip('\n\r').split('\t')
        if re.match(r'#', line[0]) or re.match(r'>', line[0]):continue
        if len(line) == 1:continue
        if '' in line:continue
        assert len(line) == 9, '\t'.join(line)
        """if line[2] == 'gene': 
            gtype, gid = '', ''
            for attb in line[-1].split(';'):
                if re.search(r'^ID=', attb):
                    gid = re.search(r'^ID=(.+)', attb).group(1)
                if re.search(r'^Note=', attb):
                    gtype = re.search(r'^Note=(.+)', attb).group(1)
            #if gtype == 'miRNA' or gtype =='other_RNA' or gtype =='pseudogene' or gtype =='snoRNA' or gtype =='snRNA' or gtype == 'tRNA':
            if gtype == 'protein_coding_gene':
                if line[0] in featuredb:
                    featuredb[line[0]].append((int(line[3]), int(line[4]), gid, line[6]))
                else:
                    featuredb[line[0]] = [(int(line[3]), int(line[4]), gid, line[6])]"""
        
        if line[2] == 'transposable_element':
            gid = ''
            for attb in line[-1].split(';'):
                if re.search(r'^ID=', attb):
                    gid = re.search(r'^ID=(.+)', attb).group(1)
            if line[0] in featuredb:
                featuredb[line[0]][(int(line[3]), int(line[4]))] = (gid, line[6])
            else:
                featuredb[line[0]] = {(int(line[3]), int(line[4])): (gid, line[6])}
        if line[2] == 'transposable_element_gene':
            gid = ''
            for attb in line[-1].split(';'):
                if re.search(r'^ID=', attb):
                    gid = re.search(r'^ID=(.+)', attb).group(1)
            if line[0] in tegenedb:
                tegenedb[line[0]][(int(line[3]), int(line[4]))] = (gid, line[6])
            else:
                tegenedb[line[0]] = {(int(line[3]), int(line[4])): (gid, line[6])}
    fh.close()
    for chrid, fte in featuredb.items():
        delete_key = []
        if chrid in tegenedb:
            for ftg in tegenedb[chrid]:
                for fs, fe in fte:
                    if ftg[0] >= fs and ftg[1] <= fe:
                        delete_key.append(ftg)
                        break
                    #break
                #break
        #break
        for delp in delete_key:
            del tegenedb[chrid][delp]

    return featuredb, tegenedb
        
def getSeq(gseq, cid, ftypes):

    fah = open(gseq, 'rU')
    for rec in SeqIO.parse(fah, 'fasta'):
        if rec.id == cid:
            for fp in ftypes:
                motif_seq = rec.seq[fp[0]-1:fp[1]]
                motif_seq = motif_seq.tomutable()
                motif_seq = 100*'N' + motif_seq + 100*'N'
                fseq = SeqRecord(motif_seq, id=fp[-1], description='A_thaliana_Other_RNA')
                #fseq = SeqRecord(motif_seq, id=ftypes[fp], description='A_thaliana_TE')
                print fseq.format("fasta")

def __main__():
    
    try:
        gff_file = sys.argv[1]
        #fasta_seq = sys.argv[2]
    except: 
        sys.stderr.write('Access denied for the input GFF file ! terminating\n')
        sys.exit(-1)

    featuredb, tegenedb = GFFParse(gff_file)
    
    for contig, feature in featuredb.items():
        for fl, fp in feature.items():
            print fp[0] + '\t' + fp[1]
            #break
    for contig, feature in tegenedb.items():
        for fl, fp in feature.items():
            print fp[0] + '\t' + fp[1]
            #break
    #for contig, feature in featuredb.items():
    #    getSeq(fasta_seq, contig, feature)

    #for contig, feature in tegenedb.items():
    #    getSeq(fasta_seq, contig, feature)

    """fah = open(fasta_raw, 'rU')
    for rec in SeqIO.parse(fah, 'fasta'):
        rna_genome = rec.seq.tomutable()
        i = 0 
        for nt in rna_genome:
            if nt == 'U':
                rna_genome[i] = 'T'
            i += 1
        rna_genome = 100*'N' + rna_genome + 100*'N' ## adding a non-nucleotide clips at both ends of the sequence.

        gseq = SeqRecord(rna_genome, id=rec.id, description='A_thaliana_ribosomal_RNA')
        print gseq.format("fasta")
        #break
    fah.close()"""

if __name__=="__main__":__main__()
