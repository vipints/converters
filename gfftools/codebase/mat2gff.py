#!/usr/bin/python
"""
Program to convert .mat file to GFF3.

Usage: python mat2gff.py in.mat > out.gff  
Required modules:
    Numpy
    Scipy 
"""
import sys
import scipy.io
import numpy as np 

def __main__():

    try:
        mat_info = scipy.io.loadmat(sys.argv[1], squeeze_me=True, struct_as_record=False)
    except:
        print __doc__
        sys.exit(-1)

    gene_details = mat_info['genes'] #TODO automatically detect the struct identifier.
    
    for each_entry in gene_details: #Iterate over the matlab struct
        if each_entry.strain.size:
            geneLine = [str(each_entry.chr), 
                    str(np.atleast_1d(each_entry.strain)[0]),
                    str(each_entry.gene_type),
                    str(each_entry.start),
                    str(each_entry.stop),
                    '.',
                    str(each_entry.strand),
                    '.',
                    'ID='+str(each_entry.name)+';Name='+str(each_entry.name)]
        else:
            geneLine = [str(each_entry.chr), 
                    '.',
                    str(each_entry.gene_type),
                    str(each_entry.start),
                    str(each_entry.stop),
                    '.',
                    str(each_entry.strand),
                    '.',
                    'ID='+str(each_entry.name)+';Name='+str(each_entry.name)]
        print '\t'.join(geneLine) ## gene line in GFF3 
        tidx = 0 
        for transcript in np.atleast_1d(each_entry.transcripts):
            try:
                start = np.atleast_1d(each_entry.exons)[tidx][0][0]
                stop = np.atleast_1d(each_entry.exons)[tidx][-1][1]
            except:
                try:
                    start = np.atleast_1d(each_entry.exons)[0][0]
                    stop = np.atleast_1d(each_entry.exons)[-1][1]
                except:
                    start = np.atleast_1d(each_entry.exons)[0]
                    stop = np.atleast_1d(each_entry.exons)[1]
            if each_entry.strain.size:
                tLine = [str(each_entry.chr),
                    str(np.atleast_1d(each_entry.strain)[0]),
                    str(np.atleast_1d(each_entry.transcript_type)[tidx]),
                    str(start),
                    str(stop),
                    '.',
                    str(each_entry.strand),
                    '.',
                    'ID='+str(transcript)+';Parent='+str(each_entry.name)]
            else:
                tLine = [str(each_entry.chr),
                    '.',
                    str(np.atleast_1d(each_entry.transcript_type)[tidx]),
                    str(start),
                    str(stop),
                    '.',
                    str(each_entry.strand),
                    '.',
                    'ID='+str(transcript)+';Parent='+str(each_entry.name)]
            print '\t'.join(tLine) ## 
            try: ## UTR5 
                for eidx in range(len(each_entry.utr5_exons)):
                    if each_entry.strain.size:
                        u5Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[tidx][eidx][0]),
                                str(each_entry.utr5_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    else:
                        u5Line = [str(each_entry.chr),
                                '.',
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[tidx][eidx][0]),
                                str(each_entry.utr5_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    print '\t'.join(u5Line)
            except:
                try:
                    for eidx in range(len(each_entry.utr5_exons)):
                        if each_entry.strain.size:
                            u5Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[eidx][0]),
                                str(each_entry.utr5_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            u5Line = [str(each_entry.chr),
                                '.',
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[eidx][0]),
                                str(each_entry.utr5_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(u5Line)
                except:
                    if each_entry.utr5_exons.size:
                        if each_entry.strain.size:
                            u5Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[0]),
                                str(each_entry.utr5_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            u5Line = [str(each_entry.chr),
                                '.',
                                'five_prime_UTR',
                                str(each_entry.utr5_exons[0]),
                                str(each_entry.utr5_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(u5Line)
            try: ## CDS 
                for eidx in range(len(each_entry.cds_exons)):
                    if each_entry.strain.size:
                        cLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'CDS',
                                str(each_entry.cds_exons[tidx][eidx][0]),
                                str(each_entry.cds_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    else:
                        cLine = [str(each_entry.chr),
                                '.',
                                'CDS',
                                str(each_entry.cds_exons[tidx][eidx][0]),
                                str(each_entry.cds_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    print '\t'.join(cLine)
            except:
                try:
                    for eidx in range(len(each_entry.cds_exons)):
                        if each_entry.strain.size:
                            cLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'CDS',
                                str(each_entry.cds_exons[eidx][0]),
                                str(each_entry.cds_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            cLine = [str(each_entry.chr),
                                '.',
                                'CDS',
                                str(each_entry.cds_exons[eidx][0]),
                                str(each_entry.cds_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(cLine)
                except:
                    if each_entry.cds_exons.size:
                        if each_entry.strain.size:
                            cLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'CDS',
                                str(each_entry.cds_exons[0]),
                                str(each_entry.cds_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            cLine = [str(each_entry.chr),
                                '.',
                                'CDS',
                                str(each_entry.cds_exons[0]),
                                str(each_entry.cds_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(cLine)
            try: ## UTR3    
                for eidx in range(len(each_entry.utr3_exons)):
                    if each_entry.strain.size:
                        u3Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[tidx][eidx][0]),
                                str(each_entry.utr3_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    else:
                        u3Line = [str(each_entry.chr),
                                '.',
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[tidx][eidx][0]),
                                str(each_entry.utr3_exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    print '\t'.join(u3Line)
            except:
                try:
                    for eidx in range(len(each_entry.utr3_exons)):
                        if each_entry.strain.size:
                            u3Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[eidx][0]),
                                str(each_entry.utr3_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            u3Line = [str(each_entry.chr),
                                '.',
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[eidx][0]),
                                str(each_entry.utr3_exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(u3Line)
                except:
                    if each_entry.utr3_exons.size:
                        if each_entry.strain.size:
                            u3Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[0]),
                                str(each_entry.utr3_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            u3Line = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'three_prime_UTR',
                                str(each_entry.utr3_exons[0]),
                                str(each_entry.utr3_exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(u3Line)
            try: ## Exons 
                for eidx in range(len(each_entry.exons)):
                    if each_entry.strain.size:
                        eLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'exon',
                                str(each_entry.exons[tidx][eidx][0]),
                                str(each_entry.exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    else:            
                        eLine = [str(each_entry.chr),
                                '.',
                                'exon',
                                str(each_entry.exons[tidx][eidx][0]),
                                str(each_entry.exons[tidx][eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                    print '\t'.join(eLine)
            except:
                try:
                    for eidx in range(len(each_entry.exons)):
                        if each_entry.strain.size:
                            eLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'exon',
                                str(each_entry.exons[eidx][0]),
                                str(each_entry.exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            eLine = [str(each_entry.chr),
                                '.',
                                'exon',
                                str(each_entry.exons[eidx][0]),
                                str(each_entry.exons[eidx][1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(eLine)
                except:
                    if each_entry.exons.size:
                        if each_entry.strain.size:
                            eLine = [str(each_entry.chr),
                                str(np.atleast_1d(each_entry.strain)[0]),
                                'exon',
                                str(each_entry.exons[0]),
                                str(each_entry.exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        else:
                            eLine = [str(each_entry.chr),
                                '.',
                                'exon',
                                str(each_entry.exons[0]),
                                str(each_entry.exons[1]),
                                '.',
                                str(each_entry.strand),
                                '.',
                                'Parent='+str(transcript)]
                        print '\t'.join(eLine)
            tidx += 1

if __name__== "__main__":
    __main__()
