#!/usr/bin/env python 
"""Create CSV format file with Genome annotation information which can be used for AnnoJ browser backend.
"""
import sys, re
from Core import GFF

def ParseGFF(fname):
    """get contents from a decent GFF file.
    """
    smap={1:'+' , -1:'-' }
    for cid in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrM', 'ChrC']: ## change according to the GFF file TODO 
        fh = open(fname, 'rU')
        limit_info = dict(gff_id=[cid], gff_source=['TAIR10']) ## change the source TODO According to the file type 
        for rec in GFF.parse(fh, limit_info=limit_info):
            for each_rec in rec.features:
                if each_rec.type=='gene':
                    for child in each_rec.sub_features:
                        print child.id, 'NULL', rec.id, smap[child.strand], child.type, child.location._start.position, child.location._end.position, each_rec.qualifiers['Note'][0]
                        cnt=1
                        for cod in child.sub_features:
                            fid=None
                            if cod.type=="five_prime_UTR":
                                fid=child.id+'_'+str(cnt)
                                cod.type='UTR5'
                            elif cod.type=="three_prime_UTR":
                                fid=child.id+'_'+str(cnt)
                                cod.type='UTR3'
                            elif cod.type=="CDS":
                                fid=child.id+'_'+str(cnt)
                            if fid:
                                print fid, child.id, rec.id, smap[cod.strand], cod.type, cod.location._start.position, cod.location._end.position,'NULL' 
                                cnt+=1
                    break
        fh.close()
                            
def __main__():
    try:
        gffname=sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    ParseGFF(gffname)

if __name__=="__main__":
    __main__()
