#!/usr/bin/python
# DESC 

import sys, re
from BCBio import GFF
    
def WriteSortedGFF(fname, limit_parse, sorted_pos):

    fgh = open(fname)
    for rec in GFF.parse(fgh, limit_info=limit_parse):
        for sort_pos in sorted_pos:
            for feature in rec.features:
                if sort_pos == (feature.location.start, feature.location.end):
                    orient = None
                    if feature.strand == 1:
                        orient = '+'
                    else:
                        orient = '-'
                    gline = [rec.id,
                            feature.qualifiers['source'][0],
                            feature.type,
                            str(feature.location._start.position + 1),
                            str(feature.location._end.position),
                            '.',
                            orient,
                            '.',
                            'ID=' + feature.qualifiers['ID'][0] + ';Name=' + feature.qualifiers['Name'][0]
                            ]
                    print '\t'.join(gline)
                    for slevel in feature.sub_features: ## second level feature 
                        orient = None
                        if slevel.strand == 1:
                            orient = '+'
                        else:
                            orient = '-'
                        tline = [rec.id,
                                slevel.qualifiers['source'][0],
                                slevel.type,
                                str(slevel.location._start.position + 1),
                                str(slevel.location._end.position),
                                str(float(slevel.qualifiers['score'][0])),
                                orient, 
                                '.',
                                'ID=' + slevel.qualifiers['ID'][0] + ';Parent=' + slevel.qualifiers['Parent'][0]
                                ]
                        print '\t'.join(tline)
                        for tlevel in slevel.sub_features:
                            orient = None
                            if tlevel.strand == 1:
                                orient = '+'
                            else:
                                orient = '-'
                            xline = []
                            if tlevel.type == 'CDS':
                                xline = [rec.id,
                                        tlevel.qualifiers['source'][0],
                                        tlevel.type,
                                        str(tlevel.location._start.position + 1),
                                        str(tlevel.location._end.position),
                                        '.',
                                        orient,
                                        tlevel.qualifiers['phase'][0],
                                        'Parent=' + tlevel.qualifiers['Parent'][0]
                                        ]
                            else:
                                xline = [rec.id,
                                        tlevel.qualifiers['source'][0],
                                        tlevel.type,
                                        str(tlevel.location._start.position + 1),
                                        str(tlevel.location._end.position),
                                        '.',
                                        orient,
                                        '.',
                                        'Parent=' + tlevel.qualifiers['Parent'][0]
                                        ]
                            print '\t'.join(xline)    
                    break
            #break ## one sorted position 
    fgh.close()

def __main__():
    
    try:
        gff_file = sys.argv[1]
        chrid = sys.argv[2]
        #source = sys.argv[3]
    except:
        sys.stderr.write('Access denied for a GFF file !\n')
        sys.exit(-1)
    
    limit_parse = dict(
                gff_id = [chrid],
                gff_source = ['rheMac2_ensGene', 'rheMac2_refGene', 'rheMac2_refSeqAnno', 'rheMac2_transMapAlnUcscGenes']
                )
    chrpos = []
    fh = open(gff_file)
    for rec in GFF.parse(fh, limit_info=limit_parse):
        for element in rec.features:
            chrpos.append((element.location.start, element.location.end))
    fh.close()
    chrpos.sort()

    WriteSortedGFF(gff_file, limit_parse, chrpos)


if __name__=="__main__":__main__()
