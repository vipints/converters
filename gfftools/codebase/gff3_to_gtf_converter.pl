#!/usr/bin/env perl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010-2011 Vipin T Sreedharan, Friedrich Miescher Laboratory of the Max Planck Society
# Copyright (C) 2010 Max Planck Society
#
# Description : Convert data in Generic Feature Format Version 3 (GFF3) to Gene Transfer Format (GTF).
 
use strict;
use warnings;
use Bio::FeatureIO;

my $usage = q(
gff3_to_gtf_converter.pl  - Program to convert a valid GFF3 format file to GTF format.
USAGE: gff3_to_gtf_converter.pl <GFF3 file name> <output file name>
);

if (scalar(@ARGV) != 2) {
    print $usage;
    exit
}
my $inFile = $ARGV[0]; 
my $outFile = $ARGV[1];
my $inGFF = Bio::FeatureIO->new( '-file' => "$inFile",
 '-format' => 'GFF',
 '-version' => 3 );
my $outGTF = Bio::FeatureIO->new( '-file' => ">$outFile",
 '-format' => 'GFF',
 '-version' => 2.5);

my ($gene, $feature_type, $exon_exon_cnt, $cds_exon_cnt) = ('', '', 0, 0);

while (my $feature = $inGFF->next_feature() ) {
    if ($feature->type->name eq 'mRNA' || $feature->type->name eq 'miRNA' || $feature->type->name eq 'ncRNA' || $feature->type->name eq 'rRNA' || $feature->type->name eq 'snoRNA' || $feature->type->name eq 'snRNA' || $feature->type->name eq 'tRNA' || $feature->type->name eq 'misc_RNA' || $feature->type->name eq 'processed_transcript' || $feature->type->name eq 'transcript' || $feature->type->name eq 'scRNA') {
        my $parent = ($feature->get_Annotations('Parent'))[0]; 
        $gene = $parent->value;
        ($cds_exon_cnt, $exon_exon_cnt) = (1, 1); 
        $feature_type = $feature->type->name;
    }
    if ($feature->type->name eq 'exon' || $feature->type->name eq 'CDS' ||$feature->type->name eq 'stop_codon' || $feature->type->name eq 'start_codon') {
        my $parent = ($feature->get_Annotations('Parent'))[0];
        my $transcript = $parent->value;
       
        my ($col_exon_number, $protein_id) = ('', '');
        if ($feature->type->name eq 'exon') {
            $col_exon_number = Bio::Annotation::SimpleValue->new( '-value' => $exon_exon_cnt, '-tagname' => 'exon_number');
            $exon_exon_cnt++;
        } elsif ($feature->type->name eq 'CDS') {
            $col_exon_number = Bio::Annotation::SimpleValue->new( '-value' => $cds_exon_cnt, '-tagname' => 'exon_number');
            $cds_exon_cnt++;
            my $pid = $transcript;
            $pid =~s/Transcript/Protein/;
            $protein_id = Bio::Annotation::SimpleValue->new( '-value' => $pid, '-tagname' => 'protein_id');
        } elsif ($feature->type->name eq 'start_codon'){
            $col_exon_number = Bio::Annotation::SimpleValue->new( '-value' => $cds_exon_cnt, '-tagname' => 'exon_number');
        } else {
            $col_exon_number = Bio::Annotation::SimpleValue->new( '-value' => $cds_exon_cnt-1, '-tagname' => 'exon_number');
        }
        my $transcript_id = Bio::Annotation::SimpleValue->new( '-value' => $transcript, '-tagname' => 'transcript_id');
        my $gene_id = Bio::Annotation::SimpleValue->new( '-value' => $gene, '-tagname' => 'gene_id');
        my $source_identifer = Bio::Annotation::SimpleValue->new( '-value' => $feature_type, '-tagname' => 'source_id');
        $feature->add_Annotation($source_identifer);
        $feature->add_Annotation($transcript_id);
        $feature->add_Annotation($gene_id);
        $feature->add_Annotation($col_exon_number);
        $feature->add_Annotation($protein_id) if ($feature->type->name eq 'CDS');
    }
    $outGTF->write_feature($feature);
}
exit;
