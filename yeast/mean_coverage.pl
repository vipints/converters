#!/usr/bin/perl -w
#Find the percentage of coverage of whole genome, here considering yeast genome.
# samtools pileup in.bam | awk '{print $4}' | perl mean_coverage.pl
use strict;
my ($num,$den)=(0,0);
while (my $_=<STDIN>){
    chomp $_;
    $num=$num+$_;
    $den++;
}
my $mean_cov=$num/$den;
print "Mean Coverage = $mean_cov\n";
