#!/usr/bin/env perl
#	Description: - read fixedStep and variableStep wiggle input data,
#		output four column bedGraph format data
#   Author: Vipin

use warnings;
use strict;

my ($position, $chr, $step, $span) = (0, "", 1, 1);

my $usage = q(
fixedStep2BED.pl  - Program to convert a valid fixedStep or variableStep wiggle input data to BED format.
USAGE: fixedStep2BED.pl <fixedStep/variableStep Wiggle format> > <output file name>

);

if (scalar(@ARGV) != 1) {
    print $usage;
    exit
}

my $inFile = $ARGV[0];
open WIG, "<$inFile" || die "Can't open the Wiggle file \n";
while (my $dataValue = <WIG>)
{
    chomp $dataValue;
    if ( $dataValue =~ m/^track /) {
        print STDERR "Skipping track line\n";
        next;
    }
    if ( $dataValue =~ m/^fixedStep/ || $dataValue =~ m/^variableStep/ ) {
        $position = 0;
        $chr = "";
        $step = 1;
        $span = 1;
        my @a = split('\t', $dataValue);
        for (my $i = 1; $i < scalar(@a); ++$i) {
            my ($ident, $value) = split('=',$a[$i]);
            if ($ident =~ m/chrom/) { $chr = $value; }
            elsif ($ident =~ m/start/) { $position = $value-1; }
            elsif ($ident =~ m/step/) { $step = $value; }
            elsif ($ident =~ m/span/) { $span = $value; }
            else {
                print STDERR "ERROR: unrecognized fixedStep line: $dataValue\n";
                exit 255;
            }
        }
    } 
    else {
        my @b = split('\s', $dataValue);
        if (scalar(@b)>1) {
            $position = $b[0];
            $dataValue = $b[1];
        }
        my $loc_pos = $position+$span;
        print "$chr\t$position\t$loc_pos\t$dataValue\n";
        $position = $position + $step;
    }
}
close WIG;
exit;
