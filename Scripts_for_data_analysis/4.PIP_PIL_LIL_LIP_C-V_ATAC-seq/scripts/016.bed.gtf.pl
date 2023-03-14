#!/bin/perl
# input narrowPeak/bed/gappedPeak
# output gtf for featureCount

my $i = 1;
while(<>){
    next if /^#/;
    chomp;
    my @F=split "\t", $_;
    my $chr=$F[0];
    my $start=$F[1];
    my $end=$F[2];
    my $peak_id="consensusPeak_$i";
    print join "\t", $chr, "PeakCaller", "peak", $start+1, $end, ".", ".", ".", "gene_id \"$peak_id\";\n";
    $i++;
}
