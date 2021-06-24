#!/usr/bin/env perl
#
#  Generates a FASTA with sequence of a known length, to test
#  `seqrequester summarize`.
#

use strict;

my $prefix = shift @ARGV;
my @lengths;
my $tLength = 0;

die "usage: $0 output-prefix\n"   if ($prefix eq "");

open(F, "> $prefix.fasta");
for (my $n=0; $n < 100000; $n++) {
    my $small = rand(1);
    my $lseed = rand(1);
    my $length;

    if    ($small < 0.0000) {
        $length = int($lseed * 128 * 1048576) + 128 * 1048576;
    }
    elsif ($small < 0.9) {
        $length = int($lseed * 30) + 1;
    }
    else {
        $length = int($lseed * 1024) + 1024;
    }

    push @lengths, $length;

    print F ">\n";
    print F "A" x $length;
    print F "\n";

    $tLength += $length;
}
close(F);

open(F, "> $prefix.lengths");
foreach my $l (@lengths) {
    print F "$l\n";
}
close(F);

print STDERR "Wrote ", scalar(@lengths), " sequences with total length $tLength.\n";

exit(0);
