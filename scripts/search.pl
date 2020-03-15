#!/usr/bin/env perl

use strict;
use Time::HiRes qw(usleep nanosleep);

my $order  = 0;
my $weight = 0;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if ($arg eq "-order") {
        $order = shift @ARGV;
    }

    if ($arg eq "-weight") {
        $weight = shift @ARGV;
    }
}
if ($order == 0) {
    die "usage: $0 -order <order>\n";
}

#  This scales at 4^k * 4^k -- (state vectors) * (cycle length)
#
#    time seqrequester shift -search -fast -tapmin 300000000000 -tapmax 300033333333 -order 12 -report 0.9
#
#  The SV MUST NOT start with 0 or 1.  We generate the new SR state as
#
#    out = rightmost 2 bits of SR
#    mul = SV * out (* in GF4, operating on 2-bit tuples)
#    SR  = (SR >> 2) ^ mul
#
#  When the SV starts with 0, the resulting SR will always start with zero -
#  it is like we simply reduce the order by one.  The longest cycle in order
#  12 here is 4194304, exactly 1/4 the maximum expected, and exactly the
#  length of the cycle for order 11.
#
#  When the SV starts with 1, we car simply moving the output symbol to the
#  input - we turn the shift register into a cyclic shift register.  The longest
#  cycle in order 12 here is 'exactly' 1/3 the maximum - 5592406 out of 16777216.
#  I do not know why.
#
#  3's and 2's seem interchangable; swapping all 3's and 2's in a vector
#  seems to generate the same cycles length.  This probably just 'reverses'
#  the sequence.  NOT RIGOROUSLY TESTED.
#
#  
#
#    valgrind --tool=callgrind --dump-instr=yes --simulate-cache=yes --collect-jumps=yes \
#      seqrequester shift -search -fast -tapmin 300000000000 -tapmax 300003333333 -order 12 -report 0.9
#    callgrind_control -h
#
#  Hardcoded job splits.
#    order <= 10 -  0 digits -       1 job
#    order == 11 -  1 digit  -       4 jobs
#    order == 12 -  2 digits -      16 jobs
#    order == 13 -  4 digits -     256 jobs
#    order == 14 -  6 digits -    4096 jobs
#    order == 15 -  8 digits -   65536 jobs
#    order == 16 - 11 digits - 4194304 jobs
#
#  Times are on Ryzen 7 3700X unless noted.
#  Times are signifcantly worse if multiple jobs run.
#
#    order == 12 - with prefix  5 took   1.43 minutes for 30000........  (actual)     -     1024 prefixes ->      24 CPU hours
#    order == 12 - with prefix  5 took   0.93 minutes for 30000........  (actual)     -     1024 prefixes -> no detect version
#    order == 12 - with prefix  4 took   5.75 minutes for 3000.........  (actual)     -      256 prefixes ->      24 CPU hours
#    order == 12 - with prefix  4 took  14.03 minutes for 3000.........  (actual)     -      256 prefixes ->      60 CPU hours (on d)
#    order == 12 - with prefix  4 took  14.03 minutes for 3000.........  (actual)     -      256 prefixes -> but 90 minutes/job is all 24 cores are used!
#    order == 12 - with prefix  3 took        minutes for 300..........  (estimated)  -       64 prefixes -> 
#    order == 12 - with prefix  2 took        minutes for 30...........  (estimated)  -       16 prefixes -> 
#
#    order == 13 - with prefix  6 took   9    minutes for 300000.......  (actual)     -     4096 prefixes ->     614 CPU hours
#    order == 13 - with prefix  5 took  36    minutes for 30000........  (estimated)  -     1024 prefixes -> 
#    order == 13 - with prefix  4 took 144    minutes for 3000.........  (estimated)  -      256 prefixes -> 
#
#    order == 14 - with prefix  9 took   8    minutes for 300000000.....  (actual)    -   262144 prefixes -> 
#    order == 14 - with prefix  8 took  32    minutes for 30000000......  (actual)    -    65536 prefixes ->  35,000 CPU hours
#    order == 14 - with prefix  7 took 128    minutes for 3000000.......  (estimated) -    16384 prefixes -> 
#
#    order == 15 - with prefix 12 took   2.96 minutes for 300000000000... (actual)    - 16777216 prefixes -> 827,675 CPU hours
#    order == 15 - with prefix 12 took   0.28 minutes for 300000000000... (actual)    - 16777216 prefixes -> no detect version
#    order == 15 - with prefix 11 took   9.13 minutes for 30000000000.... (actual)    -  4194304 prefixes -> 638,233 CPU hours
#    order == 15 - with prefix 10 took  39.38 minutes for 3000000000..... (actual)    -  1048576 prefixes -> 688,273 CPU hours (78.6 years)
#    order == 15 - with prefix 10 took   3.76 minutes for 3000000000..... (actual)    -  1048576 prefixes -> no detect version -- 
#    order == 15 - with prefix 10 took   4.11 minutes for 3000000000..... (actual)    -  1048576 prefixes -> no detect version -- with 16 copies running
#    order == 15 - with prefix  9 took        minutes for 300000000...... (estimated) -   262144 prefixes -> 
#
#  CL time
#
#  BD time
#    time seqrequester shift -search -fast -order 21 -report 0.0 -tapmin 300000000000000000000 -tapmax 300000000000000003333 -weight 3
#    66069.844u 0.283s 18:21:14.56 99.9%     177+817k 0+0io 7pf+0w
#
#  order =    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
my @plen = (  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  3,  4,  5,  9, 11, 13, 15, 17, 19,  0 );

#  prefix 5 for 15 is way too big.  this resulted in 256 jobs (1024 prefixes,
#  but only the ones starting with 3 are computed).  each job is taking a bit
#  more than 100 hours.
#    wallclock 415000s (115 hours x 256 jobs = 29440 cpu hours)
#    user      381110s (remarkably constant)
#
#  job scaling should be x4 for each search, and there are x4 more vectors to try.
#  so for constant job time, we need to increase the prefix length by 2 for each.

#  prefix 9 for 16 resulted in the last job 12550u 3:30 on bv (x 65536 jobs = 229376 cpu hours).


print STDERR "Building job names with $plen[$order] components.\n";

my @alljobs = ( "0", "1", "2", "3" );
my @subjobs;

for (my $p=1; $p < $plen[$order]; $p++) {
    my @tmpjobs = @alljobs;
    undef @alljobs;

    foreach my $e (@tmpjobs) {
        push @alljobs, "0$e";
        push @alljobs, "1$e";
        push @alljobs, "2$e";
        push @alljobs, "3$e";
    }
}

print STDERR "Filtering job names.\n";

foreach my $pp (@alljobs) {
    push @subjobs, $pp   if (($pp =~ m/^0/) && ($order < 14));   #  Ignore 0, no complete cycles in here.
    push @subjobs, $pp   if (($pp =~ m/^1/) && ($order < 14));   #  Ignore 1, no complete cycles in here.
    push @subjobs, $pp   if (($pp =~ m/^2/) && ($order < 14));   #  Ignore 2, seems to be just a dual of 3.
    push @subjobs, $pp   if (($pp =~ m/^3/) && ($order < 22));   #  Always use 3.
}

@subjobs = sort { $b cmp $a } @subjobs;

my $tJobs = scalar(@alljobs);
my $sJobs = scalar(@subjobs);

print STDERR "Created    $tJobs jobs.\n";
print STDERR "Submitting $sJobs jobs.\n";


my $bgn = "0" x ($order - $plen[$order]);
my $end = "3" x ($order - $plen[$order]);

open(F, "> o${order}w${weight}.dat");

foreach my $pp (@subjobs) {
    print F "$pp\n";
}

close(F);



open(F, "> o${order}w${weight}.sh");
print F "#!/bin/sh\n";
print F "\n";
print F "\n";
print F "if [ x\$SGE_TASK_ID = x -o x\$SGE_TASK_ID = xundefined -o x\$SGE_TASK_ID = x0 ]; then\n";
print F "  j=\$1\n";
print F "else\n";
print F "  j=\$SGE_TASK_ID\n";
print F "fi\n";
print F "if [ x\$j = x ]; then\n";
print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "\n";
print F "p=`head -n \$j o${order}w${weight}.dat | tail -n 1`";
print F "\n";
print F "\n";
print F "if [ -e \"o${order}w${weight}-out/o${order}w${weight}p\$p.out\" ]; then\n";
print F "  echo \"o${order}w${weight}p\$p.out exists, it's done!\"\n";
print F "  exit 0\n";
print F "fi\n";
print F "\n";
print F "if [ ! -d \"o${order}w${weight}-out\" ]; then\n";
print F "  mkdir -p o${order}w${weight}-out\n";
print F "fi\n";
print F "\n";
print F "/work/seqrequester/FreeBSD-amd64/bin/seqrequester \\\n";
print F "  shift -search -fast \\\n";
print F "    -order  $order \\\n";
print F "    -weight $weight \\\n";
print F "    -report 0.0 \\\n";
print F "    -tapmin \${p}$bgn \\\n";
print F "    -tapmax \${p}$end \\\n";
print F "> o${order}w${weight}-out/o${order}w${weight}p\$p.out.WORKING \\\n";
print F "&& \\\n";
print F "mv o${order}w${weight}-out/o${order}w${weight}p\$p.out.WORKING o${order}w${weight}-out/o${order}w${weight}p\$p.out\n";
print F "\n";
print F "\n";
print F "#\n";
print F "#  qsub -h -cwd -j y -o o${order}w${weight}-err/o${order}w${weight}-\\\$TASK_ID.err -l memory=1g -t 1-$sJobs -N o${order}w${weight} ./o${order}w${weight}.sh\n";
print F "#\n";
print F "\n";
print F "exit 0\n";
print F "\n";

close(F);

chmod(0755, "o${order}w${weight}.sh");

system("mkdir o${order}w${weight}-err");
system("mkdir o${order}w${weight}-out");

print  "qsub -h -cwd -j y -o o${order}w${weight}-err/o${order}w${weight}-\\\$TASK_ID.err -l memory=1g -t 1-$sJobs -N o${order}w${weight} ./o${order}w${weight}.sh\n";
system("qsub -h -cwd -j y -o o${order}w${weight}-err/o${order}w${weight}-\\\$TASK_ID.err -l memory=1g -t 1-$sJobs -N o${order}w${weight} ./o${order}w${weight}.sh")    if (1);

exit;
