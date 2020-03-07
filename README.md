# seqrequester

This is 'seqrequester', a tool for summarizing, extracting, generating and
modifying DNA sequences.

# Summarizing

The summarize mode will generate a table of Nx lengths, a lovely ASCII
plot of the histogram of sequence lengths, report GC content, and di-
and tri-nucleotide frequencies.

It can optionally split sequences at N's before computing the length of a sequence.

You can also get a simple histogram of the sequence lengths and the number of sequences
at each length, or just a simple list of all sequence lengths.

It will, of course, read FASTA and FASTQ, uncompressed or compressed with
gzip, bzip2 or xz.

Only one report is generated, regardless of how many sequence files are supplied.


```
% seqrequester summarize
usage: seqrequester [mode] [options] [sequence_file ...]

OPTIONS for summarize mode:
  -size          base size to use for N50 statistics
  -1x            limit NG table to 1x coverage

  -split-n       split sequences at N bases before computing length
  -simple        output a simple 'length numSequences' histogram
  -lengths       output a list of the sequence lengths

  -assequences   load data as complete sequences (for testing)
  -asbases       load data as blocks of bases    (for testing)
```

```
% seqrequester summarize /archive/mothra/FLX/*gz

G=6462464889                       sum of  ||               length     num
NG         length     index       lengths  ||                range    seqs
----- ------------ --------- ------------  ||  ------------------- -------
00010          652    801160    646246790  ||         42-112          4768|-
00020          582   1862887   1292493013  ||        113-183         16961|-
00030          555   3002684   1938739802  ||        184-254         89381|--
00040          538   4186751   2584986254  ||        255-325        536862|--------
00050          523   5405461   3231232945  ||        326-396       1463599|--------------------
00060          509   6657839   3877479295  ||        397-467       1960924|---------------------------
00070          488   7952426   4523725460  ||        468-538       4616863|---------------------------------------------------------------
00080          447   9329218   5169971940  ||        539-609       2858982|----------------------------------------
00090          389  10872299   5816218777  ||        610-680        625376|---------
00100           42  12803136   6462464889  ||        681-751        252454|----
001.000x            12803137   6462464889  ||        752-822        134849|--
                                           ||        823-893         78435|--
                                           ||        894-964         47976|-
                                           ||        965-1035        30852|-
                                           ||       1036-1106        21127|-
                                           ||       1107-1177        14817|-
                                           ||       1178-1248        28461|-
                                           ||       1249-1319         4930|-
                                           ||       1320-1390         3655|-
                                           ||       1391-1461         2657|-
                                           ||       1462-1532         2120|-
                                           ||       1533-1603         1597|-
                                           ||       1604-1674         1268|-
                                           ||       1675-1745          953|-
                                           ||       1746-1816          766|-
                                           ||       1817-1887          573|-
                                           ||       1888-1958          443|-
                                           ||       1959-2029          344|-
                                           ||       2030-2100         1022|-
                                           ||       2101-2171           21|-
                                           ||       2172-2242           23|-
                                           ||       2243-2313           20|-
                                           ||       2314-2384           17|-
                                           ||       2385-2455            8|-
                                           ||       2456-2526            9|-
                                           ||       2527-2597            4|-
                                           ||       2598-2668            2|-
                                           ||       2669-2739            6|-
                                           ||       2740-2810            2|-
                                           ||       2811-2881            6|-
                                           ||       2882-2952            1|-
                                           ||       2953-3023            0|
                                           ||       3024-3094            0|
                                           ||       3095-3165            1|-
                                           ||       3166-3236            0|
                                           ||       3237-3307            0|
                                           ||       3308-3378            0|
                                           ||       3379-3449            1|-
                                           ||       3450-3520            0|
                                           ||       3521-3591            1|-

--------------------- --------------------- ----------------------------------------------------------------------------------------------
       mononucleotide          dinucleotide                                                                                  trinucleotide
--------------------- --------------------- ----------------------------------------------------------------------------------------------
  1959571306 0.3032 A   665030151 0.1031 AA   237235545 0.0369 AAA    132268487 0.0205 AAC    136675399 0.0212 AAG    158473516 0.0246 AAT
  1247489432 0.1930 C   389352138 0.0604 AC   115665542 0.0180 ACA     87346626 0.0136 ACC     70986769 0.0110 ACG    114582435 0.0178 ACT
  1345011807 0.2081 G   397219280 0.0616 AG   121659180 0.0189 AGA     65811037 0.0102 AGC    102037062 0.0159 AGG    106854671 0.0166 AGT
  1910392344 0.2956 T   507072196 0.0786 AT   152454159 0.0237 ATA     89877335 0.0140 ATC    106195089 0.0165 ATG    158544503 0.0246 ATT
                        380831936 0.0590 CA   132169383 0.0205 CAA     76839888 0.0119 CAC     67197045 0.0104 CAG    104566859 0.0162 CAT
      --GC--  --AT--    281892951 0.0437 CC    86178881 0.0134 CCA     65022575 0.0101 CCC     50576089 0.0079 CCG     79660170 0.0124 CCT
      40.12%  59.88%    208535008 0.0323 CG    60164341 0.0093 CGA     27649662 0.0043 CGC     52322022 0.0081 CGG     67296554 0.0105 CGT
                        374626420 0.0581 CT    95122699 0.0148 CTA     75643338 0.0118 CTC     74304266 0.0115 CTG    129554475 0.0201 CTT
                        383528854 0.0595 GA   128291282 0.0199 GAA     70244915 0.0109 GAC     88104990 0.0137 GAG     96746854 0.0150 GAT
                        218253748 0.0338 GC    72696062 0.0113 GCA     51118632 0.0079 GCC     27797659 0.0043 GCG     66512705 0.0103 GCT
                        361154273 0.0560 GG    87662449 0.0136 GGA     51591104 0.0080 GGC    124820908 0.0194 GGG     89162545 0.0139 GGT
                        371793122 0.0576 GT   112773785 0.0175 GTA     61960408 0.0096 GTC     66540524 0.0103 GTG    130508651 0.0203 GTT
                        526153182 0.0816 TA   165978181 0.0258 TAA    109335141 0.0170 TAC    104450634 0.0162 TAG    146068463 0.0227 TAT
                        355583693 0.0551 TC   105474843 0.0164 TCA     77928750 0.0121 TCC     58773061 0.0091 TCG    113158614 0.0176 TCT
                        375551419 0.0582 TG   113251416 0.0176 TGA     72683975 0.0113 TGC     81429936 0.0127 TGG    107781308 0.0167 TGT
                        653083381 0.1013 TT   164899014 0.0256 TTA    127390884 0.0198 TTC    127703170 0.0198 TTG    233082150 0.0362 TTT
```

# Extracting

```
% seqrequester extract
usage: ./FreeBSD-amd64/bin/seqrequester [mode] [options] [sequence_file ...]

OPTIONS for extract mode:
  -bases     baselist extract bases as specified in the 'list' from each sequence
  -sequences seqlist  extract ordinal sequences as specified in the 'list'

  -reverse            reverse the bases in the sequence
  -complement         complement the bases in the sequence
  -rc                 alias for -reverse -complement

  -compress           compress homopolymer runs to one base

  -upcase
  -downcase

  -length min-max     print sequence if it is at least 'min' bases and at most 'max' bases long
  
                      a 'baselist' is a set of integers formed from any combination
                      of the following, seperated by a comma:
                           num       a single number
                           bgn-end   a range of numbers:  bgn <= end
                      bases are spaced-based; -bases 0-2,4 will print the bases between
                      the first two spaces (the first two bases) and the base after the
                      fourth space (the fifth base).
  
                      a 'seqlist' is a set of integers formed from any combination
                      of the following, seperated by a comma:
                           num       a single number
                           bgn-end   a range of numbers:  bgn <= end
                      sequences are 1-based; -sequences 1,3-5 will print the first, third,
                      fourth and fifth sequences.
```

# Sampling

```
% seqrequester sample
usage: ./FreeBSD-amd64/bin/seqrequester [mode] [options] [sequence_file ...]

OPTIONS for sample mode:
  -paired             treat inputs as paired sequences; the first two files form the
                      first pair, and so on.

  -copies C           write C different copies of the sampling (without replacement).
  -output O           write output sequences to file O.  If paired, two files must be supplied.

  -coverage C         output C coverage of sequences, based on genome size G.
  -genomesize G       

  -bases B            output B bases.

  -reads R            output R reads.
  -pairs P            output P pairs (only if -paired).

  -fraction F         output fraction F of the input bases.

```

# Generating

Undocumented.

# Simulating

```
seqrequester simulate
usage: ./FreeBSD-amd64/bin/seqrequester [mode] [options] [sequence_file ...]

OPTIONS for simulate mode:
  -genome G           sample reads from these sequences
  -circular           treat the sequences in G as circular

  -genomesize g       genome size to use for deciding coverage below
  -coverage c         generate approximately c coverage of output
  -nreads n           generate exactly n reads of output
  -nbases n           generate approximately n bases of output

  -distribution F     generate read length by sampling the distribution in file F
                        one column  - each line is the length of a sequence
                        two columns - each line has the 'length' and 'number of sequences'

                      if file F doesn't exist, use a built-in distribution
                        ultra-long-nanopore
                        pacbio
                        pacbio-hifi

  -length min[-max]   (not implemented)
  -output x.fasta     (not implemented)
```
