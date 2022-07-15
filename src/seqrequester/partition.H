
/******************************************************************************
 *
 *  This file is part of seqrequester, a tool for summarizing, extracting,
 *  generating and modifying DNA sequences.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef SEQREQUESTER_INCLUDE_H
#error  Include partition.H only from seqrequester.H
#endif


class partitionParameters {
public:
  partitionParameters() {}
  ~partitionParameters() {}

  bool      parseOption(opMode &mode, int32 &arg, int32 argc, char **argv);
  void      showUsage(opMode mode);

  bool      checkOptions(opMode mode, std::vector<char const *> &inputs, std::vector<char const *> &errors);

  bool    isPaired         = false;

  double  desiredCoverage  = 0.0;
  uint64  genomeSize       = 0;

  uint64  desiredNumReads  = uint64max;
  uint64  desiredNumBases  = uint64max;

  double  desiredFraction  = 0.0;

  bool    randomSeedValid  = false;
  uint32  randomSeed       = 0;

  bool    outputFASTA      = false;
  bool    outputFASTQ      = false;
  uint8   outputQV         = 20;

  uint32  numWriters       = 4;

  char    output1[FILENAME_MAX+1] = {0};
  char    output2[FILENAME_MAX+1] = {0};
};

