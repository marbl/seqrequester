
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
#error  Include simulate.H only from seqrequester.H
#endif


class simulateParameters {
public:
  simulateParameters() {}
  ~simulateParameters() {}

  bool      parseOption(opMode &mode, int32 &arg, int32 argc, char **argv);
  void      showUsage(opMode mode);

  bool      checkOptions(opMode mode, std::vector<char const *> &inputs, std::vector<char const *> &errors, char const *progName);

public:
  uint64  genomeSize         = 0;
  bool    circular           = false;
  bool    truncate           = false;

  bool    test               = false;

  double  desiredCoverage    = 0.0;
  uint64  desiredNumReads    = 0;
  uint64  desiredNumBases    = 0;

  uint32  desiredMinLength   = uint32max;
  uint32  desiredMaxLength   = 0;

  double  rcProb             = 0.5;

  mtRandom             mt;
  sampledDistribution  dist;

  char const    *genomeName  = nullptr;
  char const    *distribName = nullptr;
};
