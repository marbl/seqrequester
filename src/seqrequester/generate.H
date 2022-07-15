
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
#error  Include generate.H only from seqrequester.H
#endif


class generateParameters {
public:
  generateParameters() {}
  ~generateParameters() {}

  bool      parseOption(opMode &mode, int32 &arg, int32 argc, char **argv);
  void      showUsage(opMode mode);

  bool      checkOptions(opMode mode, std::vector<char const *> &inputs, std::vector<char const *> &errors);

public:
  char const  *ident = nullptr;

  uint64       minLength = 0;
  uint64       maxLength = 10000;

  uint64       nSeqs  = 0;
  uint64       nBases = 0;

  bool         useGaussian = true;
  double       gMean       = 0;
  double       gStdDev     = 0;

  bool         useMirror = false;
  char        *mirrorInput = nullptr;
  double       mirrorDistribution = 0.0;
  uint64       mirrorDistributionLen = 0;
  uint64       mirrorDistributionMax = 0;
  uint64       mirrorDistributionSum = 0;

  mtRandom     mt;

  double       aFreq = 0.25;
  double       cFreq = 0.25;
  double       gFreq = 0.25;
  double       tFreq = 0.25;
};
