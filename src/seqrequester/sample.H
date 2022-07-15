
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
#error  Include sample.H only from seqrequester.H
#endif


class sampleParameters {
public:
  sampleParameters() {
    //output1 = new char [FILENAME_MAX+1];
    //output2 = new char [FILENAME_MAX+1];
  }
  ~sampleParameters() {}

  bool      parseOption(opMode &mode, int32 &arg, int32 argc, char **argv);
  void      showUsage(opMode mode);

  bool      checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors);

public:
  bool      isPaired = false;

  uint32    numCopies = 1;

  double    desiredCoverage = 0.0;
  uint64    genomeSize = 0;

  uint64    desiredNumReads = 0;
  uint64    desiredNumBases = 0;

  double    desiredFraction = 0.0;

  mtRandom  mt;

  bool      outputFASTA = false;
  bool      outputFASTQ = false;
  uint8     outputQV = 20;

  char const *output1 = nullptr;//[FILENAME_MAX+1] = {0};
  char const *output2 = nullptr;//[FILENAME_MAX+1] = {0};
};
