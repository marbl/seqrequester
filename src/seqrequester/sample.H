
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
    isPaired        = false;

    numCopies       = 1;

    desiredCoverage = 0.0;
    genomeSize      = 0;

    desiredNumReads = 0;
    desiredNumBases = 0;

    desiredFraction = 0.0;

    randomSeedValid = false;
    randomSeed      = 0;

    outputFASTA     = false;
    outputFASTQ     = false;
    outputQV        = 20;

    memset(output1, 0, FILENAME_MAX+1);
    memset(output2, 0, FILENAME_MAX+1);
  }

  ~sampleParameters() {
  }


  void      initialize(void) {
  };


  void      finalize(void) {
  }


  bool    isPaired;

  uint32  numCopies;

  double  desiredCoverage;
  uint64  genomeSize;

  uint64  desiredNumReads;
  uint64  desiredNumBases;

  double  desiredFraction;

  bool    randomSeedValid;
  uint32  randomSeed;

  bool    outputFASTA;
  bool    outputFASTQ;
  uint8   outputQV;

  char    output1[FILENAME_MAX+1];
  char    output2[FILENAME_MAX+1];
};
