
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

#include "sampledDistribution.H"

class simulateParameters {
public:
  simulateParameters() {
    genomeSize      = 0;
    circular        = false;

    desiredCoverage = 0.0;
    desiredNumReads = 0;
    desiredNumBases = 0;

    desiredMinLength = 0;
    desiredMaxLength = UINT32_MAX;

    memset(genomeName,  0, FILENAME_MAX+1);
    memset(distribName, 0, FILENAME_MAX+1);
    memset(outputName,  0, FILENAME_MAX+1);
  }

  ~simulateParameters() {
  }


  void      initialize(void) {
  };


  void      finalize(void);


  uint64  genomeSize;
  bool    circular;

  double  desiredCoverage;
  uint64  desiredNumReads;
  uint64  desiredNumBases;

  uint32  desiredMinLength;
  uint32  desiredMaxLength;

  sampledDistribution  dist;

  char    genomeName[FILENAME_MAX+1];
  char    distribName[FILENAME_MAX+1];
  char    outputName[FILENAME_MAX+1];
};
