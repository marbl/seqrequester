
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


class summarizeParameters {
public:
  summarizeParameters() {
    genomeSize   = 0;

    limitTo1x    = false;
    breakAtN     = false;

    asSimple     = false;
    asLength     = false;

    asSequences  = true;
    asBases      = false;
  };

  ~summarizeParameters() {
  };


  void      finalize(void) {
  }


  uint64    genomeSize;

  bool      limitTo1x;
  bool      breakAtN;

  bool      asSimple;
  bool      asLength;

  bool      asSequences;
  bool      asBases;
};
