
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
#error  Include shiftregister.H only from seqrequester.H
#endif


class shiftRegisterParameters {
public:
  shiftRegisterParameters() {
    search = false;
    fast   = false;

    report = 1.0;

    order  = 0;
    length = 0;

    weight = 0;

    memset(sr,    0, sizeof(char) * 65);
    memset(svmin, 0, sizeof(char) * 65);
    memset(svmax, 0, sizeof(char) * 65);
  };
  ~shiftRegisterParameters() {
  };

  void    initialize(void);

  uint64   getEncodedSR(void);
  uint64   getCycleLen(void);
  uint64   getCycleMax(void);
  uint64   getEncodedSVmin(void);
  uint64   getEncodedSVmax(void);
  uint64   getEncodedSVmask(void);

  char     numberToBase(uint32 number) {
    number = (number << 1) + 'A';

    if (number == 'E')
      number = 'T';

    return(number);
  };

  bool    search;
  bool    fast;

  double  report;

  uint32  order;
  uint64  length;

  uint32  weight;

  char    sr[65];
  char    svmin[65];
  char    svmax[65];
};
