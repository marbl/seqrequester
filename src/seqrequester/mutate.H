
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
#error  Include mutate.H only from seqrequester.H
#endif


class mutateParameters {
public:
  mutateParameters() {
    for (uint32 ii=0; ii<256; ii++) {
      for (uint32 jj=0; jj<256; jj++)
        pS[ii][jj] = 0.0;
      pI[ii] = 0.0;
      pD[ii] = 0.0;

      pSubstitute[ii] = 0.0;
    }

    pInsert = 0.0;
  };

  ~mutateParameters() {
  };

public:
  void      setProbabilitySubstititue(double p, char baseToReplace, char replacementBase) {
    fprintf(stderr, "sub %c -> %c %f\n", baseToReplace, replacementBase, p);
    pS[baseToReplace][replacementBase] = p;
  };

  void      setProbabilityInsert(double p, char baseToInsert) {
    fprintf(stderr, "ins %c %f\n", baseToInsert, p);
    pI[baseToInsert] = p;
  };

  void      setProbabilityDelete(double p, char baseToDelete) {
    fprintf(stderr, "del %c %f\n", baseToDelete, p);
    pD[baseToDelete] = p;
  };

  void      finalize(void) {

    for (uint32 ii=0; ii<256; ii++)
      pSubstitute[ii] = 0.0;

    pInsert = 0.0;

    for (uint32 ii=0; ii<256; ii++) {
      for (uint32 jj=0; jj<256; jj++)
        pSubstitute[ii] += pS[ii][jj];

      pInsert           += pI[ii];
      pDelete           += pD[ii];
    }
  };


  //  probability of substituting base [a] with [b]    (base based)
  //  probability of inserting    an A, C, G, T, N     (space based)
  //  probability of deleting     base [a]             (base based)
  //
  //  extensions:
  //  probability of inserting ACGTN between bases a and b.
  //
  double  pS[256][256];
  double  pI[256];
  double  pD[256];

  double  pSubstitute[256];   //  Probability of substituting 'a' for anything.
  double  pInsert;            //  Probability of inserting any base.
  double  pDelete;
};
