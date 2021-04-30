
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
  generateParameters() {
    ident                 = nullptr;

    minLength             = 0;
    maxLength             = 10000;

    nSeqs                 = 0;
    nBases                = 0;

    randomSeedValid       = false;
    randomSeed            = 0;

    useGaussian           = true;
    gMean                 = 0;
    gStdDev               = 0;

    useMirror             = false;
    mirrorInput           = NULL;
    mirrorDistribution    = 0.0;
    mirrorDistributionLen = 0;
    mirrorDistributionMax = 0;
    mirrorDistributionSum = 0;

    aFreq                 = 0.25;
    cFreq                 = 0.25;
    gFreq                 = 0.25;
    tFreq                 = 0.25;
  };

  ~generateParameters() {
  };


  void      finalize(void) {

    //  Check for invalid.  If not set up, just return.

    if ((nSeqs == 0) && (nBases == 0))
      return;

    if (minLength > maxLength)
      fprintf(stderr, "ERROR:  Told to generate sequences with min length larger than max length.\n"), exit(1);

    //  Unlimit any unset limit.

    fprintf(stderr, "Generating up to " F_U64 " sequences and up to " F_U64 " bases.\n", nSeqs, nBases);

    if (nSeqs  == 0)      nSeqs = UINT64_MAX;
    if (nBases == 0)     nBases = UINT64_MAX;

    //  Set Gaussian mean and standard deviation

    gMean   = (minLength + maxLength) / 2.0;
    gStdDev = (maxLength - minLength) / 6.0;

    //  Load lengths to mirror.

    //  Set base frequencies.

    double  fSum = aFreq + cFreq + gFreq + tFreq;

    fprintf(stderr, "Using base frequencies:\n");
    fprintf(stderr, "  A = %7.5f / %7.5f = %7.5f\n", aFreq, fSum, aFreq / fSum);
    fprintf(stderr, "  C = %7.5f / %7.5f = %7.5f\n", cFreq, fSum, cFreq / fSum);
    fprintf(stderr, "  G = %7.5f / %7.5f = %7.5f\n", gFreq, fSum, gFreq / fSum);
    fprintf(stderr, "  T = %7.5f / %7.5f = %7.5f\n", tFreq, fSum, tFreq / fSum);

    aFreq /= fSum;
    cFreq /= fSum;
    gFreq /= fSum;
    tFreq /= fSum;
  };

  char const  *ident;

  uint64       minLength;
  uint64       maxLength;

  uint64       nSeqs;
  uint64       nBases;

  bool         randomSeedValid;
  uint32       randomSeed;

  bool         useGaussian;
  double       gMean;
  double       gStdDev;

  bool         useExponential;

  bool         useMirror;
  char        *mirrorInput;
  double       mirrorDistribution;
  uint64       mirrorDistributionLen;
  uint64       mirrorDistributionMax;
  uint64       mirrorDistributionSum;

  double       aFreq;
  double       cFreq;
  double       gFreq;
  double       tFreq;
};
