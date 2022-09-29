
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

#include "seqrequester.H"

bool
generateParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "generate") == 0) {
    mode = modeGenerate;
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-ident") == 0)) {
    ident = argv[++arg];
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-min") == 0)) {
    minLength = strtouint64(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-max") == 0)) {
    maxLength = strtouint64(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-sequences") == 0)) {
    nSeqs = strtouint64(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-bases") == 0)) {
    nBases = strtouint64(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-seed") == 0)) {
    mt.mtSetSeed(strtouint32(argv[++arg]));
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-gaussian") == 0)) {
    useGaussian = true;
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-mirror") == 0)) {
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-gc") == 0)) {
    double  gc = strtodouble(argv[++arg]);
    double  at = 1.0 - gc;

    gFreq = cFreq = gc / 2.0;
    aFreq = tFreq = at / 2.0;
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-at") == 0)) {
    double  at = strtodouble(argv[++arg]);
    double  gc = 1.0 - at;

    gFreq = cFreq = gc / 2.0;
    aFreq = tFreq = at / 2.0;
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-a") == 0) ){
    aFreq = strtodouble(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-c") == 0)) {
    cFreq = strtodouble(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-g") == 0)) {
    gFreq = strtodouble(argv[++arg]);
  }

  else if ((mode == modeGenerate) && (strcmp(argv[arg], "-t") == 0)) {
    tFreq = strtodouble(argv[++arg]);
  }

  else {
    return(false);
  }

  return(true);
}



void
generateParameters::showUsage(opMode mode) {

  if (mode != modeGenerate)
    return;

  fprintf(stderr, "OPTIONS for generate mode:\n");
  fprintf(stderr, "  -ident S       name reads using printf() pattern in S; default 'random%%08lu'\n");
  fprintf(stderr, "  -min M         minimum sequence length\n");
  fprintf(stderr, "  -max M         maximum sequence length\n");
  fprintf(stderr, "  -sequences N   generate N sequences\n");
  fprintf(stderr, "  -bases B       generate at least B bases, no more than B+maxLength-1 bases.\n");
  fprintf(stderr, "  -gaussian      99.73%% of the reads (3 standard deviations) will be between min and max\n");
  fprintf(stderr, "  -mirror F      (not implemented; mirror the length distrib of F)\n");
  fprintf(stderr, "  -gc bias       sets GC/AT composition (default 0.50)\n");
  fprintf(stderr, "  -at bias       sets GC/AT composition (default 0.50)\n");
  fprintf(stderr, "  -a freq        sets frequency of A bases (default 0.25)\n");
  fprintf(stderr, "  -c freq        sets frequency of C bases (default 0.25)\n");
  fprintf(stderr, "  -g freq        sets frequency of G bases (default 0.25)\n");
  fprintf(stderr, "  -t freq        sets frequency of T bases (default 0.25)\n");
  fprintf(stderr, "  -seed s        seed the pseudo-random number generate with integer 's'\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "The -gc option is a shortcut for setting all four base frequencies at once.  Order matters!\n");
  fprintf(stderr, "  -gc 0.6 -a 0.1 -t 0.3 -- sets G = C = 0.3, A = 0.1, T = 0.3\n");
  fprintf(stderr, "  -a 0.1 -t 0.3 -gc 0.6 -- sets G = C = 0.3, A = T = 0.15\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Base frequencies are scaled to sum to 1.0.\n");
  fprintf(stderr, "  -a 1.25 -- results in a sum of 2.0 (1.25 + 0.25 + 0.25 + 0.25) so final frequencies will be:\n");
  fprintf(stderr, "             A =         1.25/2 = 0.625\n");
  fprintf(stderr, "             C = G = T = 0.25/2 = 0.125.\n");
  fprintf(stderr, "  -gc 0.8 -a 1.0 -t 0.2 -- sum is also 2.0, final frequencies will be:\n");
  fprintf(stderr, "             A =         1.00/2 = 0.5\n");
  fprintf(stderr, "             C = G =     0.40/2 = 0.2\n");
  fprintf(stderr, "             T =         0.20/2 = 0.1\n");
  fprintf(stderr, "\n");
}



bool
generateParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeGenerate)
    return(false);

  if (inputs.size() > 0)
    sprintf(errors, "ERROR:  generate mode does not work with input sequence files.\n");

  //  Check for invalid.  If not set up, just return.

  if ((nSeqs == 0) && (nBases == 0))
    sprintf(errors, "ERROR:  at least one of -sequences and -bases must be supplied.\n");

  if (minLength > maxLength)
    sprintf(errors, "ERROR:  Told to generate sequences with min length larger than max length.\n");

  //  Unlimit any unset limit.

  if (nSeqs  == 0)      nSeqs = UINT64_MAX;
  if (nBases == 0)     nBases = UINT64_MAX;

  //  Set Gaussian mean and standard deviation

  gMean   = (minLength + maxLength) / 2.0;
  gStdDev = (maxLength - minLength) / 6.0;

  //  Load lengths to mirror.

  return(errors.size() > 0);
}



void
doGenerate(generateParameters &genPar) {
  uint64  nSeqs  = 0;
  uint64  nBases = 0;

  uint64  seqLen = 0;
  uint64  seqMax = 65536;
  char   *seq    = new char  [seqMax + 1];
  uint8  *qlt    = new uint8 [seqMax + 1];

  double Athresh = genPar.aFreq;
  double Cthresh = genPar.aFreq + genPar.cFreq;
  double Gthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq;
  double Tthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq + genPar.tFreq;

  fprintf(stderr, "Generating up to " F_U64 " sequences and up to " F_U64 " bases.\n", nSeqs, nBases);

  double  fSum = genPar.aFreq + genPar.cFreq + genPar.gFreq + genPar.tFreq;

  fprintf(stderr, "Using base frequencies:\n");
  fprintf(stderr, "  A = %7.5f / %7.5f = %7.5f\n", genPar.aFreq, fSum, genPar.aFreq / fSum);
  fprintf(stderr, "  C = %7.5f / %7.5f = %7.5f\n", genPar.cFreq, fSum, genPar.cFreq / fSum);
  fprintf(stderr, "  G = %7.5f / %7.5f = %7.5f\n", genPar.gFreq, fSum, genPar.gFreq / fSum);
  fprintf(stderr, "  T = %7.5f / %7.5f = %7.5f\n", genPar.tFreq, fSum, genPar.tFreq / fSum);

  genPar.aFreq /= fSum;
  genPar.cFreq /= fSum;
  genPar.gFreq /= fSum;
  genPar.tFreq /= fSum;


  if (genPar.ident == nullptr)
    genPar.ident = merylutil::duplicateString("random%08lu");

  while ((nSeqs  < genPar.nSeqs) &&
         (nBases < genPar.nBases)) {
    double   len = genPar.mt.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    while (len < -0.5)
      len = genPar.mt.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    seqLen = (uint64)round(len);

    if (seqLen+1 >= seqMax)
      merylutil::resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1);

    for (uint64 ii=0; ii<seqLen; ii++) {
      double  bp = genPar.mt.mtRandomRealOpen();

      if        (bp < Athresh) {
        seq[ii] = 'A';
        qlt[ii] = 0;

      } else if (bp < Cthresh) {
        seq[ii] = 'C';
        qlt[ii] = 0;

      } else if (bp < Gthresh) {
        seq[ii] = 'G';
        qlt[ii] = 0;

      } else {
        seq[ii] = 'T';
        qlt[ii] = 0;
      }
    }

    seq[seqLen] = 0;
    qlt[seqLen] = 0;

    merylutil::outputFASTA(stdout, seq, seqLen, 0, genPar.ident, nSeqs);

    nSeqs  += 1;
    nBases += seqLen;
  }

  delete [] seq;
  delete [] qlt;
}

