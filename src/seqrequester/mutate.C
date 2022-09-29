
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
mutateParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "mutate") == 0) {
    mode = modeMutate;
  }

  else if ((mode == modeMutate) && (strcmp(argv[arg], "-s") == 0)) {
    double  p = strtodouble(argv[arg+1]);
    char    a = argv[arg+2][0];
    char    b = argv[arg+3][0];

    arg += 3;

    setProbabilitySubstititue(p, a, b);
  }

  else if ((mode == modeMutate) && (strcmp(argv[arg], "-i") == 0)) {
    double  p = strtodouble(argv[arg+1]);
    char    a = argv[arg+2][0];

    arg += 2;

    setProbabilityInsert(p, a);
  }

  else if ((mode == modeMutate) && (strcmp(argv[arg], "-d") == 0)) {
    double  p = strtodouble(argv[arg+1]);
    char    a = argv[arg+2][0];

    arg += 2;

    setProbabilityDelete(p, a);
  }

  else if ((mode == modeMutate) && (strcmp(argv[arg], "-seed") == 0)) {
    mt.mtSetSeed(strtouint32(argv[++arg]));
  }

  else {
    return(false);
  }

  return(true);
}



void
mutateParameters::showUsage(opMode mode) {

  if (mode != modeMicroSatellite)
    return;

  fprintf(stderr, "OPTIONS for mutate mode:\n");
  fprintf(stderr, "  -s p a b   set the probability of substituting 'a' with 'b' to 'p'.\n");
  fprintf(stderr, "  -i p a     set the probability of inserting an 'a' to 'p'.\n");
  fprintf(stderr, "  -d p a     set the probability of deleting an 'a' to 'p'.\n");
}



bool
mutateParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeMutate)
    return(false);

  for (uint32 ii=0; ii<256; ii++)
    pSubstitute[ii] = 0.0;

  pInsert = 0.0;

  for (uint32 ii=0; ii<256; ii++) {
    for (uint32 jj=0; jj<256; jj++)
      pSubstitute[ii] += pS[ii][jj];

    pInsert           += pI[ii];
    pDelete           += pD[ii];
  }

  return(errors.size() > 0);
}



void
doMutate_substitute(double p,
                    char   base, uint8  qual,
                    char  *ns,   uint8 *nq,    uint32 &oo,
                    mutateParameters &mutPar) {

  for (uint32 xx=0; xx<256; xx++) {
    if (p < mutPar.pS[base][xx]) {
      //fprintf(stderr, "sub %c -> %c at pos %u\n", base, xx, oo);
      ns[oo] = xx;
      nq[oo] = 0;
      oo++;
      break;
    }

    p -= mutPar.pS[base][xx];
  }
}



void
doMutate_insert(double p,
                char   base, uint8  qual,
                char  *ns,   uint8 *nq,    uint32 &oo,
                mutateParameters &mutPar) {

  for (uint32 xx=0; xx<256; xx++) {
    if (p < mutPar.pI[xx]) {
      ns[oo] = xx;
      nq[oo] = 0;
      oo++;

      ns[oo] = base;
      nq[oo] = qual;
      oo++;

      break;
    }

    p -= mutPar.pI[xx];
  }
}



void
doMutate(vector<char const *> &inputs, mutateParameters &mutPar) {
  dnaSeq     seq;

  uint32     nMax   = 0;
  char      *nBases = nullptr;
  uint8     *nQuals = nullptr;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf     = new dnaSeqFile(inputs[ff]);
    uint64       num    = 0;

    while (sf->loadSequence(seq)) {
      uint32  iPos = 0;   //  Position in the input read
      uint32  oLen = 0;   //  Position in the output read

      //  Resize the output to fit the input string with the expected
      //  number of insert/delete made.

      uint32  expectedLen = seq.length() * (1 + 2 * mutPar.pInsert - mutPar.pDelete);

      merylutil::resizeArrayPair(nBases, nQuals, 0, nMax, expectedLen);

      //  Over every base, randomly substitue, insert or delete bases.

      for (iPos=0, oLen=0; iPos<seq.length(); iPos++) {
        char   base = seq.bases()[iPos];
        uint8  qual = seq.quals()[iPos];
        double p    = 0.0;

        //  Whoops!  Resize again?
        if (oLen + 2 > nMax)
          merylutil::resizeArrayPair(nBases, nQuals, oLen, nMax, nMax + 1000);

        //  If a chance of doing something, make a random number.
        if (mutPar.pSubstitute[base] + mutPar.pInsert + mutPar.pD[base] > 0.0)
          p = mutPar.mt.mtRandomRealClosed();

        //  If a substitution, make it.
        if (p < mutPar.pSubstitute[base]) {
          doMutate_substitute(p, base, qual, nBases, nQuals, oLen, mutPar);
          continue;
        }
        p -= mutPar.pSubstitute[base];


        //  If an insertion, make it.
        if (p < mutPar.pInsert) {
          doMutate_insert(p, base, qual, nBases, nQuals, oLen, mutPar);
          continue;
        }
        p -= mutPar.pInsert;


        //  If a deletion, make it.  Hamm, nothing to do here,
        //  just don't add the base to the output.
        if (p < mutPar.pD[base]) {
          continue;
        }
        p -= mutPar.pD[base];


        //  Otherwise, no change.  Just copy the base.
        nBases[oLen] = base;
        nQuals[oLen] = qual;
        oLen++;
      }


      //  And one more for an inserting at the end of the read.
      double p = mutPar.mt.mtRandomRealClosed();

      if (p < mutPar.pInsert)
        doMutate_insert(p, 0, 0, nBases, nQuals, oLen, mutPar);


      //  All done changing.  Output the modified read.

      merylutil::outputFASTA(stdout, nBases, oLen, 0, seq.ident());
    }

    delete sf;
  }
}


