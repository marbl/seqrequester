
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
#include "sequence.H"



uint64
doExtract_compress(char    *outputBases,
                   uint8   *outputQuals,
                   uint64   outputBasesLen) {

  if (outputBasesLen == 0)
    return(0);

  uint64  cc = 0;                  //  NOTE:  Also used in stores/sqStore.C
  uint64  rr = 1;

  while (rr < outputBasesLen) {
    if (outputBases[cc] == outputBases[rr])
      rr++;
    else {
      cc++;
      outputBases[cc] = outputBases[rr];
      outputQuals[cc] = outputBases[rr];
      rr++;
    }
  }

  cc++;
  outputBases[cc] = 0;
  outputQuals[cc] = 0;

  return(cc);

}



void
doExtract(vector<char *>    &inputs,
          extractParameters &extPar) {

  char            C[256] = {0};
  char            U[256] = {0};
  char            L[256] = {0};

  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  uint64  outputBasesLen  = 0;
  uint64  outputBasesMax  = 0;
  char   *outputBases     = NULL;
  uint8  *outputQuals     = NULL;

  //  Initialize complement. toUpper and toLower arrays.

  C['a'] = 't';  U['a'] = 'A';  L['a'] = 'a';
  C['c'] = 'g';  U['c'] = 'C';  L['c'] = 'c';
  C['g'] = 'c';  U['g'] = 'G';  L['g'] = 'g';
  C['t'] = 'a';  U['t'] = 'T';  L['t'] = 't';

  C['A'] = 'T';  U['A'] = 'A';  L['A'] = 'a';
  C['C'] = 'G';  U['C'] = 'C';  L['C'] = 'c';
  C['G'] = 'C';  U['G'] = 'G';  L['G'] = 'g';
  C['T'] = 'A';  U['T'] = 'T';  L['T'] = 't';



  for (uint32 fi=0; fi<inputs.size(); fi++) {
    dnaSeqFile  *sf   = new dnaSeqFile(inputs[fi], true);

    //  Allocate a string big enough to hold the largest output.
    //
    //  Later, maybe, we can analyze the bases to output and make this exactly the correct size.

    uint64  maxStringLength = 0;

    for (uint32 ss=0; ss<sf->numberOfSequences(); ss++)
      maxStringLength = max(maxStringLength, sf->sequenceLength(ss));

    resizeArrayPair(outputBases, outputQuals, 0, outputBasesMax, maxStringLength + 1);

    //fprintf(stderr, "seqs - length %u first %u %u\n", extPar.seqsBgn.size(), extPar.seqsBgn[0], extPar.seqsEnd[0]);

    for (uint32 si=0; si<extPar.seqsBgn.size(); si++) {
      uint64  sbgn = extPar.seqsBgn[si];
      uint64  send = extPar.seqsEnd[si];

      sbgn = min(sbgn, sf->numberOfSequences());
      send = min(send, sf->numberOfSequences());

      //fprintf(stderr, "sbgn %lu send %lu\n", sbgn, send);

      for (uint32 ss=sbgn; ss<send; ss++) {
        uint64  seqLen = sf->sequenceLength(ss);

        //fprintf(stderr, "lens - length %lu first %lu %lu\n", extPar.lensBgn.size(), extPar.lensBgn[0], extPar.lensEnd[0]);

        for (uint32 li=0; li<extPar.lensBgn.size(); li++) {
          uint64  lmin = extPar.lensBgn[li];
          uint64  lmax = extPar.lensEnd[li];

          if ((seqLen < lmin) ||
              (lmax < seqLen))
            seqLen = UINT64_MAX;
        }

        if (seqLen == UINT64_MAX)
          continue;

        if (sf->findSequence(ss) == false) {
          //fprintf(stderr, "Failed to find sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        if (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen) == false) {
          //fprintf(stderr, "Failed to load sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        //fprintf(stderr, "base - length %lu first %lu %lu\n", extPar.baseBgn.size(), extPar.baseBgn[0], extPar.baseEnd[0]);

        outputBasesLen = 0;

        for (uint32 bi=0; bi<extPar.baseBgn.size(); bi++) {
          uint64  bbgn = extPar.baseBgn[bi];
          uint64  bend = extPar.baseEnd[bi];

          bbgn = min(bbgn, seqLen);
          bend = min(bend, seqLen);

          //fprintf(stderr, "base - seqLen %lu base[%u] %lu %lu limited %lu %lu\n", seqLen, bi, extPar.baseBgn[bi], extPar.baseEnd[bi], bbgn, bend);

          if (bbgn == bend)
            continue;

          memcpy(outputBases + outputBasesLen, seq + bbgn, bend - bbgn);
          memcpy(outputQuals + outputBasesLen, qlt + bbgn, bend - bbgn);

          outputBasesLen += bend - bbgn;
        }

        outputBases[outputBasesLen] = 0;
        outputQuals[outputBasesLen] = 0;

        if (extPar.asReverse) {
          reverse(outputBases, outputBases + outputBasesLen);
          reverse(outputQuals, outputQuals + outputBasesLen);
        }

        if (extPar.asComplement)
          for (uint32 ii=0; ii<outputBasesLen; ii++)
            outputBases[ii] = C[outputBases[ii]];

        if (extPar.asUpperCase)
          for (uint32 ii=0; ii<outputBasesLen; ii++)
            outputBases[ii] = U[outputBases[ii]];

        if (extPar.asLowerCase)
          for (uint32 ii=0; ii<outputBasesLen; ii++)
            outputBases[ii] = L[outputBases[ii]];

        if (extPar.asCompressed)
          outputBasesLen = doExtract_compress(outputBases, outputQuals, outputBasesLen);

        outputSequence(stdout,
                       name, outputBases, outputQuals, outputBasesLen,
                       sf->isFASTA(),
                       sf->isFASTQ(),
                       extPar.outputFASTA,
                       extPar.outputFASTQ,
                       extPar.outputQV);
      }
    }

    //  Done with this file.  Get rid of it.

    delete sf;
  }

  //  Cleanup.

  delete [] name;
  delete [] seq;
  delete [] qlt;
  delete [] outputBases;
}





