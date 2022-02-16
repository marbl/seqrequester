
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

#include "arrays.H"
#include "sequence.H"



extractParameters::extractParameters() {
}

extractParameters::~extractParameters() {
  for (uint64 ii=0; ii<seqsName.size(); ii++)
    delete [] seqsName[ii];
}



void
extractParameters::finalize(void) {

  //  Decode base ranges.

  for (auto arg : baseArgs)
    decodeRange(arg, baseBgn, baseEnd);

  //  Decode sequence ranges.

  for (auto arg : seqsArgs) {
    fprintf(stderr, "Parse '%s'\n", arg);
    if (fileExists(arg) == true)
      seqsName.load(arg, splitWords);
    else
      decodeRange(arg, seqsBgn, seqsEnd);
  }

  //  Decode length ranges.

  for (auto arg : lensArgs)
    decodeRange(arg, lensMin, lensMax);

  //  If no base range specified, output all bases.

  if (baseBgn.size() == 0) {
    baseBgn.push_back(0);
    baseEnd.push_back(UINT64_MAX);
  }

  //  If no sequence range specified, output all sequences.

  if (seqsBgn.size() == 0) {
    seqsBgn.push_back(1);
    seqsEnd.push_back(UINT64_MAX);
  }

  //  If no length restriction, output all lengths.

  if (lensMin.size() == 0) {
    lensMin.push_back(0);
    lensMax.push_back(UINT64_MAX);
  }

  //  Check and adjust the sequence ranges.
  //
  //  To the user, sequences begin at ONE, not ZERO.
  //  To us, sequences begin at zero.

  for (uint32 si=0; si<seqsBgn.size(); si++) {
    if (seqsBgn[si] == 0) {
      fprintf(stderr, "ERROR: sequences begin at 1, not zero.\n");
      exit(1);
    }

    seqsBgn[si] -= 1;
  }

  //  Check and adjust the base ranges.  These are space based.  A quirk in the
  //  command line parsing results in bgn == end if a single number is supplied;
  //  we interpret that to mean 'output the base at space N'.

  for (uint32 bi=0; bi<baseBgn.size(); bi++) {
    if (baseBgn[bi] == baseEnd[bi])
      baseEnd[bi] += 1;

    if (baseEnd[bi] <= baseBgn[bi]) {
      fprintf(stderr, "ERROR: base range %lu-%lu is invalid, must be increasing.\n",
              baseBgn[bi], baseEnd[bi]);
      exit(1);
    }
  }

  //  Allocate a big string, but since we don't know the actual max length,
  //  we'll have to be careful to reallocate as needed.

  outputBasesLen  = 0;
  outputBasesMax  = 1048576;
  outputBases     = new char  [outputBasesMax];
  outputQuals     = new uint8 [outputBasesMax];

  //  Initialize complement. toUpper and toLower arrays.

  C['a'] = 't';  U['a'] = 'A';  L['a'] = 'a';
  C['c'] = 'g';  U['c'] = 'C';  L['c'] = 'c';
  C['g'] = 'c';  U['g'] = 'G';  L['g'] = 'g';
  C['t'] = 'a';  U['t'] = 'T';  L['t'] = 't';

  C['A'] = 'T';  U['A'] = 'A';  L['A'] = 'a';
  C['C'] = 'G';  U['C'] = 'C';  L['C'] = 'c';
  C['G'] = 'C';  U['G'] = 'G';  L['G'] = 'g';
  C['T'] = 'A';  U['T'] = 'T';  L['T'] = 't';
}




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



bool
isDesired(uint32             seqIndex,
          char const        *seqName,
          uint32             seqLength,
          extractParameters &extPar) {
  bool  desired = false;

  //  Check if this is on our desired index list.
  if (extPar.seqsBgn.size() > 0) {
    for (uint32 si=0; (desired == false) && (si < extPar.seqsBgn.size()); si++)
      if ((extPar.seqsBgn[si] <= seqIndex) &&
          (seqIndex <= extPar.seqsEnd[si]))
        desired = true;
  }

  //  Check if this is on our desired name list (with a special case for .
  if ((extPar.seqsName.size() > 0) && (seqName != nullptr)) {
    for (uint32 si=0; (desired == false) && (si < extPar.seqsName.size()); si++)
      if (strcmp(seqName, extPar.seqsName[si]) == 0)
        desired = true;
  }

  //  Check if the length of this sequence is in our desired range.
  if (extPar.lensMin.size() > 0) {
    bool  d = false;

    for (uint32 li=0; (d == false) && (li < extPar.lensMin.size()); li++)
      if ((extPar.lensMin[li] <= seqLength) &&
          (seqLength <= extPar.lensMax[li]))
        d = true;

    if (d == false)
      desired = false;
  }

  return(desired);
}



void
printSequence(dnaSeq &seq,
              dnaSeqFile *sf,
              extractParameters &extPar) {

  //  Increase space for the output copy

  resizeArrayPair(extPar.outputBases, extPar.outputQuals, 0, extPar.outputBasesMax, seq.length() + 1);

  //  Copy desired bases to the output string.

  extPar.outputBasesLen = 0;

  for (uint32 bi=0; bi<extPar.baseBgn.size(); bi++) {
    uint64  bbgn = extPar.baseBgn[bi];
    uint64  bend = extPar.baseEnd[bi];

    bbgn = min(bbgn, seq.length());
    bend = min(bend, seq.length());

    if (bbgn == bend)
      continue;

    memcpy(extPar.outputBases + extPar.outputBasesLen, seq.bases() + bbgn, bend - bbgn);
    memcpy(extPar.outputQuals + extPar.outputBasesLen, seq.quals() + bbgn, bend - bbgn);

    extPar.outputBasesLen += bend - bbgn;
  }

  extPar.outputBases[extPar.outputBasesLen] = 0;
  extPar.outputQuals[extPar.outputBasesLen] = 0;

  //  Reverse?

  if (extPar.asReverse) {
    reverse(extPar.outputBases, extPar.outputBases + extPar.outputBasesLen);
    reverse(extPar.outputQuals, extPar.outputQuals + extPar.outputBasesLen);
  }

  //  Complement?

  if (extPar.asComplement)
    for (uint32 ii=0; ii<extPar.outputBasesLen; ii++)
      extPar.outputBases[ii] = extPar.C[extPar.outputBases[ii]];

  //  Uppercase?

  if (extPar.asUpperCase)
    for (uint32 ii=0; ii<extPar.outputBasesLen; ii++)
      extPar.outputBases[ii] = extPar.U[extPar.outputBases[ii]];

  //  Lowercase?

  if (extPar.asLowerCase)
    for (uint32 ii=0; ii<extPar.outputBasesLen; ii++)
      extPar.outputBases[ii] = extPar.L[extPar.outputBases[ii]];

  //  Homopolymer compressed?

  if (extPar.asCompressed)
    extPar.outputBasesLen = doExtract_compress(extPar.outputBases, extPar.outputQuals, extPar.outputBasesLen);

  //  Output!

  outputSequence(stdout,
                 seq.ident(), extPar.outputBases, extPar.outputQuals, extPar.outputBasesLen,
                 sf->isFASTQ(),
                 extPar.outputFASTA,
                 extPar.outputFASTQ,
                 extPar.outputQV);
}



void
doExtract(vector<char *>    &inputs,
          extractParameters &extPar) {
  dnaSeq    seq;

  for (uint32 fi=0; fi<inputs.size(); fi++) {
    dnaSeqFile  *sf   = new dnaSeqFile(inputs[fi], false);

    //  If desired names exist or the file is compressed, load every sequence
    //  and output a sequence if is desired.
    //
    if ((extPar.seqsName.size() > 0) || (sf->isCompressed() == true)) {
      while (sf->loadSequence(seq) == true)
        if (isDesired(sf->seqIdx(), seq.ident(), seq.length(), extPar))
          printSequence(seq, sf, extPar);
    }

    //  But with no names and an index, we can load exactly the sequences we
    //  want to output.  isDesired() here is only checking the length;
    //  false from either findSequence() or loadSequence() should probably be
    //  reported as errors.
    //
    else {
      for (uint32 si=0; si<extPar.seqsBgn.size(); si++) {
        uint64  sbgn = min(extPar.seqsBgn[si], sf->numberOfSequences());
        uint64  send = min(extPar.seqsEnd[si], sf->numberOfSequences());

        for (uint64 ss=sbgn; ss<send; ss++)
          if (isDesired(sf->seqIdx(), nullptr, sf->sequenceLength(ss), extPar) &&
              sf->findSequence(ss) &&
              sf->loadSequence(seq))
            printSequence(seq, sf, extPar);
      }
    }

    delete sf;
  }
}
