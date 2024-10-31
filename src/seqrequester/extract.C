
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
extractParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "extract") == 0) {
    mode = modeExtract;
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-sequences") == 0)) {
    seqsArgs.push_back(argv[++arg]);
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-length") == 0)) {
    lensArgs.push_back(argv[++arg]);
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-bases") == 0)) {
    decodeRange(argv[++arg], baseBgn, baseEnd);
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-invert-seqslist") == 0)) {
    invertSeqsList = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-invert-lenslist") == 0)) {
    invertLensList = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-invert-baselist") == 0)) {
    invertBaseList = true;
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-full-header") == 0)) {
    asIdentOnly = false;
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-reverse") == 0)) {
    asReverse = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-complement") == 0)) {
    asComplement = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-rc") == 0)) {
    asReverse = true;
    asComplement = true;
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-compress") == 0)) {
    asCompressed = true;
  }

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-upcase") == 0)) {
    asUpperCase = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-downcase") == 0)) {
    asLowerCase = true;
  }

  //else if ((mode == modeExtract) && (strcmp(argv[arg], "-lowermask") == 0)) {
  //  doMasking = true;
  //  maskWithN = false;
  //}
  //else if ((mode == modeExtract) && (strcmp(argv[arg], "-nmask") == 0)) {
  //  doMasking = true;
  //  maskWithN = true;
  //}

  else if ((mode == modeExtract) && (strcmp(argv[arg], "-fasta") == 0)) {
    outputFASTA = true;
  }
  else if ((mode == modeExtract) && (strcmp(argv[arg], "-fastq") == 0)) {
    outputFASTQ = true;

    if ((arg+1 < argc) && ('0' <= argv[arg+1][0]) && (argv[arg+1][0] <= '9'))
      outputQV = strtouint32(argv[++arg]);
  }

  else {
    return(false);
  }

  return(true);
}



void
extractParameters::showUsage(opMode mode) {

  if (mode != modeExtract)
    return;

  fprintf(stderr, "OPTIONS for extract mode:\n");
  fprintf(stderr, "  SEQUENCE SELECTION:\n");
  fprintf(stderr, "    -sequences seqlist    select sequences with ordinal listed in the 'seqlist'\n");
  fprintf(stderr, "    -sequences namefile   select sequences with name listed in 'namefile'\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -length min-max       select sequence if it is at least 'min' bases and at\n");
  fprintf(stderr, "                          most 'max' bases long - that is, inclusive min-to-max\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -bases baselist       select bases listed in 'baselist' from each sequence\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -invert-seqslist      select all unselected sequences for output\n");
  fprintf(stderr, "    -invert-lenslist      select all unselected lengths for output\n");
  fprintf(stderr, "    -invert-baselist      select all unselected bases for output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  OUTPUT FORMAT:\n");
  fprintf(stderr, "    -full-header          output the full header line instead of the first word\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -reverse              reverse the bases in the sequence\n");
  fprintf(stderr, "    -complement           complement the bases in the sequence\n");
  fprintf(stderr, "    -rc                   alias for -reverse -complement\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -compress             compress homopolymer runs to one base\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -upcase\n");
  fprintf(stderr, "    -downcase\n");
  fprintf(stderr, "\n");
  //fprintf(stderr, "    -lowermask\n");
  //fprintf(stderr, "    -nmask\n");
  //fprintf(stderr, "\n");
  fprintf(stderr, "    -fasta                write output as FASTA\n");
  fprintf(stderr, "    -fastq [q]            write output as FASTQ; if no quality values, use q\n");
  fprintf(stderr, "                          (integer, 0-based) for all bases; default q=20\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    a 'seqlist' is a set of integers formed from any combination\n");
  fprintf(stderr, "    of the following, separated by commas:\n");
  fprintf(stderr, "         num       a single number\n");
  fprintf(stderr, "         bgn-end   a range of numbers:  bgn <= end\n");
  fprintf(stderr, "    sequences are 1-based; -sequences 1,3-5 will print the first, third,\n");
  fprintf(stderr, "    fourth and fifth sequences in the input.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    a 'namefile' contains the names of sequences to select, one per line,\n");
  fprintf(stderr, "    without any leading '>' or '@'.  a name is the first word on the ident line.\n");
  fprintf(stderr, "    NOTA BENE: at most one namefile can be supplied.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    a 'baselist' is a set of integer ranges, separated by commas, indicating\n");
  fprintf(stderr, "    the spaces-between-letters bounding the sequence to select.\n");
  fprintf(stderr, "    -bases 0-2,4-7 will select letters ABEFG from the following sequence:\n");
  fprintf(stderr, "       0th A 1st B 2nd C 3ed D 4th E 5th F 6th G 7th H 8th ...\n");
  fprintf(stderr, "       ---   ---   ---   ---   ---   ---   ---   ---   ---  <- spaces\n");
  fprintf(stderr, "\n");
}


template<typename TT>
void
invertBgnEnd(std::vector<TT> &oBgn, std::vector<TT> &oEnd) {
  std::vector<TT>  nBgn;
  std::vector<TT>  nEnd;

  if (oBgn.size() == 0)                      //  If no ranges, do nothing.
    return;

  if (oBgn[0] > 0) {                         //  If the first interval is not at the start, add
    nBgn.push_back(0);                       //  a new interval covering the start of the
    nEnd.push_back(oBgn[0]);                 //  sequence.
  }

  oBgn.push_back(uint64max);                 //  Make the last new interval end at infinity.

  for (uint32 ii=0; ii<oEnd.size(); ii++) {  //  Loop over each interval, adding a new
    nBgn.push_back(oEnd[ii]);                //  one covering the gap between ii and ii+1.
    nEnd.push_back(oBgn[ii+1]);              //  The last new one we add is from the end
  }                                          //  of th last interval to infinity.

  oBgn = nBgn;
  oEnd = nEnd;
}


bool
extractParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeExtract)
    return(false);

  if (inputs.size() == 0)
    sprintf(errors, "ERROR:  No input sequence files supplied.\n");

  //  BASES
  //   - space based

  for (auto arg : baseArgs)                       //  Decode base ranges.
    decodeRange(arg, baseBgn, baseEnd);

  for (uint32 bi=0; bi<baseBgn.size(); bi++)      //  Fix a quirk in the decoding;
    if (baseBgn[bi] == baseEnd[bi])               //  single numbers decode as bgn == end,
      baseEnd[bi] += 1;                           //  but we need space-based.

  for (uint32 bi=0; bi<baseBgn.size(); bi++)      //  Check for invalid ranges.
    if (baseEnd[bi] <= baseBgn[bi])
      sprintf(errors, "ERROR: base range %lu-%lu is invalid, must be increasing.\n",
              baseBgn[bi], baseEnd[bi]);

  if (invertBaseList == true)                     //  Invert ranges?
    invertBgnEnd(baseBgn, baseEnd);

  if (baseBgn.size() == 0) {                      //  Add default range of
    baseBgn.push_back(0);                         //  "everything" if no ranges
    baseEnd.push_back(UINT64_MAX);                //  supplied.
  }

#if 0
  for (uint32 ii=0; ii<baseBgn.size(); ii++)
    fprintf(stderr, "BASES: #%02u - %4lu-%-4lu\n", ii, baseBgn[ii], baseEnd[ii]);
#endif

  //  SEQUENCES
  //   - to the user, sequences begin at ONE and are inclusive.
  //   - to us, sequences are C-style.

  for (auto arg : seqsArgs)                       //  Decode sequence ranges, either
    if (merylutil::fileExists(arg) == true)       //  a file or an integer range.
      seqsName.load(arg, merylutil::splitWords);
    else
      decodeRange(arg, seqsBgn, seqsEnd);

  for (uint32 si=0; si<seqsBgn.size(); si++)      //  Check for invalid ranges.
    if (seqsBgn[si] == 0)
      sprintf(errors, "ERROR: sequences begin at 1, not zero.\n");

  for (uint32 si=0; si<seqsBgn.size(); si++)      //  Convert to C-style.
    seqsBgn[si] -= 1;

  if (invertSeqsList == true)                     //  Invert ranges?
    invertBgnEnd(seqsBgn, seqsEnd);

  if ((seqsName.size() == 0) &&                   //  Add default range of
      (seqsBgn.size() == 0)) {                    //  "everything" if no ranges
    seqsBgn.push_back(0);                         //  supplied.
    seqsEnd.push_back(UINT64_MAX);
  }

#if 0
  for (uint32 ii=0; ii<seqsBgn.size(); ii++)
    fprintf(stderr, "SEQS:  #%02u - %4lu-%-4lu\n", ii, seqsBgn[ii], seqsEnd[ii]);
#endif

  //  LENGTHS

  for (auto arg : lensArgs)                       //  Decode length ranges.
    decodeRange(arg, lensMin, lensMax);

  if (invertLensList == true)                     //  Invert ranges?
    invertBgnEnd(lensMin, lensMax);

  if (lensMin.size() == 0) {                      //  Add default range of
    lensMin.push_back(0);                         //  "everything" if no ranges
    lensMax.push_back(UINT64_MAX);                //  supplied.
  }

#if 0
  for (uint32 ii=0; ii<lensMin.size(); ii++)
    fprintf(stderr, "LENS:  #%02u - %4lu-%-4lu\n", ii, lensMin[ii], lensMax[ii]);
#endif



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

  return(errors.size() > 0);
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
extractParameters::isDesired(uint32      seqIndex,
                             char const *seqName,
                             uint32      seqLength) {
  bool  desired = false;

  //  Check if this is on our desired index list.
  if (seqsBgn.size() > 0) {
    for (uint32 si=0; (desired == false) && (si < seqsBgn.size()); si++)
      if ((seqsBgn[si] <= seqIndex) &&
          (seqIndex < seqsEnd[si]))
        desired = true;
  }

  //  Check if this is on our desired name list.
  if ((seqsName.size() > 0) && (seqName != nullptr)) {
    for (uint32 si=0; (desired == false) && (si < seqsName.size()); si++)
      if (((invertSeqsList == false) && (strcmp(seqName, seqsName[si]) == 0)) ||
          ((invertSeqsList ==  true) && (strcmp(seqName, seqsName[si]) != 0)))
        desired = true;
  }

  //  Check if the length of this sequence is in our desired range.
  if (lensMin.size() > 0) {
    bool  d = false;

    for (uint32 li=0; (d == false) && (li < lensMin.size()); li++)
      if ((lensMin[li] <= seqLength) &&
          (seqLength <= lensMax[li]))
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

  merylutil::resizeArrayPair(extPar.outputBases, extPar.outputQuals, 0, extPar.outputBasesMax, seq.length() + 1);

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

  if (extPar.outputBasesLen > 0)
    merylutil::outputSequence(stdout,
                              extPar.asIdentOnly ? seq.ident() : seq.header(),
                              extPar.outputBases, extPar.outputQuals, extPar.outputBasesLen,
                              sf->isFASTQ(),
                              extPar.outputFASTA,
                              extPar.outputFASTQ,
                              extPar.outputQV);
  //fprintf(stderr, "\033{%d;20H%9u bp - %s\n", seq.length(), seq.header());
}



void
doExtract(vector<char const *> &inputs,
          extractParameters    &extPar) {
  dnaSeq    seq;

  for (uint32 fi=0; fi<inputs.size(); fi++) {
    dnaSeqFile  *sf  = openSequenceFile(inputs[fi], false);

    //  If no desired names and an index was loaded, we can load exactly the
    //  sequences to output.

    if ((extPar.seqsName.size() == 0) &&    //  Are not doing name lookups.
        (sf->loadIndex() == true)) {        //  Do have an existing index.
      for (uint32 si=0; si<extPar.seqsBgn.size(); si++) {
        uint64  sbgn = min(extPar.seqsBgn[si], sf->numberOfSequences());
        uint64  send = min(extPar.seqsEnd[si], sf->numberOfSequences());

        for (uint64 ss=sbgn; ss<send; ss++)
          if (extPar.isDesired(ss, nullptr, sf->sequenceLength(ss)) &&
              sf->findSequence(ss) &&
              sf->loadSequence(seq))
            printSequence(seq, sf, extPar);
      }
    }

    //  Otherwise (filtering based on sequence name or no index exists), we
    //  must stream the whole file.

    else {
      while (sf->loadSequence(seq) == true)
        if (extPar.isDesired(sf->seqIdx(), seq.ident(), seq.length()))
          printSequence(seq, sf, extPar);
    }

    delete sf;
  }
}
