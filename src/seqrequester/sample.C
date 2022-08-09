
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
sampleParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "sample") == 0) {
    mode = modeSample;
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-paired") == 0)) {
    isPaired = true;
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-copies") == 0)) {
    numCopies = strtouint32(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-output") == 0)) {
    output1 = duplicateString(argv[++arg]);   //  MUST be writable; #'s in the name will
    output2 = duplicateString(argv[  arg]);   //  be replaced by '1' or '2' later.
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-fasta") == 0)) {
    outputFASTA = true;
  }
  else if ((mode == modeSample) && (strcmp(argv[arg], "-fastq") == 0)) {
    outputFASTQ = true;

    if ((arg+1 < argc) && ('0' <= argv[arg+1][0]) && (argv[arg+1][0] <= '9'))
      outputQV = strtouint32(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-coverage") == 0)) {      //  Sample reads up to some coverage C
    desiredCoverage = strtodouble(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-genomesize") == 0)) {
    genomeSize = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-bases") == 0)) {         //  Sample B bases
    desiredNumBases = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-reads") == 0)) {         //  Sample N reads
    desiredNumReads = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-pairs") == 0)) {         //  Sample N pairs of reads
    desiredNumReads = strtouint64(argv[++arg]) * 2;
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-fraction") == 0)) {      //  Sample F fraction
    desiredFraction = strtodouble(argv[++arg]);
  }

  else if ((mode == modeSample) && (strcmp(argv[arg], "-seed") == 0)) {          //  Seed for pseudo random number generator
    mt.mtSetSeed(strtouint32(argv[++arg]));
  }

  else {
    return(false);
  }

  return(true);
}



void
sampleParameters::showUsage(opMode mode) {

  if (mode != modeSample)
    return;

  fprintf(stderr, "OPTIONS for sample mode:\n");
  fprintf(stderr, "  -paired             treat inputs as paired sequences; the first two files form the\n");
  fprintf(stderr, "                      first pair, and so on.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -copies C           write C different copies of the sampling (without replacement).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -output O           write output sequences to file O.  If paired, two files must be supplied.\n");
  fprintf(stderr, "  -fasta              write output as FASTA\n");
  fprintf(stderr, "  -fastq [q]          write output as FASTQ; if no quality values, use q (integer, 0-based) for all\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -coverage C         output C coverage of sequences, based on genome size G.\n");
  fprintf(stderr, "  -genomesize G       \n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -bases B            output B bases.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -reads R            output R reads.\n");
  fprintf(stderr, "  -pairs P            output P pairs (only if -paired).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -fraction F         output fraction F of the input bases.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -seed s        seed the pseudo-random number generate with integer 's'\n");
  fprintf(stderr, "\n");
}



bool
sampleParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeSample)
    return(false);

  if (inputs.size() == 0)
    sprintf(errors, "ERROR:  No input sequence files supplied.\n");

  return(errors.size() > 0);
}



class seqEntry {
public:
  seqEntry(mtRandom &mt, uint64 pos_, uint32 len_) {
    rnd = mt.mtRandomRealOpen();
    pos = pos_;
    len = len_;
    out = UINT32_MAX;
  };

  double   rnd;    //  A random number.
  uint64   pos;    //  Position in the seqLengths vector.
  uint32   len;    //  Length of this sequence.
  uint32   out;    //  Output file index for this sequence.
};


bool seqOrderRandom(const seqEntry &a, const seqEntry &b) {
  return(a.rnd < b.rnd);
}

bool seqOrderNormal(const seqEntry &a, const seqEntry &b) {
  return(a.pos < b.pos);
}


static
void
doSample_logInto(sampleParameters &samPar) {
  char into[64];
  char emit[99];


  if (samPar.numCopies == 1)
    snprintf(into, 64, " into one output file");
  else
    snprintf(into, 64, " into each of %u output files", samPar.numCopies);


  if       (samPar.desiredNumReads > 0)
    snprintf(emit, 99, "%lu read%s", samPar.desiredNumReads, samPar.desiredNumReads == 1 ? "" : "s");

  else if ((samPar.desiredNumBases > 0) && (samPar.desiredCoverage > 0))
    snprintf(emit, 99, "%.3fx coverage (%lu base%s)", samPar.desiredCoverage, samPar.desiredNumBases, samPar.desiredNumBases == 1 ? "" : "s");

  else if ((samPar.desiredNumBases > 0) && (samPar.desiredFraction > 0))
    snprintf(emit, 99, "%.4f fraction (%lu base%s)", samPar.desiredFraction, samPar.desiredNumBases, samPar.desiredNumBases == 1 ? "" : "s");

  else if  (samPar.desiredNumBases > 0)
    snprintf(emit, 99, "%lu base%s", samPar.desiredNumBases, samPar.desiredNumBases == 1 ? "" : "s");

  else {
    strcpy(into, "");
    strcpy(emit, "nothing");
  }


  fprintf(stderr, "  %s%s.\n", emit, into);
}


void
doSample_sample(sampleParameters &samPar,
                uint64            numSeqsTotal,
                uint64            numBasesTotal,
                vector<seqEntry> &seqOrder) {

  //  Do some math to figure out what sequences to report.

  if (samPar.desiredCoverage > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredCoverage * samPar.genomeSize);
  }

  if (samPar.desiredFraction > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredFraction * numBasesTotal);
  }

  //  Randomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderRandom);

  doSample_logInto(samPar);

  //  Scan the randomized reads, assigning each to an output file,
  //  and moving to the next file when the current one is too big.

  uint64  nr  = 0;
  uint64  nbe = 0;

  uint32  of = 0;

  if (samPar.desiredNumReads > 0) {
    for (uint64 ii=0; ii<numSeqsTotal; ii++) {
      if (nr < samPar.desiredNumReads) {
        nr++;
        seqOrder[ii].out = of;
      } else {
        nr = 1;
        seqOrder[ii].out = ++of;
      }
    }
  }

  if (samPar.desiredNumBases > 0) {
    for (uint64 nbe=0, ii=0; ii<numSeqsTotal; ii++) {
      if (nbe < samPar.desiredNumBases) {
        nbe += seqOrder[ii].len;
        seqOrder[ii].out = of;
      } else {
        nbe  = seqOrder[ii].len;
        seqOrder[ii].out = ++of;
      }
    }
  }

  //  Unrandomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderNormal);

  //  Warn if we ran out of input.

  if (of <= samPar.numCopies) {
    fprintf(stderr, "\n");
    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "WARNING: ran out of input before all output files filled.\n");
    fprintf(stderr, "WARNING:\n");
  }
}



void
doSample_paired(vector<char const *> &inputs, sampleParameters &samPar) {

  vector<uint64>    numSeqsPerFile;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  //  No support for multiple copies.

  if (samPar.numCopies != 1)
    fprintf(stderr, "ERROR: No support for -copies in paried-end mode.\n"), exit(1);

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  compressedFileWriter *outFile1 = nullptr;
  compressedFileWriter *outFile2 = nullptr;

  {
    char  *a = strrchr(samPar.output1, '#');
    char  *b = strrchr(samPar.output2, '#');

    if (a == nullptr)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output1), exit(1);
    if (b == nullptr)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output2), exit(1);

    *a = '1';
    *b = '2';

    outFile1 = new compressedFileWriter(samPar.output1, 9);
    outFile2 = new compressedFileWriter(samPar.output2, 9);
  }

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;
  dnaSeq   seq2;

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      seqOrder.push_back(seqEntry(samPar.mt, numSeqsTotal, seq1.length() + seq2.length()));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length() + seq2.length();

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
    delete sf2;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      uint32 of = seqOrder[num].out;

      if (of < samPar.numCopies) {
        outputFASTA(outFile1->file(), seq1.bases(), seq1.length(), 0, seq1.ident());
        outputFASTA(outFile2->file(), seq2.bases(), seq2.length(), 0, seq2.ident());
      }

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    delete sf1;
    delete sf2;
  }

  delete outFile1;
  delete outFile2;
}



compressedFileWriter *
doSample_single_openOutput(sampleParameters &samPar, uint32 ii) {
  uint32  ap = 0;

  while ((samPar.output1[ap] != 0) &&
         (samPar.output1[ap] != '#'))
    ap++;

  //  If no #'s and multiple copies requested, fail.  But if only one copy, just open
  //  the file and return.

  if ((samPar.output1[ap] == 0) && (samPar.numCopies > 1))
    fprintf(stderr, "ERROR: Failed to find '#' in output name '%s', and asked to make multiple copies.\n", samPar.output1), exit(1);

  if  (samPar.output1[ap] == 0) {
    fprintf(stderr, "  %s\n", samPar.output1);
    return(new compressedFileWriter(samPar.output1, 9));
  }

  //  We've got #'s in the string.  We want to replace the last block of 'em
  //  with digits (ap found above is the start of the first block, sigh).
  //  Search backwards.

  while (samPar.output1[ap] != 0)                   //  Find the end of the string.
    ap++;

  while ((ap > 0) && (samPar.output1[ap] != '#'))   //  Find the last #.
    ap--;

  while ((ap > 0) && (samPar.output1[ap] == '#'))   //  Find the start of that run.
    ap--;

  if (samPar.output1[ap] != '#')                    //  Handle the stupid edge case.
    ap++;

  //  We're guaranteed to have some #'s in the name.  Count 'em.

  uint32 dp = 0;

  while (samPar.output1[ap + dp] == '#')
    dp++;

  //  Make a copy of the name, insert the digits, and return the file.
  //
  //  Suppose we have three #'s in the string; dp will be 4.
  //  We'll copy in digs[7-3 = 4]; dp will be 2 after this
  //                digs[7-2 = 5]; dp will be 1 after this
  //                digs[7-1 = 6]; dp will be 0 after this
  //                digs[7-0 = 7]; the NUL byte is not copied.

  char  name[FILENAME_MAX+1] = {0};
  char  digs[8];

  snprintf(digs, 8, "%07u", ii);

  for (uint32 ii=0; samPar.output1[ii]; ii++) {
    if ((ii < ap) || (dp == 0))
      name[ii] = samPar.output1[ii];
    else
      name[ii] = digs[7 - dp--];
  }

  fprintf(stderr, "  %s\n", name);
  return(new compressedFileWriter(name, 9));
}



void
doSample_single(vector<char const *> &inputs, sampleParameters &samPar) {
  vector<uint64>    numSeqsPerFile;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  fprintf(stderr, "Opening %u output file%s.\n", samPar.numCopies, (samPar.numCopies == 1) ? "" : "s");

  compressedFileWriter **outFiles = new compressedFileWriter * [samPar.numCopies];

  for (uint32 ii=0; ii<samPar.numCopies; ii++)
    outFiles[ii] = doSample_single_openOutput(samPar, ii);

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  fprintf(stderr, "\n");
  fprintf(stderr, "Scanning %lu input file%s.\n", inputs.size(), (inputs.size() == 1) ? "" : "s");

  dnaSeq   seq1;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      seqOrder.push_back(seqEntry(samPar.mt, numSeqsTotal, seq1.length()));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length();

      num += 1;

      if ((num % 61075) == 0)
        fprintf(stderr, "  scanned %10lu sequences in '%s'\r", num, inputs[ff]);
    }

    fprintf(stderr, "  scanned %10lu sequences in '%s'\n", num, inputs[ff]);

    numSeqsPerFile.push_back(num);

    delete sf1;
  }

  //  Figure out what to output.

  fprintf(stderr, "\n");
  fprintf(stderr, "Randomizing.\n");

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing outputs.\n");

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1) == true) {
      uint32 of = seqOrder[num].out;

      if (of < samPar.numCopies)
        outputSequence(outFiles[of]->file(),
                       seq1.ident(), seq1.bases(), seq1.quals(), seq1.length(),
                       sf1->isFASTQ(),
                       samPar.outputFASTA,
                       samPar.outputFASTQ,
                       samPar.outputQV);

      if ((++num % 61075) == 0)
        fprintf(stderr, "  wrote %10lu sequences from '%s'\r", num, inputs[ff]);
    }

    fprintf(stderr, "  wrote %10lu sequences from '%s'\n", num, inputs[ff]);

    delete sf1;
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Closing %u output file%s.\n", samPar.numCopies, (samPar.numCopies == 1) ? "" : "s");

  for (uint32 ii=0; ii<samPar.numCopies; ii++)
    delete outFiles[ii];

  fprintf(stderr, "\n");
  fprintf(stderr, "Done.\n");
}



void
doSample(vector<char const *> &inputs, sampleParameters &samPar) {

  if (samPar.isPaired == false)
    doSample_single(inputs, samPar);
  else
    doSample_paired(inputs, samPar);
}

