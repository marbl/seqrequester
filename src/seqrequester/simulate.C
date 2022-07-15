
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
simulateParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "simulate") == 0) {
    mode = modeSimulate;
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-genomesize") == 0)) {
    genomeSize = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-coverage") == 0)) {
    desiredCoverage = strtodouble(argv[++arg]);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-nreads") == 0)) {
    desiredNumReads = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-nbases") == 0)) {
    desiredNumBases = strtouint64(argv[++arg]);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-circular") == 0)) {
    circular = true;
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-truncate") == 0)) {
    truncate = true;
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-genome") == 0)) {
    genomeName = argv[++arg];
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-distribution") == 0)) {
    distribName = argv[++arg];
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-length") == 0)) {
    decodeRange(argv[++arg], desiredMinLength, desiredMaxLength);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-reverse") == 0)) {
    rcProb = strtodouble(argv[++arg]);
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-seed") == 0)) {
    mt.mtSetSeed(strtouint32(argv[++arg]));
  }

  else if ((mode == modeSimulate) && (strcmp(argv[arg], "-test") == 0)) {
    test = true;
  }

  else {
    return(false);
  }

  return(true);
}



void
simulateParameters::showUsage(opMode mode) {

  if (mode != modeSimulate)
    return;

  fprintf(stderr, "OPTIONS for simulate mode:\n");
  fprintf(stderr, "  -genome G           sample reads from these sequences\n");
  fprintf(stderr, "  -circular           treat the sequences in G as circular\n");
  fprintf(stderr, "  -truncate           sample uniformly, at expense of read length\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -genomesize g       genome size to use for deciding coverage below\n");
  fprintf(stderr, "  -coverage c         generate approximately c coverage of output\n");
  fprintf(stderr, "  -nreads n           generate exactly n reads of output\n");
  fprintf(stderr, "  -nbases n           generate approximately n bases of output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -distribution F     generate read length by sampling the distribution in file F\n");
  fprintf(stderr, "                        one column  - each line is the length of a sequence\n");
  fprintf(stderr, "                        two columns - each line has the 'length' and 'number of sequences'\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "                      if file F doesn't exist, use a built-in distribution\n");
  fprintf(stderr, "                        ultra-long-nanopore\n");
  fprintf(stderr, "                        pacbio\n");
  fprintf(stderr, "                        pacbio-hifi\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -length min[-max]   generate read length uniformly from range min-max\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -reverse p          output a reverse-complement read with probability p\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -test               generate a read at every possible start position\n");
  fprintf(stderr, "\n");
}



bool
simulateParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeSimulate)
    return(false);

  if ((genomeName == nullptr) ||
      (genomeName[0] == 0))
    sprintf(errors, "ERROR:  No reference genome sequence (-genome) supplied.\n");

  //  Load any read length distribution.

  if ((distribName == nullptr) ||
      (distribName[0] == 0)) {
    char const *path = findSharedFile("share/sequence", distribName);

    if ((path == nullptr) ||
        (fileExists(path) == false)) {
      sprintf(errors, "ERROR: File '%s' doesn't exist, and not in any data directory I know about.\n", distribName);
      return(true);
    }

    dist.loadDistribution(path);
  }

  return(errors.size() > 0);
}



class simRead {
public:
  simRead() {
  };

  ~simRead() {
    delete [] seq;
  };

  void
  setSize(uint64 l) {
    resizeArray(seq, 0, seqMax, l+1);
  };

  bool
  reverseComplement(bool flip) {
    if (flip)
      reverseComplementSequence(seq, seqLen);

    assert(seq[0] != 0);

    return(flip);
  };

  uint64   srcID  = 0;
  uint64   srcBgn = 0;
  uint64   srcEnd = 0;

  uint64   seqLen = 0;
  uint64   seqMax = 1024;
  char    *seq    = new char [seqMax];
};



void
doSimulate_loadSequences(simulateParameters  &simPar,
                         vector<dnaSeq *>    &seqs,
                         uint64              &seqsLen) {
  fprintf(stderr, "Loading sequences from '%s'\n", simPar.genomeName);

  dnaSeq           *seq    = new dnaSeq;
  dnaSeqFile       *sf     = new dnaSeqFile(simPar.genomeName);

  for (; sf->loadSequence(*seq); seq = new dnaSeq) {
    seqs.push_back(seq);

    seqsLen += seq->length();
  }

  delete seq;
  delete sf;

  fprintf(stderr, "Loaded %lu sequences.\n", seqs.size());
}



void
doSimulate_findSequenceCircular(vector<dnaSeq *> &seqs,
                                uint64            desiredStart,
                                uint64            desiredLength,
                                simRead          *read) {

  read->setSize(desiredLength);

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = seqs[ss]->length();

    if (desiredStart >= sl) {   //  Desired read is
      desiredStart -= sl;       //  after this
      continue;                 //  sequence ends.
    }

    read->srcID  = ss;
    read->srcBgn = desiredStart;
    read->srcEnd = desiredStart + desiredLength;
    read->seqLen = desiredLength;

    //  If the read is entirely within the sequence, copy consecutive bases.

    if (read->srcEnd <= sl) {
      memcpy(read->seq, seqs[ss]->bases() + read->srcBgn, sizeof(char) * read->seqLen);

      read->seq[read->seqLen] = 0;
      return;
    }

    //  Otherwise, copy bases from the end, then whole copies, then bases
    //  from the start.

    else {
      uint64  l1 = seqs[ss]->length() - read->srcBgn;
      uint64  l2 = desiredLength - l1;

      memcpy(read->seq, seqs[ss]->bases() + read->srcBgn, sizeof(char) * l1);

      while (l2 > sl) {
        memcpy(read->seq+l1, seqs[ss]->bases(), sizeof(char) * sl);
        l1 += sl;
        l2 -= sl;
      }

      memcpy(read->seq+l1, seqs[ss]->bases(), sizeof(char) * l2);

      read->srcEnd = l2;

      read->seq[read->seqLen] = 0;
      return;
    }
  }

  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



void
doSimulate_findSequenceTruncated(vector<dnaSeq *> &seqs,
                                 uint64            desiredStart,
                                 uint64            desiredLength,
                                 simRead          *read) {

  read->setSize(desiredLength);

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = (desiredLength-1) + seqs[ss]->length();

    if (desiredStart >= sl) {   //  Desired read is
      desiredStart -= sl;       //  after this
      continue;                 //  sequence ends.
    }

    read->srcID  = ss;
    read->srcBgn = desiredStart;
    read->srcEnd = desiredStart + desiredLength;
    read->seqLen = desiredLength;

    //  We've conceptually extended the sequence by desiredLength-1 bases at
    //  the start.  If the read starts in those bases, the read is
    //  'truncated' to the start of the read.  Likewise for the end of the
    //  sequence.

    //  Truncated at both ends.
    if ((desiredStart                 < desiredLength) &&
        (desiredStart + desiredLength > sl)) {
      read->srcBgn = 0;
      read->srcEnd = seqs[ss]->length();
    }

    //  Truncated at the start.
    else if (desiredStart < desiredLength) {
      read->srcBgn = 0;
      read->srcEnd = desiredStart - (desiredLength-1) + desiredLength;
    }

    //  Truncated at the end.
    else if (desiredStart + desiredLength > sl) {
      read->srcBgn = desiredStart - (desiredLength-1);
      read->srcEnd = seqs[ss]->length();
    }

    //  Full read in the middle of the sequence.
    else {
      read->srcBgn = desiredStart - (desiredLength-1);
      read->srcEnd = desiredStart - (desiredLength-1) + desiredLength;
    }

    read->seqLen = read->srcEnd - read->srcBgn;   //  Updated for truncation!

    memcpy(read->seq, seqs[ss]->bases() + read->srcBgn, sizeof(char) * read->seqLen);

    read->seq[read->seqLen] = 0;

    return;
  }


  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



void
doSimulate_findSequenceContained(vector<dnaSeq *> &seqs,
                                 uint64            desiredStart,
                                 uint64            desiredLength,
                                 simRead          *read) {

  read->setSize(desiredLength);

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = seqs[ss]->length() - (desiredLength-1);

    if (seqs[ss]->length() < desiredLength-1)   //  Completely ignore short
      continue;                                 //  sequences.

    if (desiredStart >= sl) {   //  Desired read is
      desiredStart -= sl;       //  after this
      continue;                 //  sequence ends.
    }

    read->srcID  = ss;
    read->srcBgn = desiredStart;
    read->srcEnd = desiredStart + desiredLength;
    read->seqLen = read->srcEnd - read->srcBgn;

    memcpy(read->seq, seqs[ss]->bases() + read->srcBgn, sizeof(char) * read->seqLen);

    read->seq[read->seqLen] = 0;

    return;
  }

  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



void
doSimulate_findSequence(simulateParameters &simPar,
                        vector<dnaSeq *>   &seqs,
                        uint64              desiredStart,
                        uint64              desiredLength,
                        simRead            *read,
                        uint64              readID) {

  //  Search for the sequence that has bases that start at 'seqsPos'.

  if        (simPar.circular == true) {
    doSimulate_findSequenceCircular(seqs, desiredStart, desiredLength, read);
  }

  else if (simPar.truncate == true) {
    doSimulate_findSequenceTruncated(seqs, desiredStart, desiredLength, read);
  }

  else {
    doSimulate_findSequenceContained(seqs, desiredStart, desiredLength, read);
  }

  //  Randomly reverse-complement.

  bool flipped = read->reverseComplement(simPar.mt.mtRandomRealOpen() < simPar.rcProb);

  //  And output.

  fprintf(stdout, ">read=%lu,%s,position=%lu-%lu,length=%lu,%s\n",
          readID,
          (flipped == false) ? "forward" : "reverse",
          read->srcBgn, read->srcEnd, read->seqLen,
          seqs[read->srcID]->ident());
  fprintf(stdout, "%s\n", read->seq);
}



void
doSimulate_extract(simulateParameters &simPar,
                   vector<dnaSeq *>   &seqs,
                   uint64              seqsLen,        //  Total length of all sequences
                   uint64              nReadsMax,
                   uint64              nBasesMax) {
  uint64    nReads = 0;
  uint64    nBases = 0;

  simRead  *read   = new simRead;

  while ((nReads < nReadsMax) &&
         (nBases < nBasesMax)) {

    //  Based on the input length distribution, generate a random read length.

    uint64  desiredLength;

    if (simPar.dist.empty() == true)
      desiredLength = simPar.desiredMinLength + (simPar.desiredMaxLength - simPar.desiredMinLength) * simPar.mt.mtRandomRealOpen();
    else
      desiredLength = simPar.dist.getValue(simPar.mt.mtRandomRealOpen());

    //  Compute a position in the sequences.
    //
    //  If we're circular, any position is valid.
    //
    //  If we're allowed to truncate reads, we can start a read anywhere from
    //  -readLen to the end of the sequence.
    //
    //  If neither, we cannot start a new read in the last N bases of the
    //  sequence.  If the read length is longer than the reference length, we
    //  cannot get any read out of this sequence.  Thus, we need to
    //  explicitly sum the lengths of sequences for each read length.

    uint64  effectiveLength = 0;
    uint64  desiredStart    = 0;

    if        (simPar.circular == true) {
      effectiveLength = seqsLen;
    }

    else if (simPar.truncate == true) {
      effectiveLength = seqsLen + (desiredLength-1) * seqs.size();
    }

    else {
      for (uint64 ss=0; ss < seqs.size(); ss++)
        if (desiredLength <= seqs[ss]->length())
          effectiveLength += seqs[ss]->length() - (desiredLength-1);
    }

    desiredStart = (uint64)floor(simPar.mt.mtRandomRealOpen() * effectiveLength);

    doSimulate_findSequence(simPar, seqs, desiredStart, desiredLength, read, nReads+1);

    //  Account for the read we just emitted.

    nReads += 1;
    nBases += read->seqLen;
  }

  delete read;
}



void
doSimulate_test(simulateParameters &simPar,
                vector<dnaSeq *>   &seqs,
                uint64              seqsLen,        //  Total length of all sequences
                uint64              nReadsMax,
                uint64              nBasesMax) {
  uint64  nReads = 0;
  uint64  nBases = 0;

  simRead  *read   = new simRead;

  //  Compute the length of the sequence we're sampling reads from.

  uint64  desiredLength = simPar.desiredMinLength;
  uint64  effectiveLength  = 0;

  if        (simPar.circular == true) {
    effectiveLength = seqsLen;
  }

  else if (simPar.truncate == true) {
    effectiveLength = seqsLen + (desiredLength-1) * seqs.size();
  }

  else {
    for (uint64 ss=0; ss < seqs.size(); ss++)
      if (desiredLength <= seqs[ss]->length())
        effectiveLength += seqs[ss]->length() - (desiredLength-1);
  }

  //  Iterate over every start position, extract the sequence, and output.

  for (uint64 desiredStart=0; desiredStart < effectiveLength; desiredStart++) {
    doSimulate_findSequence(simPar, seqs, desiredStart, desiredLength, read, nReads+1);

    nReads += 1;
    nBases += read->seqLen;
  }

  delete read;
}



void
doSimulate(vector<char const *> &inputs,
           simulateParameters   &simPar) {

  //  Decide how many reads or bases to make.

  uint64  nReadsMax = uint64max;
  uint64  nBasesMax = uint64max;

  if (simPar.desiredCoverage > 0)
    nBasesMax = simPar.desiredCoverage * simPar.genomeSize;

  if (simPar.desiredNumReads > 0)
    nReadsMax = simPar.desiredNumReads;

  if (simPar.desiredNumBases > 0)
    nBasesMax = simPar.desiredNumBases;

  //  Fail?

  if ((nBasesMax == uint64max) &&
      (nReadsMax == uint64max)) {
    fprintf(stderr, "ERROR: Don't know how much data to simulate.  Set -coverage, -nreads or -nbases.\n");
    exit(1);
  }

  if ((simPar.dist.empty() == true) &&
      (simPar.desiredMaxLength < simPar.desiredMinLength)) {
    fprintf(stderr, "ERROR: Don't know how long to make the reads.  Set -distribution or -length.\n");
    exit(1);
  }

  //  Load the genome sequences.

  vector<dnaSeq *>  seqs;
  uint64            seqsLen = 0;

  doSimulate_loadSequences(simPar, seqs, seqsLen);
  // Reset nBaseMax when the coverage is given while the genomeSize is not
  if ( (simPar.desiredCoverage > 0) && (simPar.genomeSize == 0) ) {
      simPar.genomeSize = seqsLen;
      nBasesMax = simPar.desiredCoverage * simPar.genomeSize;
  }

  //  Make reads!

  if (simPar.test == false)
    doSimulate_extract(simPar, seqs, seqsLen, nReadsMax, nBasesMax);
  else
    doSimulate_test(simPar, seqs, seqsLen, nReadsMax, nBasesMax);

  //  Clean up the reference sequences we loaded.

  for (uint64 ii=0; ii<seqs.size(); ii++)
    delete seqs[ii];
}


