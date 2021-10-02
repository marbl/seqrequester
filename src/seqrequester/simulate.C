
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
#include "mt19937ar.H"

#include <vector>

using namespace std;



void
simulateParameters::finalize() {

  if (distribName[0]) {
    char const *path = findSharedFile("share/sequence", distribName);

    if (path == NULL) {
      fprintf(stderr, "ERROR: File '%s' doesn't exist, and not in any data directory I know about.\n", distribName);
      exit(1);
    }

    //  If the path is a file -- which it should be -- load it.  Else, fail.

    if (fileExists(path) == true) {
      fprintf(stderr, "load '%s'\n", path);
      dist.loadDistribution(path);
    }

    else {
      fprintf(stderr, "ERROR: File '%s' doesn't exist, and not in any data directory I know about.\n", distribName);
      exit(1);
    }
  }
}



void
doSimulate_loadSequences(simulateParameters  &simPar,
                         vector<dnaSeq *>    &seqs,
                         uint64              &seqLen) {
  fprintf(stderr, "Loading sequences from '%s'\n", simPar.genomeName);

  dnaSeq           *seq    = new dnaSeq;
  dnaSeqFile       *sf     = new dnaSeqFile(simPar.genomeName);

  for (; sf->loadSequence(*seq); seq = new dnaSeq) {
    seqs.push_back(seq);

    seqLen += seq->length();
  }

  delete seq;
  delete sf;

  fprintf(stderr, "Loaded %lu sequences.\n", seqs.size());
}



uint64
doSimulate_findSequenceCircular(vector<dnaSeq *> &seqs,
                                uint64            position,
                                uint64           &bgn,
                                uint64           &end,
                                char             *r,
                                uint64            rLen) {

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = seqs[ss]->length();

    if (position < sl) {
      bgn = position;
      end = position + rLen;

      if (end <= sl) {
        memcpy(r, seqs[ss]->bases() + bgn, sizeof(char) * (end-bgn));

        r[rLen] = 0;
        return(ss);
      }
      else {
        end -= sl;

        uint64  l1 = seqs[ss]->length() - bgn;
        uint64  l2 =                      end;

        memcpy(r, seqs[ss]->bases() + bgn, sizeof(char) * l1);

        while (l2 > sl) {
          memcpy(r+l1, seqs[ss]->bases(), sizeof(char) * sl);
          l1 += sl;
          l2 -= sl;
        }

        memcpy(r+l1, seqs[ss]->bases(), sizeof(char) * l2);

        r[rLen] = 0;
        return(ss);
      }
    }

    position -= sl;
  }

  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



uint64
doSimulate_findSequenceTruncated(vector<dnaSeq *> &seqs,
                                 uint64            position,
                                 uint64           &bgn,
                                 uint64           &end,
                                 char             *r,
                                 uint64           &rLen) {

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = (rLen-1) + seqs[ss]->length();

    if (position < sl) {

      //  Truncated at the start.
      if        (position <      rLen) {
        bgn  = 0;
        end  = position+1;
        memcpy(r, seqs[ss]->bases(), sizeof(char) * (end-bgn));
      }

      //  Full read in the middle of the sequence.
      else if (position < sl - rLen) {
        bgn  = position - (rLen-1);
        end  = bgn + rLen;
        memcpy(r, seqs[ss]->bases() + position - (rLen-1), sizeof(char) * (end-bgn));
      }

      //  Truncated at the end.
      else {
        bgn  = position - (rLen-1);
        end  = seqs[ss]->length();
        memcpy(r, seqs[ss]->bases() + position - (rLen-1), sizeof(char) * (end-bgn));
      }

      rLen = end - bgn;
      r[rLen] = 0;
      return(ss);
    }

    position -= sl;
  }

  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



uint64
doSimulate_findSequence(vector<dnaSeq *> &seqs,
                        uint64            position,
                        uint64           &bgn,
                        uint64           &end,
                        char             *r,
                        uint64            rLen) {

  for (uint64 ss=0; ss < seqs.size(); ss++) {
    uint64  sl = seqs[ss]->length() - (rLen-1);

    if (seqs[ss]->length() < rLen-1)
      continue;

    if (position < sl) {
      bgn = position;
      end = bgn + rLen;

      memcpy(r, seqs[ss]->bases() + bgn, sizeof(char) * (end-bgn));

      r[end-bgn] = 0;
      return(ss);
    }

    position -= sl;
  }

  fprintf(stderr, "Not possible to make reads; desired read length too long?\n");
  exit(1);
  assert(0);
}



void
doSimulate_extract(simulateParameters &simPar,
                   vector<dnaSeq *>   &seqs,
                   uint64              seqLen,        //  Total length of all sequences
                   mtRandom           &mt,
                   uint64              nReadsMax,
                   uint64              nBasesMax) {

  uint64  nReads = 0;
  uint64  nBases = 0;

  uint64  rLen = 0;
  uint64  rMax = 1048576;
  char   *r    = new char [rMax];

  uint64  bgn = 0;
  uint64  end = 0;

  while ((nReads < nReadsMax) &&
         (nBases < nBasesMax)) {

    //  Based on the input length distribution, generate a random read length.

    if (simPar.dist.empty() == false)
      rLen = simPar.dist.getValue(mt.mtRandomRealOpen());
    else
      rLen = mt.mtRandomRealOpen() * (simPar.desiredMaxLength - simPar.desiredMinLength) + simPar.desiredMinLength;

    resizeArray(r, 0, rMax, rLen+1);

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

    uint64  sl = 0;

    if        (simPar.circular == true) {
      sl = seqLen;
    }

    else if (simPar.truncate == true) {
      sl = seqLen + (rLen-1) * seqs.size();
    }

    else {
      for (uint64 ss=0; ss < seqs.size(); ss++)
        if (rLen <= seqs[ss]->length())
          sl += seqs[ss]->length() - (rLen-1);
    }

    uint64  seqID    = 0;
    uint64  position = (uint64)floor(mt.mtRandomRealOpen() * sl);

    //  Search for the sequence that has bases that start at 'position'.

    if        (simPar.circular == true) {
      seqID = doSimulate_findSequenceCircular(seqs, position, bgn, end, r, rLen);
    }

    else if (simPar.truncate == true) {
      seqID = doSimulate_findSequenceTruncated(seqs, position, bgn, end, r, rLen);
    }

    else {
      seqID = doSimulate_findSequence(seqs, position, bgn, end, r, rLen);
    }

    //  Terminate the read and randomly flip it.

    bool flip = (mt.mtRandomRealOpen() < simPar.rcProb) ? true : false;

    if (flip)
      reverseComplementSequence(r, rLen);

    //  And output.

    fprintf(stdout, ">read=%lu,%s,position=%lu-%lu,length=%lu,%s\n",
            nReads+1,
            (flip == false) ? "forward" : "reverse",
            bgn, end, rLen,
            seqs[seqID]->ident());
    fprintf(stdout, "%s\n", r);

    //  Account for the read we just emitted.

    nReads += 1;
    nBases += rLen;
  }

  delete [] r;
}



void
doSimulate_test(simulateParameters &simPar,
                vector<dnaSeq *>   &seqs,
                uint64              seqLen,        //  Total length of all sequences
                mtRandom           &mt,
                uint64              nReadsMax,
                uint64              nBasesMax) {

  uint64  nReads = 0;
  uint64  nBases = 0;

  uint64  rLen = 50;
  uint64  rMax = 1048576;
  char   *r    = new char [rMax];

  uint64  bgn = 0;
  uint64  end = 0;

  //  Compute the length of the sequence we're sampling reads from.

  uint64  sl = 0;

  if        (simPar.circular == true) {
    sl = seqLen;
  }

  else if (simPar.truncate == true) {
    sl = seqLen + (rLen-1) * seqs.size();
  }

  else {
    for (uint64 ss=0; ss < seqs.size(); ss++)
      if (rLen <= seqs[ss]->length())
        sl += seqs[ss]->length() - (rLen-1);
  }

  //  Iterate over every start position, extract the sequence, and output.

  for (uint64 position=0; position < sl; position++) {
    uint64  seqID    = 0;

    if        (simPar.circular == true) {
      seqID = doSimulate_findSequenceCircular(seqs, position, bgn, end, r, rLen);
    }

    else if (simPar.truncate == true) {
      seqID = doSimulate_findSequenceTruncated(seqs, position, bgn, end, r, rLen);
    }

    else {
      seqID = doSimulate_findSequence(seqs, position, bgn, end, r, rLen);
    }

    bool flip = false;

    fprintf(stdout, ">read=%lu,%s,position=%lu-%lu,length=%lu,%s\n",
            nReads+1,
            (flip == false) ? "forward" : "reverse",
            bgn, end, rLen,
            seqs[seqID]->ident());
    fprintf(stdout, "%s\n", r);

    nReads += 1;
    nBases += rLen;
  }

  delete [] r;
}



void
doSimulate(vector<char *>     &inputs,
           simulateParameters &simPar) {
  mtRandom   mt;

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
  uint64            seqLen = 0;

  doSimulate_loadSequences(simPar, seqs, seqLen);

  //  Make reads!

  if (simPar.test == false)
    doSimulate_extract(simPar, seqs, seqLen, mt, nReadsMax, nBasesMax);
  else
    doSimulate_test(simPar, seqs, seqLen, mt, nReadsMax, nBasesMax);

  //  Clean up the reference sequences we loaded.

  for (uint64 ii=0; ii<seqs.size(); ii++)
    delete seqs[ii];
}
