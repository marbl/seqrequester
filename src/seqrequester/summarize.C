
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
summarizeParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if      (strcmp(argv[arg], "summarize") == 0) {
    mode = modeSummarize;
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-size") == 0)) {
    genomeSize = strtoull(argv[++arg], nullptr, 10);
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-1x") == 0)) {
    limitTo1x = true;
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-split-n") == 0)) {
    breakAtN = true;
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-simple") == 0)) {
    asSimple();
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-lengths") == 0)) {
    asLengths();
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-seqlen") == 0)) {
    asSeqLen();
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-assequences") == 0)) {
    asSequences = true;
    asBases     = false;
  }

  else if ((mode == modeSummarize) && (strcmp(argv[arg], "-asbases") == 0)) {
    asSequences = false;
    asBases     = true;
  }

  else {
    return(false);
  }

  return(true);
}



void
summarizeParameters::showUsage(opMode mode) {

  if (mode != modeSummarize)
    return;

  fprintf(stderr, "OPTIONS for summarize mode:\n");
  fprintf(stderr, "  -size          base size to use for N50 statistics\n");
  fprintf(stderr, "  -1x            limit NG table to 1x coverage\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -split-n       split sequences at N bases before computing length\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " Output format:\n");
  fprintf(stderr, "  (default)      output an N50 table, a histogram picture and mono-,\n");
  fprintf(stderr, "                   di- and tri-nucleotide frequencies\n");
  fprintf(stderr, "  -simple        output a simple 'length numSequences' histogram\n");
  fprintf(stderr, "  -lengths       output a list of 'length' for each sequence\n");
  fprintf(stderr, "  -seqlen        output a list of 'length name' for each sequence\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -assequences   load data as complete sequences (for testing)\n");
  fprintf(stderr, "  -asbases       load data as blocks of bases    (for testing)\n");
  fprintf(stderr, "\n");
}



bool
summarizeParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeSummarize)
    return(false);

  if (inputs.size() == 0)
    sprintf(errors, "ERROR:  No input sequence files supplied.\n");

  return(errors.size() > 0);
}



bool
doSummarize_loadSequence(dnaSeqFile  *sf,
                         bool         asSequences,
                         char       *&name,   uint32    &nameMax,
                         char       *&seq,
                         uint8      *&qlt,    uint64    &seqMax,
                         uint64      &seqLen) {
  uint32 error = 0;

  if (asSequences)
    return(sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen, error));

  //  Otherwise, piece it together from multiple calls to get bases.
  //  loadBases() returns true if bases were loaded, and sets endOfSeq
  //  if the block returned is to the end of a sequence.

  uint64   bufferMax = 23;
  uint64   bufferLen = 0;
  char    *buffer    = new char [bufferMax];
  bool     endOfSeq  = false;

  merylutil::resizeArray(name, 0, nameMax, (uint32)1024);
  merylutil::resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1);

  name[0] = 0;
  seq[0]  = 0;
  qlt[0]  = 0;

  seqLen = 0;

  while (sf->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
    if (seqLen + bufferLen >= seqMax)
      merylutil::resizeArrayPair(seq, qlt, seqLen, seqMax, 2 * (seqLen + bufferLen + 1));

    assert(seqLen + bufferLen + 1 < seqMax);

    memcpy(seq + seqLen, buffer, sizeof(char) * bufferLen);
    seqLen += bufferLen;

    seq[seqLen] = 0;

    if (endOfSeq)
      break;
  }

  //  We get here in two ways:
  //    loadBases() immediately hit EOF, it will return false and set endOfSeq = false.
  //    Otherwise, it found a sequence and endOfSeq (here) _must_ be true.
  //  So, we can just return endOfSeq to mean 'there might be another sequence to load'.

  delete [] buffer;

  return(endOfSeq);
}



void
doSummarize_lengthHistogramSimple(uint64           *shortLengths,
                                  uint32            shortLengthsLen,
                                  vector<uint64>   &longLengths) {

  //  Short lengths are super easy; they're already a histogram.

  for (uint64 ll=0; ll<shortLengthsLen; ll++)
    if (shortLengths[ll] > 0)
      fprintf(stdout, "%lu\t%lu\n", ll, shortLengths[ll]);

  //  Long lengths need to be sorted then counted.

  sort(longLengths.begin(), longLengths.end(), less<uint64>());

  for (uint64 bgn=0, end=1; bgn < longLengths.size(); bgn=end, end++) {
    while ((end < longLengths.size()) &&
           (longLengths[bgn] == longLengths[end]))
      end++;

    fprintf(stdout, "%lu\t%lu\n", longLengths[bgn], end - bgn);
  }
}



void
doSummarize_dumpLengths(uint64           *shortLengths,
                        uint32            shortLengthsLen,
                        vector<uint64>   &longLengths) {

  //  Dump the shorter lengths.

  for (uint64 ll=0; ll<shortLengthsLen; ll++) {
    if (shortLengths[ll] == 0)
      continue;

    for (uint64 cc=0; cc<shortLengths[ll]; cc++)
      fprintf(stdout, "%lu\n", ll);
  }

  //  Dump the longer lengths.

  sort(longLengths.begin(), longLengths.end(), less<uint64>());

  for (uint64 ii=0; ii<longLengths.size(); ii++)
    fprintf(stdout, "%lu\n", longLengths[ii]);
}



void
doSummarize_lengthHistogram(uint64          *shortLengths,
                            uint32           shortLengthsLen,
                            vector<uint64>  &longLengths,
                            uint64           genomeSize,
                            bool             limitTo1x) {

  //  NG table dimensions
  uint32   nNG      = 0;                      //  Number of lines in the NG table.

  uint64   lSum  = 0;                         //  Sum of the lengths we've encountered so far

  uint32   nStep = 10;                        //  Step of each N report.
  uint32   nVal  = nStep;                     //  Index of the threshold we're next printing.
  uint64   nThr  = genomeSize * nVal / 100;   //  Threshold lenth; if sum is bigger, emit and move to the next threshold

  //  These need to be sorted at some point, so might as well do it now.
  //  Note!  Sorted high to low.

  sort(longLengths.begin(), longLengths.end(), greater<uint64>());

  //
  //  Find the minimum and maximum lengths, and count the number of
  //  non-zero-length sequences.
  //

  uint64   minLength = uint64max;
  uint64   maxLength = uint64min;
  uint64   nSeqs     = longLengths.size();    //  Number of sequences we're summarizing.

  for (uint64 ll=0; ll<shortLengthsLen; ll++) {
    if (shortLengths[ll] > 0) {
      minLength = std::min(minLength, ll);
      maxLength = std::max(maxLength, ll);
      nSeqs    += shortLengths[ll];
    }
  }

  if (longLengths.size() > 0) {
    minLength = std::min(minLength, longLengths.back());
    maxLength = std::max(minLength, longLengths.front());
  }

  //
  //  Count the number of lines we expect to get in the NG table.
  //

  auto setStep = [&]()
                 {
                   while (lSum >= nThr) {
                     nNG++;

                     if      (nVal <    200)  nVal += nStep;
                     else if (nVal <   2000)  nVal += nStep * 10;
                     else if (nVal <  20000)  nVal += nStep * 100;
                     else if (nVal < 200000)  nVal += nStep * 1000;
                     else                     nVal += nStep * 10000;

                     nThr  = genomeSize * nVal / 100;
                   }
                 };

  for (uint64 ll=0; ll<shortLengthsLen; ll++) {
    for (uint64 cc=0; cc<shortLengths[ll]; cc++) {
      lSum += ll;
      setStep();
    }
  }

  for (uint64 ii=0; ii<longLengths.size(); ii++) {
    lSum += longLengths[ii];
    setStep();
  }

  if (nNG < 10)
    nNG = 10;

  //
  //  Decide how many rows to make in the length histogram table, and how
  //  many sequences each row will capture.
  //
  //  Set it to max(40, nNG) -- either a comfortable table or whatever the NG
  //  table is -- but reset to 1 sequence per row if it is too big.
  //

  uint32   nLHcols    = 63;                     //  Magic number to make the histogram the same width as the trinucleotide list
  uint32   nLH        = 40;                     //  Default height of the histogram; dynamically set below.

  uint32   bucketSize = (uint32)ceil((double)(maxLength - minLength) / std::max(nLH, nNG));

  if (bucketSize == 0) {
    nLH        = 1;
    bucketSize = 1;
  }

  nLH = 1 + (maxLength - minLength) / bucketSize;   //  With new bucketSize set, compute actual number of rows.

  if (lSum == 0)
    nLH = 0;

  //
  //  Compute the length histogram.
  //

  uint32  *nSeqPerLen = new uint32 [nLH];

  for (uint32 rr=0; rr<nLH; rr++)                      //  Clear the histogram.
    nSeqPerLen[rr] = 0;

  for (uint64 ll=0; ll<shortLengthsLen; ll++) {         //  Count the number of sequences per size range.
    if (shortLengths[ll] > 0) {
      uint32 r = (ll - minLength) / bucketSize;

      assert(r < nLH);
      nSeqPerLen[r] += shortLengths[ll];
    }
  }

  for (uint64 ii=0; ii<longLengths.size(); ii++) {
    uint32 r = (longLengths[ii] - minLength) / bucketSize;

    assert(r < nLH);
    nSeqPerLen[r]++;
  }

  //  Find the maximum size of any bucket.

  uint64  maxCount = 1;

  for (uint32 rr=0; rr<nLH; rr++)
    if (maxCount < nSeqPerLen[rr])
      maxCount = nSeqPerLen[rr];

  //
  //  Generate the length histogram table.
  //

  lSum  = 0;                                            //  Reset for actually generating the length histogram.
  nStep = 10;
  nVal  = nStep;
  nThr  = genomeSize * nVal / 100;

  //  For each line in the histogram table, write the size range and draw a picture.

  char **histPlot = new char * [nLH];

  for (uint32 rr=0; rr<nLH; rr++) {
    uint32  nn = (uint32)ceil(nSeqPerLen[rr] * nLHcols / (double)maxCount);
    uint64  lo = (rr+0) * bucketSize + minLength;
    uint64  hi = (rr+1) * bucketSize + minLength - 1;

    //fprintf(stdout, "rr %u bucketsize %u minLength %lu lo %lu hi %lu\n", rr, bucketSize, minLength, lo, hi);

    histPlot[rr] = new char [28 + nLHcols + 1];         //  28 = 9 + 1 + 9 + 1 + 7 + 1

    if (lo == hi)                                       //  Size range and number of sequences.
      sprintf(histPlot[rr], "%9lu           %7u|",
              lo, nSeqPerLen[rr]);
    else
      sprintf(histPlot[rr], "%9lu-%-9lu %7u|",
              lo, hi, nSeqPerLen[rr]);

    for (uint32 cc=0; cc<nn; cc++)                      //  ...histogram bars.
      histPlot[rr][28 + cc] = '-';

    histPlot[rr][28 + nn]   = 0;                        //  ...terminate the line.
    //fprintf(stdout, "histplot'%s'\n", histPlot[rr]);
  }

  //
  //  Output N table, with length histogram appended at the end of each line.
  //
  uint32  np = 0;
  uint32  hp = 0;

  fprintf(stdout, "\n");
  fprintf(stdout, "G=%-12lu"      "                     sum of  ||               length     num\n", genomeSize);
  fprintf(stdout,   "NG         length     index       lengths  ||                range    seqs\n");
  fprintf(stdout,   "----- ------------ --------- ------------  ||  ------------------- -------\n");

  //  Write lines if we're showing all data, or if we're below 1x coverage.

  auto emitLine = [&](uint64 seqlen, uint64 seqnum)
                  {
                    lSum += seqlen;

                    while (lSum >= nThr) {
                      if ((limitTo1x == false) ||
                          (nVal <= 100)) {
                        if (hp < nLH)
                          fprintf(stdout, "%05u %12lu %9lu %12lu  ||  %s\n", nVal, seqlen, seqnum, lSum, histPlot[hp++]);
                        else
                          fprintf(stdout, "%05u %12lu %9lu %12lu  ||\n",     nVal, seqlen, seqnum, lSum);
                      }

                      if      (nVal <    200)   nVal += nStep;
                      else if (nVal <   2000)   nVal += nStep * 10;
                      else if (nVal <  20000)   nVal += nStep * 100;
                      else if (nVal < 200000)   nVal += nStep * 1000;
                      else                      nVal += nStep * 10000;

                      nThr  = genomeSize * nVal / 100;
                    }
                  };

  uint64  ns = 1;   //  Sequences start at 1, not zero!

  for (uint64 ii=0; ii<longLengths.size(); ii++, ns++)
    emitLine(longLengths[ii], ns);

  for (uint64 ll=shortLengthsLen; ll-- > 0; )
    for (uint64 cc=0; cc<shortLengths[ll]; cc++, ns++)
      emitLine(ll, ns);

  assert(ns-1 == nSeqs);

  //  Output up to 1x coverage regardless of how much is there.

  while (nVal <= 100) {
    if (hp < nLH)
      fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||  %s\n",
              nVal, "-", "-", "-",
              histPlot[hp++]);
    else
      fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||\n",
              nVal, "-", "-", "-");

    nVal += nStep;
  }

  //  If we're displaying exactly 1x, write empty lines to get to there.
#if 0
  if (limitTo1x == true) {
    while (nVal <= 100) {
      if (hp < nLH)
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||  %s\n",
                nVal, "-", "-", "-",
                histPlot[hp++]);
      else
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||\n",
                nVal, "-", "-", "-");

      nVal += nStep;
    }
  }
#endif

  //  Now the summary line for the NG table.

  if (genomeSize == 0) {
    if (hp < nLH)
      fprintf(stdout, "%07.3fx           %9lu %12lu  ||  %s\n", 0.0, nSeqs, lSum, histPlot[hp++]);
    else
      fprintf(stdout, "%07.3fx           %9lu %12lu  ||\n",     0.0, nSeqs, lSum);
  }
  else {
    if (hp < nLH)
      fprintf(stdout, "%07.3fx           %9lu %12lu  ||  %s\n", (double)lSum / genomeSize, nSeqs, lSum, histPlot[hp++]);
    else
      fprintf(stdout, "%07.3fx           %9lu %12lu  ||\n",     (double)lSum / genomeSize, nSeqs, lSum);
  }


  //  And any remaining length table lines.

  while (hp < nLH)
    fprintf(stdout, "                                           ||  %s\n", histPlot[hp++]);

  fprintf(stdout, "\n");

  //  Cleanup.

  for (uint32 rr=0; rr<nLH; rr++)
    delete [] histPlot[rr];

  delete [] histPlot;
  delete [] nSeqPerLen;
}



void
doSummarize(vector<char const *> &inputs,
            summarizeParameters  &sumPar) {
  uint32          shortLengthsLen = 0;
  uint64         *shortLengths    = nullptr;
  vector<uint64>  longLengths;

  uint64          nSeqs  = 0;
  uint64          nBases = 0;

  uint32          mer = 0;

  uint64          mn[4]     = {0};
  uint64          dn[4*4]   = {0};
  uint64          tn[4*4*4] = {0};

  double          nmn = 0;
  double          ndn = 0;
  double          ntn = 0;

  uint32          nameMax = 0;
  char           *name    = nullptr;
  uint64          seqMax  = 0;
  char           *seq     = nullptr;
  uint8          *qlt     = nullptr;
  uint64          seqLen  = 0;

  if (sumPar.isSeqLen()) {
    fprintf(stdout, "--------- ------------  ----------------\n");
    fprintf(stdout, "-index          length  sequence-name\n");
    fprintf(stdout, "--------- ------------  ----------------\n");
  }

  merylutil::allocateArray(shortLengths, shortLengthsLen, 1048576);

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf = new dnaSeqFile(inputs[ff]);

    //  If sequence,
    //    Count mono-, di- and tri-nucleotides.
    //    Count number of mono-, di- and tri-nucleotides.
    //    Count number of sequences and total bases.
    //    Save the lengths of sequences.

    while (doSummarize_loadSequence(sf, sumPar.asSequences, name, nameMax, seq, qlt, seqMax, seqLen) == true) {
      uint64  pos = 0;
      uint64  bgn = 0;

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
      }

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
      }

      while (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
        tn[mer & 0x3f]++;
      }

      nmn +=                    (seqLen-0);
      ndn += (seqLen < 2) ? 0 : (seqLen-1);
      ntn += (seqLen < 3) ? 0 : (seqLen-2);

      //  Report seq-len if requensted.

      if (sumPar.isSeqLen())
        fprintf(stdout, "%-9lu %12lu  %s\n", sf->seqIdx()+1, seqLen, name);

      //  If we're NOT splitting on N, add one sequence of the given length.

      if (sumPar.breakAtN == false) {
        nSeqs  += 1;
        nBases += seqLen;

        if (seqLen < shortLengthsLen)
          shortLengths[seqLen]++;
        else
          longLengths.push_back(seqLen);

        continue;
      }

      //  But if we ARE splitting on N, add multiple sequences.

      pos = 0;
      bgn = 0;

      while (pos < seqLen) {
        while ((pos < seqLen) && ((seq[pos] == 'n') ||    //  Skip N's.
                                  (seq[pos] == 'N')))
          pos++;

        bgn = pos;                                        //  Remember our start position.

        while ((pos < seqLen) && ((seq[pos] != 'n') &&    //  Move ahead until the end of sequence or an N.
                                  (seq[pos] != 'N')))
          pos++;

        if (pos - bgn > 0) {                              //  If a non-empty sequence
          nSeqs  += 1;                                    //  summarize it.
          nBases += pos - bgn;

          if (pos - bgn < shortLengthsLen)
            shortLengths[pos - bgn]++;
          else
            longLengths.push_back(pos - bgn);
        }
      }
    }   //  Over sequences in the file.

    delete sf;
  }   //  Over input files.

  delete [] name;
  delete [] seq;
  delete [] qlt;

  if (sumPar.genomeSize == 0)      //  If no genome size supplied, set it to the sum of lengths.
    sumPar.genomeSize = nBases;

  //  If only a simple histogram of lengths is requested, dump and done.

  if (sumPar.isSimple() == true) {
    doSummarize_lengthHistogramSimple(shortLengths, shortLengthsLen, longLengths);
  }

  //  If only the read lengths are requested, dump and done.

  else if (sumPar.isLengths() == true) {
    doSummarize_dumpLengths(shortLengths, shortLengthsLen, longLengths);
  }

  //  Otherwise, generate a fancy histogram plot.
  //  And finish with the mono-, di- and tri-nucleotide frequencies.

  else if (sumPar.isComplex() == true) {
    doSummarize_lengthHistogram(shortLengths, shortLengthsLen, longLengths, sumPar.genomeSize, sumPar.limitTo1x);

    if (nmn == 0)  nmn = 1;   //  Avoid divide by zero.
    if (ndn == 0)  ndn = 1;
    if (ntn == 0)  ntn = 1;

    double gc = 100.0 * (mn[0x01] + mn[0x03]) / nmn;
    double at = 100.0 * (mn[0x00] + mn[0x02]) / nmn;

#define FMT "%12lu %6.4f"
#define GC "%05.02f%%"

    fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "       mononucleotide          dinucleotide                                                                                  trinucleotide\n");
    fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
    fprintf(stdout, ""             FMT " A" FMT " AA" FMT " AAA " FMT " AAC " FMT " AAG " FMT " AAT\n", mn[0x00], mn[0x00] / nmn, dn[0x00], dn[0x00] / ndn, tn[0x00], tn[0x00] / ntn, tn[0x01], tn[0x01] / ntn, tn[0x03], tn[0x03] / ntn, tn[0x02], tn[0x02] / ntn);
    fprintf(stdout, ""             FMT " C" FMT " AC" FMT " ACA " FMT " ACC " FMT " ACG " FMT " ACT\n", mn[0x01], mn[0x01] / nmn, dn[0x01], dn[0x01] / ndn, tn[0x04], tn[0x04] / ntn, tn[0x05], tn[0x05] / ntn, tn[0x07], tn[0x07] / ntn, tn[0x06], tn[0x06] / ntn);
    fprintf(stdout, ""             FMT " G" FMT " AG" FMT " AGA " FMT " AGC " FMT " AGG " FMT " AGT\n", mn[0x03], mn[0x03] / nmn, dn[0x03], dn[0x03] / ndn, tn[0x0c], tn[0x0c] / ntn, tn[0x0d], tn[0x0d] / ntn, tn[0x0f], tn[0x0f] / ntn, tn[0x0e], tn[0x0e] / ntn);
    fprintf(stdout, ""             FMT " T" FMT " AT" FMT " ATA " FMT " ATC " FMT " ATG " FMT " ATT\n", mn[0x02], mn[0x02] / nmn, dn[0x02], dn[0x02] / ndn, tn[0x08], tn[0x08] / ntn, tn[0x09], tn[0x09] / ntn, tn[0x0b], tn[0x0b] / ntn, tn[0x0a], tn[0x0a] / ntn);
    fprintf(stdout, "                     " FMT " CA" FMT " CAA " FMT " CAC " FMT " CAG " FMT " CAT\n",                           dn[0x04], dn[0x04] / ndn, tn[0x10], tn[0x10] / ntn, tn[0x11], tn[0x11] / ntn, tn[0x13], tn[0x13] / ntn, tn[0x12], tn[0x12] / ntn);
    fprintf(stdout, "      --GC--  --AT-- " FMT " CC" FMT " CCA " FMT " CCC " FMT " CCG " FMT " CCT\n",                           dn[0x05], dn[0x05] / ndn, tn[0x14], tn[0x14] / ntn, tn[0x15], tn[0x15] / ntn, tn[0x17], tn[0x17] / ntn, tn[0x16], tn[0x16] / ntn);
    fprintf(stdout, "      " GC "  " GC " " FMT " CG" FMT " CGA " FMT " CGC " FMT " CGG " FMT " CGT\n", gc, at,                   dn[0x07], dn[0x07] / ndn, tn[0x1c], tn[0x1c] / ntn, tn[0x1d], tn[0x1d] / ntn, tn[0x1f], tn[0x1f] / ntn, tn[0x1e], tn[0x1e] / ntn);
    fprintf(stdout, "                     " FMT " CT" FMT " CTA " FMT " CTC " FMT " CTG " FMT " CTT\n",                           dn[0x06], dn[0x06] / ndn, tn[0x18], tn[0x18] / ntn, tn[0x19], tn[0x19] / ntn, tn[0x1b], tn[0x1b] / ntn, tn[0x1a], tn[0x1a] / ntn);
    fprintf(stdout, "                     " FMT " GA" FMT " GAA " FMT " GAC " FMT " GAG " FMT " GAT\n",                           dn[0x0c], dn[0x0c] / ndn, tn[0x30], tn[0x30] / ntn, tn[0x31], tn[0x31] / ntn, tn[0x33], tn[0x33] / ntn, tn[0x32], tn[0x32] / ntn);
    fprintf(stdout, "                     " FMT " GC" FMT " GCA " FMT " GCC " FMT " GCG " FMT " GCT\n",                           dn[0x0d], dn[0x0d] / ndn, tn[0x34], tn[0x34] / ntn, tn[0x35], tn[0x35] / ntn, tn[0x37], tn[0x37] / ntn, tn[0x36], tn[0x36] / ntn);
    fprintf(stdout, "                     " FMT " GG" FMT " GGA " FMT " GGC " FMT " GGG " FMT " GGT\n",                           dn[0x0f], dn[0x0f] / ndn, tn[0x3c], tn[0x3c] / ntn, tn[0x3d], tn[0x3d] / ntn, tn[0x3f], tn[0x3f] / ntn, tn[0x3e], tn[0x3e] / ntn);
    fprintf(stdout, "                     " FMT " GT" FMT " GTA " FMT " GTC " FMT " GTG " FMT " GTT\n",                           dn[0x0e], dn[0x0e] / ndn, tn[0x38], tn[0x38] / ntn, tn[0x39], tn[0x39] / ntn, tn[0x3b], tn[0x3b] / ntn, tn[0x3a], tn[0x3a] / ntn);
    fprintf(stdout, "                     " FMT " TA" FMT " TAA " FMT " TAC " FMT " TAG " FMT " TAT\n",                           dn[0x08], dn[0x08] / ndn, tn[0x20], tn[0x20] / ntn, tn[0x21], tn[0x21] / ntn, tn[0x23], tn[0x23] / ntn, tn[0x22], tn[0x22] / ntn);
    fprintf(stdout, "                     " FMT " TC" FMT " TCA " FMT " TCC " FMT " TCG " FMT " TCT\n",                           dn[0x09], dn[0x09] / ndn, tn[0x24], tn[0x24] / ntn, tn[0x25], tn[0x25] / ntn, tn[0x27], tn[0x27] / ntn, tn[0x26], tn[0x26] / ntn);
    fprintf(stdout, "                     " FMT " TG" FMT " TGA " FMT " TGC " FMT " TGG " FMT " TGT\n",                           dn[0x0b], dn[0x0b] / ndn, tn[0x2c], tn[0x2c] / ntn, tn[0x2d], tn[0x2d] / ntn, tn[0x2f], tn[0x2f] / ntn, tn[0x2e], tn[0x2e] / ntn);
    fprintf(stdout, "                     " FMT " TT" FMT " TTA " FMT " TTC " FMT " TTG " FMT " TTT\n",                           dn[0x0a], dn[0x0a] / ndn, tn[0x28], tn[0x28] / ntn, tn[0x29], tn[0x29] / ntn, tn[0x2b], tn[0x2b] / ntn, tn[0x2a], tn[0x2a] / ntn);
    fprintf(stdout, "\n");
  }

  delete [] shortLengths;
}

