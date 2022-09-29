void
outGA(dnaSeqFile   *seqFile,
      char const   *outPrefix,
      int           window) {


  char    outName[FILENAME_MAX+1];

  // exact bed
  sprintf(outName, "%s.GA.bed", outPrefix);
  FILE *fGA = merylutil::openOutputFile(outName);
  sprintf(outName, "%s.TC.bed", outPrefix);
  FILE *fTC = merylutil::openOutputFile(outName);

  // per window count bed
  sprintf(outName, "%s.GA.%u.bed", outPrefix, window);
  FILE *fGAw = merylutil::openOutputFile(outName);
  sprintf(outName, "%s.TC.%u.bed", outPrefix, window);
  FILE *fTCw = merylutil::openOutputFile(outName);

  dnaSeq    seq;

  while (seqFile->loadSequence(seq) == true) {
    const char* bases = seq.bases();

    bool    hasG = false, hasA = false, inGA = false;
    bool    hasT = false, hasC = false, inTC = false;

    fprintf(stderr, "Processing %s\n", seq.ident());
    int     winMax = seq.length()/window + 1;
    uint32  gaWinCounts[winMax];
    uint32  tcWinCounts[winMax];

    int64     begGA = -1, endGA = -1;
    int64     begTC = -1, endTC = -1;

    for (uint32 ii=0; ii<winMax; ii++) {
      gaWinCounts[ii] = 0;
      tcWinCounts[ii] = 0;
    }

    for (uint32 ii=0; ii<seq.length(); ii++) {
      switch (bases[ii]) {
        case 'c':  //  C
        case 'C':
          hasC = true;
          inTC = true;
          if (begTC == -1)  begTC = ii;
          endTC = ii + 1;
          inGA = false;
          break;
        case 't':  //  T
        case 'T':
          hasT = true;
          inTC = true;
          if (begTC == -1)  begTC = ii;
          endTC = ii + 1;
          inGA = false;
          break;
        case 'a':
        case 'A':
          hasA = true;
          inGA = true;
          if (begGA == -1) begGA = ii;
          endGA = ii + 1;
          inTC = false;
          break;
        case 'g':
        case 'G':
          hasG = true;
          inGA = true;
          if (begGA == -1) begGA = ii;
          endGA = ii + 1;
          inTC = false;
          break;
        default :
          inTC = false;
          inGA = false;
          break;
      }
      if (!inGA) {
        // is not in GA anymore, but had G and A
        if (hasG && hasA) {
          fprintf(fGA, "%s\t%lu\t%lu\n", seq.ident(), begGA, endGA);
          for (uint32 jj = begGA; jj < endGA; jj++) {    // count for each pos
            gaWinCounts[jj/window]++;
          }
        }
        // reset
        begGA = -1;
        hasG = false;
        hasA = false;
      }

      if (!inTC) {
        // is not in TC anymore, but had T and C
        if (hasT && hasC) {
          fprintf(fTC, "%s\t%lu\t%lu\n", seq.ident(), begTC, endTC);
          for (uint32 jj = begTC; jj < endTC; jj++)
            tcWinCounts[jj/window]++;
        }
        // reset
        begTC = -1;
        hasT = false;
        hasC = false;
      }
    }

    // end of each sequence for loop
    // before going to the next chr, write down if there was a GC
    if (inGA && hasG && hasA) {
      fprintf(fGA, "%s\t%lu\t%lu\n", seq.ident(), begGA, endGA);
      for (uint32 jj = begGA; jj < endGA; jj++) {
        gaWinCounts[jj/window]++;
      }
    }

    if (inTC && hasT && hasC) {
      fprintf(fTC, "%s\t%lu\t%lu\n", seq.ident(), begTC, endTC);
      for (uint32 jj = begTC; jj < endTC; jj++)
        tcWinCounts[jj/window]++;
    }

    for (uint32 ii = 0; ii < winMax; ii++) {
      fprintf(fGAw, "%s\t%u\t%u\t%u\t%.2f\n", seq.ident(), ii*window, (ii+1) * window, gaWinCounts[ii], ((float) gaWinCounts[ii] * 100) / window);
      fprintf(fTCw, "%s\t%u\t%u\t%u\t%.2f\n", seq.ident(), ii*window, (ii+1) * window, tcWinCounts[ii], ((float) tcWinCounts[ii] * 100) / window);
    }

  }

  fclose(fGA);
  fclose(fTC);
  fclose(fGAw);
  fclose(fTCw);
}
