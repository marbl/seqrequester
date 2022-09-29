void
outAT(dnaSeqFile   *seqFile,
      char const   *outPrefix,
      int           window) {

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.AT.bed", outPrefix);
  FILE *f = merylutil::openOutputFile(outName);

  sprintf(outName, "%s.AT.%u.bed", outPrefix, window);
  FILE *fw = merylutil::openOutputFile(outName);

  dnaSeq    seq;

  while (seqFile->loadSequence(seq) == true) {
    const char* bases = seq.bases();

    bool    hasA = false, hasT = false, inAT = false;
    int64    beg = -1, end = -1;

    int     winMax = seq.length()/window + 1;
    uint32  atWinCounts[winMax];
    for (uint32 ii = 0; ii < winMax; ii++)
        atWinCounts[ii] = 0;

    fprintf(stderr, "Processing %s\n", seq.ident());

    for (uint32 ii=0; ii<seq.length(); ii++) {

      switch (bases[ii]) {
        case 'a':  //  A
        case 'A':
          hasA = true;
          inAT = true;
          if (beg == -1)  beg = ii;
          end = ii + 1;
          break;
        case 't':  //  T
        case 'T':
          hasT = true;
          inAT = true;
          if (beg == -1)  beg = ii;
          end = ii + 1;
          break;
        default :
          inAT = false;
          break;
      }

      if (!inAT) {
        // is not in GC anymore, but had C and G
        if (hasA && hasT) {
          fprintf(f, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
          for (uint32 jj = beg; jj < end; jj++)     // count for each pos
              atWinCounts[jj/window]++;
        }

        // reset
        beg = -1;
        hasA = false;
        hasT = false;
      }
    }

    //  end of sequence for loop
    //  before going to the next chr, write down if there was a GC
    if  (hasA && hasT) {
        fprintf(f, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
        for (uint32 jj = beg; jj < end; jj++)     // count for each pos
            atWinCounts[jj/window]++;
    }

    for (uint32 ii = 0; ii < winMax; ii++)
        fprintf(fw, "%s\t%u\t%u\t%u\t%.2f\n", seq.ident(), ii*window, (ii+1) * window, atWinCounts[ii], ((float) atWinCounts[ii] * 100) / window);
  }
  fclose(f);
  fclose(fw);
}
