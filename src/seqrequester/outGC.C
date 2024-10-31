void
outGC(char const   *seqFileName,
      char const   *outPrefix,
      int           window) {

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.GC.bed", outPrefix);
  FILE *f = merylutil::openOutputFile(outName);

  sprintf(outName, "%s.GC.%u.bed", outPrefix, window);
  FILE *fw = merylutil::openOutputFile(outName);

  dnaSeqFile *seqFile = openSequenceFile(seqFileName);
  dnaSeq      seq;

  while (seqFile->loadSequence(seq) == true) {
    const char* bases = seq.bases();

    bool    hasC = false, hasG = false, inGC = false;
    int64    beg = -1, end = -1;

    int     winMax = seq.length()/window + 1;
    uint32 *gcWinCounts = new uint32 [winMax];
    for (uint32 ii = 0; ii < winMax; ii++)
        gcWinCounts[ii] = 0;

    fprintf(stderr, "Processing %s\n", seq.ident());

    for (uint32 ii=0; ii<seq.length(); ii++) {

      switch (bases[ii]) {
        case 'c':  //  C
        case 'C':
          hasC = true;
          inGC = true;
          if (beg == -1)  beg = ii;
          end = ii + 1;
          break;
        case 'g':  //  G
        case 'G':
          hasG = true;
          inGC = true;
          if (beg == -1)  beg = ii;
          end = ii + 1;
          break;
        default :
          inGC = false;
          break;
      }

      if (!inGC) {
        // is not in GC anymore, but had C and G
        if (hasC && hasG) {
          fprintf(f, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
          for (uint32 jj = beg; jj < end; jj++)     // count for each pos
              gcWinCounts[jj/window]++;
        }

        // reset
        beg = -1;
        hasC = false;
        hasG = false;
      }
    }

    //  end of sequence for loop
    //  before going to the next chr, write down if there was a GC
    if  (hasC && hasG) {
        fprintf(f, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
        for (uint32 jj = beg; jj < end; jj++)     // count for each pos
            gcWinCounts[jj/window]++;
    }

    for (uint32 ii = 0; ii < winMax; ii++)
        fprintf(fw, "%s\t%u\t%u\t%u\t%.2f\n", seq.ident(), ii*window, (ii+1) * window, gcWinCounts[ii], ((float) gcWinCounts[ii] * 100) / window);

    delete [] gcWinCounts;
  }

  delete seqFile;

  fclose(f);
  fclose(fw);
}
