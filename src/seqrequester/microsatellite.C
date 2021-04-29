
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
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

#include "runtime.H"

#include "kmers.H"
#include "sequence.H"
#include "bits.H"

#define OP_NONE   0
#define OP_GA     1
#define OP_GC     2
#define OP_AT     3

void
outAT(dnaSeqFile*   seqFile,
      char*         outPrefix,
      bool          verbose,
      int           window) {

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.AT.bed", outPrefix);
  FILE *f = AS_UTL_openOutputFile(outName);

  sprintf(outName, "%s.AT.%u.bed", outPrefix, window);
  FILE *fw = AS_UTL_openOutputFile(outName);

  dnaSeq    seq;

  while (seqFile->loadSequence(seq) == true) {
    const char* bases = seq.bases();

    bool    hasA = false, hasT = false, inAT = false;
    int64    beg = -1, end = -1;

    int     winMax = seq.length()/window; 
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

void
outGC(dnaSeqFile*   seqFile,
      char*         outPrefix,
      bool          verbose,
      int           window) {

  char    outName[FILENAME_MAX+1];
  sprintf(outName, "%s.GC.bed", outPrefix);
  FILE *fGC = AS_UTL_openOutputFile(outName);

  sprintf(outName, "%s.GC.%u.bed", outPrefix, window);
  FILE *fGCw = AS_UTL_openOutputFile(outName);

  dnaSeq    seq;
  uint64    seqIdx = 0;

  while (seqFile->loadSequence(seq) == true) {
    seqIdx  = seqFile->seqIdx();
    const char* bases = seq.bases();

    bool    hasC = false, hasG = false, inGC = false;
    int64    beg = -1, end = -1;

    int     winMax = seq.length()/window;    
    uint32  gcWinCounts[winMax];
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
          fprintf(fGC, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
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
        fprintf(fGC, "%s\t%lu\t%lu\n", seq.ident(), beg, end);
        for (uint32 jj = beg; jj < end; jj++)     // count for each pos
            gcWinCounts[jj/window]++;
    }

    for (uint32 ii = 0; ii < winMax; ii++)
        fprintf(fGCw, "%s\t%u\t%u\t%u\t%.2f\n", seq.ident(), ii*window, (ii+1) * window, gcWinCounts[ii], ((float) gcWinCounts[ii] * 100) / window);
  }
  fclose(fGC);
  fclose(fGCw);
}

void
outGA(dnaSeqFile*   seqFile,
      char*         outPrefix,
      bool          verbose,
      int           window) {


  char    outName[FILENAME_MAX+1];
  
  // exact bed
  sprintf(outName, "%s.GA.bed", outPrefix);
  FILE *fGA = AS_UTL_openOutputFile(outName);  
  sprintf(outName, "%s.TC.bed", outPrefix);
  FILE *fTC = AS_UTL_openOutputFile(outName);

  // per window count bed
  sprintf(outName, "%s.GA.%u.bed", outPrefix, window);
  FILE *fGAw = AS_UTL_openOutputFile(outName);
  sprintf(outName, "%s.TC.%u.bed", outPrefix, window);
  FILE *fTCw = AS_UTL_openOutputFile(outName);

  dnaSeq    seq;
  uint64    seqIdx = 0;
  
  while (seqFile->loadSequence(seq) == true) {

    seqIdx  = seqFile->seqIdx();
    const char* bases = seq.bases();

    bool    hasG = false, hasA = false, inGA = false;
    bool    hasT = false, hasC = false, inTC = false;

    fprintf(stderr, "Processing %s\n", seq.ident());
    int     winMax = seq.length()/window;
    uint32  gaWinCounts[winMax];
    uint32  tcWinCounts[winMax];

    int64     begGA = -1, endGA = -1;
    int64     begTC = -1, endTC = -1;

    for (uint32 ii=0; ii<seq.length()/window; ii++) {
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
          for (uint32 jj = begGA; jj < endGA; jj++)     // count for each pos
            gaWinCounts[jj/window]++;

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
      for (uint32 jj = begGA; jj < endGA; jj++)
        gaWinCounts[jj/window]++;
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




int
main(int argc, char **argv) {
  char   *seqName     = NULL;
  char   *outPrefix   = NULL;
  bool    verbose     = false;
  uint32  reportType  = OP_NONE;
  int     window      = 128;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  int                        arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-seq") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-prefix") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-ga") == 0) {
      reportType = OP_GA;

    } else if (strcmp(argv[arg], "-gc") == 0) {
      reportType = OP_GC;
    
    } else if (strcmp(argv[arg], "-at") == 0) {
      reportType = OP_AT;


    } else if (strcmp(argv[arg], "-window") == 0) {
      window = atoi(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("No input sequence (-seq) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -seq <in.fasta> -prefix <prefix> (-ga | -gc | -at) \n", argv[0]);
    fprintf(stderr, "  outputs ga/tc, gc, at dimer region\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  fprintf(stderr, "Open sequence '%s'.\n", seqName);
  dnaSeqFile  *seqFile = NULL;
  seqFile = new dnaSeqFile(seqName);

  if (reportType == OP_GA)
    outGA(seqFile, outPrefix, verbose, window);

  if (reportType == OP_GC)
    outGC(seqFile, outPrefix, verbose, window);

  if (reportType == OP_AT)
    outAT(seqFile, outPrefix, verbose, window);

  fprintf(stderr, "Clean up..\n\n");

  delete seqFile;

  fprintf(stderr, "Bye!\n");

  exit(0);
}
