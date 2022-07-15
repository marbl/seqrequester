
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

#include "seqrequester.H"

#include "outAT.C"
#include "outGC.C"
#include "outGA.C"

bool
microsatelliteParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if ((strcmp(argv[arg], "microsatellite") == 0) ||
      (strcmp(argv[arg], "microsat") == 0)) {
    mode = modeMicroSatellite;
  }

  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-prefix") == 0)) {
    outPrefix = argv[++arg];
  }

  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-window") == 0)) {
    window = strtouint32(argv[++arg]);
  }

  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-ga") == 0)) {
    report_ga = true;
  }
  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-gc") == 0)) {
    report_gc = true;
  }
  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-at") == 0)) {
    report_at = true;
  }

  else if ((mode == modeMicroSatellite) && (strcmp(argv[arg], "-legacy") == 0)) {
    report_legacy = true;
  }

  else {
    return(false);
  }

  return(true);
}



void
microsatelliteParameters::showUsage(opMode mode) {

  if (mode != modeMicroSatellite)
    return;

  fprintf(stderr, "OPTIONS for microsatellite mode:\n");
  fprintf(stderr, "  -prefix P      write output to <P>.<pattern>.bed\n");
  fprintf(stderr, "  -window w      compute in windows of size w; default write output to <prefix.<pattern>.bed\n");

  fprintf(stderr, "  -ga            compute GA/TC repeat content\n");
  fprintf(stderr, "  -gc            compute GC repeat content\n");
  fprintf(stderr, "  -at            compute AT repeat content\n");
}



bool
microsatelliteParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeMicroSatellite)
    return(false);

  if (inputs.size() == 0)
    sprintf(errors, "ERROR:  No input sequence files supplied.\n");

  if ((report_ga == false) &&
      (report_gc == false) &&
      (report_at == false))
    sprintf(errors, "ERROR:  No microsatellite pattern supplied.\n");

  return(errors.size() > 0);
}



void
computeMicroSat(dnaSeq    &seq,
                char const f1, char const f2, FILE *fwdBfile, FILE *fwdWfile,
                char const r1, char const r2, FILE *revBfile, FILE *revWfile,
                uint32     window) {

  bool    hasf1 = false, hasf2 = false;
  bool    hasr1 = false, hasr2 = false;

  //  Allocate and clear output.
  //
  //  We save counts for both the forward and reverse patterns.

  char const *bases    = seq.bases();
  uint64      basesLen = seq.length();

  uint32      wcMax    = basesLen / window + 1;
  uint32     *wcf      = new uint32 [wcMax];
  uint32     *wcr      = new uint32 [wcMax];

  for (uint32 ii = 0; ii < wcMax; ii++) {
    wcf[ii] = 0;
    wcr[ii] = 0;
  }

  //  Some state variables.

  uint32   inf = 0;   //  In a forward run, and which letters we've seen.
  uint32   inr = 0;   //  In a revcomp run, and which letters we've seen.

  //  Process.
  //
  //  Iterate over all the bases PLUS ONE ADDITIONAL letter (expected to be a
  //  NUL byte) to terminate whatever runs we have open.

  uint64   begf = uint64max,  endf = 0;
  uint64   begr = uint64max,  endr = 0;

  assert(bases[basesLen] == 0);

  for (uint32 ii=0; ii <= basesLen; ii++) {
    char  bp = 0;

    if (ii < basesLen)                  //  Grab sequence if not past the end; if we
      bp = bases[ii];                   //  are past, do one more loop to close runs.

    //  If a match to the forward letters:

    if      (bp == f1) {                //  If we find the first letter, set the
      begf  = (inf == 0) ? ii : begf;   //  the begin position if this is the
      endf  =              ii;          //  first letter of the run, update the end
      inf  |= 0x01;                     //  position to the current location, and
    }                                   //  remember we've seen the first letter.

    else if (bp == f2) {
      begf  = (inf == 0) ? ii : begf;
      endf  =              ii;
      inf   |= 0x02;
    }

    else {                                     //  If no match to either the first
      if (inf == 0x03) {                       //  or second letter and we just exited
        for (uint32 jj=begf; jj<=endf; jj++)   //  a run, increment window counts.
          wcf[ jj/window ]++;                  //
        if (fwdBfile)
          fprintf(fwdBfile, "%s\t%lu\t%lu\n", seq.ident(), begf, endf+1);
        assert(endf / window < wcMax);
      }
      inf = 0;                                 //  But always reset the 'in a run' flag
    }                                          //  since we're no longer in a run.

    //  If a match to the reverse letters:

    if      (bp == r1) {
      begr  = (inr == 0) ? ii : begr;
      endr  =              ii;
      inr  |= 0x01;
    }

    else if (bp == r2) {
      begr  = (inr == 0) ? ii : begr;
      endr  =              ii;
      inr   |= 0x02;
    }

    else {
      if (inr == 0x03) {
        for (uint32 jj=begr; jj<=endr; jj++)
          wcr[ jj/window ]++;
        if (revBfile)
          fprintf(revBfile, "%s\t%lu\t%lu\n", seq.ident(), begr, endr+1);
        assert(endr / window < wcMax);
      }
      inr = 0;
    }
  }

  //  Done scanning the sequence.  Output results.  If the file isn't opened,
  //  don't output.

  if (fwdWfile) {
    for (uint32 ww=0; ww<wcMax; ww++)
      fprintf(fwdWfile, "%s\t%u\t%u\t%u\t%.2f\n",
              seq.ident(), 
              ww * window, ww * window + window,
              wcf[ww],
              wcf[ww] * 100.0 / window);
  }

  if (revWfile) {
    for (uint32 ww=0; ww<wcMax; ww++)
      fprintf(revWfile, "%s\t%u\t%u\t%u\t%.2f\n",
              seq.ident(), 
              ww * window, ww * window + window,
              wcr[ww],
              wcr[ww] * 100.0 / window);
  }
}



void
computeMicroSat(dnaSeqFile *seqFile,
                char const *outPrefix,
                char const  f1, char const  f2,
                char const  r1, char const  r2,
                uint32      window) {
  char    fwdBname[FILENAME_MAX+1], revBname[FILENAME_MAX+1];
  char    fwdWname[FILENAME_MAX+1], revWname[FILENAME_MAX+1];

  //FILE   *fwdBfile = nullptr,  *fwdWfile = nullptr;
  //FILE   *revBfile = nullptr,  *revWfile = nullptr;

  sprintf(fwdBname, "%s.%c%c.bed",    outPrefix, f1, f2);
  sprintf(fwdWname, "%s.%c%c.%u.bed", outPrefix, f1, f2, window);

  sprintf(revBname, "%s.%c%c.bed",    outPrefix, r1, r2);
  sprintf(revWname, "%s.%c%c.%u.bed", outPrefix, r1, r2, window);

  //  Open files.  If the fwd and rev letters are the same, only
  //  open one set of files.

  FILE *fwdBfile =              AS_UTL_openOutputFile(fwdBname);
  FILE *fwdWfile =              AS_UTL_openOutputFile(fwdWname);

  FILE *revBfile = (f1 != f2) ? AS_UTL_openOutputFile(revBname) : nullptr;
  FILE *revWfile = (r1 != r2) ? AS_UTL_openOutputFile(revWname) : nullptr;

  for (dnaSeq seq; (seqFile->loadSequence(seq) == true); )
    computeMicroSat(seq,
                    f1, f2, fwdBfile, fwdWfile,
                    r1, r2, revBfile, revWfile,
                    window);

  AS_UTL_closeFile(fwdBfile);
  AS_UTL_closeFile(fwdWfile);

  AS_UTL_closeFile(revBfile);
  AS_UTL_closeFile(revWfile);
}




void
doMicroSatellite(vector<char const *>     &inputs,
                 microsatelliteParameters &msPar) {

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf = new dnaSeqFile(inputs[ff]);

    if (msPar.report_legacy == false) {
      if      (msPar.report_ga == true)
        computeMicroSat(sf, msPar.outPrefix, 'G', 'A', 'T', 'C', msPar.window);
      else if (msPar.report_gc == true)
        computeMicroSat(sf, msPar.outPrefix, 'G', 'C', 'G', 'C', msPar.window);
      else if (msPar.report_at == true)
        computeMicroSat(sf, msPar.outPrefix, 'A', 'T', 'A', 'T', msPar.window);
    }

    else {
      if      (msPar.report_ga == true)
        outGA(sf, msPar.outPrefix, msPar.window);
      else if (msPar.report_gc == true)
        outGC(sf, msPar.outPrefix, msPar.window);
      else if (msPar.report_at == true)
        outAT(sf, msPar.outPrefix, msPar.window);
    }

    delete sf;
  }
}





#if 0
int
main(int argc, char **argv) {
  char   *seqName     = nullptr;
  char   *outPrefix   = nullptr;
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

  if (seqName == nullptr)
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
  dnaSeqFile  *seqFile = nullptr;
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
#endif
