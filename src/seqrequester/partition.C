
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
#include "sequence.H"
#include "mt19937ar.H"



compressedFileWriter *
openOutput(char const *pattern, uint32 index) {
  uint32  ap = 0;

  while ((pattern[ap] != 0) &&
         (pattern[ap] != '#'))
    ap++;

  //  If no #'s and multiple copies requested, fail.  But if only one copy,
  //  just open the file and return.

  if ((pattern[ap] == 0) && (index != uint32max))
    fprintf(stderr, "ERROR: Failed to find '#' in output pattern '%s', and asked to make multiple outputs.\n", pattern), exit(1);

  if  (pattern[ap] == 0)
    return(new compressedFileWriter(pattern));

  //  We've got #'s in the string.  We want to replace the last block of 'em
  //  with digits (ap found above is the start of the first block, sigh).
  //  Search backwards.

  while (pattern[ap] != 0)                   //  Find the end of the string.
    ap++;

  while ((ap > 0) && (pattern[ap] != '#'))   //  Find the last #.
    ap--;

  while ((ap > 0) && (pattern[ap] == '#'))   //  Find the start of that run.
    ap--;

  if (pattern[ap] != '#')                    //  Handle the stupid edge case.
    ap++;

  //  We're guaranteed to have some #'s in the name.  Count 'em.

  uint32 dp = 0;

  while (pattern[ap + dp] == '#')
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

  snprintf(digs, 8, "%07u", index);

  for (uint32 ii=0; pattern[ii]; ii++) {
    if ((ii < ap) || (dp == 0))
      name[ii] = pattern[ii];
    else
      name[ii] = digs[7 - dp--];
  }

  return(new compressedFileWriter(name));
}





void
doPartition_single(vector<char *> &inputs, partitionParameters &parPar) {

  parPar.initialize();

  //  Figure out how many bases/sequences to write in each output file.
  //
  //  Since we have no idea how much is in the input, we cannot compute how
  //  many output files there will be, nor can we partition into N equal
  //  pieces.  For that, the 'sample' method needs to be used.
  //
  //  If desiredCoverage is supplied, use that to compute the number of bases
  //  per file; otherwise, assume desiredNumReads and/or desiredNumBases is
  //  set.

  if (parPar.desiredCoverage > 0.0)
    parPar.desiredNumBases = (uint64)ceil(parPar.desiredCoverage * parPar.genomeSize);

  if (parPar.desiredNumBases == 0)   parPar.desiredNumBases = uint64max;
  if (parPar.desiredNumReads == 0)   parPar.desiredNumReads = uint64max;

  //  To get a little bit of parallelism, allow a few writers that will
  //  accept sequences round-robin.

  uint64                 *nBases = new uint64 [parPar.numWriters];
  uint64                 *nSeqs  = new uint64 [parPar.numWriters];
  uint32                  nFiles = 0;
  compressedFileWriter  **oFile  =  new compressedFileWriter * [parPar.numWriters];
  uint32                  oIndex = 0;

  fprintf(stderr, "Opening %u output file%s.\n", parPar.numWriters, (parPar.numWriters == 1) ? "" : "s");

  for (uint32 ii=0; ii<parPar.numWriters; ii++) {
    oFile[ii]  = openOutput(parPar.output1, nFiles++);
    nBases[ii] = 0;
    nSeqs[ii]  = 0;

    fprintf(stderr, "  %s\n", oFile[ii]->filename());
  }

  //  Read sequences from the input, copying to each output in turn.

  dnaSeq   seq1;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1  = new dnaSeqFile(inputs[ff]);
    uint64       nseq = 0;

    while (sf1->loadSequence(seq1) == true) {
      outputSequence(oFile[oIndex]->file(),
                     seq1.ident(), seq1.bases(), seq1.quals(), seq1.length(),
                     sf1->isFASTQ(),
                     parPar.outputFASTA,
                     parPar.outputFASTQ,
                     parPar.outputQV);

      nBases[oIndex] += seq1.length();
      nSeqs[oIndex]  += 1;

      if ((++nseq % 61075) == 0)
        fprintf(stderr, "  %10lu sequences in '%s'\r", nseq, inputs[ff]);

      if ((nBases[oIndex] >= parPar.desiredNumBases) ||
          (nSeqs[oIndex]  >= parPar.desiredNumReads)) {
        fprintf(stderr, "  %s closed; %lu bases %lu reads\n",
                oFile[oIndex]->filename(), nBases[oIndex], nSeqs[oIndex]);

        delete oFile[oIndex];

        oFile[oIndex]  = openOutput(parPar.output1, nFiles++);
        nBases[oIndex] = 0;
        nSeqs[oIndex]  = 0;

        fprintf(stderr, "  %s\n", oFile[oIndex]->filename());
      }

      oIndex++;
      if (oIndex >= parPar.numWriters)
        oIndex = 0;
    }

    fprintf(stderr, "  %10lu sequences in '%s'\n", nseq, inputs[ff]);

    delete sf1;
  }

  for (uint32 ii=0; ii<parPar.numWriters; ii++)
    delete oFile[ii];

  fprintf(stderr, "Done.\n");
}



void
doPartition(vector<char *> &inputs, partitionParameters &parPar) {

  if (parPar.isPaired == false)
    doPartition_single(inputs, parPar);
  //else
  //  doPartition_paired(inputs, parPar);
}

