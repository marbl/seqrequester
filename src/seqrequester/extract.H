
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

#ifndef SEQREQUESTER_INCLUDE_H
#error  Include extract.H only from seqrequester.H
#endif


class extractParameters {
public:
  extractParameters() {
    asReverse      = false;
    asComplement   = false;

    asUpperCase    = false;
    asLowerCase    = false;

    asCompressed   = false;

    doMasking      = false;

    maskWithN      = true;

    outputFASTA    = false;
    outputFASTQ    = false;
    outputQV       = 20;
  };

  ~extractParameters() {
  };


  void      finalize(void) {

    //  If no base range specified, output all bases.

    if (baseBgn.size() == 0) {
      baseBgn.push_back(0);
      baseEnd.push_back(UINT64_MAX);
    }

    //  If no sequence range specified, output all sequences.

    if (seqsBgn.size() == 0) {
      seqsBgn.push_back(1);
      seqsEnd.push_back(UINT64_MAX);
    }

    //  If no length restriction, output all lengths.

    if (lensBgn.size() == 0) {
      lensBgn.push_back(0);
      lensEnd.push_back(UINT64_MAX);
    }

    //  Check and adjust the sequence ranges.
    //
    //  To the user, sequences begin at ONE, not ZERO.
    //  To us, sequences begin at zero.

    for (uint32 si=0; si<seqsBgn.size(); si++) {
      if (seqsBgn[si] == 0) {
        fprintf(stderr, "ERROR: sequences begin at 1, not zero.\n");
        exit(1);
      }

      seqsBgn[si] -= 1;
    }

    //  Check and adjust the base ranges.  These are space based.  A quirk in the
    //  command line parsing results in bgn == end if a single number is supplied;
    //  we interpret that to mean 'output the base at space N'.

    for (uint32 bi=0; bi<baseBgn.size(); bi++) {
      if (baseBgn[bi] == baseEnd[bi])
        baseEnd[bi] += 1;

      if (baseEnd[bi] <= baseBgn[bi]) {
        fprintf(stderr, "ERROR: base range %lu-%lu is invalid, must be increasing.\n",
                baseBgn[bi], baseEnd[bi]);
        exit(1);
      }
    }
  };


  vector<uint64>  baseBgn;    //  Base ranges to print
  vector<uint64>  baseEnd;    //

  vector<uint64>  seqsBgn;    //  Sequence ranges to print
  vector<uint64>  seqsEnd;    //

  vector<uint64>  lensBgn;    //  Length ranges to print
  vector<uint64>  lensEnd;    //

  bool          asReverse;
  bool          asComplement;

  bool          asUpperCase;
  bool          asLowerCase;

  bool          asCompressed;

  bool          doMasking;    //  Mask out any base not in baseBgn/baseEnd with 'N'

  bool          maskWithN;    //  Mask with lowercase sequence instead of 'N'

  bool          outputFASTA;
  bool          outputFASTQ;
  uint8         outputQV;
};
