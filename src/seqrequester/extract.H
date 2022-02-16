
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
  extractParameters();
  ~extractParameters();

  void      finalize(void);

private:
  void      loadNames(void);

public:
  vector<char const *>  baseArgs;               //  Base ranges to print
  vector<uint64>        baseBgn;                //
  vector<uint64>        baseEnd;                //

  vector<char const *>  seqsArgs;               //  Sequence ranges to print
  vector<uint64>        seqsBgn;                //
  vector<uint64>        seqsEnd;                //
  stringList            seqsName;               //

  vector<char const *>  lensArgs;
  vector<uint64>        lensMin;                //  Length ranges to print
  vector<uint64>        lensMax;                //

  bool                  asReverse    = false;
  bool                  asComplement = false;

  bool                  asUpperCase  = false;
  bool                  asLowerCase  = false;

  bool                  asCompressed = false;

  bool                  doMasking    = false;   //  Mask out any base not in baseBgn/baseEnd with 'N'

  bool                  maskWithN    = true;    //  Mask with lowercase sequence instead of 'N'

  bool                  outputFASTA  = false;
  bool                  outputFASTQ  = false;
  uint8                 outputQV     = 20;

  //  Data for doing the extraction.

  char                  C[256] = {0};   //  Complement a base
  char                  U[256] = {0};   //  Uppercase a base
  char                  L[256] = {0};   //  Lowercase a base

  uint64                outputBasesLen  = 0;
  uint64                outputBasesMax  = 0;        //1048576;
  char                 *outputBases     = nullptr;  //new char [outputBasesMax];
  uint8                *outputQuals     = nullptr;  //new char [outputBasesMax];
};
