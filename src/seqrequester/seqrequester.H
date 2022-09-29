
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

#include <vector>
#include <set>
#include <algorithm>

#include "types.H"
#include "strings.H"
#include "arrays.H"
#include "files.H"
#include "math.H"
#include "system.H"
#include "sequence.H"

using namespace std;

using merylutil::stringList;

using merylutil::mtRandom;
using merylutil::sampledDistribution;

using merylutil::dnaSeq;
using merylutil::dnaSeqFile;

enum opMode {
  modeSummarize,       //  Summarize sequences in FASTA or FASTQ inputs
  modeExtract,         //  Extract sequences or subsequences from FASTA or FASTQ inputs
  modeGenerate,        //  Generate random sequences
  modeSimulate,        //  Simulate reads from FASTA or FASTQ inputs (ontigs/scaffolds/chromosomes)
  modePartition,       //  Rewrite inputs into fixed size outputs
  modeSample,          //  Extract random sequences from FASTA or FASTQ inputs
  modeShift,           //  Generate sequence based on a shift register
  modeMicroSatellite,  //  Find microsatellite repeats
  modeMutate,          //  Randomly mutate bases
  modeUnset            //  Cause an error
};

#define SEQREQUESTER_INCLUDE_H

#include "extract.H"
#include "generate.H"
#include "microsatellite.H"
#include "mutate.H"
#include "partition.H"
#include "sample.H"
#include "shiftregister.H"
#include "simulate.H"
#include "summarize.H"

#undef  SEQREQUESTER_INCLUDE_H

void doSummarize     (vector<char const *> &inputs, summarizeParameters       &sumPar);
void doExtract       (vector<char const *> &inputs, extractParameters         &extPar);
void doGenerate      (                              generateParameters        &genPar);
void doSimulate      (vector<char const *> &inputs, simulateParameters        &simPar);
void doPartition     (vector<char const *> &inputs, partitionParameters       &parPar);
void doSample        (vector<char const *> &inputs, sampleParameters          &samPar);
void doShiftRegister (                              shiftRegisterParameters   &srPar);
void doMicroSatellite(vector<char const *> &inputs, microsatelliteParameters  &msPar);
void doMutate        (vector<char const *> &inputs, mutateParameters          &mutPar);
