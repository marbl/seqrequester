
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

#include "runtime.H"

#include "strings.H"
#include "files.H"

#include <vector>
#include <set>
#include <algorithm>

using namespace std;


#define SEQREQUESTER_INCLUDE_H

#include "extract.H"
#include "generate.H"
#include "microsatellite.H"
#include "mutate.H"
#include "sample.H"
#include "shiftregister.H"
#include "simulate.H"
#include "summarize.H"

#undef  SEQREQUESTER_INCLUDE_H


void doSummarize     (vector<char *> &inputs, summarizeParameters       &sumPar);
void doExtract       (vector<char *> &inputs, extractParameters         &extPar);
void doGenerate      (                        generateParameters        &genPar);
void doSimulate      (vector<char *> &inputs, simulateParameters        &simPar);
void doSample        (vector<char *> &inputs, sampleParameters          &samPar);
void doShiftRegister (                        shiftRegisterParameters   &srPar);
void doMicroSatellite(vector<char *> &inputs, microsatelliteParameters  &msPar);
void doMutate        (vector<char *> &inputs, mutateParameters          &mutPar);
