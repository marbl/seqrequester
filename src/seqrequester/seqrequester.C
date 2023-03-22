
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

bool
isInput(int32 arg, int32 argc, char **argv, std::vector<char const *> &inputs) {

  if ((strcmp(argv[arg], "-") == 0) ||
      (merylutil::fileExists(argv[arg]) == true)) {
    inputs.push_back(argv[arg]);
    return(true);
  }

  return(false);
}

int
main(int argc, char **argv) {
  opMode                      mode = modeUnset;

  extractParameters           extPar;
  generateParameters          genPar;
  microsatelliteParameters    msPar;
  mutateParameters            mutPar;
  partitionParameters         parPar;
  sampleParameters            samPar;
  shiftRegisterParameters     srPar;
  simulateParameters          simPar;
  summarizeParameters         sumPar;

  std::vector<char const *>   inputs;
  std::vector<char const *>   errors;

  argc = AS_configure(argc, argv);

  //  parseOption() returns true if the option (and potentially some words
  //  following it) are consumed.

  for (int arg=1; arg < argc; arg++)
    if ((extPar.parseOption(mode, arg, argc, argv) == false) &&
        (genPar.parseOption(mode, arg, argc, argv) == false) &&
        ( msPar.parseOption(mode, arg, argc, argv) == false) &&
        (mutPar.parseOption(mode, arg, argc, argv) == false) &&
        (parPar.parseOption(mode, arg, argc, argv) == false) &&
        (samPar.parseOption(mode, arg, argc, argv) == false) &&
        ( srPar.parseOption(mode, arg, argc, argv) == false) &&
        (simPar.parseOption(mode, arg, argc, argv) == false) &&
        (sumPar.parseOption(mode, arg, argc, argv) == false) &&
        (isInput(arg, argc, argv, inputs)          == false))
      sprintf(errors, "ERROR:  Unknown parameter '%s'\n", argv[arg]);

  if (mode == modeUnset)
    sprintf(errors, "ERROR:  No mode (summarize, extract, generate, et cetera) specified.\n");

  //  Check for required options.  checkOptions() returns true if any of the
  //  supplied options generate errors; false if there are no errors or
  //  the mode is not active.

  if ((extPar.checkOptions(mode, inputs, errors)          == true) ||
      (genPar.checkOptions(mode, inputs, errors)          == true) ||
      ( msPar.checkOptions(mode, inputs, errors)          == true) ||
      (mutPar.checkOptions(mode, inputs, errors)          == true) ||
      (parPar.checkOptions(mode, inputs, errors)          == true) ||
      (samPar.checkOptions(mode, inputs, errors)          == true) ||
      ( srPar.checkOptions(mode, inputs, errors)          == true) ||
      (simPar.checkOptions(mode, inputs, errors, argv[0]) == true) ||
      (sumPar.checkOptions(mode, inputs, errors)          == true) ||
      (mode == modeUnset))
    sprintf(errors, "ERROR:  Supplied options don't make sense.\n");

  //  If errors, report usage.

  if (errors.size() > 0) {
    fprintf(stderr, "usage: %s [mode] [options] [sequence_file ...]\n", argv[0]);
    fprintf(stderr, "\n");

    if (mode == modeUnset) {
      fprintf(stderr, "MODES:\n");
      fprintf(stderr, "  summarize      report N50, length histogram, mono-, di- and tri-nucleotide frequencies\n");
      fprintf(stderr, "  extract        extract the specified sequences\n");
      fprintf(stderr, "  sample         emit existing sequences randomly\n");
      fprintf(stderr, "  generate       generate random sequences\n");
      fprintf(stderr, "  microsatellite compute microsatellite percent per window size\n");
      fprintf(stderr, "  simulate       errors in existing sequences\n");
      fprintf(stderr, "\n");
    }

    extPar.showUsage(mode);
    genPar.showUsage(mode);
    msPar.showUsage(mode);
    mutPar.showUsage(mode);
    parPar.showUsage(mode);
    samPar.showUsage(mode);
    srPar.showUsage(mode);
    simPar.showUsage(mode);
    sumPar.showUsage(mode);

    for (auto error : errors)
      if (error)
        fputs(error, stderr);

    return(1);
  }

  //  Now run!

  switch (mode) {
    case modeExtract:        doExtract       (inputs, extPar); break;
    case modeGenerate:       doGenerate      (        genPar); break;
    case modeMicroSatellite: doMicroSatellite(inputs,  msPar); break;
    case modeMutate:         doMutate        (inputs, mutPar); break;
    case modePartition:      doPartition     (inputs, parPar); break;
    case modeSample:         doSample        (inputs, samPar); break;
    case modeShift:          doShiftRegister (         srPar); break;
    case modeSimulate:       doSimulate      (inputs, simPar); break;
    case modeSummarize:      doSummarize     (inputs, sumPar); break;
    default:
      break;
  }

  return(0);
}
