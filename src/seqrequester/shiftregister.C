
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

#include "bits.H"

bool
shiftRegisterParameters::parseOption(opMode &mode, int32 &arg, int32 argc, char **argv) {

  if (strcmp(argv[arg], "shift") == 0) {
    mode = modeShift;
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-order") == 0)) {
    order  = strtouint32(argv[++arg]);
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-weight") == 0)) {
    weight  = strtouint32(argv[++arg]);
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-emit") == 0)) {
    search = false;
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-length") == 0)) {
    length = strtouint64(argv[++arg]);
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-search") == 0)) {
    search = true;
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-report") == 0)) {
    report = strtodouble(argv[++arg]);
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-fast") == 0)) {
    fast   = true;
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-state") == 0)) {   //  Initial sequence
    strcpy(sr, argv[++arg]);                                              //  ACGTGGTAA
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-tap") == 0)) {     //  SR control bits
    strcpy(svmin, argv[++arg]);                                           //  011010011
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-tapmin") == 0)) {  //  SR control bits
    strcpy(svmin, argv[++arg]);                                           //  011010011
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "-tapmax") == 0)) {  //  SR control bits
    strcpy(svmax, argv[++arg]);                                           //  011010011
  }

  else if ((mode == modeShift) && (strcmp(argv[arg], "") == 0)) {
  }

  else {
    return(false);
  }

  return(true);
}



void
shiftRegisterParameters::showUsage(opMode mode) {

  if (mode != modeShift)
    return;

  fprintf(stderr, "No help for mode shift.\n");
}



bool
shiftRegisterParameters::checkOptions(opMode mode, vector<char const *> &inputs, vector<char const *> &errors) {

  if (mode != modeShift)
    return(false);

  uint32 srLen = strlen(sr);
  uint32 snLen = strlen(svmin);
  uint32 sxLen = strlen(svmax);
  bool   fail  = false;

  fail |=  (order == 0);

  fail |= ((srLen > 0) && (snLen > 0) && (srLen != snLen));
  fail |= ((srLen > 0) && (sxLen > 0) && (srLen != sxLen));
  fail |= ((snLen > 0) && (sxLen > 0) && (snLen != sxLen));

  fail |= ((srLen > 0) && (srLen != order));
  fail |= ((snLen > 0) && (snLen != order));
  fail |= ((sxLen > 0) && (snLen != order));

  if (fail == true) {
    sprintf(errors, "ERROR: order %u\n", order);
    sprintf(errors, "ERROR: sr    %s len %u\n", sr,    srLen);
    sprintf(errors, "ERROR: svmin %s len %u\n", svmin, snLen);
    sprintf(errors, "ERROR: svmax %s len %u\n", svmax, sxLen);
  }
  assert(fail == false);

  if (srLen > 0)   std::reverse(sr,    sr    + srLen);
  if (snLen > 0)   std::reverse(svmin, svmin + snLen);
  if (sxLen > 0)   std::reverse(svmax, svmax + sxLen);

  return(errors.size() > 0);
}



uint64
shiftRegisterParameters::getEncodedSR(void) {
  uint64  r = 0llu;

  if (sr[0] == 0) {
    r = 1llu;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= sr[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getCycleLen(void) {
  return(0llu);
}



uint64
shiftRegisterParameters::getCycleMax(void) {
  return(1llu << (2 * order));
}



uint64
shiftRegisterParameters::getEncodedSVmin(void) {
  uint64  r = 0llu;

  if (svmin[0] == 0) {
    r   = 1llu;
    r <<= 2 * order - 2;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= svmin[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getEncodedSVmax(void) {
  uint64  r = 0llu;

  if (svmax[0] == 0) {
    r   = 1llu;
    r <<= 2 * order;
    r  -= 1;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= svmax[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getEncodedSVmask(void) {
  return((1llu << (2 * order)) - 1);
}




void  searchShiftRegisterFast(shiftRegisterParameters &srPar);
void  searchShiftRegisterSlow(shiftRegisterParameters &srPar);
void  emitShiftRegisterFast(shiftRegisterParameters &srPar);

void  emitShiftRegisterSlow(shiftRegisterParameters &srPar) {
}



void
doShiftRegister(shiftRegisterParameters &srPar) {

  fprintf(stderr, "VERSION 7\n");

  if      ((srPar.search == true) && (srPar.fast == true))
    searchShiftRegisterFast(srPar);

  else if ((srPar.search == true) && (srPar.fast == false))
    searchShiftRegisterSlow(srPar);

  else if ((srPar.search == false) && (srPar.fast == true))
    emitShiftRegisterFast(srPar);

  else if ((srPar.search == false) && (srPar.fast == false))
    emitShiftRegisterSlow(srPar);
}

