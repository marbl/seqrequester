
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
#include "shiftregister-gf4.H"
#include "bits.H"


void
searchShiftRegisterFast(shiftRegisterParameters &srPar) {
 //  Allocate space for the loop detector, set local variables.

  uint64  srinit = srPar.getEncodedSR();      //  The shift register state
  uint64  sr     = srPar.getEncodedSR();      //  The shift register state
  uint64  cyclen = srPar.getCycleLen();       //
  uint64  cycmax = srPar.getCycleMax();       //  The length of the maximum cycle
  uint64  svmin  = srPar.getEncodedSVmin();   //  The tap vector
  uint64  svmax  = srPar.getEncodedSVmax();   //  The first vector we don't want to examine
  uint64  svmask = srPar.getEncodedSVmask();
  uint64  gf4widemult[4];

  //  Log.

  fprintf(stderr, "Finding cycles for length %u bases.\n", srPar.order);
  fprintf(stderr, "  srinit %0*lo\n",  srPar.order, expandTo3(srinit));
  fprintf(stderr, "  cyclen %lu\n",    cyclen);
  fprintf(stderr, "  cycmax %lu\n",    cycmax);
  fprintf(stderr, "  svmin  %0*lo\n", srPar.order, expandTo3(svmin));
  fprintf(stderr, "  svmax  %0*lo\n", srPar.order, expandTo3(svmax));
  fprintf(stderr, "  svmask %0*lo\n", srPar.order, expandTo3(svmask));
  fprintf(stderr, "\n");

  //  We can precompute the result of gf4mult[sv[kk]][out] used in the addition.
  //  There are four possible 'out' values, and sv is fixed for each cycle.
  //
  //  The addition itself needs an extra bit for overflow, which
  //  does make our detect array 50% larger.

#ifdef DETECT
  bitArray  *detect = new bitArray(cycmax);
#endif

  for (uint64 sv = svmin; sv <= svmax; sv++) {

    if ((sv % 99999999) == 0)
      fprintf(stderr, "%0*lo cycle %14lu / %14lu - %6.2f%% - %0*lo\r",
              srPar.order, expandTo3(sv),
              cyclen,
              cycmax,
              100.0 * cyclen / cycmax,
              srPar.order, expandTo3(sr));


    //  Count the number of taps.  If this is different than the requested
    //  weight, skip it.


    //  EXTREMELY SLOW at high order

    if (srPar.weight > 0) {
      uint64 svwl  = (sv & 0x5555555555555555llu);   //  Low  order bits of each tap.
      uint64 svwr  = (sv & 0xaaaaaaaaaaaaaaaallu);   //  High order bits of each tap.
      uint64 svw   = (svwr >> 1) | (svwl);           //  One bit set for every positive tap.

      uint32 w = countNumberOfSetBits64(svw);

      if (w != srPar.weight) {
        //fprintf(stderr, "SKIP sv %0*lo\n", srPar.order, expandTo3(sv));
        continue;
      }
    }

    //  Build the multiplication table.
    //
    //  Including all these for small sizes (e.g., size 9) doesn't seem to
    //  result in any slow down.

    for (uint64 out=0; out<4; out++) {
      uint64 mult = 0;

      //lt |= gf4mult[(sv >> 42) & 0x03][out];  mult <<= 2;   // 22 - too big for expandTo3()
      mult |= gf4mult[(sv >> 40) & 0x03][out];  mult <<= 2;   // 21
      mult |= gf4mult[(sv >> 38) & 0x03][out];  mult <<= 2;   // 20
      mult |= gf4mult[(sv >> 36) & 0x03][out];  mult <<= 2;   // 19
      mult |= gf4mult[(sv >> 34) & 0x03][out];  mult <<= 2;   // 18
      mult |= gf4mult[(sv >> 32) & 0x03][out];  mult <<= 2;   // 17
      mult |= gf4mult[(sv >> 30) & 0x03][out];  mult <<= 2;   // 16
      mult |= gf4mult[(sv >> 28) & 0x03][out];  mult <<= 2;   // 15
      mult |= gf4mult[(sv >> 26) & 0x03][out];  mult <<= 2;   // 14
      mult |= gf4mult[(sv >> 24) & 0x03][out];  mult <<= 2;   // 13
      mult |= gf4mult[(sv >> 22) & 0x03][out];  mult <<= 2;   // 12
      mult |= gf4mult[(sv >> 20) & 0x03][out];  mult <<= 2;   // 11
      mult |= gf4mult[(sv >> 18) & 0x03][out];  mult <<= 2;   // 10
      mult |= gf4mult[(sv >> 16) & 0x03][out];  mult <<= 2;   //  9
      mult |= gf4mult[(sv >> 14) & 0x03][out];  mult <<= 2;   //  8
      mult |= gf4mult[(sv >> 12) & 0x03][out];  mult <<= 2;   //  7
      mult |= gf4mult[(sv >> 10) & 0x03][out];  mult <<= 2;   //  6
      mult |= gf4mult[(sv >>  8) & 0x03][out];  mult <<= 2;   //  5
      mult |= gf4mult[(sv >>  6) & 0x03][out];  mult <<= 2;   //  4
      mult |= gf4mult[(sv >>  4) & 0x03][out];  mult <<= 2;   //  3
      mult |= gf4mult[(sv >>  2) & 0x03][out];  mult <<= 2;   //  2
      mult |= gf4mult[(sv >>  0) & 0x03][out];                //  1

      mult &= svmask;

#ifdef DEBUG
      fprintf(stderr, "widemult[%0*lo][%02x] %0*lo\n", srPar.order, expandTo3(sv), out, srPar.order, expandTo3(mult));
#endif

      gf4widemult[out] = mult;
    }

    //  Loop until we hit a cycle.
    //
    //  Note that cyclen will be one less than the max possible, because we
    //  cannot ever get the all zero state.  This is nice for the
    //  implementation, because we can then just compare against svmax, the
    //  last possible state vector, to decide if we cycled through all
    //  possible 'svmax + 1' states (including the all zero state we can't
    //  ever get).

#ifdef DETECT
    detect->clear();                     //  Reset the cycle detector.
#endif

    sr     = srinit;
    cyclen = 1;

#ifdef DETECT
    for (cyclen=1; (detect->flipBit(sr) == false); cyclen++) {
#else
    do {
#endif
      uint64  out = sr & 0x03;           //  Save the output value
      uint64  mul = gf4widemult[out];    //  Compute the multiplier

      if ((cyclen & 0x1fffffff) == 0)
        fprintf(stderr, "%0*lo cycle %14lu / %14lu - %6.2f%% - %0*lo\r",
                srPar.order, expandTo3(sv),
                cyclen,
                cycmax,
                100.0 * cyclen / cycmax,
                srPar.order, expandTo3(sr));

#ifdef DEBUG
      fprintf(stderr, "cycle %14lu %6.2f%% out %02lx sr %0*lo\n", cyclen, 100.0 * cyclen / cycmax, out, srPar.order, expandTo3(sr >> 2));
      fprintf(stderr, "                                   add %0*lo\n", srPar.order, expandTo3(mul));
      fprintf(stderr, "                                 final %0*lo\n", srPar.order, expandTo3((sr >> 2) ^ mul));
#endif

      sr = (sr >> 2) ^ mul;
#ifdef DETECT
    }
#else
    cyclen++;
    } while ((sr != srinit) && (sr != 0));
#endif

    assert((sr == 0) || (sr == srinit));

    //  Report the cycle if it's sufficiently large.

    if (cyclen + 1 >= srPar.report * cycmax)
      fprintf(stdout, "%14lu/%14lu %7.3f%% for vector %0*lo sr %0*lo\n",
              cyclen,
              cycmax,
              100.0 * cyclen / cycmax,
              srPar.order, expandTo3(sv),
              srPar.order, expandTo3(sr));
  }

#ifdef DETECT
  delete detect;
#endif
}


