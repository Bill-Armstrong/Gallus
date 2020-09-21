// ALN Library
// Copyright (C) 2018 William W. Armstrong.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// Version 3 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// 
// For further information contact 
// William W. Armstrong
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4

// alnrand.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
/*
///////////////////////////////////////////////////////////////////////////////
// faster rand()

//static int g_nFastRandSeed = 0x38549391;		// initial value

inline unsigned long DoFastRand()
{
  // Using the current seed, generate a new random value and seed and
	// return it.
  
  // This algorithm is similar to one proposed in Numerical Recipes in C,
  // which is in turn based on a Donald Knuth algorithm.
  
  return (g_nFastRandSeed = 1664525L * g_nFastRandSeed + 1013904223L);
}

// seeding for ALN internal pseudo-random number generator
ALNIMP void ALNAPI ALNSRand(unsigned int nSeed)
{
  g_nFastRandSeed = nSeed;
}

ALNIMP unsigned long ALNAPI ALNRand()
{
	return DoFastRand();
}

ALNIMP float ALNAPI ALNRandFloat() 
{
  
}
*/

#include <random>

unsigned long DoFastRand();
ALNIMP unsigned long ALNAPI ALNRand();
ALNIMP float ALNAPI ALNRandFloat();
ALNIMP void ALNAPI ALNSRand(unsigned int nSeed);

std::mt19937_64 generator;
std::uniform_int_distribution<int> distribution1(0, 1474836UL);
std::uniform_real_distribution<float> distribution2(0.0,1.0);


unsigned long DoFastRand()
{
    return distribution1(generator);
}
  
ALNIMP unsigned long ALNAPI ALNRand()
{
	return DoFastRand(); // this does a random int
}

ALNIMP float ALNAPI ALNRandFloat() 
{
   return distribution2(generator);
}


ALNIMP void ALNAPI ALNSRand(unsigned int nSeed)
{
	return; // not used here
}
