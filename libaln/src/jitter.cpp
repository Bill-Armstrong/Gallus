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

// jitter.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

inline
float ALNAPI Noise()
{
//   This generates a random value with a triangular distribution from -1 to 1
//   This by subtracting two random uniformly distributed values in (0, 1).
	return ALNRandFloat() - ALNRandFloat();
}

void ALNAPI Jitter(ALN* pALN, float* afltX)
{
  ASSERT(pALN);
  ASSERT(afltX);

  int nDim = pALN->nDim;
  int nOutput = pALN->nOutput;
  
  // save output value
  float fltOutput = afltX[nOutput];

  for (int i = 0; i < nDim; i++)
  {
#ifdef _DEBUG
    float flt = afltX[i];
#endif
		
    afltX[i] += Noise() * pALN->aRegions[0].aConstr[i].fltEpsilon;

#ifdef _DEBUG
    float fltEps = pALN->aRegions[0].aConstr[i].fltEpsilon;
    ASSERT(afltX[i] >= flt - fltEps && 
           afltX[i] <= flt + fltEps);
#endif
  }

  // restore output value
  afltX[nOutput] = fltOutput;
}