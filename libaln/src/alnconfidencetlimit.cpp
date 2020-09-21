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
// alnconfidencetlimit.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include "boost\math\special_functions\beta.hpp"
using namespace boost::math;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static int ALNAPI ValidateALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
                                              float fltInterval,
                                              float* pfltTLimit);

#ifdef _DEBUG
static void DebugValidateALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
                                             float fltInterval,
                                             float* pfltTLimit);
#endif

// much of this is derived from theory documented in Master95 p302-323
// and Press et al p228-229

ALNIMP int ALNAPI ALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
                                      float fltInterval,
                                      float* pfltTLimit)
{
	int nReturn = ValidateALNConfidenceTLimit(pConfidence, fltInterval, pfltTLimit);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  #ifdef _DEBUG
    DebugValidateALNConfidenceTLimit(pConfidence, fltInterval, pfltTLimit);
  #endif

  // calc number of samples in tail
  int nTailSamples = (int)floor((float)pConfidence->nSamples * pConfidence->fltP - 1);
  
  // calculate probablity of exceeding desired interval
  *pfltTLimit = (float)1.0 - ibeta((float)(pConfidence->nSamples - 2 * nTailSamples + 1), // need incomp beta fn
                             (float)(2 * nTailSamples), 
                             fltInterval);

	return nReturn;
}


// validate params
static int ALNAPI ValidateALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
                                              float fltInterval,
                                              float* pfltTLimit)
{
  if (pConfidence == NULL || pfltTLimit == NULL)
    return ALN_GENERIC;

  return ALN_NOERROR;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateALNConfidenceTLimit(const ALNCONFIDENCE* pConfidence,
                                             float fltInterval,
                                             float* pfltTLimit)
{
  ASSERT(fltInterval >= 0.0 && fltInterval <= 1.0);
}
#endif
