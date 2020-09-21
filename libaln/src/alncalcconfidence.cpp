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

// alncalcconfidence.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static int ALNAPI ValidateALNCalcConfidence(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      ALNCONFIDENCE* pConfidence);

#ifdef _DEBUG
static void DebugValidateALNCalcConfidence(const ALN* pALN,
                                     ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     ALNCONFIDENCE* pConfidence);
#endif

static int __cdecl CompareErrors(const void* pElem1, const void* pElem2)
{
  float flt1 = *(float*)pElem1;
  float flt2 = *(float*)pElem2;

  if (flt1 < flt2)
    return -1;
  else if (flt1 > flt2)
    return 1;

  return 0;
}

// set ALN confidence intervals based on a data set
// much of this is derived from theory documented in Master95 p302-323
// and Press et al p228-229

ALNIMP int ALNAPI ALNCalcConfidence(const ALN* pALN,
                                    ALNDATAINFO* pDataInfo,
                                    const ALNCALLBACKINFO* pCallbackInfo,
                                    ALNCONFIDENCE* pConfidence)
{
	int nReturn = ValidateALNCalcConfidence(pALN, pDataInfo, pCallbackInfo, pConfidence);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  #ifdef _DEBUG
    DebugValidateALNCalcConfidence(pALN, pDataInfo, pCallbackInfo, pConfidence);
  #endif

  // result array
  float* afltResult = NULL;

  try
  {
    // see how many points there are
    long nTRcurrSamples = pDataInfo->nTRcurrSamples;

    // allocate results array
    afltResult = new float[nTRcurrSamples];
    
    // evaluate on data
    int nStart, nEnd;
    nReturn = EvalTree(pALN->pTree, pALN, pDataInfo, pCallbackInfo, 
                       afltResult, &nStart, &nEnd, TRUE);
    if (nReturn != ALN_NOERROR)
    {
      ThrowALNException();
    }

    // set the number of errors and error vector
    int nErr = (nEnd - nStart + 1);
    float* afltErr = afltResult + nStart;
    if (nErr <= 0)
    {
      ThrowALNException();  // bad error count
    }

    // EvalTree returned the errors in afltResult... now we sort them!
    qsort(afltErr, nErr, sizeof(float), CompareErrors);

    // calculate upper an lower bound indexes by discarding np-1 from each end
    // (conservative approach.. see Masters95 p305)
    int nDiscard = (int)floor((float)nErr * pConfidence->fltP - 1);
    
    // if nErr * pConfidence->fltP is less than 1, then nDiscard will be less than 0
    if (nDiscard < 0)
      nDiscard = 0;

    ASSERT(nDiscard >= 0 && nDiscard < nErr / 2);
    
    // set number of samples used 
    pConfidence->nSamples = nErr;

    // set lower bound
    int nLower = nDiscard;
    ASSERT(nLower >= 0 && nLower < nErr);
    pConfidence->fltLowerBound = afltErr[nLower];
    
    // set upper bound
    int nUpper = nErr - nDiscard - 1;
    ASSERT(nUpper >= 0 && nUpper < nErr);
    pConfidence->fltUpperBound = afltErr[nUpper];
  }
  catch (CALNUserException* e)	  // user abort exception
	{
		nReturn = ALN_USERABORT;
    e->Delete();
	}
	catch (CALNMemoryException* e)	// memory specific exceptions
	{
		nReturn = ALN_OUTOFMEM;
    e->Delete();
	}
	catch (CALNException* e)	      // anything other exception we recognize
	{
		nReturn = ALN_GENERIC;
    e->Delete();
	}
	catch (...)		                  // anything else, including FP errs
	{
		nReturn = ALN_GENERIC;
	}

  // deallocate mem
  delete[] afltResult;
  
	return nReturn;
}


// validate params
static int ALNAPI ValidateALNCalcConfidence(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      ALNCONFIDENCE* pConfidence)
{
  int nReturn = ValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  if (pConfidence->fltP <= 0.0 || pConfidence->fltP >= 0.5)
    return ALN_GENERIC;

  return ALN_NOERROR;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateALNCalcConfidence(const ALN* pALN,
                                     ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     ALNCONFIDENCE* pConfidence)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  ASSERT(pConfidence->fltP > 0.0 && pConfidence->fltP < 0.5);
}
#endif
