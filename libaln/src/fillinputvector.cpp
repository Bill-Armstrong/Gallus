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

// file  fillinputvector.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern float* h_B;
extern "C"
void execGPU(void);

void ALNAPI FillInputVector(const ALN* pALN,
                            float* afltX, 
                            long nSample,
                            int nStart,
                            ALNDATAINFO* pDataInfo,
                            const ALNCALLBACKINFO* pCallbackInfo)
{
  ASSERT(afltX);
  ASSERT(pALN);

  int nDim = pALN->nDim;  // register statement removed
  int nCols = pDataInfo->nTRcols; // This is 2 * nDim + 1
  const float* afltTRdata = pDataInfo->afltTRdata;
  const VARINFO* aVarInfo = pDataInfo->aVarInfo; // Not used

  // fill input vector
  VECTORINFO vectorinfo;        
  vectorinfo.bNeedData = FALSE;
  if (afltTRdata != NULL)
  {
    memcpy(afltX, afltTRdata + (nSample * nCols), nDim * sizeof(float));
	memcpy(h_B + 1, afltX, (nDim - 1) * sizeof(float)); // Just copy the domain values shifted by one float to put in a 1 in the 0th place to match with afltW[0]
	h_B[0] = 1.0; // This is the part of h_B which adds in the bias term during scalar products
	// We can compute the ALN value and the active LFN after using the GPU to compute scalar products
	execGPU();
  }
   // send vector info message
	if (pCallbackInfo && CanCallback(AN_VECTORINFO, pCallbackInfo->pfnNotifyProc,
									pCallbackInfo->nNotifyMask))
	{
	vectorinfo.nSample = nSample + nStart; // nStart based
	vectorinfo.aVarInfo = aVarInfo;
	vectorinfo.afltX = afltX;
	Callback(pALN, AN_VECTORINFO, &vectorinfo, pCallbackInfo->pfnNotifyProc,
				pCallbackInfo->pvData);
	
	}
}