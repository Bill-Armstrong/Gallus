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

// evaltree.cpp
// eval support routines


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifdef _DEBUG
static void DebugValidateEvalTreeInfo(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      float* afltResult,
                                      int* pnStart, int* pnEnd,
                                      ALNNODE** apActiveLFNs,
                                      float* afltInput,
                                      float* afltOutput);
#endif

// evaluation of ALN on data

int ALNAPI EvalTree(const ALNNODE* pNode, 
                    const ALN* pALN,
                    ALNDATAINFO* pDataInfo,
                    const ALNCALLBACKINFO* pCallbackInfo,
                    float* afltResult,
                    int* pnStart, int* pnEnd,
                    BOOL bErrorResults /*= FALSE*/,
                    ALNNODE** apActiveLFNs /*= NULL*/,
                    float* afltInput /*= NULL*/,
                    float* afltOutput /*= NULL*/)
{
  ASSERT(pNode);
#ifdef _DEBUG
  DebugValidateEvalTreeInfo(pALN, pDataInfo, pCallbackInfo,
                            afltResult, pnStart, pnEnd, 
                            apActiveLFNs, afltInput, afltOutput);
#endif
  
  int nDim = pALN->nDim;
  int nDimt2p1 = 2 * nDim + 1;
  int nDimt2p1ti;
  long nTRcurrSamples = pDataInfo->nTRcurrSamples;

  // calc start and end points
  long nStart, nEnd; 
  //CalcDataEndPoints(nStart, nEnd, pALN, pDataInfo);
  nStart = 0;
  nEnd = nTRcurrSamples - 1;

  if (pnStart != NULL)
    *pnStart = nStart;
  if (pnEnd != NULL)
    *pnEnd = nEnd;

  // evaluation loop
  int nReturn = ALN_NOERROR;        // assume OK
  float* afltX = NULL;             // eval vector
  ALNNODE* pTree = pALN->pTree;		  // on stack for quicker access
  CCutoffInfo* aCutoffInfo = NULL;  
 
	try
 	{
    // reset active lfn array
    if (apActiveLFNs != NULL)
    {
      memset(apActiveLFNs, 0, pDataInfo->nTRcurrSamples * sizeof(ALNNODE*));
    }


  	// allocate input vector     
   	afltX = new float[nDim];   
    if (!afltX) ThrowALNMemoryException();
    memset(afltX, 0, sizeof(float) * nDim);

    // main loop
    ALNNODE* pActiveLFN = NULL;
    for (int i = nStart; i <= nEnd; i++)
    {
		nDimt2p1ti = nDimt2p1 * i;
      // fill input vector
      FillInputVector(pALN, afltX, i - nStart, nStart, pDataInfo, pCallbackInfo);

      // copy input vector?
      if (afltInput)
      {
        // get the input row
        float* afltRow = afltInput + nDimt2p1ti; // Changed this to give the correct buffer width.

        // copy values
        memcpy(afltRow, afltX, pALN->nDim * sizeof(float));
        
        // set the bias value in the output var spot
        afltRow[pALN->nOutput] = 1.0; // should we change to - 1.0 from 1.0 ?
      }

      // copy desired output 
      // ... do this before setting output value in input vector to zero below
      if (afltOutput)
      {
        afltOutput[i] = afltX[pALN->nOutput];
      }

      // CutoffEval returns distance from surface to point in the direction of
	    // the output variable, so we need to add that to the existing output value
      // to get the actual surface value
      if (!bErrorResults)
      {
        afltX[pALN->nOutput] = 0; // set output value to zero...

        // ... since output value is zero, the distance CutoffEval returns
        // is the value of the function surface
      }

      // get the distance from the point to the surface defined by the ALN
      afltResult[i] = CutoffEval(pTree, pALN, afltX, CEvalCutoff(),
                                 &pActiveLFN);
      
      // save the active LFN
      if (apActiveLFNs != NULL)
      {
        apActiveLFNs[i] = pActiveLFN;
      }
    }
  }
  catch (CALNUserException* e)
  {
  	nReturn = ALN_USERABORT;
    e->Delete();
  }
  catch (CALNMemoryException* e)
  {
  	nReturn = ALN_OUTOFMEM;
    e->Delete();
  }
  catch (CALNException* e)
  {
  	nReturn = ALN_GENERIC;
    e->Delete();
  }
  catch (...)
  {
  	nReturn = ALN_GENERIC;
  }

  // clear memory	
	delete[] afltX;
  return nReturn;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateEvalTreeInfo(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      float* afltResult,
                                      int* pnStart, int* pnEnd,
                                      ALNNODE** apActiveLFNs,
                                      float* afltInput,
                                      float* afltOutput)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  
  ASSERT(afltResult != NULL);
}
#endif
