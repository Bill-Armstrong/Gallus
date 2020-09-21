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

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

int ALNAPI ValidateALNDataInfo(const ALN* pALN,
                               ALNDATAINFO* pDataInfo,
                               const ALNCALLBACKINFO* pCallbackInfo)
{
  // parameter variance
  if (pDataInfo == NULL)
  {
    return ALN_GENERIC;
  }

  if (pALN == NULL)
  {
    return ALN_GENERIC;
  }
  /*
  // must have at least one training point or not have a training set
	if ((pDataInfo->nTRcurrSamples <= 0) && (pDataInfo->fltMSEorF <= 0))
  {
    return ALN_GENERIC;
  }

  // need proc if no data
  if(pDataInfo->afltTRdata == NULL && 
     (pCallbackInfo == NULL || pCallbackInfo->pfnNotifyProc == NULL)   )
  {
    return ALN_GENERIC;
  }

  // need AN_VECTORINFO if no data
  if(pDataInfo->afltTRdata == NULL && 
     (pCallbackInfo == NULL || !(pCallbackInfo->nNotifyMask & AN_VECTORINFO)))
  {
    return ALN_GENERIC;
  }

  // make sure columns valid
  if (pDataInfo->aVarInfo == NULL && 
      pDataInfo->afltTRdata != NULL && 
      pDataInfo->nTRcols < pALN->nDim)
  {
    return ALN_GENERIC;
  }
  
  // check var info column, delta validity
  // const VARINFO* aVarInfo = pDataInfo->aVarInfo; this is not used
  */
  return ALN_NOERROR;
}

#ifdef _DEBUG
void ALNAPI DebugValidateALNDataInfo(const ALN* pALN,
                                     ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo)
{
  ASSERT(pDataInfo != NULL);

  // valid aln pointer
	ASSERT(pALN != NULL);

  // valid number of points
  ASSERT(pDataInfo->nTRcurrSamples > 0);

  // valid data cols
  ASSERT((pDataInfo->afltTRdata != NULL && pDataInfo->nTRcols > 0) || 
         pDataInfo->afltTRdata == NULL);

  // valid notify proc
  ASSERT(pDataInfo->afltTRdata != NULL || 
         (pCallbackInfo != NULL && 
          pCallbackInfo->pfnNotifyProc != NULL && 
          (pCallbackInfo->nNotifyMask & AN_VECTORINFO)));
  
  // valid varinfo
  ASSERT(pDataInfo->aVarInfo != NULL || pDataInfo->nTRcols >= (2 * pALN->nDim + 1));
  if (pDataInfo->aVarInfo != NULL)
  {
    for (int i = 0; i < pALN->nDim; i++)
    {
      // column validity
      ASSERT(pDataInfo->afltTRdata);
    }
  }
}
#endif
