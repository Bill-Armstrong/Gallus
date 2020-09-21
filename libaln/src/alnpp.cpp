// ALN Library C++ wrapper
// ALNfit Learning Engine for approximation of functions defined by samples.
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

// alnpp.cpp

#ifdef __GNUC__
#include <typeinfo>
#endif

#include <stdlib.h>
#include <memory.h>
#include <limits>
#include "alnpp.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifndef ASSERT
#define ASSERT ALNASSERT
#endif

/////////////////////////////////////////////////////////////////////////////
// class CAln

CAln::CAln()
{
  m_pALN = NULL;
  memset(&m_datainfo, 0, sizeof(m_datainfo));
  m_nLastError = ALN_GENERIC; // no ALN pointer yet!
}

CAln::~CAln()
{
  Destroy();
  ASSERT(m_pALN == NULL);
}

void ALNAPI CAln::addTRsample(float* afltX, const int nDim)
{
	float sum;
	int nDimm1 = nDim - 1;
	int nDimt2 = 2 * nDim;
	int nDimt2m1 = nDimt2 - 1;
	int nDimt2p1 = nDimt2 + 1;
	// nDimt2p1 is the number of columns in the data buffer
	// nDimt2p1 * i is the index of first domain component of a sample vector in the i-th row.
	// In order to get the index of
	// 1. the desired output value of a sample: add nDimm1,
	// 2. the first element of the difference vector domain components: add nDim,
	// 3. the difference of desired output values: add  nDimt2m1;
	// 4. the squared distance between two closest samples: add nDimt2

	ALNDATAINFO* thisDataInfo = this->GetDataInfo();
	// Put some items on the stack
	long nTRmaxSamples = thisDataInfo->nTRmaxSamples;
	long nTRcurrSamples = thisDataInfo->nTRcurrSamples;
	int nTRcols = thisDataInfo->nTRcols;
	long nTRinsert = thisDataInfo->nTRinsert;
	float fltMSEorF = thisDataInfo->fltMSEorF;
	float* afltTRdata; // This is a pointer to the buffer on the stack
	ASSERT(nTRcols == nDimt2p1); // Check

	// The new sample afltX will be placed here
	long nDimt2p1tTRinsert = nDimt2p1 * nTRinsert;

	// Initialize the data buffer afltTRdata when the first sample is added
	if (nTRcurrSamples == 0) // when adding the first sample
	{
		long bufferSize = nTRmaxSamples * nTRcols;
		//Allocate the buffer on first use 
		thisDataInfo->afltTRdata = (float*)malloc( bufferSize * sizeof(float));
		// make the stack pointer point to the allocated space
		afltTRdata = thisDataInfo->afltTRdata;
		// First we fill TRbuff with 0's and an FLT_MAX
		memset(afltTRdata, 0, bufferSize * sizeof(float));
		/*
		for (long i = 0; i < nTRmaxSamples; i++)
		{
			int nDimt2p1ti = nDimt2p1 * i;
			
			for (int j = 0; j < nDimt2; j++)
			{
				afltTRdata[nDimt2p1ti + j] = 0; // We manipulate the buffer using the local pointer
			}
		}
		*/
		afltTRdata[nDimt2] = FLT_MAX; // the value is huge to make sure a new difference vector is inserted at the second insertion
		for (int j = 0; j < nDim; j++)
		{
			afltTRdata[j] = afltX[j];
		}
		// update the buffer values in the ALN
		thisDataInfo->nTRcurrSamples = 1;
		thisDataInfo->nTRinsert++;
		return;
	}// End of initializing the data buffer

	// Now comes the processing after initialization of the buffer
	// There are two major cases depending on whether or not fltMSEorF is < 0.
	// A sample can always be added; the sample at TRinsert has to be replaced if the buffer is already at max size
	// The following is done in the case fltMSEorF < 0 when the F-test is called for. Otherwise, finding a closest sample etc. is not done.
	// Of course in that case we can't change later to use the F-test.
	// Start by comparing the new sample to all in the buffer to find
	// 1. what other samples it is closest to
	// 2. what other sample is closest to it.	(N.B. "closest" means * among * the closest if it is not unique)
	
	afltTRdata = thisDataInfo->afltTRdata; // restore the stack-based pointer.
	if (fltMSEorF < 0)
	{
		float* afltYtemp = (float*)malloc((nDim + 1) * sizeof(float)); // stores the difference vector and square distance 
					// of the sample which is currently the closest to afltX which is being inserted;
		float* afltYdiff = (float*)malloc(nDim * sizeof(float)); // stores the differences of afltX to the current buffer sample
		afltYtemp[nDim] = FLT_MAX; // This forces a new sample to initially be closest to the first when there is no valid afltYtemp vector 
		for (long i = 0; i < nTRcurrSamples; i++)
		{
			int nDimt2p1ti = nDimt2p1 * i; // this is the starting field of the index i sample in afltTRdata
			// Get the current square distance of the new sample afltX from sample i in afltTRdata (there is at least one already)
			sum = 0;
			for (int j = 0; j < nDim; j++)
			{
				afltYdiff[j] = afltX[j] - afltTRdata[nDimt2p1ti + j];
				if (j == nDim - 1) break; // use domain components only to calculate square distance
				sum += afltYdiff[j] * afltYdiff[j]; // TBD later: different axes will contribute with different weights to noise variance
			}
			// Compare the sum to the squared distance of the closest to the i-th sample
			if (sum < afltTRdata[nDimt2p1ti + nDimt2])
			{
				// If the sum is less than for the existing closest
				// replace the closest to sample i with the new differences and update the square distance
				for (int j = 0; j < nDim; j++)
				{

					afltTRdata[nDimt2p1ti + nDim + j] = afltYdiff[j];
				}
				afltTRdata[nDimt2p1ti + nDimt2] = sum;
			}
			// If the ith sample in afltTRdata is closer to afltX than anything was earlier, then update the temporary values
			if (sum < afltYtemp[nDim])
			{
				// Store sample i in the temporarily closest sample to afltX with its square distance.
				for (int j = 0; j < nDim; j++)
				{
					afltYtemp[j] = afltTRdata[nDimt2p1ti + j] - afltX[j];
				}
				afltYtemp[nDim] = sum;
			}
		}
		// Insert the new sample and its difference vector and closest distance at nTRinsert
		for (int j = 0; j < nDim; j++)
		{
			afltTRdata[nDimt2p1tTRinsert + j] = afltX[j];
			afltTRdata[nDimt2p1tTRinsert + nDim + j] = afltYtemp[j];
		}
		afltTRdata[nDimt2p1tTRinsert + nDimt2] = afltYtemp[nDim];
		free(afltYdiff);
		free(afltYtemp);
	} // end of insertion loop for F-test part
	else
	{
		// Insert the new sample, its zeroed difference vector and closest distance at nTRinsert
		// Note: the two lines involving differences would be useful if we changed from doing an F-test to not doing it during a run
		for (int j = 0; j < nDim; j++)
		{
			afltTRdata[nDimt2p1tTRinsert + j] = afltX[j];
			afltTRdata[nDimt2p1tTRinsert + nDim + j] = 0; // zero out the noise variance part (e.g.from a previous sample with fltMSEorF < 0)
		}
		afltTRdata[nDimt2p1tTRinsert + nDimt2] = 0; // Once we have decided not to use the F-test, we have to stick with that decision
	}
	// The buffer has been updated, now we see where the next insertion can be done
	if (nTRcurrSamples < nTRmaxSamples)
	{
		nTRcurrSamples++; // This will stay at the max if it gets to it
	}
	thisDataInfo->nTRcurrSamples = nTRcurrSamples;
	if(++nTRinsert == nTRmaxSamples) nTRinsert = 0; // Too TRICKY!!! See if it works!
	thisDataInfo->nTRinsert = nTRinsert; // Pass the information back to the ALN.
	
	return;
}

void ALNAPI CAln::reduceNoiseVariance()
{
	// This routine should only be used when there are many samples in afltTRdata since
	// the closest other sample to a sample should have a close value of the ideal function
	// In that case we replace all sample values by the average of two and get 1/2 the noise variance.
	ALNDATAINFO* thisDataInfo = this->GetDataInfo();
	float* afltTRdata = thisDataInfo->afltTRdata;
	//int nTRmaxSamples = thisDataInfo->nTRmaxSamples;
	long nTRcurrSamples = thisDataInfo->nTRcurrSamples;
	int nTRcols = thisDataInfo->nTRcols;
	long nTRinsert = thisDataInfo->nTRinsert;
	//float fltMSEorF = thisDataInfo->fltMSEorF;
	if (!afltTRdata) return; // this avoids a crash
	if (thisDataInfo->fltMSEorF > 0) return; // we are not using noise variance
	int nDim = (nTRcols - 1)/2;
	int nDimm1 = nDim - 1;
	int nDimt2m1 = nDim * 2 - 1;
	float value;
	for (long i = 0; i < nTRcurrSamples; i++)
	{
		value = afltTRdata[nTRcols * i + nDimm1 ];
		value -= 0.5 * afltTRdata[nTRcols * i + nDimt2m1];
		afltTRdata[nTRcols * i + nDimm1] = value;
	}
}

ALNREGION* CAln::GetRegion(int nRegion)
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  return m_pALN->aRegions + nRegion;
}

const ALNREGION* CAln::GetRegion(int nRegion) const
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  return m_pALN->aRegions + nRegion;
}

ALNCONSTRAINT* CAln::GetConstraint(int nVar, int nRegion)
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  ASSERT(nVar >= 0 && nVar < m_pALN->nDim);
  return m_pALN->aRegions[nRegion].aConstr + nVar;
}

const ALNCONSTRAINT* CAln::GetConstraint(int nVar, int nRegion) const
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  ASSERT(nVar >= 0 && nVar < m_pALN->nDim);
  return m_pALN->aRegions[nRegion].aConstr + nVar;
}

float CAln::GetEpsilon(int nVar, int nRegion) const
{
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->fltEpsilon; 
}

void CAln::SetEpsilon(float fltEpsilon, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  pConstr->fltEpsilon = fltEpsilon;
	pConstr->fltSqEpsilon = fltEpsilon * fltEpsilon;
}

float CAln::GetWeightMin(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->fltWMin; 
}

void CAln::SetWeightMin(float fltWMin, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  ASSERT(nVar != m_pALN->nOutput);  
    // can't change output var weight!
  
  if (nVar != m_pALN->nOutput)
    pConstr->fltWMin = fltWMin; 
}

float CAln::GetWeightMax(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->fltWMax; 
}

void CAln::SetWeightMax(float fltWMax, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  ASSERT(nVar != m_pALN->nOutput);  
    // can't change output var weight!
  
  if (nVar != m_pALN->nOutput)
    pConstr->fltWMax = fltWMax; 
}  

float CAln::GetMin(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->fltMin; 
}

void CAln::SetMin(float fltMin, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  
  pConstr->fltMin = fltMin; 
}

float CAln::GetMax(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->fltMax; 
}

void CAln::SetMax(float fltMax, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  
  pConstr->fltMax = fltMax; 
}  

ALNNODE* CAln::GetTree()
{
  if (m_pALN == NULL)
    return NULL;

  return m_pALN->pTree;
}

const ALNNODE* CAln::GetTree() const
{
  if (m_pALN == NULL)
    return NULL;

  return m_pALN->pTree;
}

void CAln::SetDataInfo(float* afltTRdata, const long nTRmaxSamples, long nTRcurrSamples, int nTRcols, long nTRinsert,	const float fltMSEorF)
{
  m_datainfo.afltTRdata = afltTRdata;
  m_datainfo.nTRmaxSamples = nTRmaxSamples;
  m_datainfo.nTRcurrSamples = nTRcurrSamples;
  m_datainfo.nTRcols = nTRcols;
  m_datainfo.nTRinsert = nTRinsert;
  m_datainfo.fltMSEorF = fltMSEorF;
}

BOOL CAln::Create(int nDim, int nOutput)
{
  Destroy();
  ASSERT(m_pALN == NULL);

  m_pALN = ALNCreateALN(nDim, nOutput);
  if (m_pALN != NULL)
  {
    m_nLastError = ALN_NOERROR;
  }
  else
  {
    m_nLastError = ALN_GENERIC;
  }

  return m_pALN != NULL;
}

#ifdef ENABLE_REGIONS
int CAln::AddRegion(int nParentRegion, float fltLearnFactor, 
                    int nConstr, int* anConstr)
{
  m_nLastError = ALN_NOERROR;
  return ALNAddRegion(m_pALN, nParentRegion, fltLearnFactor,
                      nConstr, anConstr);
}
#endif

BOOL CAln::AddLFNs(ALNNODE* pParent, int nParentMinMaxType, 
                   int nLFNs, ALNNODE** apLFNs /*= NULL*/)
{
  m_nLastError = ALNAddLFNs(m_pALN, pParent, nParentMinMaxType, nLFNs, apLFNs);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::AddMultiLayer(ALNNODE* pParent, int nParentMinMaxType, 
                         int nLayers, int nMinMaxFanin, int nLFNFanin, 
                         int nFlags)
{
  m_nLastError = ALNAddMultiLayer(m_pALN, pParent, nParentMinMaxType, nLayers, 
                                  nMinMaxFanin, nLFNFanin, nFlags);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::AddTreeString(ALNNODE* pParent, const char* pszTreeString, 
                         int& nParsed)
{
  m_nLastError = ALNAddTreeString(m_pALN, pParent, pszTreeString, &nParsed);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::SetGrowable(ALNNODE* pNode)
{
  m_nLastError = ALN_NOERROR;
  return ALNSetGrowable(m_pALN, pNode);
}

BOOL CAln::Destroy()
{
  m_nLastError = ALN_NOERROR;
  BOOL b = ALNDestroyALN(m_pALN);
  if (b) m_pALN = NULL;
  return b;
}

// private callback data struct
struct CALLBACKDATA
{
  CAln* pALN;
  void* pvData;
};

// training - uses internal ALNDATAINFO if pData == NULL
BOOL CAln::Train(int nMaxEpochs, float fltMinRMSErr, float fltLearnRate,
                 BOOL bJitter, int nNotifyMask /*= AN_NONE*/, 
                 ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNTrain(m_pALN, pData, &callback, nMaxEpochs, fltMinRMSErr, fltLearnRate, bJitter);

	return (m_nLastError == ALN_NOERROR || m_nLastError == ALN_USERABORT);
}

float CAln::CalcRMSError(int nNotifyMask /*= AN_NONE*/, 
                          ALNDATAINFO* pData /*= NULL*/, 
                          void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;
  
  float flt;
  m_nLastError = ::ALNCalcRMSError(m_pALN, pData, &callback, &flt);

  if (m_nLastError != ALN_NOERROR)
    flt = -1.0;

  return flt;
}

// eval
BOOL CAln::Eval(float* afltResult, int* pnStart /*=NULL*/, 
                int* pnEnd /*= NULL*/, int nNotifyMask /*= AN_NONE*/,
                ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNEval(m_pALN, pData, &callback, afltResult, pnStart, pnEnd);

  return m_nLastError == ALN_NOERROR;
}

// quick eval
float CAln::QuickEval(const float* afltX, ALNNODE** ppActiveLFN /*= NULL*/)
{
  m_nLastError = ALN_NOERROR;
  return ALNQuickEval(m_pALN, afltX, ppActiveLFN);
}

// get variable monotonicicty, returns -1 on failure
int CAln::VarMono(int nVar)
{
  int nMono;
  m_nLastError = ALNVarMono(m_pALN, nVar, &nMono);
  return nMono;
}

// invert the aln to get an aln for a different output variable WWA
BOOL CAln::Invert(int nVar)
{
	m_nLastError = ALNInvert(m_pALN, nVar);
	return (m_nLastError == ALN_NOERROR);
}

// save ALN to disk file
BOOL CAln::Write(const char* pszFileName)
{
  m_nLastError = ALNWrite(m_pALN, pszFileName);
  return m_nLastError == ALN_NOERROR;
}

// read ALN from disk file... destroys any existing ALN
BOOL CAln::Read(const char* pszFileName)
{
  Destroy();
  ASSERT(m_pALN == NULL);

  m_nLastError = ALNRead(pszFileName, &m_pALN);
  return m_nLastError == ALN_NOERROR;
}

// conversion to dtree
DTREE* CAln::ConvertDtree(int nMaxDepth)
{
  DTREE* pDtree = NULL;
	DTREE** ppDtree = &pDtree;
  m_nLastError = ALNConvertDtree(m_pALN, nMaxDepth, ppDtree);
	pDtree = *ppDtree;
  return pDtree;
}

// confidence intervals
BOOL CAln::CalcConfidence(ALNCONFIDENCE* pConfidence, 
                          int nNotifyMask /*= AN_NONE*/, 
                          ALNDATAINFO* pData /*= NULL*/, 
                          void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNCalcConfidence(m_pALN, pData, &callback, pConfidence);

  return m_nLastError == ALN_NOERROR;
}

float CAln::ConfidencePLimit(const ALNCONFIDENCE* pConfidence, 
                              float fltSignificance)
{
  float fltPLimit;
  
  if (ALNConfidencePLimit(pConfidence, fltSignificance, &fltPLimit) != ALN_NOERROR)
    return -1.0;
  
  return fltPLimit;
}

float CAln::ConfidenceTLimit(const ALNCONFIDENCE* pConfidence, 
                              float fltInterval)
{
  float fltTLimit;
  
  if (ALNConfidenceTLimit(pConfidence, fltInterval, &fltTLimit) != ALN_NOERROR)
    return -1.0;
  
  return fltTLimit;
}

// lfn analysis
BOOL CAln::LFNAnalysis(void*& pvAnalysis,
                       int& nLFNStats,
                       int nNotifyMask /*= AN_NONE*/, 
                       ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo; // In new format
	CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNLFNAnalysis(m_pALN, pData, &callback,
                                pvAnalysis, nLFNStats);

  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::LFNFreeAnalysis(void* pvAnalysis)
{
  return ALNLFNFreeAnalysis(pvAnalysis) == ALN_NOERROR;
}

BOOL CAln::LFNStats(void* pvAnalysis, 
                    int nLFNStat,
                    LFNSTATS& LFNStats,
                    int& nWeights,
                    ALNNODE*& pLFN)
{
  return ALNLFNStats(pvAnalysis, nLFNStat, &LFNStats, 
                     &nWeights, &pLFN) == ALN_NOERROR;
}

BOOL CAln::LFNWeightStats(void* pvAnalysis,
                          int nLFNStat,
                          int nWeightStat,
                          LFNWEIGHTSTATS& LFNWeightStats)
{
  return ALNLFNWeightStats(pvAnalysis, nLFNStat, nWeightStat, 
                           &LFNWeightStats) == ALN_NOERROR;
}

int ALNAPI CAln::ALNNotifyProc(const ALN* pALN, int nCode, void* pParam, 
                               void* pvData)
{
  CALLBACKDATA* pData = (CALLBACKDATA*)pvData;
  CAln* pALNObj = pData->pALN;
  
  ASSERT(pALNObj != NULL && pALN == pALNObj->m_pALN);

  BOOL bContinue = FALSE;

  switch (nCode)
	{
		case AN_TRAINSTART:
			bContinue = pALNObj->OnTrainStart((TRAININFO*)pParam, pData->pvData);
			break;

		case AN_TRAINEND:
		  bContinue = pALNObj->OnTrainEnd((TRAININFO*)pParam, pData->pvData);
			break;

		case AN_EPOCHSTART:
		  bContinue = pALNObj->OnEpochStart((EPOCHINFO*)pParam, pData->pvData);
			break;

		case AN_EPOCHEND:
      bContinue = pALNObj->OnEpochEnd((EPOCHINFO*)pParam, pData->pvData);
			break;

		case AN_ADAPTSTART:
		  bContinue = pALNObj->OnAdaptStart((ADAPTINFO*)pParam, pData->pvData);
      break;

		case AN_ADAPTEND:
		  bContinue = pALNObj->OnAdaptEnd((ADAPTINFO*)pParam, pData->pvData);
      break;

    case AN_LFNADAPTSTART:
      bContinue = pALNObj->OnLFNAdaptStart((LFNADAPTINFO*)pParam, pData->pvData);
      break;
    
    case AN_LFNADAPTEND:
      bContinue = pALNObj->OnLFNAdaptEnd((LFNADAPTINFO*)pParam, pData->pvData);
      break;
    
    case AN_VECTORINFO:
      bContinue = pALNObj->OnVectorInfo((VECTORINFO*)pParam, pData->pvData);
      break;
	}

	return bContinue;
}
