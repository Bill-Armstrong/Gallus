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

// alnlfnanalysis.cpp


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

static int ALNAPI ValidateALNLFNAnalysis(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      void*& pvAnalysis,
                                      int& nLFNStats);

#ifdef _DEBUG
static void DebugValidateALNLFNAnalysis(const ALN* pALN,
                                     ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     void*& pvAnalysis,
                                     int& nLFNStats);
#endif


// helper to sort data set by active LFN
struct LFNSORT
{
  ALNNODE* pActiveLFN;
  float* afltInputRow;
  float* afltOutputRow;
	float* afltResultRow; //Added August 28, 1999.
};
static int __cdecl CompareLFNs(const void* pElem1, const void* pElem2);

// helpers to calc T and F probability functions
static float CalcProbF(float fltF, float fltV1, float fltV2);
static float CalcProbT(float fltT, float fltV);

// helpers to calc error, regression, and total sum of squares
static void CalcSquares(float* afltDesired, float* afltResult,
                        int nStart, int nEnd,
                        float& fltESS, float& fltRSS, float& fltTSS);

// helpers to allocate / deallocate analysis storage
struct LFNINFO
{
  ALNNODE* pLFN;
  LFNSTATS LFNStats;
  int nWStats;
  LFNWEIGHTSTATS aWStats[1];   // actually, nWStats elements
};

struct LFNANALYSIS
{
  int nLFNInfoSize;     // size of each LFNINFO element
  int nLFNInfo;         // nuumber of LFNINFO elements
  LFNINFO aLFNInfo[1];  // actually, nLFNInfo elements
};

static LFNANALYSIS* AllocLFNAnalysis(const ALN* pALN);
static void FreeLFNAnalysis(LFNANALYSIS* pLFNAnalysis);

inline
LFNINFO* GetLFNInfo(LFNANALYSIS* pLFNAnalysis, int nLFN)
{
  return (LFNINFO*)(((char*)(pLFNAnalysis->aLFNInfo)) + (nLFN * pLFNAnalysis->nLFNInfoSize));
}


// analyze LFNs
ALNIMP int ALNAPI ALNLFNAnalysis(const ALN* pALN,
                                 ALNDATAINFO* pDataInfo,
                                 const ALNCALLBACKINFO* pCallbackInfo,
                                 void*& pvAnalysis,
                                 int& nLFNStats)
{
	
  int nReturn = ValidateALNLFNAnalysis(pALN, pDataInfo, pCallbackInfo,
                                       pvAnalysis, nLFNStats);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  #ifdef _DEBUG
    DebugValidateALNLFNAnalysis(pALN, pDataInfo, pCallbackInfo,
                                pvAnalysis, nLFNStats);
  #endif

  // arrays

  // evaltree results
  float* afltResult = NULL;
  ALNNODE** apActiveLFNs = NULL;
  float* afltInput = NULL;
  float* afltOutput = NULL;
  LFNSORT* aLFNSort = NULL;

  // covariance calculation (on single block of inputs and outputs for one LFN)
  float* afltX = NULL;      // RMA based ALN input vectors (nRows * nCols); nCols = nDim
  float* afltY = NULL;      // RMA based desired column vector (nRows); nRows = #points this LFN
  float* afltC = NULL;      // RMA covariance matrix (nCols * nCols)
  float* h_A = NULL;      // RMA fitted parameter vector (nCols); h_A can be computed by SVD
  float* afltS = NULL;      // RMA std dev vector (nRows); output of SVD using afltX
  float* afltU = NULL;      // RMA U matrix (nRows * nCols); output of SVD using afltX
  float* afltV = NULL;      // RMA V matrix (nCols * nCols); output of SVD using afltX
  float* afltW = NULL;      // RMA W matrix (nCols); singular values from the SVD of matrix afltX 
  float* afltALN = NULL;    // RMA based ALN result column vector (nRows); the given fit to the data

	// idea: Use  std::vector< std::vector<Indices> > indices(layers, std::vector<Indices>(corrections));

  // LFN analysis
  LFNANALYSIS* pLFNAnalysis = NULL;

  try
  {
    // allocate analysis struct
    pLFNAnalysis = AllocLFNAnalysis(pALN);
    if (pLFNAnalysis == NULL)
      ThrowALNMemoryException();

    // see how many points there are
    long nTRcurrSamples = pDataInfo->nTRcurrSamples;
    
    // get dim... we use all the input variables of the ALN,
    // plus an explicit bias variable always equal to one,
    int nDim = pALN->nDim;
		int nCols = nDim - 1;//??
		int nOutput = pALN->nOutput;

    // allocate arrays
    afltResult = new float[nTRcurrSamples];
    apActiveLFNs = new ALNNODE*[nTRcurrSamples];
    afltInput = new float[nTRcurrSamples * nDim];
    afltOutput = new float[nTRcurrSamples];
    aLFNSort = new LFNSORT[nTRcurrSamples];
    afltX = new float[nTRcurrSamples * nCols];
		afltY = new float[nTRcurrSamples];
    afltC = new float[nCols * nCols];
    h_A = new float[nCols];
    afltS = new float[nTRcurrSamples];
    afltU = new float[nTRcurrSamples * nTRcurrSamples];
    afltV = new float[nCols * nCols];
    afltW = new float[nCols];
    afltALN = new float[nTRcurrSamples];


    if (!(afltResult && apActiveLFNs && afltInput && afltOutput && 
          aLFNSort && afltX && afltALN && afltC && h_A && afltS &&
          afltU && afltV && afltW && afltY))
    {
      ThrowALNMemoryException();
    }

    // init S
    for (long i = 0; i < nTRcurrSamples; i++)
    {
      afltS[i] = 1.0; // start off with std dev 1.0 in each axis
    }
    
    // evaluate on data
    int nStart, nEnd;
    nReturn = EvalTree(pALN->pTree, pALN, pDataInfo, pCallbackInfo, 
                       afltResult, &nStart, &nEnd, FALSE,
                       apActiveLFNs, afltInput, afltOutput);
    if (nReturn != ALN_NOERROR)
    {
      ThrowALNException();
    }

    // sort input/output arrays by responsible LFN
		// pointers in LFNSORT entries refer to  *unsorted* arrays
		// In the afltInput array, the nOutput components have been set to 1.0 by EvalTree.
    for (int i = nStart; i <= nEnd; i++)
    {
      LFNSORT& lfnsort = aLFNSort[i];
      lfnsort.pActiveLFN = apActiveLFNs[i];
      lfnsort.afltInputRow = afltInput + (i * nDim);
			// WWA I think it would be better to use afltOutputCol and afltResultCol as names:
			// The way the SVD is presented as
			// Matrix of input rows x parameter column = Result column.
			// The Output column is originally in the matrix of input rows
			// as the column of desired values (later replaced during evaluation of the DTREE)
      lfnsort.afltOutputRow = afltOutput+i; // desired
			afltOutput[i];
			lfnsort.afltResultRow = afltResult+i; // ALN output
    }
    // sort the LFN index structure
    qsort(aLFNSort + nStart, nEnd - nStart + 1, sizeof(LFNSORT), CompareLFNs);
  
    // calc covariance matrix for each block of data points
    int nBlockStart = nStart; // Monroe had 0, changed by wwa Aug 21 1999
    ALNNODE* pBlockLFN = aLFNSort[nStart].pActiveLFN;
    nLFNStats = 0;
    for (int i = nStart; i <= nEnd; i++)
    {
      if (i == nEnd || aLFNSort[i + 1].pActiveLFN != pBlockLFN)
      {
        // get pointer to current LFN info
        LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStats);

        // increment LFN stat count
        nLFNStats += 1;

        // end of LFN block... copy data to scratch input and output
        int nBlockPoints = i - nBlockStart + 1;
        for (int j = 0; j < nBlockPoints; j++)
        {
          // input row (with 1.0 in nOutput component)  WWA we don't copy the 1.0 into X
          memcpy(afltX + j * nCols,                             
                 aLFNSort[nBlockStart + j].afltInputRow, 
                 nCols * sizeof(float)); 
          // ALN output
          afltALN[j] = *(aLFNSort[nBlockStart + j].afltResultRow); //contains ALN output

				  // desired output
					afltY[j] = *(aLFNSort[nBlockStart + j].afltOutputRow); // contains desired
        }
				// The result of t4 above is correct
        // calc covariance

				// NB The actual ALN outputs play no role in the covariance calculation
				// Y below is the desired output
				// The ALN outputs could be compared to the outputs using parameter vector A found below
        float fltChiSq = 0;
        BOOL bSuccess = CalcCovariance(nCols,
                                       nBlockPoints, 
                                       afltX,
                                       afltY,
                                       afltC,
                                       h_A,
                                       afltS,
                                       afltU,
                                       afltV,
                                       afltW,
                                       fltChiSq);
        // redo the Y copy,just in case it changed WWA because of a memcpy that had nDim instead of nCols in
				// some places, the Y values for j = 0 to j = 4 were wrong, both inside and leaving the calc routine.
				// After correction, this part is no longer necessary.
        // for (int j = 0; j < nBlockPoints; j++)
				//           
				  // desired 
					// 	afltY[j] = *(aLFNSort[nBlockStart + j].afltOutputRow); 

          // ALN output
					//    afltALN[j] = *(aLFNSort[nBlockStart + j].afltResultRow); 
					// }

        // calc RSS, ESS, TSS
        float fltRSS, fltESS, fltTSS;
        CalcSquares(afltY, afltALN, 0, nBlockPoints -1, fltESS, fltRSS, fltTSS);
        
        pLFNInfo->LFNStats.fltRSS = fltRSS;
        pLFNInfo->LFNStats.fltESS = fltESS;
        
        // calc R2 for the LFN
        pLFNInfo->LFNStats.fltR2 = fltRSS / fltTSS;

        // calc degrees of freedom
        pLFNInfo->LFNStats.fltDF = nBlockPoints - nDim; 

				pLFNInfo->LFNStats.fltSEE = NAN; // set to NAN if we don't have fltDF>0
        pLFNInfo->LFNStats.fltF   = NAN;
        pLFNInfo->LFNStats.fltFp  = NAN;
        if (pLFNInfo->LFNStats.fltDF > 0)
        {
          ASSERT(nDim >= 1);
          float fltV1 = nDim - 1;
          float fltV2 = pLFNInfo->LFNStats.fltDF;

          // calc SEE "standard error of estimate" for the LFN
					pLFNInfo->LFNStats.fltSEE = sqrt(fltESS / fltV2);

					// calc F and corresponding p value for the LFN
          pLFNInfo->LFNStats.fltF = (fltRSS / fltV1) / (fltESS / fltV2);
          pLFNInfo->LFNStats.fltFp = CalcProbF(pLFNInfo->LFNStats.fltF, fltV1, fltV2);
        }
              
        // set the node
        pLFNInfo->pLFN = pBlockLFN;

        // set number of weight stats
        pLFNInfo->nWStats = nDim;

        // calc T, and corresponding p values for each LFN weight
        const float* afltLFNW = LFN_W(pBlockLFN);
        float fltBias = *afltLFNW++;  // skip past bias
        float fltV = pLFNInfo->LFNStats.fltDF;
				float fltLFNweight; // used for ALN weights to be compared to SVD result
        for (int k = 0; k < nDim; k++)
				{
          // get LFNWEIGHTSTATS pointer
          LFNWEIGHTSTATS* pWStat = &(pLFNInfo->aWStats[k]);

          // set weights from LFN!!!
          if (k == pALN->nOutput)
          fltLFNweight = fltBias;		// bias replaces output weight for stats
          else
          fltLFNweight = afltLFNW[k];	// kth input weight

				 // the weights are now obtained from SVD and are in h_A
				 // (as Monroe suggested doing). 
				 // fltBias should approximate h_A[k]if k = nOutput 
				 pWStat->fltW = h_A[k];

          // set T stat
          pWStat->fltSEw = 0;
          pWStat->fltT = NAN;
          pWStat->fltTp = 1.0;
          if (pLFNInfo->LFNStats.fltDF > 0)
          {
            // calc standard error
            float fltCkk = afltC[k*nDim + k];  // diagonal element Ckk 

            pWStat->fltSEw = sqrt(fltCkk) * pLFNInfo->LFNStats.fltSEE;

            if (pWStat->fltSEw != 0.0)
            {
              // calc T, Tp
							// WWA the following line didn't include the term afltLFNW[k],
							// thus making it refer to the k-th weight, not to the
							// difference of ALN computed and SVD computed weights.
							// In the original interpretation, that would decide whether
							// the weight could just as well be zero, ie not significant.
							// Now we want to use it to say if the difference of ALN
							// weight from the SVD computed one is significant.  If ALN
							// does well, then all of the T's should be less than 2.0 or 3.0 
							// and the probabilities Tp should be not to far from 1.0.
							// Is this OK?
              pWStat->fltT = (fltLFNweight - pWStat->fltW) / pWStat->fltSEw;
              pWStat->fltTp = CalcProbT(fabs(pWStat->fltT), fltV);
            }
          }
        }
        
        // start of new LFN block
        if (i < nEnd)
        {
          nBlockStart = i + 1;  
          pBlockLFN = aLFNSort[nBlockStart].pActiveLFN;
        }
      }   
    }  // end of loop start to end        
  } // end of try block

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
	delete[] afltALN;
	delete[] afltW;
	delete[] afltV;
	delete[] afltU; // this is where the heap corruption happens
	delete[] afltS;
	delete[] h_A;
	delete[] afltC;
	delete[] afltY;
	//delete[] afltX;
	delete[] aLFNSort;
	delete[] afltOutput;
	delete[] afltInput;
	delete[] apActiveLFNs;
	delete[] afltResult;

  // set the analysis value
  if (nReturn != ALN_NOERROR && pLFNAnalysis != NULL)
  {
    FreeLFNAnalysis(pLFNAnalysis);
  }
  else
  {
    pvAnalysis = (void*)pLFNAnalysis;
  }
  
  return nReturn;
}

// get LFN stats
ALNIMP int ALNAPI ALNLFNStats(void* pvAnalysis, 
                              int nLFNStat,
                              LFNSTATS* pLFNStats,
                              int* pnWeightStats,
                              ALNNODE** ppLFN)
{
  if (pvAnalysis == NULL)
    return ALN_GENERIC;

  // make sure LFN index in range
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)pvAnalysis;
  if (nLFNStat < 0 || nLFNStat >= pLFNAnalysis->nLFNInfo)
    return ALN_GENERIC;

  // get LFNINFO
  LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStat);
  
  // copy stat values
  if (pLFNStats != NULL)
    memcpy(pLFNStats, &(pLFNInfo->LFNStats), sizeof(LFNSTATS));
  
  if (pnWeightStats != NULL)
    *pnWeightStats = pLFNInfo->nWStats;

  if (ppLFN != NULL)
    *ppLFN = pLFNInfo->pLFN;

  return ALN_NOERROR;
}

// get LFN weight stats
ALNIMP int ALNAPI ALNLFNWeightStats(void* pvAnalysis, 
                                    int nLFNStat,
                                    int nWeightStat,
                                    LFNWEIGHTSTATS* pLFNWeightStats)
{
  if (pvAnalysis == NULL || pLFNWeightStats == NULL)
    return ALN_GENERIC;

  // make sure LFN index in range
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)pvAnalysis;
  if (nLFNStat < 0 || nLFNStat >= pLFNAnalysis->nLFNInfo)
    return ALN_GENERIC;

  // get LFNINFO
  LFNINFO* pLFNInfo = GetLFNInfo(pLFNAnalysis, nLFNStat);

  // make sure W index is in range
  if (nWeightStat < 0 || nWeightStat >= pLFNInfo->nWStats)
    return ALN_GENERIC;

  // copy stat values
  memcpy(pLFNWeightStats, 
         &(pLFNInfo->aWStats[nWeightStat]), 
         sizeof(LFNWEIGHTSTATS));

  return ALN_NOERROR;
}

// free analysis memory
ALNIMP int ALNAPI ALNLFNFreeAnalysis(void* pvAnalysis)
{
  if (pvAnalysis != NULL)
    FreeLFNAnalysis((LFNANALYSIS*)pvAnalysis);
    
  return ALN_NOERROR;
}


// validate params
static int ALNAPI ValidateALNLFNAnalysis(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      void*& pvAnalysis,
                                      int& nLFNStats)
{
  int nReturn = ValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  if (nReturn != ALN_NOERROR)
    return nReturn;
  
  //if (pvAnalysis == NULL || nLFNStats != 0)
  //  return ALN_GENERIC;

  return ALN_NOERROR;
}

// debug version ASSERTS if bad params
#ifdef _DEBUG
static void DebugValidateALNLFNAnalysis(const ALN* pALN,
                                     ALNDATAINFO* pDataInfo,
                                     const ALNCALLBACKINFO* pCallbackInfo,
                                     void*& pvAnalysis,
                                     int& nLFNStats)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
}
#endif

static int __cdecl CompareLFNs(const void* pElem1, const void* pElem2)
{
  LFNSORT& lfnsort1 = *(LFNSORT*)pElem1;
  LFNSORT& lfnsort2 = *(LFNSORT*)pElem2;

  if (lfnsort1.pActiveLFN < lfnsort2.pActiveLFN)
    return -1;
  else if (lfnsort1.pActiveLFN > lfnsort2.pActiveLFN)
    return 1;

  return 0;
}

// helpers to calc T and F probability functions
static float CalcProbF(float fltF, float fltV1, float fltV2)
{
  // F-Distribution probability function - not two tailed as in
  // Press et al p229, p619.  Here it is used in regression and is 1-tailed.
	// It is testing H0: all weights = 0 (except for the bias weight).
  
	return ibeta(fltV2 / 2.0, fltV1 / 2.0, fltV2 / (fltV2 + (fltV1 * fltF))); //ibeta ??

//
}

static float CalcProbT(float fltT, float fltV)
{
  // T-Distribution probability function - single tailed
  // see Press et al p228, p616
  return ibeta(fltV / 2.0, 0.5, fltV / (fltV + (fltT * fltT)));  //ibetac??
}

// helpers to calc error, regression, and total sum of squares
static void CalcSquares(float* afltDesired, float* afltResult,
                        int nStart, int nEnd,
                        float& fltESS, float& fltRSS, float& fltTSS)
{
  ASSERT(afltDesired && afltResult);
  ASSERT(nStart <= nEnd);

  fltESS = 0;
  fltTSS = 0;
  fltRSS = 0;

  float fltMeanDes = 0;
  int nElem = nEnd - nStart + 1;
  ASSERT(nElem >= 1);

  // calc error sum of squares and average of desired
  for (int i = nStart; i <= nEnd ;i++)
  {
    fltMeanDes += afltDesired[i];
    float fltError = afltResult[i] - afltDesired[i];
    fltESS += fltError * fltError;
  }
  fltMeanDes /= nElem;

  // calc total sum of squares
  for (int i = nStart; i <= nEnd; i++)
  {
    float fltError = afltDesired[i] - fltMeanDes; // Monroe had afltResult[i] - fltMeanDes;
    fltTSS += fltError * fltError;
  }

  // calc regression sum of squares
  fltRSS = fltTSS - fltESS;
}

static LFNANALYSIS* AllocLFNAnalysis(const ALN* pALN)
{
  // calc number of LFNs
  int nTotalLFNs = 0;
  int nAdaptedLFNs = 0;
  CountLFNs(pALN->pTree, nTotalLFNs, nAdaptedLFNs);

  if (nTotalLFNs == 0)
    return NULL;

  // calc number of bytes needed for entire struct
  int nWStatBytes = pALN->nDim * sizeof(LFNWEIGHTSTATS);
  int nLFNInfoBytes = sizeof(LFNINFO) - sizeof(LFNWEIGHTSTATS) + nWStatBytes;
  int nTotalBytes = sizeof(LFNANALYSIS) - sizeof(LFNINFO) + (nTotalLFNs * nLFNInfoBytes);

  // alloc block
  LFNANALYSIS* pLFNAnalysis = (LFNANALYSIS*)calloc(1, nTotalBytes);               
  if (pLFNAnalysis == NULL)
    return NULL;

  // init block
  pLFNAnalysis->nLFNInfoSize = nLFNInfoBytes;
  pLFNAnalysis->nLFNInfo = nTotalLFNs;
  
  return pLFNAnalysis;
}

static void FreeLFNAnalysis(LFNANALYSIS* pLFNAnalysis)
{
  if (pLFNAnalysis != NULL)
    free(pLFNAnalysis);
}

