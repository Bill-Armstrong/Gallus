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

// alntrain.cpp
// training support routines

// libaln/src/alntrain.cpp
// Revision date: Dec 26, 2019


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include "alnpp.h"
#include ".\cmyaln.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Helper declarations relating to ALN tree growth
void splitControl(ALN*, ALNDATAINFO*); // This does a test to see if a piece fits well or must be split.
extern BOOL bALNgrowable = TRUE; //If FALSE, no splitting happens, e.g. for linear regression.
BOOL bStopTraining = FALSE; // This causes training to stop when all leaf nodes have stopped splitting. This means all linear regression resultss will not change.
extern BOOL bDistanceOptimization; // This prevents considering leaf nodes so far away from the current input sample that they must be cut off in the max-min structure.
extern BOOL bClassify2; // We don't want to do weight decay if this is not a classification problem
extern float WeightBound; // This is a maximum and minimum bound on weights.
extern float WeightDecay; //  This is a factor close to 1.0. It is used during classification into two classes when lower weights give better generalization.
// GPU transfer items
//extern float* h_A, , h_B, d_A, d_B;
//extern int  DATA_SZ, ELEMENT_SZ:

// Train calls ALNTrain, which expects data in a monolithic array, row major order, ie,
//   row 0 col 0, row 0 col 1, ..., row 0 col n,
//   row 1 col 0, row 1 col 1, ..., row 1 col n,
//   ...,
//   row m col 0, row m col 1, ..., row m col n,
// If MSEorF is positive, this is just a level of training error above which pieces can split, then there must be 1-to-1 nDim entries in the data array.
// If MSEorF is negative, the data array also includes, for each sample, nDim entries for the closest other sample minus this one and a float for the distance between them *\

// nNotifyMask is the bitwise OR of all the notifications that should
// be sent during training callbacks. Special cases: AN_ALL, AN_NONE.
// ALN train returns an ALN_* error code, (ALN_NOERROR on success).

// Helper declarations
int ALNAPI ValidateALNTrainInfo(const ALN* pALN,
                                       ALNDATAINFO* pDataInfo,
                                       const ALNCALLBACKINFO* pCallbackInfo,
                                       int nMaxEpochs,
                                       float fltMinRMSErr,
                                       float fltLearnRate);

#ifdef _DEBUG
void DebugValidateALNTrainInfo(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      int nMaxEpochs,
                                      float fltMinRMSErr,
                                      float fltLearnRate);
#endif


static int ALNAPI DoTrainALN(ALN* pALN,
                             ALNDATAINFO* pDataInfo,
                             const ALNCALLBACKINFO* pCallbackInfo,
                             int nMaxEpochs,
                             float fltMinRMSErr,
                             float fltLearnRate,
                             BOOL bJitter);



ALNIMP int ALNAPI ALNTrain(ALN* pALN,
                           ALNDATAINFO* pDataInfo,
                           const ALNCALLBACKINFO* pCallbackInfo,
                           int nMaxEpochs,
                           float fltMinRMSErr,
                           float fltLearnRate,
                           BOOL bJitter)
{
	int nReturn; // = ValidateALNTrainInfo(pALN, pDataInfo, pCallbackInfo,
    //                             nMaxEpochs, fltMinRMSErr, fltLearnRate); MYTEST stopping validation for training
 
	//if (nReturn == ALN_NOERROR)
  {
    // train if the ALN is successfully prepped
    if (PrepALN(pALN))
    {
      nReturn = DoTrainALN(pALN, pDataInfo, pCallbackInfo,
                           nMaxEpochs, fltMinRMSErr, fltLearnRate,
                           bJitter);
    }
    else
    {
      nReturn = ALN_GENERIC;
    }
  }
  return nReturn;
}


///////////////////////////////////////////////////////////////////////////////
// workhorse of TrainALN 
static int ALNAPI DoTrainALN(ALN* pALN,
                             ALNDATAINFO* pDataInfo,
                             const ALNCALLBACKINFO* pCallbackInfo,
                             int nMaxEpochs,
                             float fltMinRMSErr,
                             float fltLearnRate,
                             BOOL bJitter)
{
	//#ifdef _DEBUG
	// DebugValidateALNTrainInfo(pALN, pDataInfo, pCallbackInfo, nMaxEpochs, 
	//                           fltMinRMSErr, fltLearnRate);
	//#endif

	int nReturn = ALN_NOERROR;		    // assume success
	int nDim = pALN->nDim;
	long nTRcurrSamples = pDataInfo->nTRcurrSamples;
	ALNNODE* pTree = pALN->pTree;
	float* afltX;                    // input vector
	long* anShuffle = NULL;				    // point index shuffle array
	CCutoffInfo* aCutoffInfo = NULL;  // eval cutoff speedup

	TRAININFO traininfo;					    // training info
	EPOCHINFO epochinfo;					    // epoch info
	TRAINDATA traindata;              // training data

	int nNotifyMask = (pCallbackInfo == NULL) ? AN_NONE : pCallbackInfo->nNotifyMask;
	void* pvData = (pCallbackInfo == NULL) ? NULL : pCallbackInfo->pvData;
	ALNNOTIFYPROC pfnNotifyProc = (pCallbackInfo == NULL) ? NULL : pCallbackInfo->pfnNotifyProc;

	// init traindata
	traindata.fltLearnRate = fltLearnRate; // an epoch is one  pass through the training data
	traindata.nNotifyMask = nNotifyMask;
	traindata.pvData = pvData;
	traindata.pfnNotifyProc = pfnNotifyProc;

	// calc start and end points of training
	long nStart, nEnd;
	//CalcDataEndPoints(nStart, nEnd, pALN, pDataInfo); This is very simple since there are no lags as in time series.
	nStart = 0;
	nEnd = pDataInfo->nTRcurrSamples - 1;

	try	// main processing block
	{
		// allocate input vector
			afltX = new float[nDim];
		if (!afltX) ThrowALNMemoryException();
		memset(afltX, 0, sizeof(float) * nDim); // this has space for all the inputs and the output value

		// allocate and init shuffle array
		anShuffle = new long[nEnd - nStart + 1L];
		if (!anShuffle) ThrowALNMemoryException();
		for (long i = nStart; i <= nEnd; i++)
			anShuffle[i - nStart] = i - nStart;

		// allocate and init cutoff info array
		// pLFN will contain a pointer to the active LFN of a piece
		// when the input is on that piece.  It will speed up cutoffs in evaluation.
		aCutoffInfo = new CCutoffInfo[nEnd - nStart + 1L];
		if (!aCutoffInfo) ThrowALNMemoryException();
		for (int i = nStart; i <= nEnd; i++)
			aCutoffInfo[i - nStart].pLFN = NULL;
		// count total number of LFNs in ALN
		int nLFNs = 0;
		int nAdaptedLFNs = 0;
		CountLFNs(pALN->pTree, nLFNs, nAdaptedLFNs);

		// init callback info 
		traininfo.nEpochs = nMaxEpochs;
		traininfo.nLFNs = nLFNs;
		traininfo.nActiveLFNs = 0;
		traininfo.fltRMSErr = 0.0;
	
		epochinfo.nEpoch = 0;
		epochinfo.nLFNs = nLFNs;
		epochinfo.nActiveLFNs = nAdaptedLFNs;
		epochinfo.fltEstRMSErr = 0.0;
    
		// notify beginning of training
		if (CanCallback(AN_TRAINSTART, pfnNotifyProc, nNotifyMask))
		{
			TRAININFO ti(traininfo);  // make copy to send!
			Callback(pALN, AN_TRAINSTART, &ti, pfnNotifyProc, pvData);
		}

		///// begin epoch loop
		// We reset counters for splitting when adaptation has had a chance to adjust pieces very
		// closely to the training samples, e.g. the limited number of pieces fits well.
		// We do under nMaxEpochs/2 training, then allow splitting
		// We must set nMaxEpochs considering epochsize, learning rate,  prescribed RMS error, etc.
		for (int nEpoch = 0; nEpoch < nMaxEpochs; nEpoch++)
		{
			if (nEpoch == nMaxEpochs / 2) bDistanceOptimization = FALSE; // MYTEST Turn optimization off until the end of the splitting epoch; see if necessary later MYTEST
			int nCutoffs = 0;
			// notify beginning of epoch
			epochinfo.nEpoch = nEpoch;
			if (CanCallback(AN_EPOCHSTART, pfnNotifyProc, nNotifyMask))
			{
				EPOCHINFO ei(epochinfo);  // make copy to send!
				Callback(pALN, AN_EPOCHSTART, &ei, pfnNotifyProc, pvData);
			}

			// track squared error
			float fltSqErrorSum = 0;

			// We prepare a random reordering of the training data for the next epoch
			Shuffle(nStart, nEnd, anShuffle);

			long nSample; // The number of training samples may be huge.
				// this does all the samples in an epoch in a randomized order.

			for (nSample = nStart; nSample <= nEnd; nSample++)
			{
				long nTrainSample = anShuffle[nSample - nStart]; // A sample is picked for training
				ASSERT((nTrainSample + nStart) <= nEnd);
				

				// fill input vector
				FillInputVector(pALN, afltX, nTrainSample, nStart, pDataInfo, pCallbackInfo);

				// if we have our first data point, init LFNs on first pass
				if (nEpoch == 0 && nSample == nStart)
					InitLFNs(pTree, pALN, afltX);


				// jitter the data point
				if (bJitter) Jitter(pALN, afltX);

				// do an adapt eval to get active LFN and distance, and to prepare
				// tree for adaptation
				ALNNODE* pActiveLFN = NULL;
				CCutoffInfo& cutoffinfo = aCutoffInfo[nSample];

				// START MYTEST insert something for speed
				// Get the cutoff info if possible from a previous epoch
				// pActiveLFN = cutoffinfo.pLFN;
				//if ((nEpoch > 1) && !LFN_CANSPLIT(pActiveLFN) && ( ALNRandFloat() > 0.5))	continue;
				// END MYTEST

				float flt = AdaptEval(pTree, pALN, afltX, &cutoffinfo, &pActiveLFN);

				// track squared error before adapt, since adapt routines
				// do not relcalculate value of adapted surface
				fltSqErrorSum += flt * flt;

				// notify start of adapt
				if (CanCallback(AN_ADAPTSTART, pfnNotifyProc, nNotifyMask))
				{
					ADAPTINFO adaptinfo;
					adaptinfo.nAdapt = nSample - nStart;
					adaptinfo.afltX = afltX;
					adaptinfo.fltErr = flt;
					Callback(pALN, AN_ADAPTSTART, &adaptinfo, pfnNotifyProc, pvData);
				}

				// do a useful adapt to correct any error
				traindata.fltGlobalError = flt;
				if(nEpoch > 0) Adapt(pTree, pALN, afltX, 1.0, TRUE, &traindata);// we should not adapt in the epoch when counting hits!!
				// notify end of adapt
				if (CanCallback(AN_ADAPTEND, pfnNotifyProc, nNotifyMask))
				{
					ADAPTINFO adaptinfo;
					adaptinfo.nAdapt = nSample - nStart;
					adaptinfo.afltX = afltX;
					adaptinfo.fltErr = flt;
					Callback(pALN, AN_ADAPTEND, &adaptinfo, pfnNotifyProc, pvData);
				}
			}	// end for each point in data set
			if (bClassify2 && WeightDecay != 1.0F)
			{
				DecayWeights(pTree, pALN, WeightBound, WeightDecay); //We decay weights only when there are a few samples on the piece at the end of training
			}

			// estimate RMS error on training set for this epoch
			epochinfo.fltEstRMSErr = sqrt(fltSqErrorSum / nTRcurrSamples);

			// calc true RMS if estimate below min, or if last epoch, or every 10 epochs when jittering
			if (epochinfo.fltEstRMSErr <= fltMinRMSErr || nEpoch == nMaxEpochs)
			{
				epochinfo.fltEstRMSErr = DoCalcRMSError(pALN, pDataInfo, pCallbackInfo);
			}

			// notify end of epoch
			nLFNs = nAdaptedLFNs = 0;
			CountLFNs(pALN->pTree, nLFNs, nAdaptedLFNs);
			epochinfo.nLFNs = nLFNs;
			epochinfo.nActiveLFNs = nAdaptedLFNs;

			// update train info, too
			traininfo.nEpochs = nEpoch;
			traininfo.nLFNs = epochinfo.nLFNs;
			traininfo.nActiveLFNs = epochinfo.nActiveLFNs;
			traininfo.fltRMSErr = epochinfo.fltEstRMSErr;	// used to terminate epoch loop
			
			if (nEpoch == nMaxEpochs - 1 && CanCallback(AN_EPOCHEND, pfnNotifyProc, nNotifyMask))
			{
				EPOCHINFO ei(epochinfo);  // make copy to send!
				Callback(pALN, AN_EPOCHEND, &ei, pfnNotifyProc, pvData);
			}

			// Split candidate LFNs in a middle epoch in this call to ALNTrain.
			if (nEpoch == nMaxEpochs/2)
			{
				bStopTraining = TRUE;  // this will be set to FALSE by any leaf node needing further training after splitControl()
				splitControl(pALN, pDataInfo);  // This leads to leaf nodes splitting
				bDistanceOptimization = TRUE;
			}
		} // end epoch loop

		// notify end of training
		
		if (CanCallback(AN_TRAINEND, pfnNotifyProc, nNotifyMask))
		{
			// don't bother copying traininfo, since this is the last message sent
			// we don't care if user changes it!
			Callback(pALN, AN_TRAINEND, &traininfo, pfnNotifyProc, pvData);
		}
		delete[] afltX;
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

	delete[] anShuffle;
	delete[] aCutoffInfo;
	return nReturn;
}

// validate ALNTRAININFO struct
static int ALNAPI ValidateALNTrainInfo(const ALN* pALN,
                                       ALNDATAINFO* pDataInfo,
                                       const ALNCALLBACKINFO* pCallbackInfo,
                                       int nMaxEpochs,
                                       float fltMinRMSErr,
                                       float fltLearnRate)
{
  int nReturn = ValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  if (nReturn != ALN_NOERROR)
    return nReturn;

  // need at least one epoch
  if (nMaxEpochs <= 0)
  {
    return ALN_GENERIC;
  }

  // need non-negative error
  if (fltMinRMSErr < 0)
  {
    return ALN_GENERIC;
  }

  // learnrate should be positive and probably no more than than 0.5
	if (fltLearnRate < 0.0 || fltLearnRate > 0.5)
	{
		return ALN_GENERIC;
	}

  return ALN_NOERROR;
}

#ifdef _DEBUG
static void DebugValidateALNTrainInfo(const ALN* pALN,
                                      ALNDATAINFO* pDataInfo,
                                      const ALNCALLBACKINFO* pCallbackInfo,
                                      int nMaxEpochs,
                                      float fltMinRMSErr,
                                      float fltLearnRate)
{
  DebugValidateALNDataInfo(pALN, pDataInfo, pCallbackInfo);
  ASSERT(nMaxEpochs > 0 && fltMinRMSErr >= 0 && fltLearnRate > 0.0 && fltLearnRate <= 0.5);
}
#endif
