// ALN Library (libaln)
// file split_ops.cpp
// 
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

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// include files
#include <stdafx.h>
#include <alnpp.h>

// include classes
#include ".\cmyaln.h"
#include "aln.h"

// We use fltRespTotal in two ways and the following definition helps.
#define NOISEVARIANCE fltRespTotal

extern BOOL bStopTraining; // This becomes TRUE and stops training when pieces are no longer splitting.
void setSplitAlpha(ALNDATAINFO* pDataInfo);
void splitControl(ALN* pALN, ALNDATAINFO* pDataInfo);
void zeroSplitValues(ALN* pALN, ALNNODE* pNode);
void splitUpdateValues(ALN * pALN, ALNDATAINFO* pDataInfo);
void doSplits(ALN* pALN, ALNNODE* pNode, float fltLimit);
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);
extern float WeightDecay;
extern float WeightBound;
extern BOOL bClassify2;
extern BOOL bConvex;
extern int SplitsAllowed;
extern int SplitCount;


// We use the first three fields in ALNLFNSPLIT (declared in aln.h)
// in two different ways: for training and between training intervals.
// The following routines use the SPLIT typedef between trainings, near the end of alntrain.cpp.
// During nMaxEpochs epochs of training, where the hyperplanes in leaf nodes change weights,
// the pieces adapt to the data the best they can. Following that, 
// a decision must be made for each leaf node whether to split or not.
// To do that, the differences of pairs of close sample values in the data
// which are stored in aNoiseSampleTool are adjusted for the current weights of the piece
// to compensate for the samples being in slightly different domain locations.
// Then splitcontrol takes the total square training errors of
// each piece. This is compared to the total of the noise variance samples.
// over the same training points.  If the average square training error is 
// greater than fltLimit times the average of the noise variance samples on the piece, then using
// an F-test, the piece is split because it does not yet fit within the limits of noise.


// Explanation of fltLimit
// fltLimit = 2.59 says that splitting of a linear piece is prevented when the mean square
// training error of a piece becomes less than 2.59 times the average of the noise variance
// samples on it. This value comes from tables of the F-test for d.o.f. > 7 and probability 90%.
// For 90% with 3 d.o.f the value is 5.39, i.e. with fewer d.o.f. training stops sooner
// and the training error will generally be larger than with a lower F-value.
// The values below for afltFconstant35 are interpolated and may not be accurate.
// The values for 25 are the reciprocals of those for 75.

static const float afltFconstant90[13]{ 9.00, 5.39, 4.11, 3.45, 3.05, 2.78, 2.59, 2.44, 2.32, 1.79, 1.61, 1.51, 1.40 };
static const float afltFconstant75[13]{ 3.00, 2.36, 2.06, 1.89, 1.78, 1.70, 1.64, 1.59, 1.55, 1.36, 1.28, 1.24, 1.19 };
static const float afltFconstant50[13]{ 1,1,1,1,1,1,1,1,1,1,1,1,1 };
// The following two have not had any beneficial effect. Who knows when they might be useful?
// static const float afltFconstant35[13]{ 0.58, 0.65, 0.70, 0.73, 0.75, 0.77, 0.78, 0.79, 0.80, 0.86, 0.88, 0.90, 0.92 };
static const float afltFconstant25[13]{ 0.333, 0.424, 0.485, 0.529, 0.562, 0.588, 0.610, 0.629, 0.645, 0.735, 0.781, 0.806, 0.840 };
static const float afltFconstant10[13]{ 0.111, 0.185, 0.243, 0.290, 0.327, 0.359, 0.386, 0.410, 0.431, 0.558, 0.621, 0.662, 0.714 };
static float afltF_Alpha[13]{ 1,1,1,1,1,1,1,1,1,1,1,1,1 };

void setSplitAlpha(ALNDATAINFO* pDataInfo)
{
	// If fltMSEorF < 0, this routine is called once before any training.
	float fltLimit = pDataInfo->fltMSEorF;
	if (fltLimit >= 0) return;
	if(-fltLimit == 50)
	{
		for (int i = 0; i < 13; i++) // We are doing an F-test
		{
			afltF_Alpha[i] = afltFconstant50[i];
		}
	}
	else  // this is according to tables for 25, 50, 75 and rest is approximate
	{
		for (int i = 0; i < 13; i++) // We are doing an F-test
		{
			afltF_Alpha[i] = pow(afltFconstant75[i], (-fltLimit -50)/25.0);
		}
	}
}

void splitControl(ALN* pALN, ALNDATAINFO* pDataInfo)  // routine
{
	if (SplitCount >= SplitsAllowed) return; // 
	float fltLimit = pDataInfo->fltMSEorF;
	ASSERT(pALN);
	ASSERT(pALN->pTree);
	// initialize all the SPLIT values to zero
	zeroSplitValues(pALN, pALN->pTree);
	// get square errors of pieces on training set and the noise variance estimates
	splitUpdateValues(pALN, pDataInfo);
	// With the above statistics, doSplits recursively determines splits of eligible pieces.
    doSplits(pALN, pALN->pTree, fltLimit);
    // Resetting the SPLIT components to zero by zeroSplitValues is done in alntrain.
}

// Routines that set some fields to zero

void zeroSplitValues(ALN* pALN, ALNNODE* pNode) // routine
{
	// initializes split counters of all leaf nodes before the next training period
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		zeroSplitValues(pALN, MINMAX_LEFT(pNode));
		zeroSplitValues(pALN, MINMAX_RIGHT(pNode));
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (LFN_CANSPLIT(pNode))
		{
			(pNode->DATA.LFN.pSplit)->nCount = 0;
			(pNode->DATA.LFN.pSplit)->fltSqError = 0;
			(pNode->DATA.LFN.pSplit)->NOISEVARIANCE = 0;
		}
	}
}

// Routines that get the training errors and noise variance values.

void splitUpdateValues(ALN * pALN, ALNDATAINFO* pDataInfo) // routine
{
	// Assign the square errors on the training set and the noise variance
	// sample values to the leaf nodes of the ALN.

	float desired = 0;
	float alnval = 0;
	int nDim = pALN->nDim;
	int nDimm1 = nDim - 1;
	int nDimt2 = nDim * 2;
	int nDimt2m1 = nDimt2 - 1;

	int nDimt2p1 = nDimt2 + 1;
	int nDimt2p1ti;
	float* afltX = (float*)malloc(nDim * sizeof(float));
	ALNNODE* pActiveLFN;
	float* afltTRdata = pDataInfo->afltTRdata;
	long nrows = pDataInfo->nTRcurrSamples;
	float fltMSEorF = pDataInfo->fltMSEorF;
	for (long i = 0; i < nrows; i++)
	{
		nDimt2p1ti = nDimt2p1 * i;
		for (int j = 0; j < nDimm1; j++) // just the domain values of the sample
		{
			afltX[j] = afltTRdata[nDimt2p1ti + j];
		}
		//afltX[nDimm1] = 0; // hide the desired value (not used anyway)
		alnval = ALNQuickEval(pALN, afltX, &pActiveLFN); // the current ALN value
		if (LFN_CANSPLIT(pActiveLFN)) // Skip this leaf node if it can't split anyway.//READ ACCESS VIOLATION pActiveLFN was 0x4E210
		{
			float noiseSampleTemp;
			desired = afltTRdata[nDimt2p1ti + nDimm1];
			(pActiveLFN->DATA.LFN.pSplit)->nCount++;
			(pActiveLFN->DATA.LFN.pSplit)->fltSqError += (alnval - desired) * (alnval - desired);
			if (fltMSEorF <= 0)
			{
				noiseSampleTemp = afltTRdata[nDimt2p1ti + nDimt2m1]; // Get the difference of desired sample values in the tool
				// This has to be corrected for the slopes of the LFN
				for (int kk = 0; kk < nDim - 1; kk++) // Just do the domain dimensions.
				{
					// get the weights for the LFN and correct the sample for slope
					// Adding 1 in kk + 1 skips the bias weight.
					noiseSampleTemp -= LFN_W(pActiveLFN)[kk + 1] * afltTRdata[nDimt2p1ti + nDim + kk];
				}
				// The following should be a sample for the noise variance of the piece
				// which is paired with data sample i in case fltMSEorF is negative and we are doing an F-test.
				(pActiveLFN->DATA.LFN.pSplit)->NOISEVARIANCE += 0.5 * noiseSampleTemp * noiseSampleTemp;
			}
		}
	} // end loop over both files
	free(afltX);
} // END of splitUpdateValues

void doSplits(ALN* pALN, ALNNODE* pNode, float fltMSEorF) // routine
{
	int nDim = pALN->nDim;
	int  CanSplitAbove = bClassify2 ? 1 : (int)(1.2F * (float)nDim + 1.0F);
	// This routine visits all the leaf nodes and determines whether or not to split.
	// If fltLimit < 0, it uses an F test with d.o.f. based on the number of samples counted,
	// but if fltLimit >= 0 it uses the actual fltLimit value to compare to the square training error.

	long Count;

	ASSERT(pNode);
	ALNNODE* pParent = NODE_PARENT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		// The following is the WWA optimization, which tries to avoid evaluating nodes not in the neighborhood of sample afltX
		// In preparation, we zero the count, centroid and sigma of this node, which will be changed in the children
		MINMAX_COUNT(pNode) = 0;
		for (int i = 0; i < nDim - 1; i++)
		{
			MINMAX_CENTROID(pNode)[i] = 0; // zero the centroid of the subtree on the depth-first path down the tree
			MINMAX_SIGMA(pNode)[i] = 0; // zero the half-width of the subtree on the way down
		}
		// Now do the children
		doSplits(pALN, MINMAX_LEFT(pNode), fltMSEorF); // N.B. left and right subtrees have centroids whose differences of components can be positive or negative
		doSplits(pALN, MINMAX_RIGHT(pNode), fltMSEorF);// i.e. the subtrees have no particular order in the nDim -1 coordinates of the domain
		// and pop back up here
		for (int i = 0; i < nDim - 1; i++)
		{
			MINMAX_CENTROID(pNode)[i] /= MINMAX_COUNT(pNode);
			// This should give the weighted centroid of this node between its two child centroids
		}
		if (pParent)
		{
			// N.B. this code does not know if it is the left or right child of the parent, but either way this works
			Count = MINMAX_COUNT(pParent) += MINMAX_COUNT(pNode);
			for (int i = 0; i < nDim - 1; i++)
			{
				MINMAX_CENTROID(pParent)[i] += MINMAX_COUNT(pNode) * MINMAX_CENTROID(pNode)[i]; // error corrected Feb 7
				// This is the weighted contribution of this leaf to the centroid of the parent
				// It has to be divided by the Count of the parent on the way up
			}
		}
	
		ALNNODE* pRightChild = MINMAX_RIGHT(pNode);
		ALNNODE* pLeftChild = MINMAX_LEFT(pNode);
		// Now get the left and right domain centroids of the children
		float* afltCL = (float*)malloc((nDim - 1) * sizeof(float));
		float* afltCR = (float*)malloc((nDim - 1) * sizeof(float));
		float* afltH  = (float*)malloc((nDim - 1) * sizeof(float)); // This is a point on the hyperplane roughly dividing the points belonging to the two branches
		if (NODE_ISMINMAX(pRightChild))
		{
			for (int i = 0; i < nDim - 1; i++)
			{
				afltCR[i] = MINMAX_CENTROID(pRightChild)[i];
			}
		}
		else
		{
			ASSERT(NODE_ISLFN(pRightChild));
			for (int i = 0; i < nDim - 1; i++)
			{
				afltCR[i] = LFN_C(pRightChild)[i];
			}
		}
		
		if (NODE_ISMINMAX(pLeftChild))
		{
			for (int i = 0; i < nDim - 1; i++)
			{
				afltCL[i] = MINMAX_CENTROID(pLeftChild)[i];
			}
		}
		else
		{
			ASSERT(NODE_ISLFN(pLeftChild));
			for (int i = 0; i < nDim - 1; i++)
			{
				afltCL[i] = LFN_C(pLeftChild)[i];
			}
		}
		
		//Now we use the left and right *domain* centroids of the children of this node to generate other needed items
		MINMAX_THRESHOLD(pNode) = 0;
		for (int i = 0; i < nDim - 1; i++)
		{
			MINMAX_NORMAL(pNode)[i] = MINMAX_CENTROID(pNode)[i] - afltCL[i];
			afltH[i] = afltCR[i] - MINMAX_NORMAL(pNode)[i]; // changed the sign of H, Jan 28
			MINMAX_THRESHOLD(pNode) += -afltH[i] * MINMAX_NORMAL(pNode)[i];
			MINMAX_SIGMA(pNode)[i] += 0.5 * fabs(afltCR[i] - afltCL[i]); // there are already two contributions from children as below
			if (pParent) MINMAX_SIGMA(pParent)[i] += 0.5 * MINMAX_SIGMA(pNode)[i]; // we add half of one of two child values to the parent sigma
		}
		free(afltCR);
		free(afltCL);
		free(afltH);
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		Count = (pNode->DATA.LFN.pSplit)->nCount;
		ALNNODE* pParent = NODE_PARENT(pNode);
		if (pParent)
		{
			MINMAX_COUNT(pParent) += Count;
			for (int i = 0; i < nDim - 1; i++)
			{
				MINMAX_CENTROID(pParent)[i] += Count * LFN_C(pNode)[i]; // This is the weighted contribution of this leaf to the centroid of the parent
																// It has to be divided by the Count of the parent on the way up
				MINMAX_SIGMA(pParent)[i] += 4.0 * sqrt(LFN_D(pNode)[i]); // We have to add the absolute value of the difference of centroid components on the way up
				// one-sided probabilities at various sigmas: sigma = 2.0 0.0227; 2.5 0.0062; 3.0 0.0013; 3.5 0.0002; 3.88 0.0001 Seelalso alnmem.cpp line 497
			}
		}
		if (Count > CanSplitAbove) // added July 29,2020 to try for better classification, adapted for regression use August 10, 2020
		{
			if (LFN_CANSPLIT(pNode))
			{
				float fltSplitLimit;
				float fltPieceSquareTrainError = (pNode->DATA.LFN.pSplit)->fltSqError; // total square error on the piece
				float fltPieceNoiseVariance = (float)Count; // Used when there is no F-test, that is when fltLimit is <= 0, otherwise it tests training MSE < fltLimit
				if (fltMSEorF < 0) // if this is TRUE, we do the F test.
				{
					fltPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->NOISEVARIANCE; // total of noise variance samples
					//the average of noise variance samples estimates the actual noise variance.
					int dofIndex; // get the fltSplitLimit corresponding to the degrees of freedom of the F test
					dofIndex = Count - 2;
					if (Count > 10) dofIndex = 8;
					if (Count > 20) dofIndex = 9;
					if (Count > 30) dofIndex = 10;
					if (Count > 40) dofIndex = 11;
					if (Count > 60) dofIndex = 12; // MYTEST  encourage splitting worked, now it's too much
					fltSplitLimit = afltF_Alpha[dofIndex]; // One can reject the H0 of a good fit with various percentages
					// 90, 75, 50, 35, 25. E.g. 90% says that if the training error is greater than the fltSplitLimit prescribes
					// it is 90% sure that the fit is bad.  A higher percentage needs less training time.
					// Note that when there are few hits on the piece, the fltSplitLimit is larger and 
					// the criterion for fitting well enough is easier to satisfy.
				}
				else
				{
					fltSplitLimit = fltMSEorF;
				}

				if (fltPieceSquareTrainError > fltPieceNoiseVariance * fltSplitLimit)
				{
					// if it has enough samples the piece needs to split
					// The piece must split if it is overdetermined (Count > nDim) and we want to fit all samples perfectly. Yet
					// if it splits and points are not shared, there will be an underdetermined piece.
					// Does this mean sharing samples is obligatory, hence smoothing must be used?

					if (Count > 2) // A very easy to satisfy criterion for splitting any over-determined LFN // CHANGED for classification MYTEST was nDim
					{
						// The piece doesn't fit and needs to split; then training must continue.
						// We want to choose the same way of splitting, max or min, as the parent. This may need some experimentation
						//if (fabs(LFN_SPLIT_T(pNode)) * Count * 20 < fltPieceSquareTrainError) LFN_SPLIT_T(pNode) = 0; MYTEST fix this!!!!!!!!!!!!!!!!!
						if (SplitCount < SplitsAllowed)
						{
							SplitLFN(pALN, pNode); // We split *every* leaf node that reaches this point.
							SplitCount++;
						}
						// We start an epoch with bStopTraining == TRUE, but if any leaf node might still split,
						bStopTraining = FALSE; //  we set it to FALSE and continue to another epoch of training.
					}
					// else the piece doesn't have enough samples to split, but it doesn't fit well -- continue training
				}
				else
				{
					// The piece fits well enough and doesn't need to split or train
					// MYTEST does the following cause more problems than I think???????
					//Change Aug 7, 2020 LFN_FLAGS(pNode) &= ~LF_SPLIT;  // this flag setting prevents further splitting It may be good in high dimensions.
					//LFN_FLAGS(pNode) |= NF_CONSTANT; // Don't train the node any more MYTEST For MNIST
					// Adjoining pieces become responsible for the rest of the fit.
					// This is a queston of speed, so we can test it without stopping sample presentation
				}
			} // end of if LFN_CANSPLIT(pNode)
		} 
	}
}

int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode)
{
	if (bClassify2  && bConvex)
	{
		// This is for convex classification
		return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL);
	}
	else
	{
		// This is splitting for function approximation, we only get here if LFN_CANSPLIT(pNode) is TRUE
		// We split the same as the parent if LFN_SPLIT_T(pNode) == 0
		if ((LFN_SPLIT_T(pNode) == 0) && (NODE_PARENT(pNode) != NULL))
		{

			if (MINMAX_ISMAX(NODE_PARENT(pNode)))
			{
				return ALNAddLFNs(pALN, pNode, GF_MAX, 2, NULL);
					// A max is convex down:  \/,  \_/ etc.
			}
			else
			{
				return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL);
					// A min of several LFNs is like a dome.
			}

		}

		if (LFN_SPLIT_T(pNode) > 0) // This ">" is TRUE if the values of the samples
			// in the training data are higher than the LFN surface some distance from the centroid.
			// This causes the LFN to split into a MAX of two LFNs.
		{
			return ALNAddLFNs(pALN, pNode, GF_MAX, 2, NULL);
			// A max is convex down:  \/,  \_/ etc.
		}
		else
		{
			return ALNAddLFNs(pALN, pNode, GF_MIN, 2, NULL);
			// A min of several LFNs is like a dome.
		}
	}
}
