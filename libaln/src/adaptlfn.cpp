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

// adaptlfn.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern BOOL bClassify2;
// LFN specific adapt

void ALNAPI AdaptLFN(ALNNODE* pNode, ALN* pALN, const float* afltX, 
                     float fltResponse, BOOL bUsefulAdapt, const TRAINDATA* ptdata)
{
	ASSERT(NODE_ISLFN(pNode));
	ASSERT(LFN_ISINIT(pNode));
	ASSERT(LFN_VARMAP(pNode) == NULL);      // var map not yet supported
	ASSERT(LFN_VDIM(pNode) == pALN->nDim);  // no different sized vectors yet
	ASSERT(NODE_ISEVAL(pNode));
	// constraining region
	ASSERT(NODE_REGION(pNode) >= 0 && NODE_REGION(pNode) < pALN->nRegions);
	ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];
	// get dimension
	int nDim = pALN->nDim;
	// count adapts
	NODE_RESPCOUNT(pNode)++; // this may not be used
	if (NODE_ISCONSTANT(pNode))	return; // no adaptation of this constant active LFN
	// This procedure adapts the nDim centroid values and the nDim - 1 weight values so that the changes
	// are likely to be individually small but together correct about fraction fltLearnrate of the error
	// of the surface w. r. t. the training point. The error (fltError) is the distance in the output axis 
	// from the sample desired value to the (sometimes smoothed) ALN function surface
	// (flt Error is positive when the ALN function is greater than the sample value).
	// We have to think of the error as that of the *surface*, not the data point!
	float fltError = ptdata->fltGlobalError; //This is the error just for this training point
	// notify begining of LFN adapt
	if (CanCallback(AN_LFNADAPTSTART, ptdata->pfnNotifyProc, ptdata->nNotifyMask))
	{
		LFNADAPTINFO lai;
		lai.afltX = afltX;
		lai.pLFN = pNode;
		lai.fltError = fltError;
		lai.fltResponse = fltResponse;
		Callback(pALN, AN_LFNADAPTSTART, &lai, ptdata->pfnNotifyProc, ptdata->pvData);
	}

	ASSERT(LFN_SPLIT(pNode) != NULL);
	LFN_SPLIT_COUNT(pNode)++;
	LFN_SPLIT_SQERR(pNode) += fltError * fltError * fltResponse;
	LFN_SPLIT_RESPTOTAL(pNode) += fltResponse;

	// copy LFN vector pointers onto the stack for faster access
	float* afltW = LFN_W(pNode);	    // weight vector( must be shifted to use the same index as afltC, afltX 
	afltW++; // afltW[0] contains a value that, for the sake of speed, avoids implicating the centroids in evaluation
	float* afltC = LFN_C(pNode);			// centroid vector
	float* afltD = LFN_D(pNode);			// average square dist of sample from centroid vector -- a sample variance
											// This is used to adapt weights and also for optimization

	// set up the learning rates
	float learnBoost = 1.0F;
	if (bClassify2 && fltError > 0.05) learnBoost = 0.1111111F;   // This attenuates the downward pull by the 9 non-target classes
	// If response is < 1 because of a smoothed minmax node, then the adaptation has less effect.
	// We need two learning rates.  The first is for the centroids and weights, which divides by 2*nDim -1, thus
	// putting them on an equal footing with respect to correcting a share of the error.
	float fltLearnRate = ptdata->fltLearnRate;
	float fltLearnRespParam = learnBoost * fltLearnRate * region.fltLearnFactor / (double)(2 * nDim - 1); // The region factor should be optimized out

	// ADAPT THE CENTROID FOR OUTPUT BY EXPONENTIAL SMOOTHING
	// L is the value of the affine function of the linear piece at the input components of X
	// V is the value of the sample X[nDim - 1]
	// fltError = L - V,   N.B. if L is greater than the sample the error is positive
	// We neglect the fillet if any and make the average V the target for afltC[nDim - 1]
	afltC[nDim -1] += (afltX[nDim -1] - afltC[nDim - 1]) * fltLearnRespParam;

	// ADAPT CENTROID AND WEIGHT FOR EACH INPUT VARIABLE  
	float fltXmC = 0;
	float fltBend = 0;
	for (int i = 0; i < nDim -1; i++) //Skip the output centroid and weight at nDim - 1. (The output weight is always -1)
	{
		// get pointer to variable constraints
		ALNCONSTRAINT* pConstr = GetVarConstraint(NODE_REGION(pNode), pALN, i);
		// skip any variables of constant monotonicity; W is constant and X is irrelevant
		if (pConstr->fltWMax == pConstr->fltWMin) continue;
		// Compute the distance of X from the old centroid in axis i
		fltXmC = afltX[i] - afltC[i];
		// UPDATE VARIANCE BY EXPONENTIAL SMOOTHING
		// We adapt this first so the adaptation of the centroid to this input will not affect it
		ASSERT(afltD[i] >= 0);
		afltD[i] += (fltXmC * fltXmC - afltD[i]) * fltLearnRate; // This learning rate is not involved in correcting fltError
		// afltD[i] is not allowed to go to 0.
		if (afltD[i] < 0.0000001)afltD[i] = 0.0000001; // It could happen for an under-determined LFN with one sample on it.

		// UPDATE THE CENTROID BY EXPONENTIAL SMOOTHING
		afltC[i] += fltXmC * fltLearnRespParam; // the centroid is moved part way to X
		// ADAPT WEIGHTS
		afltW[i] -= fltError * fltLearnRespParam * fltXmC / afltD[i]; // Note how this preserves units for weight: output unit/input unit of this axis
		// Bound the weight
		afltW[i] = max(min(pConstr->fltWMax, afltW[i]), pConstr->fltWMin);

		// COLLECT DATA FOR LATER SPLITTING THIS PIECE:
		// We analyze the errors of sample value minus ALN value V - L = -fltError (N.B. minus) on the piece which are
		// further from and closer to the centroid than the stdev of the points on the piece along the current axis.
		// If the V - L  is positive (negative) away from the centre compared to the error closer to the centre,
		// then we need a split of the LFN into a MAX (MIN) node.
		// Even if the fit is extremely bad, we need to know the convexity to split the piece in the right direction, MIN or MAX.
		// To capture the information for all the samples on the piece and all the nDim -1 domain directions, we use exponential smoothing.
		if (LFN_CANSPLIT(pNode))
		{
			fltBend = (fltXmC * fltXmC - afltD[i] > 0) ? -(fltXmC * fltXmC - afltD[i]) * fltError : (fltXmC * fltXmC - afltD[i]) * fltError;
			LFN_SPLIT_T(pNode) += (fltBend - LFN_SPLIT_T(pNode)) * fltLearnRate;
		}
	} // end loop over all nDim - 1 domain dimensions
	// compress the weighted centroid info into W[0]       
	float * const pfltW0 = LFN_W(pNode);
	*pfltW0 = afltC[nDim - 1];
	for (int i = 0; i < nDim - 1; i++)
	{
		*pfltW0 -= afltW[i] * afltC[i]; // here the afltW pointer is still shifted up by one float
	}
	// notify end of LFN adapt
	if (CanCallback(AN_LFNADAPTEND, ptdata->pfnNotifyProc, ptdata->nNotifyMask))
	{
		LFNADAPTINFO lai;
		lai.afltX = afltX;
		lai.pLFN = pNode;
		lai.fltError = fltError;
		lai.fltResponse = fltResponse;
		Callback(pALN, AN_LFNADAPTEND, &lai, ptdata->pfnNotifyProc, ptdata->pvData);
	}
}
