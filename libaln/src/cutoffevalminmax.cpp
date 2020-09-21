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

// cutoffevalminmax.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Switches for turning on/off optimizations
extern BOOL bAlphaBeta;
extern BOOL bDistanceOptimization;


///////////////////////////////////////////////////////////////////////////////
// minmax node specific eval - returns distance to surface
//  - non-destructive, ie, does not change ALN structure
// NOTE: cutoff always passed on stack!

float ALNAPI CutoffEvalMinMax(const ALNNODE* pNode, const ALN* pALN,
	const float* afltX, CEvalCutoff cutoff,
	ALNNODE** ppActiveLFN)
{
	ASSERT(NODE_ISMINMAX(pNode));

	// We use the sample counts and centroids of the child nodes to generate a hyperplane H roughly separating the child samples.
	// The branch to take first during an evaluation is the one representing the side of H afltX lies on.
	const ALNNODE* pChild0;
	const ALNNODE* pChild1;
	int nDim = pALN->nDim;
	// We take the dot product of the normal vector, in direction left centroid to right centroid, with afltX - H to find the branch which goes first.
	// Note that this handles the case where the normal is zero length as after a split and child centroids are equal.
	float dotproduct = 0;
	for (int i = 0; i < nDim - 1; i++)
	{
		dotproduct += MINMAX_NORMAL(pNode)[i] * afltX[i];
	}
	dotproduct += MINMAX_THRESHOLD(pNode); // add constant stored for speed; see split_ops.cpp
	if (dotproduct > 0)
	{
		pChild0 = MINMAX_RIGHT(pNode);
		pChild1 = MINMAX_LEFT(pNode);
	}
	else
	{
		pChild1 = MINMAX_RIGHT(pNode);
		pChild0 = MINMAX_LEFT(pNode);
	}

	/*
	// set first child -- Monroe's code prioritizes the past active node, see also buildcutoffroute.cpp
	if (MINMAX_EVAL(pNode))
		pChild0 = MINMAX_EVAL(pNode);
	else
		pChild0 = MINMAX_LEFT(pNode);

	// set next child
	if (pChild0 == MINMAX_LEFT(pNode))
		pChild1 = MINMAX_RIGHT(pNode);
	else
		pChild1 = MINMAX_LEFT(pNode);
	*/
	// get reference to region for this node
	const ALNREGION& region = pALN->aRegions[NODE_REGION(pNode)];
	float fltDist;

	// eval first child
	ALNNODE* pActiveLFN0;
	float flt0 = CutoffEval(pChild0, pALN, afltX, cutoff, &pActiveLFN0);

	// see if we can cutoff...
	if (bAlphaBeta && Cutoff(flt0, pNode, cutoff))
	{
		*ppActiveLFN = pActiveLFN0;
		return flt0;
	}
	// Check if second child is too far from sibling in any axis j; in which case
	// the tree of the second child is cut off and we return the value from the first child

	ASSERT(pChild1);
	if (bDistanceOptimization && MINMAX_SIGMA(pChild1))
	{
		for (int j = 0; j < nDim - 1; j++)
		{
			if (fabs(afltX[j] - MINMAX_CENTROID(pChild1)[j]) >  MINMAX_SIGMA(pChild1)[j]) // Do we need a safety factor?
			{
				return flt0;
			}
		}
	}

	// eval second child
	ALNNODE* pActiveLFN1;
	float flt1 = CutoffEval(pChild1, pALN, afltX, cutoff, &pActiveLFN1);

	// calc active child and distance without using CalcActiveChild()
	if((MINMAX_ISMAX(pNode) > 0) == (flt1 > flt0)) // int MINMAX_ISMAX is used as a bit-vector!
	{
		*ppActiveLFN = pActiveLFN1;
		fltDist = flt1;
	}
	else
	{
		*ppActiveLFN = pActiveLFN0;
		fltDist = flt0;
	}
  
  return fltDist;
}

