// decayweights.cpp

#include "aln.h"
#include "alnpriv.h"
extern BOOL bClassify2;

ALNIMP void ALNAPI DecayWeights(const ALNNODE* pNode, const ALN* pALN, float WeightBound, float WeightDecay)
{
	// THis routine should not be called except for two-class classification
	if (NODE_ISMINMAX(pNode))
	{
		DecayWeights(MINMAX_RIGHT(pNode), pALN, WeightBound, WeightDecay);
		DecayWeights(MINMAX_LEFT(pNode), pALN, WeightBound, WeightDecay);
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (NODE_ISCONSTANT(pNode))return;
		float* pW = LFN_W(pNode);
		int nDim = pALN->nDim;
		if (bClassify2)
		{
			// we rotate the LFN if the output centroid is not too close to samples.
			float C_output = LFN_C(pNode)[nDim - 1];
			if (C_output < -0.9 || C_output > 0.9) return;
		}
		for (int i = 1; i < nDim -1; i++)
		{
			pW[i] *= WeightDecay;
			pW[i] = max(min( WeightBound, pW[i]), -WeightBound);
		}

	}
}
