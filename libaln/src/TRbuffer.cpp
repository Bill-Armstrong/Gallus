// ALN Library example
// file buffer.cpp
// Copyright (C) 2019 William W. Armstrong.
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

#include "aln.h"
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
// include files
#include <alnpp.h>
#include <malloc.h>
#include <array>
#include "cmyaln.h"
//#include "alnextern.h"
//#include "alnintern.h"

long nCurrSize = 0;
long TRindex = 0;
const long nMaxBuff =20000;
double* adblYdiff;
double* adblYtemp;

static void ALNAPI addTRsample(CAln* pALN, double* adblX, const int nDim)
{
	// A sample can always be added; the sample at TRindex has to go if the buffer is at max size
	// Start by comparing the new sample to all in the buffer to see (N.B. "closest" means among the closest if it is not unique) 
	// 1. what other samples it is closest to
	// 2. what other sample is closest to it.
	// A TRbuff entry has the sample (nDim entries) followed starting at offset nDim by the differences
	// of the closest to that sample minus the sample itself (nDim entries) followed by the squared distance
	// between the two samples in the domain (one entry).
	// We need to keep track in the ALN of the current buffer size, nCurrSize,
	// the maximum allowed size, nMaxBuff, and the current insertion point , nTRinsert, for the next sample
	// as well as the pointer to the buffer itself (double*) aTRbuff.
	double value, sum;
	int nDimt2 = 2 * nDim;
	int add4output = nDim - 1;
	int add4diffstart = nDim;
	int add4outputdiff = add4diffstart + nDim - 1;
	int add4squaredist = add4outputdiff + 1;
	ALNDATAINFO* thisDataInfo = pALN->GetDataInfo();
	double* adblTRdata = thisDataInfo->adblTRdata;
	int nTRmaxSamples = thisDataInfo->nTRmaxSamples;
	int nTRcurrSamples = thisDataInfo->nTRcurrSamples;
	int nTRcols = thisDataInfo->nTRcols;
	int nTRinsert = thisDataInfo->nTRinsert;
	double dblMSEorF = thisDataInfo->dblMSEorF;
	if (nTRcurrSamples == 0) // for empty or one sample buffer only
	{

		// This sets up the training buffer TRbuffer
		adblYdiff = (double*)malloc(nDim * sizeof(double)); // stores the differences of adblX to the current buffer sample
		adblYtemp = (double*)malloc(nDim * sizeof(double)); // stores the buffer sample differences which are currently the closest to adblX
		// First we fill TRbuff with 0's
		for (long i = 0; i < nMaxBuff; i++)
		{
			int nDimt2p1ti = (nDimt2 + 1) * i; // this is the base of the i-th sample in TRbuff
			for (int j = 0; j < nDimt2 - 2; j++)
			{
				aTRbuff[nDimt2p1ti +j] = 0;
			}
            aTRbuff[nDimt2p1ti + nDimt2] = std::numeric_limits<double>::max(); // the value of the sample
		}
		// insert the sample at TRindex 0 and store a fictitious closest far away to be replaced at the next insertion
		for (int j = 0; j < nDim - 1; j++) // domain components only
		{
			aTRbuff[j] = adblX[j];
			aTRbuff[nDim + j] = 0;
		}
		aTRbuff[nDim - 1] = adblX[nDim - 1];
		aTRbuff[nDimt2 -1] = std::numeric_limits<double>::max();
		nTRcurrSamples++;
		// update the buffer in the ALN
		pALN->SetDataInfo(adblTRdata, nTRmaxSamples, nTRcurrSamples, nTRcols, nTRinsert, dblMSEorF);
		nCurrSize = 1;
		return;
	}

	for (long i = 0; i < nCurrSize; i++)
	{
		int nDimt2ti = nDimt2 * i; // this is the base of the i-th sample in TRbuff
		// Get the current square distance of the new sample from sample i in TRbuff
		sum = 0;
		for (int j = 0; j < nDim - 1; j++) // domain components only
		{
			value = aTRbuff[i * nDimt2 + j];
			adblYdiff[j] = adblX[j] - value;
			sum += adblYdiff[j] * adblYdiff[j]; // TBD later: different axes will contribute different weights
		}
		// Get the stored closest square distance from the i-th sample
		value = aTRbuff[nDimt2ti + nDimt2 - 1];
		if (sum < value)
		{
			// If the sum is less than for the existing closest
			// replace the closest to sample i with the new differences and update the distance
			for (int j = 0; j < nDim - 1; j++)
			{
				aTRbuff[nDimt2ti + nDim + j] = adblYdiff[j];
			}
			aTRbuff[nDimt2ti + nDimt2 - 1] = sum;
		}
		// If the ith sample in TRfile is closer to adblX than anything was earlier then update the temporary values
		value = adblYtemp[nDim - 1]; 
		if (sum < value)
		{
			// Store sample i in the temporarily closest sample for adblX.
			for (int j = 0; j < nDim - 1; j++)
			{
				adblYtemp[j] = adblX[j];
			}
			adblYtemp[nDim - 1] = sum;
		}
	}
	// Insert the new sample and its difference vector and closest distanc at TRindex
	long nDimt2tTRindex = nDimt2 * TRindex;
	for (int j = 0; j < nDim - 2; j++) // domain components only
	{
		aTRbuff[nDimt2tTRindex +j ] = adblX[j];
		aTRbuff[nDimt2tTRindex + nDim + j] = adblYtemp[j];
	}
	aTRbuff[nDimt2tTRindex + nDim -1] = adblX[nDim -1];
	aTRbuff[nDimt2tTRindex + nDimt2 -1] = sum;
	if (nCurrSize < nMaxBuff)
	{
		nCurrSize++;
	}
	TRindex++;
	TRindex %= nTRcurrSamples;
	pALN->SetDataInfo(adblTRdata, nTRmaxSamples, nTRcurrSamples, nTRcols, nTRinsert, dblMSEorF);
	return;
}



/*
void SetTrainParams(int nDim)
{
	// Set up the approximation ALN
	pBaseNeuron = new CMyAln; // NULL initialized ALN
	if (!pBaseNeuron->Create(nDim, nDim - 1))
	{
		fprintf(fpProtocol, "ALN creation failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	// Now make the tree growable
	if (!pBaseNeuron->SetGrowable(pBaseNeuron->GetTree()))
	{
		fprintf(fpProtocol, "Setting ALN growable failed!\n");
		fflush(fpProtocol);
		exit(0);
	}

	nMaxEpochs = 20; // The number of passes through the data without splitting LFNs.
	dblMinRMSE = 1e-20; // Stops training when the error is tiny.
	dblLearnRate = 0.2;


	// Set constraints on variables for ALN
	// NB The following loop excludes the output variable of the ALN, index nDim -1.
	for (int m = 0; m < nDim - 1; m++)
	{
		pBaseNeuron->SetEpsilon(adblEpsilon[m], m);
		if (adblEpsilon[m] == 0)
		{
			fprintf(fpProtocol, "Stopping: Variable %d appears to be constant. Try removing it.\n", m);
			fflush(fpProtocol);
			exit(0);
		}
		// The minimum value of the domain is a bit smaller than the min of the data points
		// in TVfile, and the maximum is a bit larger.
		pBaseNeuron->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m], m);
		pBaseNeuron->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m], m);

		// A rough value for the range of output (for a uniform dist.) divided by the
		// likely distance between samples in axis m.
		pBaseNeuron->SetWeightMin(-pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);
		pBaseNeuron->SetWeightMax(pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);

		// Impose the a priori bounds on weights which have been set by the user
		if (dblMinWeight[m] > pBaseNeuron->GetWeightMin(m))
		{
			pBaseNeuron->SetWeightMin(dblMinWeight[m], m);
		}
		if (dblMaxWeight[m] < pBaseNeuron->GetWeightMax(m))
		{
			pBaseNeuron->SetWeightMax(dblMaxWeight[m], m);
		}
	}
	
	(pBaseNeuron->GetRegion(0))->dblSmoothEpsilon = 0;
	fprintf(fpProtocol, "The smoothing for training approximation is %f\n", 0.0);
}
*/

/*
void ALNAPI Train2Split(CMyAln pBaseNeuron) // routine
{
	

	BOOL bStopTraining = FALSE; // Set to TRUE at the start of each epoch in alntrain.cpp. 
	// Set FALSE by any piece needing more training. 
	// nNumberLFNs = 1;  // initialize at 1 MYTEST where should this be??
	  // Tell the training algorithm the way to access the data using fillvector
	const double* adblData = TRfile.GetDataPtr(); // This is where training gets samples.
	// The advantage of giving the pointer adblData instead of NULL is that training permutes
	// the order of the samples and goes through all samples exactly once per epoch.
	dblLimit = -1.0; // Negative to split leaf nodes according to an F test;
	// positive to split if training MSE > dblLimit.
	pBaseNeuron.SetDataInfo(nCurrSize, nDim, adblData, NULL, dblLimit);
	// Train the ALN nMaxEpochs, after the last of which pieces are allowed to break
	if (!pBaseNeuron.Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
	{
		fprintf(fpProtocol, "Training failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	if (bStopTraining == TRUE)
	{
		fprintf(fpProtocol, "All leaf nodes have stopped changing!\n");
		fflush(fpProtocol);
	}
}



void ALNAPI cleanup() // routine
{
	// cleanup what was allocated for training in approximation
	pBaseNeuron->Destroy();
	TRfile.Destroy();
	TRdiff.Destroy();

	// clean up the rest
  fclose(fpOutput);
  fclose(fpProtocol);
}
*/