// NANO program
// This program controls the libaln NANO library.  The user has to input the name of the data file and
// the columns to be used as inputs and the (single) output. The fltMSEorF value if negative, uses an F-test
// to decide whether to break a linear piece (a multi-dimensional hyperplane in general) into two separate
// hyperplanes connected by a max or min operator. If an F-test is not used, a positive level of Mean
// Square Training Error, below which the pieces will not break, can be set.

#ifdef __GNUC__
#include <typeinfo>
#endif
#include "aln.h"
#include "alnpp.h"
#include "datafile.h"
#include "cmyaln.h"
#include "alnpriv.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>  // for high_resolution_clock

static char szInfo[] = "NANO (Noise-Attenuating Neuron Online) program\n"
"Copyright (C)  2019 William W. Armstrong\n"
"Licensed under LGPL\n\n";
// SET TARGET DIGIT HERE !!!!!
float targetDigit =2; //  Here we specify the digit to be recognized. The other digits are all in the second class

BOOL bClassify2 = TRUE; // FALSE produces the usual function learning; TRUE is for two-class classification with a target class
BOOL bConvex = FALSE;  // Used when bClassify2 is TRUE. If bConvex is TRUE, then we do convex classification, i.e. all but one split involves minima.
extern long CountLeafevals; // Global to test optimization
extern BOOL bStopTraining;
// Switches for turning on/off optimizations
BOOL bAlphaBeta = FALSE;
BOOL bDistanceOptimization = FALSE;
float fltTrainErr;
int nMaxEpochs = 20;
int nNumberLFNs = 0;
void setSplitAlpha(ALNDATAINFO* pdata); // This works in the split routine to implement an F-test.
float WeightDecay = 1.0F; //  This is a factor like 0.7599. It is used during classification into two classes when lower weights give better generalization.
float WeightBound = FLT_MAX;
float WeightBoundIncrease = 0.000003; // increase every iteration, e.g. in 200 iterations it can go from 0.0 to 0.006, in 1800 from 0.0006 to .006
int iterations =1900; // Initial iterations
int SplitsAllowed = 75; // These are the two splits when the ALN is set up below, change to 4 if the two extra maxes are used.
int SplitCount = 0;
float* h_A;
float* h_B;
float* h_C;
int MAXLFNS = 512;
const int ELEMENT_N = 197;
BOOL bUseGPU = TRUE; // To be used to allow choice of using GPU or not.
void callCUDA();
void initGPU(void);
void  CheckResult(void);

int main(int argc, char* argv[])
{
	// The first four arguments are the data file name, nDim (number of ALNinputs including one for the output),
	// nMaxEpochs value (training epochs including one epoch of splitting), the fltMSEorF value for stopping splitting, and the weight bound,
	// Example input: think(VeryNoisySine.txt, 2, 20, 0.1, 12) or think(NoisySinCos20000.txt, 3, 3, -90, 0.5).
	// When the program is executed from Visual Studio, the above inputs can be specified in the properties of the think project, Debugging>>
	// Command Arguments with quotes around the file name and spaces instead of commas.
	// The working directory can also be defined there.
	// This version of the program doesn't allow skipping input file columns because it has to handle a large number of ALN inputs.
	std::cout << argc << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6]<< std::endl;
	int argCount = argc;
	int nDim = atoi(argv[2]);

	if (argc != 7 ) // We expect arguments as listed here (argc is not counted among them)
	{
		std::cout << "Bad argument list!\n" << "Usage: " << "Data_file_name nDim nMaxEpochs fltRMSEorF WeightBound Downshift  " << std::endl;
		return 1;
	}
	CDataFile file;
	// INSERT INPUT FILE PATH AND NAME
	if (file.Read(argv[1])) // Important: this file can't have headers; it is all floats
	{
		std::cout << "Reading file " << argv[1] << " Succeeded!" << std::endl;
	}
	else
	{
		std::cout << "Reading file failed!" << std::endl;
		return 1;
	}
	std::cout << "The target for this run is digit " << targetDigit << endl;
	long nTRmaxSamples = file.RowCount(); 
	int nCols = file.ColumnCount();
	nMaxEpochs = atoi(argv[3]);
	// fltRMSEorF (fltMSEorF) can be a positive (root) mean square error limit on training, below which pieces won't split,
	// or a negative number, whose negative is the probability (e.g. -50) as a percent for an F-test.
	// The latter is a local stopping criterion comparable to the global stopping criteria for other neural nets using a validation set.
	// Special cases are -25 -50 -75 where the values come from tables; the other values are approximate, calculated by a simple formula.
	// F-tests for -75 and -90 are likely to stop training before a perfect fit, and -35 and -25 may prolong training and tend to overfit.
	// A negative value causes *much* slower loading of the training buffer.
	float fltRMSEorF = (float) atof(argv[4]); // a negative value indicates use of an F-test to stop splitting, intuition understands fltRMSEorF, but if >0 we use the square
	float fltMSEorF = fltRMSEorF > 0 ? pow(fltRMSEorF, 2) : fltRMSEorF;
	WeightBound = (float) atof(argv[5]); // Plus or minus this value bounds the weights from above and below in cases where all inputs have the same characeristics
	int* ColumnNumber = (int*)malloc(nDim * sizeof(int));
	for (int ALNinput = 0; ALNinput < nDim; ALNinput++)
	{
		ColumnNumber[ALNinput] = ALNinput;
	}
	int nTRcols = 2 * nDim + 1; // The nDim columns of sample data in afltTRdata are extended 
								// to 2 * nDim + 1 total columns for the noise-attenuation tool.
	WeightDecay = atof(argv[6]);

	// Set up the weight array so ALN creation can begin
	h_A = (float*)malloc(nDim * MAXLFNS * sizeof(float));
	h_B = (float*)malloc(nDim * sizeof(float));
	callCUDA();

	initGPU();

	// create ALN 
	std::cout << "Creating ALN";
	CMyAln aln;
	if (!aln.Create(nDim, nDim - 1))
	{
		std::cout << "failed!" << std::endl;
		return 1;
	}
	else
	{
		std::cout << "succeeded!" << std::endl;
	}
	CMyAln* pALN = &aln;
	ALNNODE* pTree = pALN->GetTree(); // The tree is initially just one leaf node
	pALN->SetGrowable(pTree);

	if (bClassify2)
	{
		for (int i = 1; i < nDim; i++)
		{
			LFN_W(pTree)[i] = -WeightBound;
		}
	}





	if (bClassify2)
	{
		const float ConstLevel = 0.95F;  // This must be > 0, e.g. 0.95, to create an interval around 0 where ALN values will lie.
		// The following sets up the special ALN structure for pattern classification into two classes denoted by 1.0 and -1.0
		// Split the root
		ALNAddLFNs(aln, pTree, GF_MIN, 2, NULL);
		ALNNODE* pChildR = MINMAX_RIGHT(pTree);
		ALNNODE* pChildL = MINMAX_LEFT(pTree);
		ASSERT(NODE_ISLFN(pChildR));
		ASSERT(NODE_ISLFN(pChildL));
		// Now split the left child
		ALNAddLFNs(aln, pChildL, GF_MAX, 2, NULL);
		ALNNODE* pGChildR = MINMAX_RIGHT(pChildL);
		ALNNODE* pGChildL = MINMAX_LEFT(pChildL);
		ASSERT(NODE_ISLFN(pGChildR));
		ASSERT(NODE_ISLFN(pGChildL));
		// The next few lines set up a maximum with 0.95 and a minimum with -0.95 for the classification problems
		// This assures that all ALN outputs are in the interval [-0.95, 0.95]
		// A constant input to a minimum node means that only values less than or equal to that constant get through, e.g. at 0.95
		// A constant input to a maximum node means that only values greater than or equal to that constant get through, e.g. at -0.95
		// Set up pChildR to cut off the ALN values above 0.95 using the maximum		
		for (int i = 0; i < nDim; i++)
		{
			LFN_W(pChildR)[i] = 0;
			LFN_C(pChildR)[i] = 0;
			LFN_D(pChildR)[i] = 0.001; // just not zero and not enough to destroy optimization
		}
		LFN_W(pChildR)[0] = ConstLevel;
		LFN_W(pChildR)[nDim] = -1.0; // this is the weight for the output(0 is for the bias weight, there is a shift by one unit)
		LFN_C(pChildR)[nDim - 1] = ConstLevel;
		LFN_FLAGS(pChildR) |= NF_CONSTANT; // Don't allow the LFN to adapt
		LFN_FLAGS(pChildR) &= ~LF_SPLIT; //Don't allow the new right leaf to split
		// Set up the right grandchild pGChild to cut off the ALN values below -0.95 using the minimum
		for (int i = 0; i < nDim; i++)
		{
			LFN_W(pGChildR)[i] = 0;
			LFN_C(pGChildR)[i] = 0;
			LFN_D(pGChildR)[i] = 0.001; // just not zero and not enough to destroy optimization
		}
		LFN_W(pGChildR)[0] = -1.0 * ConstLevel;
		LFN_W(pGChildR)[nDim] = -1.0; // this is the weight for the output(0 is for the bias weight, there is a shift by one unit)
		LFN_C(pGChildR)[nDim - 1] = -1.0 * ConstLevel;
		LFN_FLAGS(pGChildR) |= NF_CONSTANT;
		LFN_FLAGS(pGChildR) &= ~LF_SPLIT; //Don't allow the new right leaf to split
		// The left grandchild should be growable, non-constant, splittable, we set it to be flat at 0.0
		for (int i = 0; i < nDim; i++)
		{
			LFN_W(pGChildL)[i] = 0;
			LFN_C(pGChildL)[i] = 0;
			LFN_D(pGChildL)[i] = 0.001; // just not zero and not enough to destroy optimization
		}
		LFN_W(pGChildL)[0] = 0.0;
		LFN_W(pGChildL)[nDim] = -1.0; // this is the weight for the output(0 is for the bias weight, there is a shift by one unit)
		LFN_C(pGChildL)[nDim - 1] = 0.0;
		/*
		// Option to add one or more maximum nodes to simplify recognition hardware if bConvex is TRUE.
		// It uses one or several domes each doing convex classification to solve a non-convex problem.
		// This shouldn't hurt anything and additional maxes can be added following the recipe below
		ALNAddLFNs(aln, pGChildL, GF_MAX, 2, NULL);
		ALNNODE* pGGChildL = MINMAX_LEFT(pGChildL);
		ASSERT(NODE_ISLFN(pGGChildL));
		ALNAddLFNs(aln, pGGChildL, GF_MAX, 2, NULL);
		*/
	} // End of special code for two-class pattern classification

	ALNREGION* pRegion = pALN->GetRegion(0);
	pRegion->fltSmoothEpsilon = 0;
	// Restrictions on weights (0 is the region -- the region concept is not fully implemented)
	for (int m = 0; m < nDim-1; m++)
	{
		pALN->SetWeightMin(-WeightBound, m, 0);
		pALN->SetWeightMax(WeightBound, m, 0); // MYTEST  try all weights negative
	}
	// This sets up the training buffer of floats afltTRbuffer. The F-test is specified in split_ops.cpp
	ALNDATAINFO* pdata = pALN->GetDataInfo();
	pdata->nTRmaxSamples = nTRmaxSamples;
	pdata->nTRcurrSamples = 0;
	pdata->nTRcols = nTRcols;
	pdata->nTRinsert = 0;
	pdata->fltMSEorF = fltMSEorF;
	// The following sets the alpha for the F-test
	if (fltMSEorF < 0) setSplitAlpha(pdata);
	std::cout << "Loading the data buffer ... please wait" << std::endl;
	// Load the buffer; as the buffer gets each new sample, it is compared to existing samples to create a noise variance tool.
	// During training, once the weights of a piece are known, the noise variance can be estimated.
	float* afltX = (float*)malloc(nDim * sizeof(float));
	int colno;
	long samplesAdded = 0;
	float temp;
	for (long i = 0; i < nTRmaxSamples; i++)
	{
		// get the sample
		for (int j = 0; j < nDim; j++)
		{
			colno = ColumnNumber[j];
			afltX[j] = file.GetAt(i, colno, 0);
		}
		if (bClassify2)
		{
			temp = afltX[nDim - 1]; // Replace the desired value by +1.0 for the target, -1.0 for the others.
			afltX[nDim - 1] = (fabs(temp - targetDigit) < 0.1) ? 1.0 : -1.0;
		}
		pALN->addTRsample(afltX, nDim);  
		samplesAdded++;
	}
	ASSERT(pdata->nTRcurrSamples == samplesAdded);
	std::cout << "ALNDATAINFO: " << "TRmaxSamples = " << pdata->nTRmaxSamples << "  "
		<< "TRcurrSamples = " << pdata->nTRcurrSamples << "  "	<< "TRcols  = " << pdata->nTRcols << "  " << "TRinsert = " << pdata->nTRinsert << std::endl;
	int nDimt2 = nDim *2;
	int nDimt2p1 = nDimt2 + 1;
	float averageNVsample = 0;

	/*
	// Print part of the buffer contents to check
	int nBegin = 0;
	int nLinesPrinted = 15;
	for (int ii = nBegin; ii < nBegin + nLinesPrinted; ii++) // Just a sampling of the file
	{
		for (int jj = 0; jj < nDimt2p1; jj++)
		{
			if (jj == nDim) std::cout << "closest is ";
			if (jj == nDimt2) std::cout << " square dist = ";
			std::cout << pdata->afltTRdata[nDimt2p1 * ii + jj] << " ";
		}
		std::cout << std::endl;
		averageNVsample += 0.5 * pow(pdata->afltTRdata[nDimt2p1 * ii + nDimt2], 2); // MYTEST is this the right entry?
	}
	averageNVsample /= nLinesPrinted;
	std::cout << "\n\n Average noise difference = " << averageNVsample << std::endl;
	// For 6000 MNIST training samples, the average noise distance was on the order of millions.
	*/

	BOOL bJitter = FALSE;
	bStopTraining = FALSE;
	bAlphaBeta = FALSE;
	bDistanceOptimization = FALSE;
	BOOL bFirstTime = TRUE;
	float fltLearnRate = 0.2F;
	float fltMinRMSE = 0.00000001F;// This is set small and not very useful.  fltRMSEorF is used now to stop training.
	int nNotifyMask = AN_TRAIN; // required callbacks for information or insertion of data. You can OR them together with |
	ALNNODE** ppActiveLFN = NULL;
	CDataFile ExtendTR, NormalReplaceTR, PurgeReplaceTR;
	long countReplace = 0;
	int colsOut = (bClassify2 ? nCols + 2 : nCols + 1);
	ExtendTR.Create(file.RowCount(), colsOut); // Create a data file to hold the results
	NormalReplaceTR.Create(file.RowCount(), file.ColumnCount());
	long rowsReplace = file.RowCount();
	std::cout << std::endl << "Starting training " << std::endl;
	// Record start time
	auto start_training = std::chrono::high_resolution_clock::now();
	do
	{
		for (int iteration = 1; iteration <= iterations; iteration++)
		{
			if (iteration == 1 || iteration % 5 == 0) std::cout << "\nIteration " << iteration << " (of " << iterations << " )   ";
			flush(std::cout);
			auto start_iteration = std::chrono::high_resolution_clock::now();
			// VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  Training!

			//pTree = pALN->GetTree();
			if (!pALN->Train(nMaxEpochs, fltMinRMSE, fltLearnRate, bJitter, nNotifyMask))
			{
				std::cout << " Training failed!" << std::endl;
				flush(std::cout);
				return 1;
			}

			// Every iteration increases the WeightBound
			for (int m = 0; m < nDim - 1; m++) // The slopes of the LFNs discriminating classes are not known in advance, so this step reduces the effect of the initial choice.
			{
				pALN->SetWeightMin(-WeightBound, m, 0);
				pALN->SetWeightMax(WeightBound, m, 0);
			}
			if(WeightBound < 0.05) WeightBound += WeightBoundIncrease; // After this, it's manual.
			if(SplitsAllowed < nDim || iteration%10 == 1 ) SplitsAllowed++;
			std::cout << "/" << SplitsAllowed << " ";

			if ( bFirstTime && iteration == 290 ) // Optimization is not needed when there are few leaf nodes, e.g. < 256
			{
				bAlphaBeta = TRUE; // Once these optimizatons are switched on, they stay on. (This is after the first split)  MYTEST
				bDistanceOptimization = TRUE;
				// fltRMSEorF = 10.0F; // When this is applied, there should be no more splitting!
				//pdata->fltMSEorF = fltMSEorF = fltRMSEorF > 0 ? pow(fltRMSEorF, 2) : fltRMSEorF;
				//WeightBound = 0.006;
				//WeightDecay = 1.02;
				bFirstTime = FALSE;
			}
		
			// MYTEST

			//if (iteration > 200) bConvex = FALSE;


			auto finish_iteration = std::chrono::high_resolution_clock::now();
			std::chrono::duration<float> elapsed0 = finish_iteration - start_training;
			std::chrono::duration<float> elapsed1 = finish_iteration - start_iteration;
			if (iteration == 1 || iteration % 5  == 0)
			{
				std::cout << " Duration " << elapsed1.count() << " Elapsed time " << elapsed0.count() << " seconds" << std::endl;
				flush(std::cout);
			}
		}
		std::cout << "RMSEorF = " << fltRMSEorF << " Weight bound = " << WeightBound << " Weight decay  = " << WeightDecay  << endl;

		float entry;
		long countcorrect = 0;
		long counterrors = 0;
		countReplace = 0;
		float desired;
		ALNDATAINFO* pdata2 = pALN->GetDataInfo();
		ofstream BadTrainImages("BadTrainImages.txt", ios_base::trunc);
		if (!BadTrainImages.good())
		{
			cerr << "File BadTrainImages.txt could not be opened." << endl;
		}
		BadTrainImages << "The target digit for the following images was " << targetDigit << endl;
		for (long i = 0; i < nTRmaxSamples; i++)
		{
			float totalIntensity = 0;;
			// Copy the row of the input data file
			for (int j = 0; j < nDim; j++) // we include all columns
			{
				entry = file.GetAt(i, j, 0);
				ExtendTR.SetAt(i, j, entry, 0);
				totalIntensity += entry;
			}
			if (totalIntensity < 1) totalIntensity = 1;
			NormalReplaceTR.SetAt(i, nDim - 1, -1024, 0); // Every row in the replacement file is first invalidated in the desired output column, some to be over-written later by correct samples
			// Evaluate the ALN on this row
			for (int k = 0; k < nDim - 1; k++) // Get the domain coordinates
			{
				afltX[k] = file.GetAt(i, k, 0);
				NormalReplaceTR.SetAt(countReplace, k, afltX[k] *98000.0F / (totalIntensity + .001), 0); // This copies the file with normalized intensities, but bad rows will be omitted.
			}
			afltX[nDim - 1] = file.GetAt(i, nDim - 1, 0);
			NormalReplaceTR.SetAt(countReplace, nDim -1, afltX[nDim - 1], 0); //fix the label entry
			desired = afltX[nDim - 1]; // Store the desired output according to the training file
			afltX[nDim - 1] = -250; // Be sure the desired output will not help in the computation
			entry = pALN->QuickEval(afltX, ppActiveLFN); // get the ALN-computed value for entry into the output file
			ExtendTR.SetAt(i, nDim, entry, 0); // the first additional column on the right is the ALN output
			float correctClass;
			if(bClassify2)
			{
				correctClass = (abs(desired - targetDigit) < 0.1) ? 1 : -1;
				if (entry * correctClass > 0)
				{
					countcorrect++; // This is for the recognition of a single digit e.g. "9" vs all the others.
					ExtendTR.SetAt(i, nCols + 1, 0, 0); // the second additional column on the right points out mistakes
					countReplace++; // If this is not incremented, the row may be over-written by a good sample later.
				}
				else
				{
					counterrors++;
					ExtendTR.SetAt(i, nCols + 1, 9999.9F, 0); // the second additional column on the right points out mistakes
					for (int m = 0; m < nDim -1; m++)
					{
						entry = ExtendTR.GetAt(i, m, 0); // get the domain values
						BadTrainImages << ( entry > 50.0F ? 'm' : ' ');
						if (m % 14 == 13) BadTrainImages << endl;
					}
					BadTrainImages << "The label on the above digit was " << (int)desired << endl;
					countReplace++; // Incrementing this gives the normalized version of the whole training file
				}
			}
		}
		// OUTPUT OF RESULTS
		std::cout << std::endl << "Closing BadTrainImages file with " << counterrors << " images " << std::endl;
		BadTrainImages.close();
		std::cout << std::endl << "Writing ExtendedTR.txt file " << std::endl;
		if (!ExtendTR.Write("ExtendedTR.txt"))
		{
			std::cout << "ExtendedTR.txt was locked, look at ExtendedTR2.txt instead." << std::endl;
			ExtendTR.Write("ExtendedTR2.txt");
		}
		std::cout << std::endl << "Writing NormalReplaceTR.txt file " << std::endl;
		if (!NormalReplaceTR.Write("NormalReplaceTR.txt"))
		{
			std::cout << "NormalReplaceTR.txt was locked, look at NormalReplaceTR2.txt instead." << std::endl;
			ExtendTR.Write("NormalReplaceTR2.txt");
		}
		//std::cout << "Writing trained aln file... please wait" << std::endl;
		//pALN->Write("NANOoutput.aln");
		if (bClassify2)
		{
			std::cout << "The number of correct classifications was " << countcorrect << " or " << 100.0 * countcorrect / file.RowCount() << " percent " << std::endl;
			std::cout << "The number of errors was " << counterrors << std::endl;
		}
		std::cout << "Continue training?? Enter number of additional iterations (multiple of 10), 0 to quit, 1 to modify." << 
			"\nModifications: iterations, RMSEorF, WeightBound, WeightDecay" << std::endl;
		std::cin >> iterations;
		if (iterations == 1)
		{
			std::cin >> iterations >> fltRMSEorF >> WeightBound >> WeightDecay;
			float fltMSEorF = fltRMSEorF > 0 ? pow(fltRMSEorF, 2) : fltRMSEorF;
			pdata->fltMSEorF = fltMSEorF;
			for (int m = 0; m < nDim - 1; m++)
			{
				pALN->SetWeightMin(-WeightBound, m, 0);
				pALN->SetWeightMax(WeightBound, m, 0);
			}
		}
	} while (iterations > 0);
	std::cout << "Do you want to replace incorrect or incorrectly labeled images? (y or n)" << std::endl;
	char dummy;
	std::cin >> dummy;
	float entry;
	long replaceCount = 0;
	float totalIntensity = 0;
	// The PurgeReplacementTR file will not be normalized, just purged of dubious training samples.
	PurgeReplaceTR.Create(file.RowCount(), file.ColumnCount());
	if (dummy == 'y'|| dummy == 'Y')
	{
		for (long i = 0; i < nTRmaxSamples; i++)
		{
			if (ExtendTR.GetAt(i, nDim + 1, 0) < 9999.9F)
			{
				for (int m = 0; m < nDim; m++)
				{
					entry = ExtendTR.GetAt(i, m, 0);
					PurgeReplaceTR.SetAt(replaceCount, m, entry, 0);
				}
				replaceCount++;
			}
			else
			{
				totalIntensity = 0;
				for (int m = 0; m < nDim; m++)
				{
					entry = ExtendTR.GetAt(i, m, 0);
					totalIntensity += entry;
					std::cout << (entry > 50.0F ? 'm' : ' ');
					if (m % 14 == 13) std::cout << endl;
				}
				std::cout << "The label on the above digit was " << (int) entry << " Average intensity = " << totalIntensity/196.0 << std::endl;
				std::cout << "Should this image and label be kept?  (y or n)" << std::endl;
				std::cin >> dummy;
				if (dummy == 'y' || dummy == 'Y')
				{
					for (int m = 0; m < nDim; m++)
					{
						entry = ExtendTR.GetAt(i, m, 0);
						PurgeReplaceTR.SetAt(replaceCount, m, entry, 0);
					}
					replaceCount++;
				}
			}
		}
		std::cout << "The number of rows in the replacement file is " << replaceCount << endl;
		if (!PurgeReplaceTR.Write("PurgeReplaceTR.txt"))
		{
			std::cout << "PurgeReplaceTR.txt was locked, look at PurgeReplaceTR2.txt instead." << std::endl;
			ExtendTR.Write("PurgeReplaceTR2.txt");
		}

	}
	/*
	// Speed test
	std::cout << std::endl << "Starting speed test... please wait" << std::endl;
	CountLeafevals = 0;
	// Record start time
	auto start = std::chrono::high_resolution_clock::now();
	for (long i = 0; i < nTRmaxSamples; i++)
	{
		for (int j = 0; j < nDim; j++)
		{
			afltX[j] = trainfile.GetAt(i, j, 0);
		}
		afltX[nDim - 1] = 0; // This desired output value is thus not input to QuickEval
		ALNoutput = pALN->QuickEval(afltX, ppActiveLFN);
	}
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
	// Now do the same loop without doing the QuickEval
	// Record start time
	start = std::chrono::high_resolution_clock::now();
	for (long i = 0; i < nTRmaxSamples; i++)
	{
		for (int j = 0; j < nDim; j++)
		{
			afltX[j] = trainfile.GetAt(i, j, 0);
		}
		afltX[nDim - 1] = 0; // This desired output value is thus not input to QuickEval
		// ALNoutput = pALN->QuickEval(afltX, ppActiveLFN);
	}
	// Record end time
	finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed2 = finish - start;
	std::cout << std::endl << "Elapsed time for " << nTRmaxSamples <<  " evaluations: " << elapsed.count() - elapsed2.count() << " seconds" << std::endl;
	std::cout << "Leaf evaluations during testing " << CountLeafevals << std::endl;
	*/

	//*********************************************************************************
	// For MNIST
	// Now go on to the recognition phase for testing
	std::cout << "Starting evaluation on the test file ... please wait." << std::endl;
	ofstream BadTestImages("BadTestImages.txt", ios_base::trunc);
	if (!BadTestImages.good())
	{
		cerr << "File BadTestImages.txt could not be opened." << endl;
	}
	BadTestImages << "The target digit for the following images was " << targetDigit << endl;
	CDataFile testfile;
	if (testfile.Read("Normalized500MNIST_NANO_TestFile.txt"))
	{
		std::cout << "Reading file succeeded!" << std::endl;
	}
	else
	{
		std::cout << "Reading file failed!" << std::endl;
	}
	int nCorrect = 0;
	int nWrong = 0;
	std::cout << std::endl << "Starting evaluation" << std::endl;
	bAlphaBeta = FALSE; // Optimizations are not done in the interest of accuracy.
	bDistanceOptimization = FALSE;
	float DesiredALNOutput, ALNoutput, label;
	for (long i = 0; i < 10000; i++)
	{
		for (int j = 0; j < nDim; j++) // Get the domain values
		{
			afltX[j] = testfile.GetAt(i, j, 0);
		}
		DesiredALNOutput = (fabs(afltX[nDim - 1] - targetDigit) < 0.1 ? 1 : -1);
		label = afltX[nDim - 1];
		afltX[nDim - 1] = -250.0; // Put in an incorrect value here (to show there is no cheating).
		ALNoutput = pALN->QuickEval(afltX, ppActiveLFN);
		if (ALNoutput * DesiredALNOutput > 0) 
		{
			nCorrect++;
		}
		else
		{
			nWrong++;
			for (int j = 0; j < nDim -1; j++) // Get the domain values
			{
				BadTestImages << ((afltX[j] > 50.0F) ? 'M' : ' ');
				if (j % 14 == 13) BadTestImages << endl;
			}
			BadTestImages << "The label on the above digit was " << (int)label << endl;
		}
	}
	std::cout << "Correct: " << nCorrect << " Wrong: " << nWrong << endl;
	BadTestImages.close();
	free(afltX);
	free(pdata->afltTRdata);
	pALN->Destroy();
	//aln.Destroy(); which one to use?
}

