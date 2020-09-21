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

// alnpp.h

#ifndef __ALNPP_H__
#define __ALNPP_H__

#ifndef __ALN_H__
#include "aln.h"
#endif

#include <assert.h>

/////////////////////////////////////////////////////////////////////////////
// library files (Microsoft compilers only)

#ifndef ALNPP_NOFORCE_LIBS

#endif  // ALNPP_NOFORCE_LIBS

/////////////////////////////////////////////////////////////////////////////
// data type definitions

typedef int BOOL;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

/////////////////////////////////////////////////////////////////////////////
// class CAln

class CAln
{
// Construction
public:
  CAln();
  BOOL Create(int nDim, int nOutput);

// Attributes
public:

  // last error
  int GetLastError() const { return m_nLastError; }

  // pointer to aln
  ALN* GetALN() { return m_pALN; }
  
  // tree
  ALNNODE* GetTree();
  const ALNNODE* GetTree() const;

  // variable epsilon
  float GetEpsilon(int nVar, int nRegion = 0) const;
  void SetEpsilon(float fltEpsilon, int nVar, int nRegion = 0);

  // variable weight max/min
  float GetWeightMin(int nVar, int nRegion = 0) const;
  void SetWeightMin(float fltWMin, int nVar, int nRegion = 0);
  float GetWeightMax(int nVar, int nRegion = 0) const;
  void SetWeightMax(float fltWMax, int nVar, int nRegion = 0);

  // variable max/min
  float GetMin(int nVar, int nRegion = 0) const;
  void SetMin(float fltMin, int nVar, int nRegion = 0);
  float GetMax(int nVar, int nRegion = 0) const;
  void SetMax(float fltMax, int nVar, int nRegion = 0);

	// aln data 
  ALNDATAINFO* GetDataInfo() { return &m_datainfo; }
  void SetDataInfo(float* afltTRdata, const long nTRmaxSamples, long nTRcurrSamples, int nTRcols, long nTRinsert, const float fltMSEorF);
                  
  // region (nRegion must be 0)
  ALNREGION* GetRegion(int nRegion = 0);
  const ALNREGION* GetRegion(int nRegion = 0) const;

  // constraint (nRegion must be 0)
  ALNCONSTRAINT* GetConstraint(int nVar, int nRegion = 0);
  const ALNCONSTRAINT* GetConstraint(int nVar, int nRegion = 0) const;
  
// Operations
public:

  // region, tree growth
#ifdef ENABLE_REGIONS
  int AddRegion(int nParentRegion, float fltLearnFactor, 
                int nConstr, int* anConstr);
#endif
  BOOL AddLFNs(ALNNODE* pParent, int nParentMinMaxType, 
               int nLFNs, ALNNODE** apLFNs = NULL);
  BOOL AddMultiLayer(ALNNODE* pParent, int nParentMinMaxType, 
                     int nLayers, int nMinMaxFanin, int nLFNFanin, 
                     int nFlags);
  BOOL AddTreeString(ALNNODE* pParent, const char* pszTreeString, 
                     int& nParsed);
  BOOL SetGrowable(ALNNODE* pNode);
  
  BOOL Destroy();

  // training
  BOOL Train(int nMaxEpochs, float fltMinRMSErr, float fltLearnRate,
             BOOL bJitter, int nNotifyMask = AN_NONE, 
             ALNDATAINFO* pData = NULL, void* pvData = NULL);

  
  float CalcRMSError(int nNotifyMask = AN_NONE, ALNDATAINFO* pData = NULL,
                      void* pvData = NULL);

  // eval
  BOOL Eval(float* afltResult, int* pnStart = NULL, int* pnEnd = NULL, 
            int nNotifyMask = AN_NONE, ALNDATAINFO* pData = NULL, 
            void* pvData = NULL);

  // quick eval
  float QuickEval(const float* afltX, ALNNODE** ppActiveLFN = NULL);
  
  // get variable monotonicicty, returns -1 on failure
  int VarMono(int nVar);

  // inversion, nVar specifies 0 based index of new output variable WWA
  BOOL Invert(int nVar);

  // save ALN to disk file
  BOOL Write(const char* pszFileName);

  // read ALN from disk file... destroys any existing ALN
  BOOL Read(const char* pszFileName);

  // conversion to dtree
  DTREE* ConvertDtree(int nMaxDepth);

  // confidence intervals
  BOOL CalcConfidence(ALNCONFIDENCE* pConfidence, int nNotifyMask = AN_NONE, 
                      ALNDATAINFO* pData = NULL, void* pvData = NULL);

  static float ConfidencePLimit(const ALNCONFIDENCE* pConfidence, 
                                 float fltSignificance);

  static float ConfidenceTLimit(const ALNCONFIDENCE* pConfidence, 
                                 float fltInterval);

  // lfn analysis
  BOOL LFNAnalysis(void*& pvAnalysis,
                   int& nLFNStats,
                   int nNotifyMask = AN_NONE, 
                   ALNDATAINFO* pData = NULL, void* pvData = NULL);
  static BOOL LFNFreeAnalysis(void* pvAnalysis);
  static BOOL LFNStats(void* pvAnalysis, 
                       int nLFNStat,
                       LFNSTATS& LFNStats,
                       int& nWeightStats,
                       ALNNODE*& pLFN);

  // if nWeight == the output index of the ALN, 
  // then the weight is the bias weight of the LFN
  static BOOL LFNWeightStats(void* pvAnalysis, 
                             int nLFNStat,
                             int nWeightStat,
                             LFNWEIGHTSTATS& LFNWeightStats);

  // ALN internal pseudo-random number generator
  static void SRand(unsigned int nSeed) { ::ALNSRand(nSeed); }
  static unsigned long Rand() { return ::ALNRand(); }
  static float RandFloat() { return ::ALNRandFloat(); }

 

// Overrideables
public:

  // training and evaluation notifications... 
  // return FALSE to cancel training or evaluation
  virtual BOOL OnVectorInfo(VECTORINFO* pVectorInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnTrainStart(TRAININFO* pTrainInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnTrainEnd(TRAININFO* pTrainInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnEpochStart(EPOCHINFO* pEpochInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnEpochEnd(EPOCHINFO* pEpochInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnAdaptStart(ADAPTINFO* pAdaptInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnAdaptEnd(ADAPTINFO* pAdaptInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnLFNAdaptStart(LFNADAPTINFO* pLFNAdaptInfo, void* pvData) 
    { return TRUE; }
  virtual BOOL OnLFNAdaptEnd(LFNADAPTINFO* pLFNAdaptInfo, void* pvData) 
    { return TRUE; }

// Implementation
public:
  ~CAln();
  operator ALN*() { return m_pALN; } 

  // TRdata buffer
  void ALNAPI addTRsample(float* afltX, const int nDim);
  void ALNAPI reduceNoiseVariance();

protected:
  ALN* m_pALN;
  ALNDATAINFO m_datainfo;
  int m_nLastError;

  static int ALNAPI ALNNotifyProc(const ALN* pALN, int nCode, void* pParam, 
                                  void* pvData);


private:
  CAln(const CAln&) {};         // disable copying
  CAln& operator =(const CAln&) { return *this; }
};

/////////////////////////////////////////////////////////////////////////////

#endif // __ALNPP_H__

