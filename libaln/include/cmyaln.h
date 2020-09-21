//file cmyaln.h

#ifndef _CMYALN_H
#define _CMYALN_H

#include <iostream>

extern "C" int _kbhit();
extern "C" int _getch();


class CMyAln;

class CMyAln : public CAln
{
  // notification overrides
  public:
  virtual BOOL OnTrainStart(TRAININFO* pTrainInfo, void* pvData)
  { 
    //cerr << "Training starts..." << std::endl;
    return TRUE; 
  }
  
	virtual BOOL OnTrainEnd(TRAININFO* pTrainInfo, void* pvData) 
  {	
	std::cerr << "    RMSE " << pTrainInfo->fltRMSErr << " LFNs " << pTrainInfo->nLFNs ;
	return TRUE;
  }

  virtual BOOL OnEpochStart(EPOCHINFO* pEpochInfo, void* pvData) 
  { 
      return TRUE; 
  }

  virtual BOOL OnEpochEnd(EPOCHINFO* pEpochInfo, void* pvData) 
  {
	 //std::cerr << " Leaf nodes " << pEpochInfo->nLFNs << " Active leaf nodes " << pEpochInfo->nActiveLFNs;
	 return TRUE;
  }

  virtual BOOL OnVectorInfo(VECTORINFO* pVectorInfo, void* pvData) 
  {
	// getvector may be used for reinforcement learning
    if(_kbhit())
    {
      //cerr << endl << "Checking keyboard input (type 's' to stop training)" << endl;
      return _getch()!='s';
    }
    else
    {
      return TRUE;
    }
  }
};

#endif
