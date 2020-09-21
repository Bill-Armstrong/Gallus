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

// calcactivechild.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


int ALNAPI CalcActiveChild(float& fltRespActive, float& fltDistance, 
                           float flt0, float flt1, const ALNNODE* pNode)                          

{
  int nActive = -1;

  // MAX node handling
	if (MINMAX_ISMAX(pNode)) 
	{
		if(flt1 > flt0)	//  this puts child 1 100% active
		{
		  nActive = 1;
		  fltRespActive = 1.0;
		  fltDistance = flt1;
		}
		else	// child 0 is 100% active
		{
			nActive = 0;
		  fltRespActive = 1.0;
		  fltDistance = flt0;
		}	
	}
	else // this is a MIN node
	{    
  		if(flt1 < flt0) //  this puts child 1 100% active
		{
			nActive = 1;
			fltRespActive = 1.0;
			fltDistance = flt1;
  		}
		else	 // child 0 is 100% active
		{
			nActive = 0;
			fltRespActive = 1.0;
			fltDistance = flt0;
		}
	}

  ASSERT(nActive != -1);
  return nActive;
}
