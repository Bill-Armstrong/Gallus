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

// plimit.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include "boost\math\special_functions\beta.hpp"

using namespace boost::math;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// calculate probability p of an event occuring, such that
// the probablity of m or less such events occuring in n trials
// is x; currently limited to accuracy of 1.e-7

// implementation based on Ridders method documented in Press et al

// method: advance p until cumulative binomial dist 0 to m events in n 
//         trials drops to x

float ALNAPI PLimit(int n, int m, float fltX)
{
  static const float fltInc = 0.1;     // coarse increment
  static const float fltAcc = 1.0e-7;  // maximum accuracy

  if (fltX < 0.0 || fltX > 1.0 || n < 0)
  {
    return NAN;
  }
  
  // known cases
  if (m < 0)
    return 0.0;
  if (m >= n)
    return 1.0;

 
  // P is desired probability, Y is difference between desired area
  // under tail and the area under the tail given P

  // therfore, Y goes to 0 as P approaches desired value

  // lower bound
  float fltP1 = 0.0;
  float fltY1 = fltX - 1.0;
  ASSERT(fltY1 <= 0.0);

  // begin coarse approximation of P
  
  // scan upward until Y3 is +ve
  float fltP3, fltY3;
  for (fltP3 = fltInc; fltP3 < 1.0; fltP3 += fltInc)
  {
    fltY3 = fltX - (1.0 - ibeta(m + 1, n - m, fltP3));  //ibet??
                   // cumulative binomial dist 0 to m events in n trials, 
                   // see Press et al p229

    // convergence test (unlikely at this point)
    if (fabs(fltY3) < fltAcc)
      return fltP3;

    // check for sign change
    if (fltY3 > 0.0)
      break;  // we've bracketed desired value

    // else, new lower bound
    fltP1 = fltP3;
    fltY1 = fltY3;
  }

  // P1 and P3 bracket desired value... refine using ridders method
  const int nMaxIt = 100;
  
  for (int i = 0; i < nMaxIt; i++)
  {
    // get mid-values
    float fltP2 = 0.5 * (fltP1 + fltP3);
    if ((fltP3 - fltP1) < fltAcc)   // convergence test
      return fltP2;
    
    float fltY2 = fltX - (1.0 - ibeta(m + 1, n - m, fltP2)); //ibeta??

    // convergence test
    if (fabs(fltY2) < fltAcc)
      return fltP2;

    float fltDenom = sqrt(fltY2 * fltY2 - fltY1 * fltY3);  // y1, y3 opposite sign
    float fltTrial = fltP2 + (fltP1 - fltP2) * fltY2 / fltDenom;

    float fltY = fltX - (1.0 - ibeta(m + 1, n - m, fltTrial));  //ibeta

    // convergence test
    if (fabs(fltY) < fltAcc)
      return fltTrial;

    // between mid and test point?
    if ((fltY2 < 0.0) && (fltY > 0.0))
    {
      fltP1 = fltP2;    // new lower bound
      fltY1 = fltY2;
      fltP3 = fltTrial; // new upper bound
      fltY3 = fltY;
    }
    else if ((fltY < 0.0) && (fltY2 > 0.0))
    {
      fltP1 = fltTrial; // new lower bound
      fltY1 = fltY;
      fltP3 = fltP2;    // new upper bound
      fltY3 = fltY2;
    }
    else if (fltY < 0.0)  // both negative
    {
      fltP1 = fltTrial;
      fltY1 = fltY;
    }
    else  // both positive
    {
      fltP3 = fltTrial;
      fltY3 = fltY;
    }
  }

  // convergence failed... return best guess?
  return 0.5 * (fltP1 + fltP3); 
}