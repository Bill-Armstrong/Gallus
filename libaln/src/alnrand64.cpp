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

// alnrand.cpp


#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

///////////////////////////////////////////////////////////////////////////////
// faster rand()

#include "aln.h"
#include "alnpriv.h"
#include <random>
#include <iostream>

unsigned long ALNRand();
float ALNRandFloat();
void ALNSRand(unsigned int nSeed);

std::mt19937_64 generator;
std::uniform_int_distribution<int> distribution1(0, 1474836UL);
std::uniform_real_distribution<float> distribution2(0.0,1.0);

/*
unsigned long ulnumalnrand;
float flnumalnrand;
int main()
{
	for (int i = 0; i < 20; i++)
	{
		ulnumalnrand = distribution1(generator);
		flnumalnrand = distribution2(generator);
		std::cout << ulnumalnrand << " \n ";
		std::cout << flnumalnrand << " \n";
	}
}
*/

unsigned long DoFastRand()
{
    return distribution1(generator);
}
  
unsigned long ALNRand()
{
	return distribution1(generator);
}

float ALNRandFloat() 
{
   return distribution2(generator);
}


void ALNSRand(unsigned int throwAway)
{
  for (int i = 0; i < throwAway; i++)
	{
		generator();
	}
}
