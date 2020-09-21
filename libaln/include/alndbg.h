/*
// ALNfit Learning Engine for approximation of functions defined by samples.
// Copyright (C) 2018 William W. Armstrong.

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// Version 3 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

// For further information contact 
// William W. Armstrong
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4
*/

/*
// alndbg.h
*/

#ifndef __ALNDBG_H__
#define __ALNDBG_H__

#ifndef __ALNCFG_H__
  #include <alncfg.h>
#endif

/*
/////////////////////////////////////////////////////////////////////////////
// Special ALNDebugBreak: used to break into debugger at critical times
// Non-Microsoft, non-Intel compilers will need to implement DebugBreak()
// This was removed by WWA on July 22, 2019.
*/

#ifndef ALNDebugBreak
#define ALNDebugBreak()
#endif


/*
///////////////////////////////////////////////////////////////////////////////
// diagnostic support
*/

#ifdef _DEBUG

ALNIMP void __cdecl ALNTrace(const char* pszFormat, ...);

#ifndef ALNTRACE
#define ALNTRACE ALNTrace
#endif  /* ALNTRACE */

#ifndef assert
#include <assert.h>
#endif

#ifndef ALNASSERT
#define ALNASSERT(f) assert(f)
#endif  /* ALNASSERT */

#ifndef ALNVERIFY
#define ALNVERIFY(f) ALNASSERT(f)
#endif  /* ALNVERIFY */

#else   /* _DEBUG */

inline void __cdecl ALNTrace(const char* pszFmt, ...) { }

#ifndef ALNTRACE
#define ALNTRACE              1 ? (void)0 : ALNTrace
#endif  /* ALNTRACE */

#ifndef ALNASSERT
#define ALNASSERT(f) ((void)0)
#endif  /* ALNASSERT */

#ifndef ALNVERIFY
#define ALNVERIFY(f) ((void)(f))
#endif  /* ALNVERIFY */

#endif  /* !_DEBUG */


/*
///////////////////////////////////////////////////////////////////////////////
*/

#endif  /* __ALNDBG_H__ */
