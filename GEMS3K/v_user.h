//-------------------------------------------------------------------
// $Id: v_user.h 1373 2009-07-22 12:25:22Z gems $
//
// Declaration of platform-specific utility functions and classes
//
// Copyright (C) 1996,2001 A.Rysin, S.Dmytriyeva
// Uses  gstring class (C) A.Rysin 1999
//
// This file is part of the GEM-Vizor library and the GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _v_user_h_
#define _v_user_h_

#include <algorithm>
#include <iostream>
using namespace std;

#include "string.h"
#include "verror.h"

#ifdef __APPLE__

#ifndef __unix
#define __unix
#endif
#ifndef __FreeBSD
#define __FreeBSD
#endif

typedef unsigned int uint;
#endif

const int MAXKEYWD = 6+1;

#ifndef  __unix

typedef unsigned int uint;

#endif //  __noborl

void Gcvt(double number, size_t ndigit, char *buf);
double NormDoubleRound(double aVal, int digits);
void NormDoubleRound(double *aArr, int size, int digits);
void NormFloatRound(float *aArr, int size, int digits);

inline
int ROUND(double x )
{
    return int((x)+.5);
}

template <class T>
inline
void fillValue( T* arr, T value, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = value;
}

template <class T>
inline
void copyValues( T* arr, T* data, int size )
{
  if( !arr || !data )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = data[ii];
}

inline
void copyValues( double* arr, float* data, int size )
{
  if( !arr || !data )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = (double)data[ii];
}

inline
void copyValues( float* arr, double* data, int size )
{
  if( !arr || !data )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = (float)data[ii];
}

inline
void copyValues( long int* arr, short* data, int size )
{
  if( !arr || !data )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = (long int)data[ii];
}


//#define max( a, b )  ( (a) >( b) ? (a) : (b) )
//#define min( a, b )  ( (a) <( b) ? (a) : (b) )



// Combines path, directory, name and extension to full pathname
gstring
u_makepath(const gstring& dir,
           const gstring& name, const gstring& ext);

// Splits full pathname to path, directory, name and extension
void
u_splitpath(const gstring& Path, gstring& dir,
            gstring& name, gstring& ext);

#define fileNameLength 64
// Get Path of file and Reading list of file names from it, return number of files 
char  (* f_getfiles(const char *f_name, char *Path, 
		long int& nElem, char delim ))[fileNameLength];



#endif
// _v_user_h_

