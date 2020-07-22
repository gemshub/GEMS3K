//-------------------------------------------------------------------
// $Id$
/// \file v_detail.h
/// Declaration of platform-specific utility functions and classes
//
// Copyright (C) 1996,2001,2020 A.Rysin, S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#ifndef V_DETAIL_H
#define V_DETAIL_H

#include <limits>
#include <cmath>
//#include <algorithm>
#include "verror.h"


void strip( std::string& str);
void replace( std::string& str, const char* old_part, const char* new_part);
void replaceall( std::string& str, const char* old_part, const char* new_part);
void replaceall(std::string& str, char ch1, char ch2);
/// read string as: "<characters>"
std::istream& f_getline(std::istream& is, std::string& str, char delim);
/// Combines path, directory, name and extension to full pathname
std::string u_makepath(const std::string& dir,
           const std::string& name, const std::string& ext);

/// Splits full pathname to path, directory, name and extension
void u_splitpath(const std::string& Path, std::string& dir,
            std::string& name, std::string& ext);

#define FLOAT_EMPTY	          1.17549435e-38F
#define DOUBLE_EMPTY         2.2250738585072014e-308
#define CHAR_EMPTY   	     '`'

inline bool IsFloatEmpty( const float v )
{
    return ( v>0. && v <= FLOAT_EMPTY);
}
inline bool IsDoubleEmpty( const double v )
{
    return ( v>0. && v <= DOUBLE_EMPTY);
}

template <class T>
inline void fillValue( T* arr, T value, int size )
{
  if( !arr )
    return;
  for(int ii=0; ii<size; ii++)
    arr[ii] = value;
}

template <class T, class IT>
inline void copyValues( T* arr, T* data, IT size )
{
  if( !arr || !data )
    return;
  for(IT ii=0; ii<size; ii++)
    arr[ii] = data[ii];
}

template <class IT>
inline void copyValues( double* arr, float* data, IT size )
{
  if( !arr || !data )
    return;
  for(IT ii=0; ii<size; ii++)
    arr[ii] = static_cast<double>(data[ii]);
}

template <class IT>
inline void copyValues( float* arr, double* data, IT size )
{
  if( !arr || !data )
    return;
  for(IT ii=0; ii<size; ii++)
    arr[ii] = static_cast<float>(data[ii]);
}

template <class IT>
inline void copyValues( long int* arr, short* data, IT size )
{
  if( !arr || !data )
    return;
  for(IT ii=0; ii<size; ii++)
    arr[ii] = static_cast<long int>(data[ii]);
}

template<typename T>
bool approximatelyEqual( const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

template<typename T>
bool essentiallyEqual( const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

template<typename T>
bool definitelyGreaterThan( const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

template<typename T>
bool definitelyLessThan( const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

template<typename T>
bool approximatelyZero( const T& a, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return fabs(a) <=  epsilon;
}

template<typename T>
bool noZero( const T& a, const T& epsilon = std::numeric_limits<T>::epsilon() )
{
    return fabs(a) >  epsilon;
}

#endif // V_DETAIL_H
