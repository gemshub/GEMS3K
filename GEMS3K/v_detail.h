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
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "verror.h"


void strip( std::string& str);
void replace( std::string& str, const char* old_part, const char* new_part );
void replaceall( std::string& str, const std::string& old_part, const std::string& new_part );
void replaceall(std::string& str, char ch1, char ch2);

//// Extract the string value from data.
std::string regexp_extract_string( const std::string& regstr, const std::string& data );
/// Extract the string value by key from json string
std::string extract_string_json( const std::string& key, const std::string& jsondata );
/// Extract the int value by key from json string
int extract_int_json( const std::string& key, const std::string& jsondata );


/// read string as: "<characters>"
std::istream& f_getline(std::istream& is, std::string& str, char delim);
/// Combines path, directory, name and extension to full pathname
std::string u_makepath(const std::string& dir,
           const std::string& name, const std::string& ext);

/// Replace all characters to character in string (in place).
inline void replace_all(std::string &s, const std::string &characters, char to_character )
{
    std::replace_if( s.begin(), s.end(), [=](char ch) {
        return characters.find_first_of(ch)!=std::string::npos;
    }, to_character );
}

/// Splitting full pathname to path, directory, name and extension
void u_splitpath(const std::string& pathname, std::string& dir,
            std::string& name, std::string& ext);
/// Get directory from full pathname
std::string u_getpath( const std::string& pathname );

inline int ROUND(double x )
{
    return int((x)+.5);
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

template <class T, class IT>
bool isAllZero( T* arr, IT size )
{
  if( !arr  )
    return true;
  for(IT ii=0; ii<size; ii++)
    if(  !approximatelyZero (arr[ii]) )
       return false;

  return true;
}

/// Trim all whitespace characters from start (in place).
inline void ltrim(std::string &s )
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !isspace(ch);
    }));
}

/// Trim all whitespace characters from end (in place).
inline void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !isspace(ch);
    }).base(), s.end());
}

/// Trim all whitespace characters from both ends (in place).
inline void trim(std::string &s )
{
    ltrim(s);
    rtrim(s);
}

/// Trim characters from start (in place).
inline void ltrim(std::string &s, const std::string &characters )
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [=](char ch) {
        return characters.find_first_of(ch)==std::string::npos;
    }));
}

/// Trim characters from end (in place).
inline void rtrim(std::string &s, const std::string &characters)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [=](char ch) {
        return characters.find_first_of(ch)==std::string::npos;
    }).base(), s.end());
}

/// Trim characters from both ends (in place).
inline void trim(std::string &s, const std::string &characters )
{
    ltrim(s, characters);
    rtrim(s, characters);
}

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

const float FLOAT_INFPLUS = 3.4028230E+38F;
const float FLOAT_INFMINUS = -3.4028230E+38F;
const float FLOAT_NAN = 3.402800E+38F;
const double DOUBLE_INFPLUS = 1.79769313486230000E+308;
const double DOUBLE_INFMINUS = -1.79769313486230000E+308;
const double DOUBLE_NAN = 1.797693000000000E+308;
const std::string DOUBLE_INFPLUS_STR = std::to_string(DOUBLE_INFPLUS);
const std::string DOUBLE_INFMINUS_STR = std::to_string(DOUBLE_INFMINUS);
const std::string DOUBLE_NAN_STR = std::to_string(DOUBLE_NAN);

template < typename T >
T InfMinus()
{
  return std::numeric_limits<T>::min();
}
template < typename T >
T InfPlus()
{
  return std::numeric_limits<T>::max();
}
template < typename T >
T Nan()
{
  return std::numeric_limits<T>::min();
}

template < typename T >
inline bool is_minusinf(const T& value)
{
  return approximatelyEqual( value, InfMinus<T>() );
}

template < typename T >
inline bool is_plusinf(const T& value)
{
  return std::isinf(value) || approximatelyEqual( value, InfPlus<T>() );
}

template < typename T >
inline bool is_nan(const T& value)
{
  return std::isnan(value) || approximatelyEqual( value, Nan<T>() );
}

template <class T,
          typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
std::string write_floating_point( const T& value )
{
    std::string value_str;
    if( is_minusinf(value) )
        value_str = "-inf";
    else if( is_plusinf(value) )
        value_str = "inf";
    else if( is_nan(value) )
        value_str = "nan";
    else
    {
        std::ostringstream str_stream;
        str_stream << std::setprecision(15) << value;
        value_str = str_stream.str();
    }
    return value_str;
}

template < typename T >
T read_floating_point( const std::string& value_str )
{
    T value;
    if( value_str == "-inf" )
        value = InfMinus<T>();
    else if( value_str == "inf" )
        value = InfPlus<T>() ;
    else if( value_str == "nan" )
        value = Nan<T>();
    else
        value = std::stod( value_str.c_str() );

    return value;
}

template <> double InfMinus();
template <> double InfPlus();
template <> double Nan();
template <> float InfMinus();
template <> float InfPlus();
template <> float Nan();


#endif // V_DETAIL_H
