//-------------------------------------------------------------------
// $Id$
/// \file v_detail.h
/// Declaration of platform-specific utility functions and classes
//
// Copyright (C) 1996,2001,2021 A.Rysin, S.Dmytriyeva
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
#include <sstream>
#include <iomanip>
#include "verror.h"

template <typename T>
bool is( T& x, const std::string& s)
{
  std::istringstream iss(s);
  return iss >> x && !iss.ignore();
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

//#define FLOAT_EMPTY	          1.17549435e-38F
//#define DOUBLE_EMPTY         2.2250738585072014e-308
const float FLOAT_EMPTY = std::numeric_limits<float>::min();
const double DOUBLE_EMPTY = std::numeric_limits<double>::min();
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
  return std::numeric_limits<T>::max();
}

template < typename T >
inline bool is_minusinf(const T& value)
{
  return value <= InfMinus<T>();
}

template < typename T >
inline bool is_plusinf(const T& value)
{
  return std::isinf(value) || value >= InfPlus<T>();
}

template < typename T >
inline bool is_nan(const T& value)
{
  return std::isnan(value) || value >= Nan<T>();
}

template <class T,
          typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
std::string floating_point_to_string( const T& value, int precision = 15 )
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
        str_stream << std::setprecision(precision) << value;
        value_str = str_stream.str();
    }
    return value_str;
}

template <class T,
          typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
T string_to_floating_point( const std::string& value_str )
{
    T value;
    if( value_str == "-inf" )
        value = InfMinus<T>();
    else if( value_str == "inf" )
        value = InfPlus<T>();
    else if( value_str == "nan" )
        value = Nan<T>();
    else
        value = static_cast<T>(std::stod(value_str.c_str()));

    return value;
}

/// Read value from string.
template <class T,
          typename std::enable_if<!std::is_floating_point<T>::value,T>::type* = nullptr>
bool string2value( T& x, const std::string& s)
{
    if(s.empty())
        return false;
    std::istringstream iss(s);
    return iss >> x && !iss.ignore();
}

template <class T,
          typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
bool string2value( T& x, const std::string& s)
{
    if(s.empty())
        return false;
    x = string_to_floating_point<T>(s);
    return true;
}

/// Serializations a numeric value to a string.
template <class T,
          typename std::enable_if<!std::is_floating_point<T>::value,T>::type* = nullptr>
std::string value2string( const T& value, int )
{
    std::ostringstream os;
    os << value;
    return os.str();
}

template <class T,
          typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
std::string value2string( const T& value, int precision = 15 )
{
    return floating_point_to_string( value, precision );
}

template <> double InfMinus();
template <> double InfPlus();
template <> double Nan();
template <> float InfMinus();
template <> float InfPlus();
template <> float Nan();


#endif // V_DETAIL_H
