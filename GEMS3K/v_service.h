//-------------------------------------------------------------------
// $Id$
/// \file v_service.h
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

#ifndef V_SERVICE_H
#define V_SERVICE_H

#include <string>
#include <vector>
#include <algorithm>

std::string char_array_to_string(const char* data_ptr, size_t max_size);
void strip(std::string& str);
void replace(std::string& str, const char* old_part, const char* new_part);
void replaceall(std::string& str, const std::string& old_part, const std::string& new_part);
void replaceall(std::string& str, char ch1, char ch2);
std::vector<std::string> split(const std::string& str, const std::string& delimiters);

//// Extract the string value from data.
std::string regexp_extract_string( std::string regstr, std::string data );
/// Extract the string value by key from json string
std::string extract_string_json( std::string key, std::string jsondata );
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


template < typename T, typename LT >
std::string to_string( T* arr, LT size )
{
    std::string logs;
    for( LT ii=0; ii<size; ++ii) {
        logs += std::to_string(arr[ii])+" ";
    }
    return logs;
}

#endif // V_SERVICE_H
