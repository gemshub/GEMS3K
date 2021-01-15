//-------------------------------------------------------------------
// $Id$
/// \file io_simdjson.h
/// Various service functions for writing/reading arrays in files
//
// Copyright (C) 2020 S.Dmytriyeva
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

#pragma once

#include <fstream>
#include <vector>
#include <memory>
#include "verror.h"

namespace  io_formats {

/// Print fields to structure outField
class SimdJsonWrite
{

public:

    /// Constructor
    SimdJsonWrite( std::iostream& ff, const std::string& test_set_name, bool not_brief ):
        fout(ff), dense(not_brief), current_set_name(test_set_name)
    {}

    const std::string& set_name() const
    {
      return current_set_name;
    }
    void put_head( const std::string &key_name, const std::string &field_name );
    void dump( bool  );

    void write_comment( const std::string&  ) {}


    /// Writes integral field to a json.
    template <class T>
    void write_key_value( const std::string& field_name, const T& value  )
    {
        add_key( field_name );
        add_value( value );
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - interpret the flag with_comments=true as "pretty JSON" and
    ///                                   with_comments=false as "condensed JSON"
    /// \param brief_mode - Do not write data items that contain only default values
    void write_array( const std::string& field_name, const std::vector<double>& arr, long int l_size );

    /// Writes T array to a text file.
    template < typename T, typename LT >
    void write_array( const std::string& field_name, T* arr, LT size, LT l_size )
    {
        LT jj=0, sz = ( l_size > 0 ? l_size: values_in_line );
        add_key( key( field_name ) );
        fout  << ( dense ? "[\n        " : "[" );

        for( LT ii=0; ii<size; ii++, jj++ )
        {
            add_next( ii, jj, sz );
            add_value( arr[ii] );
        }
        fout << ( dense ? "\n    ]" : "]" );
    }

    /// Writes char array to a json file.
    template < typename T=char, typename LT >
    void write_array( const std::string& field_name, char*  arr, LT size, LT arr_size )
    {
        if( field_name.empty() )  // comment
            return;

        LT jj=0, sz = values_in_line;
        add_key( key( field_name ) );
        fout  << ( dense ? "[\n        " : "[" );

        for( LT ii=0; ii<size; ii++, jj++ )
        {
            add_next( ii, jj, sz );
            add_value( std::string( arr +(ii*arr_size), 0, arr_size ) );
        }

        fout << ( dense ? "\n    ]" : "]" );
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void write_sel_array( const std::string& field_name, T* arr, LT size, long int* sel_arr, LT ncolumns, LT l_size )
    {
        LT jj=0, kk=0, sz = ( l_size > 0 ? l_size: values_in_line );
        add_key( key( field_name ) );
        fout  << ( dense ? "[\n        " : "[" );

        for( LT ii=0; ii<size; ii++ )
        {
            for(LT cc=0; cc<ncolumns; cc++ )
            {
              add_next( kk++, jj++, sz );
              add_value( arr[sel_arr[ii]*ncolumns+cc] );
            }
        }

        fout << ( dense ? "\n    ]" : "]" );
    }


private:

    /// Internal structure of file data
    std::iostream& fout;
    bool dense = false;
    bool first = false;
    const long values_in_line = 10;
    std::string current_set_name;

    template <class T>
    void add_value( const T& value  )
    {
        fout << value;
    }


    void add_key( const std::string& field_name )
    {
        if( first )
        {
            first = false;
            fout << ( dense ? "\n" : "" );
        }
        else
        {
            fout << ( dense ? ",\n" : "," );
        }
        fout << ( dense ? "    \"" : "\"" );
        fout << field_name << ( dense ? "\": " : "\":" );
    }

    void add_next( long ndx, long& jj, long line_size )
    {
        if( ndx )
            fout << ( dense ? ", " : "," );
        if( jj == line_size )
        {
            jj=0;
            fout << ( dense ? "\n        " : "" );
        }
    }

    std::string key( const std::string& name ) const
    {
        return  name;
    }


};

class SimdJsonImpl;

/// Read fields of structure
class SimdJsonRead
{

public:

    /// Constructor
    SimdJsonRead( std::iostream& ff, const std::string& test_set_name,  const std::string& field_name );

    /// Reset json loop
    void reset();

    /// Read next name from file
    bool  has_next( std::string& next_field_name );

    /// Read next label from file ( must be exist, otherwise error )
    bool test_next_is( const std::string& label );

    /// Reads array from a TIO format.
    template <class T,
              typename std::enable_if<std::is_integral<T>::value,T>::type* = nullptr>
    void read_array(  const std::string& field_name, T* arr, long int size )
    {
        std::string msg;
        std::vector<int64_t> js_arr;
        read_array( field_name, js_arr );

        ErrorIf( static_cast<size_t>(size) > js_arr.size(), std::string("SimdJson read error :"),
                 "illegal array size "+ field_name );

        for( long int ii=0; ii<size; ++ii )
        {
            arr[ii] = static_cast<T>(js_arr[ii]);
        }
    }

    /// Reads array from a TIO format.
    template <class T,
              typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
    void read_array(  const std::string& field_name, T* arr, long int size )
    {
        std::string msg;
        std::vector<double> js_arr;
        read_array( field_name, js_arr );

        ErrorIf( static_cast<size_t>(size) > js_arr.size(), std::string("SimdJson read error :"),
                 "illegal array size "+ field_name );

        for( long int ii=0; ii<size; ++ii )
        {
            arr[ii] = static_cast<T>(js_arr[ii]);
        }
    }

    /// Reads strings array from a text file.
    void read_strings_array( const std::string& field_name, char* arr, long int size, long int el_size );
    /// Reads double vector from a text file.
    void read_array( const std::string& name, std::vector<double>& arr );
    /// Reads int vector from a text file.
    void read_array(const std::string &field_name, std::vector<int64_t>& arr);

    /// Empty function
    bool skip_line() { return false; }

protected:

    // Internal structure of file data
    std::shared_ptr<SimdJsonImpl> impl;

};


template <> void SimdJsonWrite::add_value( const double& );
template <> void SimdJsonWrite::add_value( const float& );
template <> void SimdJsonWrite::add_value( const char& value );
template <> void SimdJsonWrite::add_value( const std::string& value );

}  // io_formats

