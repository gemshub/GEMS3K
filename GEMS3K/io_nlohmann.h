//-------------------------------------------------------------------
// $Id$
/// \file io_nlohmann.h
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

#ifdef USE_NLOHMANNJSON

#include <nlohmann/json.hpp>
#include <fstream>
#include "verror.h"
#include "v_service.h"

namespace  io_formats {

/// Print fields to structure outField
class NlohmannJsonWrite
{

    /// Internal structure of file data
    nlohmann::json json_data;
    nlohmann::json top_data;
    std::iostream& fout;
    std::string current_set_name;
    std::string top_field_name;

    template <class T>
    void add_value( const T& value, nlohmann::json& json_arr  )
    {
        json_arr.push_back(value);
    }

    std::string key( const std::string& name ) const
    {
        return  name;
    }

public:

    /// Constructor
    NlohmannJsonWrite( std::iostream& ff, const std::string& test_set_name ):
        json_data(), top_data(), fout(ff), current_set_name(test_set_name)
    {}
    const std::string& set_name() const
    {
      return current_set_name;
    }

    void put_head( const std::string &key_name, const std::string &the_field_name )
    {
      top_field_name = the_field_name;
      top_data["set"] = current_set_name;
      top_data["_key"] = key_name;
    }

    void dump( bool not_brief )
    {
        top_data[top_field_name] = json_data;
        auto json_array =  nlohmann::json::array();
        json_array.push_back(top_data);
        fout << json_array.dump(( not_brief ? 4 : 0 ));
    }

    void write_comment( const std::string&  ) {}

    template < typename T >
    void writeValue( const T& value )
    {
        if( json_data.is_array() )
            json_data.push_back(value);
    }

    /// Writes integral field to a json.
    template <class T>
    void write_key_value( const std::string& field_name, const T& value  )
    {
        json_data[ key( field_name ) ] = value;
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void write_array( const std::string& field_name, const std::vector<double>& arr, long int  )
    {
        json_data[ key(field_name) ] = arr;
    }

    /// Writes T array to a text file.
    template < typename T, typename LT >
    void write_array( const std::string& field_name, T* arr, LT size, LT  )
    {
        auto arr_key = key( field_name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( int ii=0; ii<size; ii++ )
        {
            add_value( arr[ii], json_data[ arr_key ]);
        }
    }

    /// Writes char array to a json file.
    template < typename T=char, typename LT >
    void write_array( const std::string& field_name, char*  arr, LT size, LT arr_size )
    {
        if( field_name.empty() )  // comment
            return;
        auto arr_key = key( field_name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( int ii=0, jj=0; ii<size; ii++, jj++  )
        {
            std::string str = char_array_to_string( arr +(ii*arr_size), arr_size );
            add_value( str, json_data[ arr_key ] );
        }
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void write_sel_array( const std::string& field_name, T* arr, LT size, long int* sel_arr, LT ncolumns, LT  )
    {
        auto arr_key = key( field_name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( long int ii=0; ii<size; ii++  )
        {
            for(long int cc=0; cc<ncolumns; cc++ )
            {
                add_value( arr[sel_arr[ii]*ncolumns+cc], json_data[ arr_key ] );
            }
        }
    }

};


/// Read fields of structure
class NlohmannJsonRead
{

    /// Internal structure of file data
    nlohmann::json json_data;
    nlohmann::json::iterator json_it;

    std::string key( const std::string& name ) const
    {
        return name;
    }

public:

    /// Constructor
    NlohmannJsonRead( std::iostream& ff, const std::string& test_set_name,  const std::string& field_name );

    /// Reset json loop
    void reset()
    {
      json_it = json_data.begin();
    }

    /// Read next name from file
    bool  has_next( std::string& next_field_name );

    /// Read next label from file ( must be exist, otherwise error )
    bool test_next_is( const std::string& label )
    {
        return  json_data.find( key(label) ) != json_data.end();
    }

    /// Reads array from a TIO format.
    template <class T>
    void read_array(  const std::string& field_name, T* arr, long int size )
    {
        std::string msg;
        std::string jkey = key( field_name );
        auto json_arr_it = json_data.find(jkey);
        if( json_arr_it == json_data.end() )
        {
            msg = field_name;
            msg += " - No data where expected.\n";
            Error( "Json read error 01", msg );
        }

        if( !json_arr_it->is_structured() &&  size==1 )
        {
            json_arr_it->get_to<T>(*arr);
        }
        else
        {
            if( json_arr_it->size() != static_cast<size_t>(size) )
            {
                msg = field_name;
                msg += " - No size (";
                msg += std::to_string(size).c_str();
                msg += ") as expected ";
                msg += std::to_string(json_arr_it->size()).c_str();
                Error( "Json read error 02", msg );
            }
            for( long int ii=0; ii<size; ii++  )
            {
                json_arr_it->at(ii).get_to<T>( arr[ii] );
            }
        }
    }

    /// Reads strings array from a text file.
    void read_strings_array( const std::string& field_name, char* arr, long int size, long int el_size );

    /// Reads double vector from a text file.
    void read_array( const std::string& name, std::vector<double> arr );

    /// Empty function
    bool skip_line() { return false; }
};


template <> void NlohmannJsonWrite::add_value( const char& value, nlohmann::json& json_arr );
template <> void NlohmannJsonWrite::add_value( const std::string& value, nlohmann::json& json_arr );
template <> void NlohmannJsonWrite::write_key_value( const std::string& field_name, const char& value );
template <> void NlohmannJsonWrite::write_key_value( const std::string& field_name, const std::string& value );

}  // io_formats

#endif
