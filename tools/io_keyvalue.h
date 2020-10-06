//-------------------------------------------------------------------
// $Id$
/// \file io_keyvalue.h
/// Various service functions for writing/reading arrays in files
//
// Copyright (C) 2006-2012 S.Dmytriyeva
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

namespace  io_formats {

/// Print fields to structure outField
class KeyValueWrite
{

    /// Internal structure of file data
    std::iostream& fout;

public:

    /// Constructor
    KeyValueWrite( std::iostream& ff ): fout(ff) {}

    void dump( bool brief )
    {
       if( !brief )
         fout << "\n# End of file";
       fout << std::endl;
    }

    void write_comment( const std::string& comment )
    {
        fout << std::endl << comment;
    }

    template < typename T >
    void write_value( const T& value )
    {
        fout << value;
    }

    /// Writes integral field to a json.
    template <class T>
    void write_key_value( const std::string& field_name, const T& value  )
    {
        fout << std::endl << "<" << field_name << ">  ";
        write_value(value);
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void write_array( const std::string& field_name, const std::vector<double>& arr, long int  );

    /// Writes T array to a text file.
    template < typename T, typename LT >
    void write_array( const std::string& field_name, T* arr, LT size, LT l_size )
    {
        int sz = ( l_size > 0 ? l_size: 40 );
        fout << std::endl << "<" << field_name << ">" << std::endl;

        for( int ii=0, jj=0; ii<size; ii++, jj++  )
        {
            if( jj == sz )
            {
                jj=0;
                fout << std::endl;
            }
            write_value(arr[ii]);
            fout << " ";
        }
    }

    /// Writes char array to a json file.
    template < typename T=char, typename LT >
    void write_array( const std::string& field_name, char*  arr, LT size, LT arr_size )
    {
        bool isComment = false;

        if( field_name.empty() )
            fout << std::endl << "<" << field_name << ">" << std::endl;
        else
        {
            fout << std::endl << "#  ";
            isComment = true;
        }
        for( long int ii=0, jj=0; ii<size; ii++, jj++  )
        {
            if( jj == 40 )
            {
                jj=0;
                fout << std::endl;
                if(isComment)
                    fout << "#  ";
            }
            std::string str = std::string( arr +(ii*arr_size), 0, arr_size );
            write_value( str );
            fout << " ";
        }
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void write_sel_array( const std::string& field_name, T* arr, LT size, long int* sel_arr, LT ncolumns, LT l_size )
    {
        int sz = ( l_size > 0 ? l_size: 40 );

        fout << std::endl << "<" << field_name << ">" << std::endl;
        for( long int ii=0, jj=0; ii<size; ii++  )
        {
            for(long int cc=0; cc<ncolumns; cc++ )
            {
                if(jj == sz)
                {
                    jj=0;
                    fout << std::endl;
                }
                write_value( arr[sel_arr[ii]*ncolumns+cc] );
                fout << " ";
                jj++;
            }
        }
    }

};


/// Read fields of structure
class KeyValueRead
{

    /// Internal structure of file data
    std::iostream& fin;

    /// Reads value from a text file.
    void read_value( float& val );
    /// Reads value from a text file.
    void read_value( double& val );

    void  skip_space();

public:

    /// Constructor
    KeyValueRead( std::iostream& ff ): fin( ff )  {}

    /// Reset internal data
    void reset()  {}

    /// Read next name from file
    bool has_next( std::string& next_field_name );

    /// Read next label from file ( must be exist, otherwise error )
    bool test_next_is( const std::string& label );

    /// Reads array from a text file.
    template <class T,
              typename std::enable_if<std::is_integral<T>::value,T>::type* = nullptr>
    void read_array( const char *, T* arr, long int size )
    {
        for( long int ii=0; ii<size; ii++  )
        {
            skip_space();
            fin >> arr[ii];
        }
    }

    /// Reads array from a text file.
    template <class T,
              typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
    void read_array( const char *, T* arr, long int size )
    {
        for( long int ii=0; ii<size; ii++  )
            read_value(arr[ii]);
    }


    /// Reads strings array from a text file.
    void read_strings_array( const std::string& field_name, char* arr, long int size, long int el_size );

    /// Reads double vector from a text file.
    void read_array( const std::string& name, std::vector<double> arr );

};


template <> void KeyValueWrite::write_value( const double& );
template <> void KeyValueWrite::write_value( const float& );
template <> void KeyValueWrite::write_value( const char& value );
template <> void KeyValueWrite::write_value( const std::string& value );

}  // io_formats

