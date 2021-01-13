//-------------------------------------------------------------------
// $Id$
/// \file io_template.h
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

#include <vector>
#include "verror.h"

namespace  io_formats {

/// Internal descriptions of fields
struct outField
{
    /// Name of field in structure
    std::string name;
    /// 1 - must be read, 0 - default values can be used
    long int alws;
    /// 0; set to 1 after reading the field from input file
    long int readed;
    /// 1 - static object; 0 - undefined; <0 type of indexation, >1 number of elements in array
    long int indexation;
    /// Field description
    std::string comment;
};

/// Basic class for read/write fields of structure
class TRWArrays
{

public:

    /// Constructor
    TRWArrays( short aNumFlds, outField* aFlds ):
        num_flds(aNumFlds), flds(aFlds)
    {}

    virtual ~TRWArrays()
    {}

    /// Find field by name
    virtual  long int findFld( const std::string& name ) const;

    /// Set the data object can be skipped from the file
    /// and default value(s) can be used
    /// \param ii index in array flds
    void  setNoAlws( long int ii )
    {  flds[ii].alws = 0; }

    /// Set the data object can be skipped from the file
    /// and default value(s) can be used
    /// \param Name of field in array flds
    void  setNoAlws( const std::string& name )
    {
        long int ii = findFld( name );
        if( ii >=0 )
            setNoAlws(ii);
    }

    /// Set the data object  must be always present in the file
    /// \param ii index in array flds
    void  setAlws( long int ii )
    {  flds[ii].alws = 1; }

    /// Set the data object name must be always present in the file
    /// \param Name of field in array flds
    void  setAlws( const std::string& name )
    {
        long int ii = findFld( name );
        if( ii >=0 )
            setAlws(ii);
    }

    /// Test the data object  must be always present in the file
    /// \param ii index in array flds
    bool  getAlws( long int ii ) const
    {  return (flds[ii].alws == 1); }

    /// Test the data object name must be always present in the file
    /// \param Name of field in array flds
    bool  getAlws(  const std::string& name  ) const
    {
        long int ii = findFld( name );
        if( ii >=0 )
            return getAlws(ii);
        else
            return false;
    }

    /// Clear readed information
    void reset()
    {
        for( long int ii=0; ii < num_flds; ii++ )
            flds[ii].readed = 0;
    }

    /// Test for reading all fields must be always present in the file
    std::string testRead() const;

protected:

    /// Size of array flds
    long int num_flds;
    /// Array of permissible fields
    outField* flds;

};

/// Print fields of structure outField
template < class TIO >
class TPrintArrays: public  TRWArrays
{
    TIO& out_format;

public:

    /// Constructor
    TPrintArrays( int  aNumFlds, outField* aFlds, TIO& fout ):
        TRWArrays( aNumFlds, aFlds), out_format( fout )
    {}

    void writeComment( bool with_comments, const std::string& line )
    {
        if( with_comments )
            out_format.write_comment( line );
    }

    /// Writes long field to a text file.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    /// \param brief_mode - Do not write data items that contain only default values
    template < typename T >
    void writeField( long f_num, const T& value, bool with_comments, bool brief_mode  )
    {
        if( !brief_mode || getAlws( f_num ))
        {
            if( with_comments )
                out_format.write_comment( flds[f_num].comment );
            out_format.write_key_value( flds[f_num].name, value );
        }
    }

    /// Writes field to a text file.
    template < typename T >
    void addField( const std::string& field_name , const T& value )
    {
            out_format.write_key_value( field_name, value );
    }

    /// Writes field to a text file.
    void addField( const std::string& field_name , const char *value )
    {
            out_format.write_key_value( field_name, std::string(value) );
    }

    /// Writes array to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    /// \param brief_mode - Do not write data items that contain only default values
    template < typename T >
    void writeArray( long f_num,  T* arr,  long int size, long int l_size,
                     bool with_comments = false, bool brief_mode = false )
    {
        if( !brief_mode || getAlws(f_num ))
        {
            if( with_comments )
                out_format.write_comment( flds[f_num].comment );
            writeArray( flds[f_num].name,  arr, size, l_size );
        }
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num, const std::vector<double>& arr, long int l_size=0,
                     bool with_comments = false, bool brief_mode = false)
    {
        if( !brief_mode || getAlws(f_num) )
        {
            if( with_comments )
                out_format.write_comment( flds[f_num].comment );
            out_format.write_array( flds[f_num].name, arr, l_size );
        }
    }

    /// Writes char array to a text file.
    /// <flds[f_num].name> "arr[0]" ... "arr[size-1]"
    /// \param l_size - Setup number of characters in one element
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArrayF( long f_num, char* arr,  long int size, long int l_size,
                      bool with_comments = false, bool brief_mode = false )
    {
        if( !brief_mode || getAlws(f_num ) )
        {
            if( with_comments )
                out_format.write_comment( flds[f_num].comment.c_str() );
            writeArray( flds[f_num].name.c_str(), arr, size, l_size );
        }
    }

    /// Writes T array to a text file.
    template < typename T, typename LT >
    void writeArray(  const std::string& name, T*  arr, LT size, LT arr_size=static_cast<LT>(-1) )
    {
        out_format.write_array( name, arr, size, arr_size );
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void writeArray(  const std::string& name, T* arr, LT size, long int* sel_arr,
                     LT ncolumns=static_cast<LT>(1), LT l_size=static_cast<LT>(-1) )
    {
        if( !arr || !sel_arr )
            return;
        out_format.write_sel_array( name, arr, size, sel_arr, ncolumns, l_size );
    }

};


/// Read fields of structure
template < class TIO >
class TReadArrays : public  TRWArrays
{
    TIO& in_format;

public:

    /// Constructor
    TReadArrays( int aNumFlds, outField* aFlds, TIO& fin ):
        TRWArrays( aNumFlds, aFlds), in_format( fin ), current_readed("")
    {
       in_format.reset();
    }

    /// Read next name from file and find in fields list
    long int findNext()
    {
        long int ii;
        std::string name_key;
        while( in_format.has_next( name_key) )
        {
            ii = findFld( name_key );
            if( ii >= 0 )
            {
                flds[ii].readed = 1;
                return ii;
            }
        }
        return -3;
    }

    /// Read next label from file ( must be exist, otherwise error )
    void  readNext( const std::string& label)
    {
        if( !in_format.test_next_is( label )  )
        {
            std::string msg = label;
            msg += " - No data where expected.\n";
            msg += current_readed;
            Error( "TReadArrays read error 01", msg );
        }
    }

    /// Reads array from a TIO format.
    template <class T>
    void readArray(  const std::string& name, T* arr, long int size )
    {
        set_current_name( name, size );
        in_format.read_array( name, arr, size );
    }

    /// Reads strings array from a text file.
    void readArray( const std::string& name, char* arr, long int size, long int el_size )
    {
        set_current_name( name, size );
        in_format.read_strings_array( name, arr, size, el_size );
    }


    /// Reads double vector from a text file.
    void readArray( const std::string& name, std::vector<double> arr )
    {
        set_current_name( name, arr.size() );
        in_format.read_array( name, arr );
    }


protected:

    std::string current_readed;

    void set_current_name(  const std::string& name, long int size )
    {
       current_readed = "After successfully read <"+name+"> "+std::to_string(size)+" data items";
    }

};

}  // io_formats
