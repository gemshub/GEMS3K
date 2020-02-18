//-------------------------------------------------------------------
// $Id$
/// \file io_arrays.h
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


#include <nlohmann/json.hpp>
#include "io_arrays.h"



/// Print fields of structure outField
class TPrintJson: public  TRWArrays
{
    fstream ff;

    /// Internal structure of file data
    nlohmann::json& json_data;

    template <class T>
    void writeValue( const T& value, nlohmann::json& json_arr  )
    {
       json_arr.push_back(value);
    }

    std::string key( const char * name ) const
    {
      return std::string("<") + name + ">";
    }

public:

    /// Constructor
    TPrintJson( short  aNumFlds, outField* aFlds, nlohmann::json& json_out ):
        TRWArrays( aNumFlds, aFlds), json_data( json_out )
    {}

    /// Writes integral field to a json.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    template <class T>
    void writeField( long f_num, const T& value, bool /*with_comments*/, bool brief_mode  )
    {
        if( !brief_mode || getAlws( f_num ))
        {
            //if( with_comments && flds[f_num].comment.length()>1)
            //    ff << endl << flds[f_num].comment.c_str();
            json_data[ key( flds[f_num].name.c_str() ) ] = value;
        }
    }

    /// Writes array to a json.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    template < typename T >
    void writeArray( long f_num,  T* arr,  long int size, long int l_size,
                    bool /*with_comments*/ = false, bool brief_mode = false )
    {
      if(!brief_mode || getAlws(f_num ))
      {
        //if( with_comments )
        //     ff <<  endl << flds[f_num].comment.c_str();
        writeArray( flds[f_num].name.c_str(),  arr, size, l_size );
      }
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num, const vector<double>& arr, long int l_size=0,
                     bool with_comments = false, bool brief_mode = false);

    /// Writes char array to a text file.
    /// <flds[f_num].name> "arr[0]" ... "arr[size-1]"
    /// \param l_size - Setup number of characters in one element
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArrayF( long f_num, char* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false );

    /// Writes T array to a text file.
    template < typename T, typename LT >
    void writeArray( const char *name, T*   arr, LT size, LT =static_cast<LT>(-1) )
    {
        auto arr_key = key( name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( int ii=0; ii<size; ii++ )
        {
            writeValue(arr[ii], json_data[ arr_key ]);
        }
    }

    /// Writes char array to a json file.
    template < typename T=char, typename LT >
    void writeArray( const char *name, char*  arr, LT size, LT arr_size=static_cast<LT>(-1) )
    {
        if( !name )  // comment
          return;
        auto arr_key = key( name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( int ii=0, jj=0; ii<size; ii++, jj++  )
        {
            gstring str = gstring( arr +(ii*arr_size), 0, arr_size );
            writeValue(str, json_data[ arr_key ]);
        }
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void writeArray( const char *name, T* arr, LT size, long int* selArr,
                     LT nColumns=static_cast<LT>(1), LT =static_cast<LT>(-1) )
    {
        if( !arr )
            return;

        auto arr_key = key( name );
        json_data[ arr_key ] = nlohmann::json::array();
        for( long int ii=0; ii<size; ii++  )
        {
            for(long int cc=0; cc<nColumns; cc++ )
            {
                writeValue(arr[selArr[ii]*nColumns+cc], json_data[ arr_key ]);
            }
        }
    }

    void writeArrayS( const char *name, char* arr, long int size, long int arr_siz );

};


 class TReadJson : public  TRWArrays /// Read fields of structure
 {
    fstream& ff;
    gstring curArray;

 protected:
    /// Reads value from a text file.
    void readValue(float& val);
    /// Reads value from a text file.
    void readValue(double& val);
    /// Reads format value from a text file.
    long int readFormatValue(double& val, gstring& format);
    bool  readFormat( gstring& format );

    inline void setCurrentArray( const char* name, long int size );
 
 public:

    /// Constructor
    TReadJson( short aNumFlds, outField* aFlds, fstream& fin ):
        TRWArrays( aNumFlds, aFlds), ff( fin ), curArray("")
    {}

    void  skipSpace();
    void reset();  ///< Reset to 0 all flags (readed)

    long int findFld( const char *Name ); ///< Find field by name
    long int findNext();  ///< Read next name from file and find in fields list
    long int findNextNotAll();  ///< Read next name from file and find in fields list (if doesnot find read next name)
    void  readNext( const char* label); ///< Read next name from file

    gstring testRead();   ///< Test for reading all fields must be always present in the file

    /// Reads array from a text file.
    template <class T,
             typename std::enable_if<std::is_integral<T>::value,T>::type* = nullptr>
    void readArray( const char *name, T* arr, long int size )
    {
        setCurrentArray( name, size);
        for( long int ii=0; ii<size; ii++  )
        {
            skipSpace();
            ff >> arr[ii];
        }
    }

    /// Reads array from a text file.
    template <class T,
             typename std::enable_if<std::is_floating_point<T>::value,T>::type* = nullptr>
    void readArray( const char *name, T* arr, long int size )
    {
     setCurrentArray( name, size);
     for( long int ii=0; ii<size; ii++  )
        readValue(arr[ii]);
    }

    /// Reads array from a text file.
    void readArray( const char *name, char* arr, long int size, long int el_size );

    /// Reads string from a text file.
    void readArray( const char* name, gstring &arr, long int el_size=198 );
    /// Reads double vector from a text file.
    void readArray( const char* name, vector<double> arr );

    void readFormatArray( const char* name, double* arr,
        long int size, vector<IOJFormat>& vFormats );
};

template <> void TPrintJson::writeValue( const char& value, nlohmann::json& json_arr );
template <> void TPrintJson::writeValue( const gstring& value, nlohmann::json& json_arr );
template <> void TPrintJson::writeField( long f_num, const char& value, bool /*with_comments*/, bool brief_mode  );
template <> void TPrintJson::writeField( long f_num, const gstring& value, bool /*with_comments*/, bool brief_mode  );


 //=============================================================================
