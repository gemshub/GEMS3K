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

#pragma once

#include <iostream>
#include <vector>
//#include <cmath>
#include "verror.h"

struct outField /// Internal descriptions of fields
 {
   std::string name; ///< name of field in structure
   long int alws;    ///< 1 - must be read, 0 - default values can be used
   long int readed;  ///< 0; set to 1 after reading the field from input file
   long int indexation;  ///< 1 - static object; 0 - undefined; <0 type of indexation, >1 number of elements in array
   std::string comment;

};

enum FormatType {
        ft_Value=0,   // value
        ft_F,   // command F
        ft_L,   // command L
        ft_R,   // command R
        ft_Internal
};

struct IOJFormat /// Internal descriptions of output/input formats with JSON notation
 {
   long int index;    ///< index formatted value into reading array
   long int type;  ///< type of formatted value { F, L, R, ...}
   std::string format; ///< string with formatted data for different type

   IOJFormat( long int aType, long int aIndex, std::string aFormat ):
               index(aIndex), type(aType), format(aFormat)
       {}

   IOJFormat( const IOJFormat& data ):
       index(data.index), type(data.type),  format(data.format)
       { }
};

class TRWArrays  /// Basic class for red/write fields of structure
 {
 protected:
    long int numFlds; ///< Size of array flds
    outField* flds;   ///< Array of permissible fields

 public:

    /// Constructor
     TRWArrays( short aNumFlds, outField* aFlds ):
         numFlds(aNumFlds), flds(aFlds)
    {}

     virtual ~TRWArrays()
    {}


    /// Find field by name
    virtual  long int findFld( const char *Name );

     /// Set the data object can be skipped from the file
     /// and default value(s) can be used
     /// \param ii index in array flds
    void  setNoAlws( long int ii )
    {  flds[ii].alws = 0; }

    /// Set the data object can be skipped from the file
    /// and default value(s) can be used
    /// \param Name of field in array flds
    void  setNoAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setNoAlws(ii);
    }

    /// Set the data object  must be always present in the file
    /// \param ii index in array flds
    void  setAlws( long int ii )
    {  flds[ii].alws = 1; }

    /// Set the data object Name must be always present in the file
    /// \param Name of field in array flds
    void  setAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
            setAlws(ii);
    }

    /// Test the data object  must be always present in the file
    /// \param ii index in array flds
    bool  getAlws( long int ii )
    {  return (flds[ii].alws == 1); }

    /// Test the data object Name must be always present in the file
    /// \param Name of field in array flds
    bool  getAlws( const char *Name )
    {
    	long int ii = findFld( Name );
         if( ii >=0 )
           return getAlws(ii);
         else
        return false;	 
    }

};

/// Print fields of structure outField
class TPrintArrays: public  TRWArrays
{
    std::iostream& ff;

public:

    /// Constructor
    TPrintArrays( short  aNumFlds, outField* aFlds, std::iostream& fout ):
        TRWArrays( aNumFlds, aFlds), ff( fout )
    {}

    template < typename T >
    void writeValue( const T& value )
    {
        ff << value;// << " ";
    }

    /// Writes long field to a text file.
    /// <flds[f_num].name> value
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    template < typename T >
    void writeField( long f_num, const T& value, bool with_comments, bool brief_mode  )
    {
        if(!brief_mode || getAlws( f_num ))
        {
            if( with_comments && flds[f_num].comment.length()>1)
                ff << std::endl << flds[f_num].comment.c_str();
            ff << std::endl << "<" << flds[f_num].name.c_str() << ">  ";
            writeValue(value);
        }
    }

    /// Writes array to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    template < typename T >
    void writeArray( long f_num,  T* arr,  long int size, long int l_size,
                    bool with_comments = false, bool brief_mode = false )
    {
      if(!brief_mode || getAlws(f_num ))
      {
        if( with_comments )
             ff <<  std::endl << flds[f_num].comment.c_str();
        writeArray( flds[f_num].name.c_str(),  arr, size, l_size );
      }
    }

    /// Writes double vector to a text file.
    /// <flds[f_num].name> arr[0] ... arr[size-1]
    /// \param l_size - Setup number of elements in line
    /// \param with_comments - Write files with comments for all data entries
    /// \param brief_mode - Do not write data items that contain only default values
    void writeArray( long f_num, const std::vector<double>& arr, long int l_size=0,
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
    void writeArray( const char *name, T*   arr, LT size, LT arr_size=static_cast<LT>(-1) )
    {
        int sz = 40;
        if( arr_size > 0 )
            sz = arr_size;

        ff << std::endl << "<" << name << ">" << std::endl;
        for( int ii=0, jj=0; ii<size; ii++, jj++  )
        {
            if(jj == sz) {
                jj=0;  ff << std::endl;
            }
            writeValue(arr[ii]);
            ff << " ";
        }
    }

    /// Writes char array to a text file.
    template < typename T=char, typename LT >
    void writeArray( const char *name, char*  arr, LT size, LT arr_size=static_cast<LT>(-1) )
    {
     bool isComment = false;

     if( name )
         ff << std::endl << "<" << name << ">" << std::endl;
     else
     {
         ff << std::endl << "#  ";
         isComment = true;
     }
     for( long int ii=0, jj=0; ii<size; ii++, jj++  )
     {
        if(jj == 40 )
        {
          jj=0;  ff << std::endl;
          if(isComment)
              ff << "#  ";
        }
        std::string str = std::string( arr +(ii*arr_size), 0, arr_size );
        writeValue(str);
        ff << " ";
     }
    }

    /// Writes selected elements from float array to a text file.
    template < typename T, typename LT >
    void writeArray( const char *name, T* arr, LT size, long int* selArr,
                     LT nColumns=static_cast<LT>(1), LT l_size=static_cast<LT>(-1) )
    {
        if(!arr)
            return;

        long int sz = 40;
        if( l_size > 0 )
            sz = l_size;

        ff << std::endl << "<" << name << ">" << std::endl;
        for( long int ii=0, jj=0; ii<size; ii++  )
        {
            for(long int cc=0; cc<nColumns; cc++ )
            {
                if(jj == sz){
                    jj=0;  ff << std::endl;
                }
                writeValue(arr[selArr[ii]*nColumns+cc]);
                ff << " ";
                jj++;
            }
        }
    }

    void writeArrayS( const char *name, char* arr, long int size, long int arr_siz );

};


 class TReadArrays : public  TRWArrays /// Read fields of structure
 {
     std::iostream& ff;
     std::string curArray;

 protected:
    /// Reads value from a text file.
    void readValue(float& val);
    /// Reads value from a text file.
    void readValue(double& val);
    /// Reads format value from a text file.
    long int readFormatValue(double& val, std::string& format);
    bool  readFormat( std::string& format );

    void setCurrentArray( const char* name, long int size );
 
 public:

    /// Constructor
    TReadArrays( short aNumFlds, outField* aFlds, std::iostream& fin ):
        TRWArrays( aNumFlds, aFlds), ff( fin ), curArray("")
    {}

    void  skipSpace();
    void reset();  ///< Reset to 0 all flags (readed)

    long int findFld( const char *Name ); ///< Find field by name
    long int findNext();  ///< Read next name from file and find in fields list
    long int findNextNotAll();  ///< Read next name from file and find in fields list (if doesnot find read next name)
    void  readNext( const char* label); ///< Read next name from file

    std::string testRead();   ///< Test for reading all fields must be always present in the file

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
    void readArray( const char* name, std::string &arr, long int el_size=198 );
    /// Reads double vector from a text file.
    void readArray( const char* name, std::vector<double> arr );

    void readFormatArray( const char* name, double* arr,
        long int size, std::vector<IOJFormat>& vFormats );
};

template <> void TPrintArrays::writeValue( const double& );
template <> void TPrintArrays::writeValue( const float& );
template <> void TPrintArrays::writeValue( const char& value );
template <> void TPrintArrays::writeValue( const std::string& value );


 //=============================================================================
