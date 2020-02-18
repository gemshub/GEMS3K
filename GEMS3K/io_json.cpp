//-------------------------------------------------------------------
// $Id$
//
/// \file io_arrays.cpp
/// Implementation of service functions for writing/reading arrays in files
//
// Copyright (c) 2006-2012 S.Dmytriyeva
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

#include <iomanip>
#include  <iostream>

#include "io_json.h"
#include "v_user.h"

#ifdef IPMGEMPLUGIN

istream& f_getline(istream& is, gstring& str, char delim);

//    This constants should be 'defined' to satisfy all compilers
#define SHORT_EMPTY 	   -32768
#define LONG_EMPTY             -2147483648L
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

#else

#include "v_vals.h"

#endif


/// Write double value to file
template <> void TPrintJson::writeValue( const char& value, nlohmann::json& json_arr )
{
    json_arr.push_back( std::string(1,value) );
}

/// Write double value to file
template <> void TPrintJson::writeValue( const gstring& value, nlohmann::json& json_arr )
{
    auto val = value;
#ifdef IPMGEMPLUGIN
    strip(val);
#else
    val.strip();
#endif
    json_arr.push_back( std::string(val.c_str()) );
}

template <> void TPrintJson::writeField( long f_num, const char& value, bool /*with_comments*/, bool brief_mode  )
{
    if( !brief_mode || getAlws( f_num ))
    {
        //if( with_comments && flds[f_num].comment.length()>1)
        //    ff << endl << flds[f_num].comment.c_str();
        json_data[ key( flds[f_num].name.c_str() ) ] = std::string(1,value);
    }
}

template <> void TPrintJson::writeField( long f_num, const gstring& value, bool /*with_comments*/, bool brief_mode  )
{
    if( !brief_mode || getAlws( f_num ))
    {
        //if( with_comments && flds[f_num].comment.length()>1)
        //    ff << endl << flds[f_num].comment.c_str();
        json_data[ key( flds[f_num].name.c_str() ) ] = std::string(value.c_str());
    }
}

void TPrintJson::writeArray( long f_num, const vector<double>& arr, long int /*l_size*/,
                             bool /*with_comments*/, bool brief_mode )
{
    if(!brief_mode || getAlws(f_num ))
    {
        //if( with_comments )
        //    ff <<  endl << flds[f_num].comment.c_str();
        json_data[ key(flds[f_num].name.c_str())] = arr;
    }
}

void TPrintJson::writeArrayF( long f_num, char* arr, long int size, long int l_size,
                              bool /*with_comments*/, bool brief_mode  )
{
    if(!brief_mode || getAlws(f_num ))
    {
        //if( with_comments )
        //    ff <<  endl << flds[f_num].comment.c_str();
        writeArrayS( flds[f_num].name.c_str(),  arr,size, l_size);
    }
}

void TPrintJson::writeArrayS( const char *name, char* arr,
                              long int size, long int arr_siz )
{
    writeArray( name, arr, size, arr_siz );
}

//------------------------------------------------------------------

inline void TReadJson::setCurrentArray( const char* name, long int size )
{
    char buf[200];
    sprintf( buf, "After successfully read <%s> %ld data items", name, size);
    curArray = buf;
}

void TReadJson::readValue(float& val)
{
    char input;
    skipSpace();

    ff.get( input );
    if( input == CHAR_EMPTY )
        val = FLOAT_EMPTY;
    else
    {
        ff.putback(input);
        ff >> val;
    }
}

void TReadJson::readValue(double& val)
{
    char input;
    skipSpace();

    ff.get( input );
    if( input == CHAR_EMPTY )
        val = DOUBLE_EMPTY;
    else
    {
        ff.putback(input);
        ff >> val;
    }
}

long int TReadJson::readFormatValue(double& val, gstring& format)
{
    char input;
    format = "";
    long int type = ft_Value;

    skipSpace();

    if( ff.eof() )
        return ft_Internal;

    ff.get( input );
    switch( input )
    {
    case CHAR_EMPTY: val = DOUBLE_EMPTY;
        return type;
    case '<': ff.putback(input);
        return ft_Internal;
    case 'F': type = ft_F;
        if(!readFormat( format ))
            ff >> val;
        break;
    case 'L': type = ft_L;
        if(!readFormat( format ))
            ff >> val;
        break;
    case 'R': type = ft_R;
        if(!readFormat( format ))
            ff >> val;
        break;
    default:
    {   ff.putback(input);
        ff >> val;
        break;
    }
    }
    return type;
}

// skip  ' ',  '\n', '\t' and comments (from '#' to end of line)
void  TReadJson::skipSpace()
{
    char input;
    if( ff.eof() )
        return;
    ff.get( input );
    while( input == '#' || input == ' ' ||
           input == '\n' || input == '\t')
    {
        if( input == '#' )
            do{
            ff.get( input );
        }while( input != '\n' && input != '\0' && !ff.eof());
        if( input == '\0' || ff.eof() )
            return;
        ff.get( input );
    }
    ff.putback(input);
    // cout << ff << endl; // comented out DM 03.05.2013
}

// Read format string
bool  TReadJson::readFormat( gstring& format )
{
    char input;
    format = "";
    long int count1=0;  // counters of {}
    long int count2=0;  // counters of []

    skipSpace();
    ff.get( input );
    if( input != '{' && input != '[')
    {
        ff.putback(input);
        return false;  // no format string (only value)
    }
    do
    {
        format += input;
        if( input == '{')
            count1++;
        if( input == '}')
            count1--;
        if( input == '[')
            count2++;
        if( input == ']')
            count2--;
        ff.get( input );
        if( input == '\0' || ff.eof() )
            break;
    } while( count1 != 0 || count2 !=0 );

    return true;
}

void TReadJson::reset()
{
    for(long int ii=0; ii < numFlds; ii++ )
        flds[ii].readed = 0;
}

long int TReadJson::findFld( const char *Name )
{
    long int ii;
    gstring str = Name;
    size_t len = str.find('>');
    str = str.substr(0, len );

    for( ii=0; ii < numFlds; ii++ )
        if( !( strcmp( flds[ii].name.c_str(), str.c_str() ) ))
            return ii;
    return -1;
}

long int TReadJson::findNext()
{
    char buf[200];
    skipSpace();

    if( ff.eof() )
        return -3;

    ff >> buf;

    if( !( memcmp( "END_DIM", buf+1, 7 )) )
        return -2;

    long int ii = findFld( buf+1 );
    if(  ii < 0 )
    {
        gstring msg = buf;
        msg += " - Invalid label of data.\n";
        msg += curArray;
        Error( "Formatted read error 01", msg );
    }

    flds[ii].readed = 1;
    return ii;
}


void TReadJson::readNext( const char* label)
{
    char buf[200];
    gstring msg;
    skipSpace();

    if( ff.eof() )
    {
        msg = label;
        msg += " - No data where expected.\n";
        msg += curArray;
        Error( "Formatted read error 02", msg );
    }

    ff >> buf;
    gstring str = buf+1;
    size_t len = str.find('>');
    str = str.substr(0, len );

    if( !( strcmp( label, str.c_str() ) ))
        return;

    msg = buf;
    msg += " - Invalid label of data.\n";
    msg += curArray;
    Error( "Formatted read error 03", msg );

}

long int TReadJson::findNextNotAll()
{
    char bufx[200];
    char input;

    skipSpace();

    if( ff.eof() )
        return -3;

again:

    ff >> bufx;

    if( !( memcmp( "END_DIM", bufx+1, 7 )) )
        return -2;

    long int ii = findFld( bufx+1 );
    if(  ii < 0 )
    {
        do{
            ff.get( input );
            if( input == '#' )
            { ff.putback(input);
                skipSpace();
            }
        }while( input != '<' && input != '\0' && !ff.eof());
        if( input == '\0' || ff.eof() )
            return -3;
        ff.putback(input);
        goto again;
    }

    flds[ii].readed = 1;
    //cout << flds[ii].name << endl;
    return ii;
}

void TReadJson::readFormatArray( const char* name, double* arr,
                                   long int size, vector<IOJFormat>& vFormats )
{
    gstring format;
    long int type;

    setCurrentArray( name, size);
    // vFormats.clear();

    //ff << setprecision(15);
    for( long int ii=0; ii<size; ii++  )
    {
        type = readFormatValue(arr[ii], format);
        if(type > ft_Value && type< ft_Internal )
        {
            vFormats.push_back( IOJFormat(type, ii, format));
        }
    }
}

void TReadJson::readArray( const char* name, char* arr, long int size, long int el_size )
{
    char ch;
    char buf[200];

    setCurrentArray( name, size);

    for( long int ii=0; ii<size; ii++  )
    {
        skipSpace();
        ff.get(ch);
        //   while( ff.good() && ch != '\'' )
        //       ff.get(ch);
        ff.getline( buf, el_size+1, '\'');
        copyValues( arr +(ii*el_size), buf, el_size );
    }

}

void TReadJson::readArray( const char* name, vector<double> arr )
{
    int retSimb= 0; // next field is only value
    double value;
    gstring str;

    setCurrentArray( name, 0);
    //ff << setprecision(15);

    do{
        retSimb = readFormatValue(value, str);
        if( retSimb == ft_Value )
            arr.push_back( value );
    } while(retSimb == ft_Value );
}

// DM corrected added gstring& arr instead of gstring arr 18.04.2013
void TReadJson::readArray( const char* name, gstring& arr, long int el_size )
{
    char ch;
    char buf[10000]; // DM changed form 400 to be able to read long character sections like DataSelect

    setCurrentArray( name, 1);
    skipSpace();
    ff.get(ch);
    ff.getline( buf, el_size+1, '\'');
    arr = buf;
}

gstring TReadJson::testRead()
{
    gstring ret = "";
    for(long int ii=0; ii < numFlds; ii++ )
        if( flds[ii].alws==1 && flds[ii].readed != 1 )
        {  if( !ret.empty() )
                ret += ", ";
            ret += flds[ii].name;
        }
    return ret;
}

//=============================================================================
// io_arrays.cpp
