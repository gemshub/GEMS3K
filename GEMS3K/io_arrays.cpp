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

#include "io_arrays.h"
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

long int TRWArrays::findFld( const char *Name )
{
    long int ii;

    for( ii=0; ii < numFlds; ii++ )
        if( !( strcmp( flds[ii].name.c_str(), Name ) ))
            return ii;
    return -1;
}

//---------------------------------------------------------//
// print Arrays ( fields of structure )


/// Write float value to file
template <> void TPrintArrays::writeValue( const float& val )
{
    if( IsFloatEmpty( val ))
        ff << CHAR_EMPTY << " ";
    else
        //    ff << setprecision(10) << scientific << arr[ii] << " ";
        ff << setprecision(15) << val;// << " ";
}

/// Write double value to file
template <> void TPrintArrays::writeValue( const double& val )
{
    if( IsDoubleEmpty( val ))
        ff << CHAR_EMPTY << " ";
    else
        //    ff << setprecision(18) << scientific << arr[ii] << " ";
        ff << setprecision(15) << val;// << " ";
}

/// Write double value to file
template <> void TPrintArrays::writeValue( const char& value )
{
    ff << "\'" << value << "\'";
}

/// Write double value to file
template <> void TPrintArrays::writeValue( const gstring& value )
{
    auto val = value;
#ifdef IPMGEMPLUGIN // 24/08/2010
    strip(val);
#else
    val.strip();
#endif
    ff  << "\'" << val.c_str() << "\'" /*<< " "*/;
    // commented out (space after text conflicts with gemsfit2 read-in) DM 16.07.2013
}

void TPrintArrays::writeArray( long f_num, const vector<double>& arr, long int l_size,
                               bool with_comments, bool brief_mode )
{
    long int jj;
    if(!brief_mode || getAlws(f_num ))
    {
        if( with_comments )
            ff <<  endl << flds[f_num].comment.c_str();

        int sz = 40;
        if( l_size > 0 )
            sz = l_size;

        ff << endl << "<" << flds[f_num].name.c_str() << ">" << endl;
        jj=0;
        for( size_t ii=0; ii<arr.size(); ii++, jj++  )
        {
            if(jj == sz) {
                jj=0;  ff << endl;
            }
            writeValue(arr[ii]);
            ff << " ";
        }
    }
}

void TPrintArrays::writeArrayF( long f_num, char* arr,
                                long int size, long int l_size, bool with_comments, bool brief_mode  )
{
    if(!brief_mode || getAlws(f_num ))
    {
        if( with_comments )
            ff <<  endl << flds[f_num].comment.c_str();
        writeArrayS( flds[f_num].name.c_str(),  arr,size, l_size);
    }
}


void TPrintArrays::writeArrayS( const char *name, char* arr,
                                long int size, long int arr_siz )
{
    writeArray( name, arr, size, arr_siz );
}

//------------------------------------------------------------------

inline void TReadArrays::setCurrentArray( const char* name, long int size )
{
    char buf[200];
    sprintf( buf, "After successfully read <%s> %ld data items", name, size);
    curArray = buf;
}

void TReadArrays::readValue(float& val)
{
    char input;
    skipSpace();

    ff.get( input );
    if( input == CHAR_EMPTY )
        val = FLOAT_EMPTY;
    else if( input =='i') //inf
    {
        ff.get( input );
        ff.get( input );
        val = FLOAT_EMPTY;
    }
    else
    {
        ff.putback(input);
        ff >> val;
    }
}

void TReadArrays::readValue(double& val)
{
    //char input;
    skipSpace();

    std::string buf;
    ff >> buf;
    if( buf == "`" || buf == "inf" || buf == "-inf" )
        val = DOUBLE_EMPTY;
    else
        val = std::atof( buf.c_str());

    /*ff.get( input );
    if( input == CHAR_EMPTY )
        val = DOUBLE_EMPTY;
    else if( input =='i') //inf
    {
        ff.get( input );
        ff.get( input );
        val = DOUBLE_EMPTY;
    }
    else
    {
        ff.putback(input);
        ff >> val;
    }*/
}

long int TReadArrays::readFormatValue(double& val, gstring& format)
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
void  TReadArrays::skipSpace()
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
bool  TReadArrays::readFormat( gstring& format )
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

void TReadArrays::reset()
{
    for(long int ii=0; ii < numFlds; ii++ )
        flds[ii].readed = 0;
}

long int TReadArrays::findFld( const char *Name )
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

long int TReadArrays::findNext()
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


void TReadArrays::readNext( const char* label)
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

long int TReadArrays::findNextNotAll()
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

void TReadArrays::readFormatArray( const char* name, double* arr,
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

void TReadArrays::readArray( const char* name, char* arr, long int size, long int el_size )
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

void TReadArrays::readArray( const char* name, vector<double> arr )
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
void TReadArrays::readArray( const char* name, gstring& arr, long int el_size )
{
    char ch;
    char buf[10000]; // DM changed form 400 to be able to read long character sections like DataSelect

    setCurrentArray( name, 1);
    skipSpace();
    ff.get(ch);
    ff.getline( buf, el_size+1, '\'');
    arr = buf;
}

gstring TReadArrays::testRead()
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
