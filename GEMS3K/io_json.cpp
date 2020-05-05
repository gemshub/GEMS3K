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

/// Write double value to file
template <> void TPrintJson::writeValue( const char& value, nlohmann::json& json_arr )
{
    json_arr.push_back( std::string(1,value) );
}

/// Write double value to file
template <> void TPrintJson::writeValue( const gstring& value, nlohmann::json& json_arr )
{
    auto val = value;
    strip(val);
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

void TReadJson::setCurrentArray( const char* name, long int size )
{
    char buf[200];
    sprintf( buf, "After successfully read <%s> %ld data items", name, size);
    curArray = buf;
}

void TReadJson::reset()
{
    for(long int ii=0; ii < numFlds; ii++ )
        flds[ii].readed = 0;
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
    long int ii;
    std::string jkey;
    while( json_it != json_data.end() )
    {
        jkey = json_it.key();
        json_it++;
        ii = findFld( jkey.c_str()+1 );
        if( ii >= 0 )
        {
            flds[ii].readed = 1;
            return ii;
        }
    }
    return -3;
}


void TReadJson::readNext( const char* label )
{
    std::string jkey = key( label );
    gstring msg;
    if( json_data.find(jkey) == json_data.end() )
    {
        msg = label;
        msg += " - No data where expected.\n";
        msg += curArray;
        Error( "Json read error 01", msg );
    }
}


void TReadJson::readArray( const char* name, char* arr, long int size, long int el_size )
{
    setCurrentArray( name, size);

    std::string jkey = key( name );
    gstring msg;
    std::string val;
    auto json_arr_it = json_data.find(jkey);
    if( json_arr_it == json_data.end() )
    {
        msg = name;
        msg += " - No data where expected.\n";
        Error( "Json read error 01", msg );
    }

    if( !json_arr_it->is_structured() &&  size==1 )
    {
        json_arr_it->get_to(val);
        memcpy( arr, val.c_str(), el_size );
    }
    else
    {
        if( json_arr_it->size() != static_cast<size_t>(size) )
        {
            msg = name;
            msg += " - No size (";
            msg += std::to_string(size).c_str();
            msg += ") as expected ";
            msg += std::to_string(json_arr_it->size()).c_str();
            Error( "Json read error 04", msg );
        }
        for( long int ii=0; ii<size; ii++  )
        {
            json_arr_it->at(ii).get_to( val );
            memcpy( arr +(ii*el_size), val.c_str(), el_size );
        }
    }
}


//=============================================================================
// io_json.cpp
