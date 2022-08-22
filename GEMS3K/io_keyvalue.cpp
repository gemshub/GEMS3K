//-------------------------------------------------------------------
// $Id$
//
/// \file io_keyvalue.cpp
/// Implementation of service functions for writing/reading arrays in files
//
// Copyright (c) 2020 S.Dmytriyeva
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

#include "io_keyvalue.h"
#include "v_detail.h"
#include "v_service.h"

namespace  io_formats {

/// Write float value to file
template <> void KeyValueWrite::writeValue( const float& val )
{
    if( IsFloatEmpty( val ))
        fout << CHAR_EMPTY << " ";
    else
        fout << floating_point_to_string( val );
}

/// Write double value to file
template <> void KeyValueWrite::writeValue( const double& val )
{
    if( IsDoubleEmpty( val ))
        fout << CHAR_EMPTY << " ";
    else
        fout << floating_point_to_string( val );
}

/// Write double value to file
template <> void KeyValueWrite::writeValue( const char& value )
{
    fout << "\'" << value << "\'";
}

/// Write double value to file
template <> void KeyValueWrite::writeValue( const std::string& value )
{
    auto val = value;
    strip(val);
    fout  << "\'" << val << "\'";
}

template <> void KeyValueWrite::write_key_value( const std::string& field_name, const std::string& value )
{
    fout << std::endl << "<" << field_name << ">  ";
    fout  << "\"" << value << "\"";
}

void KeyValueWrite::write_array(const std::string &field_name, const std::vector<double> &arr, long l_size )
{
    int sz = ( l_size > 0 ? l_size: 40 );
    int jj = 0;
    fout << std::endl << "<" << field_name << ">" << std::endl;

    for( size_t ii=0; ii<arr.size(); ii++, jj++  )
    {
        if(jj == sz)
        {
            jj=0;
            fout << std::endl;
        }
        writeValue( arr[ii] );
        fout << " ";
    }
}


//------------------------------------------------------------------


void KeyValueRead::read_value(float& val)
{
    skip_space();
    std::string buf;
    fin >> buf;
    if( buf == "`" )
        val = FLOAT_EMPTY;
    else
        val = string_to_floating_point<float>( buf );
}

void KeyValueRead::read_value(double& val)
{
    skip_space();
    std::string buf;
    fin >> buf;
    if( buf == "`" )
        val = DOUBLE_EMPTY;
    else
        val = string_to_floating_point<double>( buf );
}


// skip  ' ',  '\n', '\t' and comments (from '#' to end of line)
void  KeyValueRead::skip_space()
{
    char input;
    if( fin.eof() )
        return;
    fin.get( input );
    while( input == '#' || isspace(input) )
    {
        if( input == '#' )
        {
            do{
                fin.get( input );
            }
            while( input != '\n' && input != '\0' && !fin.eof());
        }
        if( input == '\0' || fin.eof() )
            return;
        fin.get( input );
    }
    fin.putback(input);
}


bool KeyValueRead::has_next(std::string &next_field_name)
{
    skip_space();
    if( fin.eof() )
        return false;

    //    Try skip not readed
    //    char input;
    //    do{
    //        fin.get( input );
    //        if( input == '#' )
    //        {
    //            fin.putback(input);
    //            skip_space();
    //        }
    //    }
    //    while( input != '<' && input != '\0' && !fin.eof() );

    //    if( input == '\0' || fin.eof() )
    //        return false;

    //    fin.putback(input);

    fin >> next_field_name;
    trim(next_field_name, "<>");
    return( next_field_name != "END_DIM" );
}

bool KeyValueRead::test_next_is(const std::string &label)
{
    skip_space();
    if( fin.eof() )
        return false;

    std::string next_field_name;
    fin >> next_field_name;
    trim(next_field_name, "<>");

    return label==next_field_name;
}

void KeyValueRead::read_strings_array( const std::string&, char *arr, long size, long el_size )
{
    char ch, buf[200];
    for( long int ii=0; ii<size; ii++  )
    {
        skip_space();
        fin.get(ch);
        //   while( fin.good() && ch != '\'' && ch != '\"' )
        //       fin.get(ch);
        fin.getline( buf, el_size+1, ch /*'\''*/);
        copyValues( arr +(ii*el_size), buf, el_size );
    }
}

void KeyValueRead::read_array( const std::string &, std::vector<double> arr )
{
    double dval;
    char input;

    if( fin.eof() )
        return;
    fin.get( input );

    while( input != '<' &&  input != '\0')
    {
        fin.putback(input);
        read_value( dval );
        arr.push_back( dval );
        fin.get( input );
        if(  fin.eof() )
            break;
    }
}

bool KeyValueRead::skip_line()
{
    char ch;
    skip_space();
    fin.get(ch);
    if( ch ==  '<' )
    {
        fin.putback(ch);
        return false;
    }
    char buf[300];
    fin.getline( buf, 300, '\n');
    return true;
}

}
