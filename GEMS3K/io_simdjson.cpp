//-------------------------------------------------------------------
// $Id$
//
/// \file io_json.cpp
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

#ifndef USE_OLD_KV_IO_FILES

#include "io_simdjson.h"
//#include "simdjson.cpp"
#include "v_detail.h"


namespace  io_formats {

/// Write char value to file
template <> void SimdJsonWrite::add_value( const char& value, nlohmann::json& json_arr )
{
    json_arr.push_back( std::string(1,value) );
}

/// Write string value to file
template <> void SimdJsonWrite::add_value( const std::string& value, nlohmann::json& json_arr )
{
    auto val = value;
    strip(val);
    json_arr.push_back( val );
}

template <> void SimdJsonWrite::write_key_value( const std::string& field_name, const char& value )
{
    json_data[ key( field_name ) ] = std::string(1,value);
}

template <> void SimdJsonWrite::write_key_value( const std::string& field_name, const std::string& value )
{
    auto val = value;
    strip(val);
    json_data[ key( field_name ) ] = std::string(value);
}

//------------------------------------------------------------------------------------------

io_formats::SimdJsonRead::SimdJsonRead(std::iostream &ff): json_data()
{
    std::stringstream buffer;
    buffer << ff.rdbuf();
    auto input_str = buffer.str();

    simdjson::dom::parser parser;
    auto error = parser.parse(input_str).get(json_data); // do the parsing
    test_simdjson_error( error );
    json_it = json_data.begin();
}

bool SimdJsonRead::has_next(std::string &next_field_name)
{
    next_field_name.clear();
    if( json_it != json_data.end() )
    {
        next_field_name = json_it.key();
        trim(next_field_name, "<>");
        json_it++;
        return true;
    }
    return false;
}

void SimdJsonRead::read_strings_array(const std::string &field_name, char *arr, long size, long el_size)
{
        std::string msg;
        std::string_view val;
        std::string jkey = key( field_name );

        simdjson::dom::element json_arr;
        auto error = json_data.at_key(jkey).get(json_arr);
        test_simdjson_error( error );

        if( json_arr.type() != simdjson::dom::element_type::ARRAY &&  size==1 )
        {
            error = json_arr.get( val );
            test_simdjson_error( error );
            memcpy( arr, val.data(), el_size );
        }
        else
        {
            for( long int ii=0; ii<size; ++ii )
            {
                error = json_arr.at(ii).get(val);
                // error if different size
                test_simdjson_error( error );
                memcpy( arr +(ii*el_size), val.data(), el_size );
            }
        }
}

void SimdJsonRead::read_array(const std::string &field_name, std::vector<double> arr)
{
    double value;
    std::string jkey = key( field_name );
    arr.clear();

    simdjson::dom::element json_arr;
    auto error = json_data.at_key(jkey).get(json_arr);
    test_simdjson_error( error );

    if(  json_arr.type() == simdjson::dom::element_type::ARRAY  )
    {
        for (simdjson::dom::element arr_element : json_arr)
        {
            error = arr_element.get(value);
            test_simdjson_error( error );
            arr.push_back(value);
        }
    }
}


}  // io_formats

#endif
