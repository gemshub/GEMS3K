//-------------------------------------------------------------------
// $Id$
//
/// \file io_json.cpp
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

#ifdef USE_NLOHMANNJSON

#include <sstream>
#include "io_nlohmann.h"
#include "v_detail.h"
#include "v_service.h"

namespace  io_formats {

/// Write char value to file
template <> void NlohmannJsonWrite::add_value( const char& value, nlohmann::json& json_arr )
{
    json_arr.push_back( std::string(1,value) );
}

/// Write string value to file
template <> void NlohmannJsonWrite::add_value( const std::string& value, nlohmann::json& json_arr )
{
    auto val = value;
    strip(val);
    json_arr.push_back( val );
}

template <> void NlohmannJsonWrite::write_key_value( const std::string& field_name, const char& value )
{
    json_data[ key( field_name ) ] = std::string(1,value);
}

template <> void NlohmannJsonWrite::write_key_value( const std::string& field_name, const std::string& value )
{
    auto val = value;
    strip(val);
    json_data[ key( field_name ) ] = std::string(value);
}

NlohmannJsonRead::NlohmannJsonRead(std::iostream &ff, const std::string &test_set_name, const std::string &field_name):
    json_data()
{
    std::stringstream buffer;
    buffer << ff.rdbuf();
    auto input_string = buffer.str();
    replaceall( input_string, "-inf", DOUBLE_INFMINUS_STR);
    replaceall( input_string, "inf", DOUBLE_INFPLUS_STR);
    replaceall( input_string, "nan", DOUBLE_NAN_STR);
    nlohmann::json json_arr = nlohmann::json::parse(input_string);

    json_data = json_arr[0];
    if( !test_set_name.empty() )
    {
        auto  json_set = json_data["set"].get<std::string>();
        if( json_set.find(test_set_name) == std::string::npos ) {
            gems_logger->warn(" Read the document from another set: {} , current set {}", json_set, test_set_name);
        }
    }
    json_data =  json_data[field_name];
    json_it = json_data.begin();
}

bool NlohmannJsonRead::has_next(std::string &next_field_name)
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

void NlohmannJsonRead::read_strings_array(const std::string &field_name, char *arr, long size, long el_size)
{
    std::string msg, val;
    std::string jkey = key( field_name );

    auto json_arr_it = json_data.find(jkey);
    if( json_arr_it == json_data.end() )
    {
        msg = field_name;
        msg += " - No data where expected.\n";
        Error( "Json read error 03", msg );
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
            msg = field_name;
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

void NlohmannJsonRead::read_array(const std::string &field_name, std::vector<double> arr)
{
    std::string jkey = key( field_name );

    auto json_arr_it = json_data.find(jkey);
    if( json_arr_it == json_data.end() )
    {
        std::string msg;
        msg = field_name;
        msg += " - No data where expected.\n";
        Error( "Json read error 05", msg );
    }

    if( !json_arr_it->is_array()  )
        json_arr_it->get_to(arr);
}

}  // io_formats

#endif
