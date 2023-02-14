//-------------------------------------------------------------------
// $Id$
//
/// \file io_simdjson.cpp
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

#ifndef USE_NLOHMANNJSON

#include <cmath>
#include <string_view>
#include "simdjson/simdjson.h"
#include "simdjson/simdjson.cpp"
#include "io_simdjson.h"
#include "v_detail.h"
#include "v_service.h"

namespace  io_formats {


/// Read fields of structure
class SimdJsonImpl
{

public:

    /// Constructor
    explicit SimdJsonImpl( const std::string& json_string, const std::string& test_set_name,  const std::string& field_name ): json_data()
    {
        try {
            auto input_string = json_string;
            trim(input_string);
            trim(input_string, "[]");
            replaceall( input_string, "-inf", DOUBLE_INFMINUS_STR);
            replaceall( input_string, "inf", DOUBLE_INFPLUS_STR);
            replaceall( input_string, "nan", DOUBLE_NAN_STR);
            json_data = parser.parse( input_string );

            if( !test_set_name.empty() )
            {
                std::string_view json_set = json_data["set"];
                if( json_set.compare( test_set_name.c_str() )) {
                    gems_logger->warn(" Read the document from another set: {} , current set {}", json_set, test_set_name);
                }
            }
            json_data =  json_data[field_name];
            json_it = json_data.begin();
            //gems_logger->trace("{}", json_data);
        }
        catch( simdjson::simdjson_error& err )
        {
            gems_logger->error("SimdJson read error : {}.{}", std::to_string(err.error()), err.what());
            Error( std::string("SimdJson read error :") + std::to_string(err.error()), err.what() );
        }

    }

    /// Reset json loop
    void reset()
    {
        json_it = json_data.begin();
    }

    /// Read next name from file
    bool  has_next( std::string& next_field_name )
    {
        next_field_name.clear();
        if( json_it != json_data.end() )
        {
            auto key = json_it.key();
            next_field_name = std::string( key.begin(), key.end() );
            json_it++;
            return true;
        }
        return false;
    }

    /// Read next label from file ( must be exist, otherwise error )
    bool test_next_is( const std::string& label )
    {
        return  !json_data.at_key( label ).error();
    }

    /// Reads strings array from a text file.
    void read_strings_array( const std::string& field_name, char* arr, long int size, long int el_size );
    /// Reads double vector from a text file.
    void read_double_array( const std::string& name, std::vector<double>& arr );
    /// Reads int vector from a text file.
    void read_int_array(const std::string &field_name, std::vector<int64_t>& arr);

protected:

    // Internal structure of file data
    simdjson::dom::parser parser;
    simdjson::dom::object json_data;
    simdjson::dom::object::iterator json_it;

};

void SimdJsonImpl::read_strings_array(const std::string &field_name, char *arr, long size, long el_size)
{
    try {
        std::string_view val;
        simdjson::dom::element json_arr = json_data[field_name];

        if( json_arr.type() != simdjson::dom::element_type::ARRAY &&  size==1 )
        {
            val = json_arr;
            memcpy( arr, val.data(), el_size );
        }
        else
        {
            for( long int ii=0; ii<size; ++ii )
            {
                val = json_arr.at(ii);
                memcpy( arr +(ii*el_size), val.data(), el_size );
            }
        }
    }
    catch( simdjson::simdjson_error& err )
    {
        gems_logger->error("SimdJson read error : {}.{}", std::to_string(err.error()), err.what());
        Error( std::string("SimdJson read error :") + std::to_string(err.error()), err.what() );
    }
}

void SimdJsonImpl::read_double_array(const std::string &field_name, std::vector<double>& arr)
{
    try {
        arr.clear();
        simdjson::dom::element json_arr = json_data[field_name];

        if(  json_arr.type() == simdjson::dom::element_type::ARRAY  )
        {
            for (double value : json_arr)
                arr.push_back(value);
        }
        else
        {
            double value = json_arr;
            arr.push_back(value);
        }
    }
    catch( simdjson::simdjson_error& err )
    {
        gems_logger->error("SimdJson read error : {}.{}", std::to_string(err.error()), err.what());
        Error( std::string("SimdJson read error :") + std::to_string(err.error()), err.what() );
    }
}

void SimdJsonImpl::read_int_array(const std::string &field_name, std::vector<int64_t>& arr)
{
    try {
        arr.clear();
        simdjson::dom::element json_arr = json_data[field_name];

        if(  json_arr.type() == simdjson::dom::element_type::ARRAY  )
        {
            for (int64_t value : json_arr)
                arr.push_back(value);
        }
        else
        {
            int64_t value = json_arr;
            arr.push_back(value);
        }
    }
    catch( simdjson::simdjson_error& err )
    {
        gems_logger->error("SimdJson read error : {}.{}", std::to_string(err.error()), err.what());
        Error( std::string("SimdJson read error :") + std::to_string(err.error()), err.what() );
    }
}



//------------------------------------------------------------------------------------------


SimdJsonRead::SimdJsonRead(std::iostream &ff, const std::string& test_set_name,  const std::string& field_name): impl()
{
    std::stringstream buffer;
    buffer << ff.rdbuf();
    auto input_str = buffer.str();
    impl = std::make_shared<SimdJsonImpl>(input_str, test_set_name, field_name );
}

void SimdJsonRead::reset()
{
    impl->reset();
}

bool SimdJsonRead::has_next(std::string &next_field_name)
{
    return impl->has_next(next_field_name);
}

bool SimdJsonRead::test_next_is(const std::string &label)
{
    return  impl->test_next_is(label);
}

void SimdJsonRead::read_strings_array(const std::string &field_name, char *arr, long size, long el_size)
{
    impl->read_strings_array(field_name, arr, size, el_size );
}

void SimdJsonRead::read_double_array(const std::string &field_name, std::vector<double>& arr)
{
    impl->read_double_array( field_name, arr);
}

void SimdJsonRead::read_int_array(const std::string &field_name, std::vector<int64_t>& arr)
{
    impl->read_int_array( field_name, arr);
}

template <> float SimdJsonRead::internal_cast( double value )
{
    float cast_value;
    if( is_minusinf(value) )
        cast_value = InfMinus<float>();
    else if( is_plusinf(value) )
        cast_value = InfPlus<float>();
    else if( is_nan(value) )
        cast_value = Nan<float>();
    else
        cast_value = static_cast<float>(value);
    return cast_value;
}

//-------------------------------------------------------------------------------------


void SimdJsonWrite::put_head(const std::string &key_name, const std::string &field_name)
{
    first = true;
    fout << ( dense ? "[\n{" : "[{" );
    fout << ( dense ? "\n  " : "" );
    fout << "\"_key\"" << ( dense ? " : \"" : ":\"" ) << key_name << "\",";
    fout << ( dense ? "\n  " : "" );
    fout << "\"set\"" << ( dense ? " : \"" : ":\"" ) << current_set_name << "\",";
    fout << ( dense ? "\n  " : "" );
    fout << "\"" << field_name << ( dense ? "\" : " : "\":" );
    fout << "{";
}


void SimdJsonWrite::dump(bool)
{
    fout << ( dense ? "\n  " : "" );
    fout << ( dense ? "}\n}\n]" : "}}]" );

}

/// Write float value to file
template <> void SimdJsonWrite::add_value( const float& val )
{
    fout << floating_point_to_string( val );
}

/// Write double value to file
template <> void SimdJsonWrite::add_value( const double& val )
{
    fout << floating_point_to_string( val );
}

/// Write double value to file
template <> void SimdJsonWrite::add_value( const char& value )
{
    fout << "\"" << std::string(1,value) << "\"";
}

/// Write double value to file
template <> void SimdJsonWrite::add_value( const std::string& value )
{
    auto val = value;
    strip(val);
    fout  << "\"" << val << "\"";
}

void SimdJsonWrite::write_array(const std::string &field_name, const std::vector<double> &arr, long l_size)
{
    long jj=0, sz = ( l_size > 0 ? l_size: values_in_line );
    add_key( key( field_name ) );
    fout  << ( dense ? "[\n        " : "[\n" );

    for( size_t ii=0; ii<arr.size(); ii++, jj++ )
    {
        add_next( ii, jj, sz );
        add_value( arr[ii] );
    }
    fout << ( dense ? "\n    ]" : "\n]" );
}

}  // io_formats

// https://stackoverflow.com/questions/8610571/what-is-rvalue-reference-for-this/8610714#8610714

#endif
