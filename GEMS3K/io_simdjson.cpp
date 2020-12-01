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

#include <cmath>
#include "io_simdjson.h"
#include "simdjson/simdjson.h"
#include "simdjson/simdjson.cpp"
#include "v_detail.h"

namespace  io_formats {


/// Read fields of structure
class SimdJsonImpl
{

public:

    /// Constructor
    explicit SimdJsonImpl( const std::string& json_string ): json_data()
    {
        auto input_string = json_string;
        replaceall( input_string, "inf", "0");
        std::cout << input_string << std::endl;
        auto error = parser.parse(input_string).get(json_data); // do the parsing
        test_simdjson_error( error );
        json_it = json_data.begin();
        //    std::cout <<  json_data << std::endl;
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
    void read_array( const std::string& name, std::vector<double>& arr );
    /// Reads int vector from a text file.
    void read_array(const std::string &field_name, std::vector<int64_t>& arr);

protected:

    // Internal structure of file data
    simdjson::dom::parser parser;
    simdjson::dom::object json_data;
    simdjson::dom::object::iterator json_it;

    void test_simdjson_error( simdjson::error_code  error ) const
    {
        ErrorIf( error, std::string("SimdJson read error :") + std::to_string(error) ,
                 simdjson::error_message(error) );
    }

};

void SimdJsonImpl::read_strings_array(const std::string &field_name, char *arr, long size, long el_size)
{
        std::string msg;
        std::string_view val;

        simdjson::dom::element json_arr;
        auto error = json_data.at_key(field_name).get(json_arr);
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

void SimdJsonImpl::read_array(const std::string &field_name, std::vector<double>& arr)
{
    arr.clear();

    simdjson::dom::element json_arr;
    auto error = json_data.at_key(field_name).get(json_arr);
    test_simdjson_error( error );

    if(  json_arr.type() == simdjson::dom::element_type::ARRAY  )
    {
        for (double value : json_arr)
        {
            arr.push_back(value);
        }
    }
    else
    {
        double value;
        error = json_arr.get(value);
        test_simdjson_error( error );
        arr.push_back(value);
    }
}

void SimdJsonImpl::read_array(const std::string &field_name, std::vector<int64_t>& arr)
{
    arr.clear();

    simdjson::dom::element json_arr;
    auto error = json_data.at_key(field_name).get(json_arr);
    test_simdjson_error( error );

    if(  json_arr.type() == simdjson::dom::element_type::ARRAY  )
    {
        for (int64_t value : json_arr)
        {
            arr.push_back(value);
        }
    }
    else
    {
        int64_t value;
        error = json_arr.get(value);
        test_simdjson_error( error );
        arr.push_back(value);
    }
}


//------------------------------------------------------------------------------------------


SimdJsonRead::SimdJsonRead(std::iostream &ff): impl()
{
    std::stringstream buffer;
    buffer << ff.rdbuf();
    auto input_str = buffer.str();
    impl = std::make_shared<SimdJsonImpl>(input_str);
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

void SimdJsonRead::read_array(const std::string &field_name, std::vector<double>& arr)
{
    impl->read_array( field_name, arr);
}

void SimdJsonRead::read_array(const std::string &field_name, std::vector<int64_t>& arr)
{
    impl->read_array( field_name, arr);
}


//-------------------------------------------------------------------------------------


/// Write float value to file
template <> void SimdJsonWrite::add_value( const float& val )
{
  fout << std::setprecision(15) << val;
}

/// Write double value to file
template <> void SimdJsonWrite::add_value( const double& val )
{
   fout << std::setprecision(15) << val;
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

