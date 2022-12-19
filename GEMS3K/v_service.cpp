#include <istream>
#include <regex>
#include <cstring>
#include "v_service.h"


std::string char_array_to_string(const char* data_ptr, size_t max_size)
{
    std::string data_str;
    if(!data_ptr)
        return data_str;
    for(size_t pos=0; pos<max_size; ++pos)
    {
        if( data_ptr[pos] == '\0')
            break;
        data_str += data_ptr[pos];
    }
    return data_str;
}

void strip(std::string& str)
{
    trim(str);
}

void replace( std::string& str, const char* old_part, const char* new_part)
{
    size_t pos = str.find( old_part ); //rfind( old_part );
    if( pos != std::string::npos )
    {
        std::string res(str.substr(0, pos));
        res += new_part;
        res += str.substr( pos+strlen(old_part));
        str = res;
    }
}

void replaceall( std::string& str, const std::string& old_part, const std::string& new_part)
{
    size_t posb=0, pos = str.find( old_part ); //rfind( old_part );
    while( pos != std::string::npos )
    {
        std::string res(str.substr(0, pos));
        res += new_part;
        res += str.substr( pos+old_part.length() );
        str = res;
        posb = pos + new_part.length();
        pos = str.find( old_part, posb );
    }
}

void replaceall(std::string& str, char ch1, char ch2)
{
  for(size_t ii=0; ii<str.length(); ii++ )
   if( str[ii] == ch1 )
            str[ii] = ch2;
}

/// read string as: "<characters>"
std::istream& f_getline(std::istream& is, std::string& str, char delim)
{
    char ch;
    is.get(ch);
    str="";

    while( is.good() && ( ch==' ' || ch=='\n' || ch== '\t') )
        is.get(ch);
    if(ch == '\"')
        is.get(ch);
    while( is.good() &&  ch!=delim && ch!= '\"' )
    {
        str += ch;
        is.get(ch);
    }
    while( is.good() &&  ch!=delim )
            is.get(ch);

   return is;
}

std::string
u_makepath(const std::string& dir,  const std::string& name, const std::string& ext)
{
    std::string Path(dir);
    if( dir != "")
      Path += "/";
    Path += name;
    Path += ".";
    Path += ext;

    return Path;
}

std::string u_getpath( const std::string& file_path )
{
    std::size_t pos = file_path.find_last_of("/\\");
    if( pos != std::string::npos )
        return file_path.substr(0, pos);
    return "";
}


void u_splitpath(const std::string& Path, std::string& dir,
            std::string& name, std::string& ext)
{
    // Get path
    std::size_t pos = Path.find_last_of("/\\");
    if( pos != std::string::npos )
    {
        dir = Path.substr(0, pos);
        pos++;
    }
    else
    {
        dir = "";
        pos = 0;
    }
    name = Path.substr(pos);
    size_t pose = name.rfind(".");
    if( pose != std::string::npos )
    {
        ext = name.substr( pose+1, std::string::npos );
        name = name.substr(0, pose);
    }
    else
    {
        ext = "";
    }
}

std::string regexp_extract_string( std::string regstr, std::string data )
{
    std::string token = "";
    std::regex re( regstr );
    std::smatch match;

    if( std::regex_search( data, match, re ))
    {
        if (match.ready())
            token = match[1];
    }
    return token;
}

// Extract the string value by key from jsonstring
std::string extract_string_json( std::string key, std::string jsondata )
{
    size_t key_size = key.length()+2;
    std::string key_find = "\"" + key + "\"";
    std::string field_value;

    auto pos_set = jsondata.find( key_find, 0);
    while( pos_set != std::string::npos )
    {
       pos_set += key_size;
       while( isspace( jsondata[pos_set] ) )
           ++pos_set;
       if( jsondata[pos_set] == ':')
       {
          ++pos_set;
          while( isspace( jsondata[pos_set] ) )
               ++pos_set;
          if( jsondata[pos_set] == '\"')
          {
            ++pos_set;
            auto pos_end = jsondata.find_first_of( "\"", pos_set);
            field_value =  jsondata.substr( pos_set, pos_end-pos_set);
            break;
          }
       }
       pos_set = jsondata.find( key_find, pos_set);
    }
    return field_value;
}

// Extract the string value by key from jsonstring
std::string extract_string_json_old( std::string key, std::string jsondata )
{
    std::string data = jsondata;
    replace_all( data, "\'", '\"');
    std::string regstr =  std::string(".*\"")+key+"\"\\s*:\\s*\"([^\"]*)\".*";
    return regexp_extract_string( regstr, data );
}

// Extract the string value by key from query
int extract_int_json( const std::string& key, const std::string& jsondata )
{
    std::string data = jsondata;
    replace_all( data, "\'", '\"');
    std::string regstr =  std::string(".*\"")+key+"\"\\s*:\\s*([+-]?[1-9]\\d*|0).*";
    auto token = regexp_extract_string( regstr, data );
    if( token.empty() )
        return 0;
    return stoi(token);
}

std::vector<std::string> split(const std::string &str, const std::string &delimiters)
{
    std::vector<std::string> v;
    std::string vv;

    if( str.empty() )
        return v;

    std::string::size_type start = 0;
    auto pos = str.find_first_of(delimiters.c_str(), start);
    while(pos != std::string::npos)
    {
        vv = std::string(str, start, pos - start);
        strip( vv );
        v.push_back( vv );
        start = pos + 1;
        pos = str.find_first_of(delimiters.c_str(), start);
    }

    vv = std::string (str, start, str.length() - start);
    strip( vv );
    if( !vv.empty() )
        v.push_back( vv );
    return v;
}
