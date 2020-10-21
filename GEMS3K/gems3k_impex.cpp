#include "fstream"
#include "sstream"
#include "gems3k_impex.h"
#include "v_detail.h"


GEMS3KGenerator::IOModes GEMS3KGenerator::default_type_f =
        GEMS3KGenerator::f_json;


std::string GEMS3KGenerator::get_dbr_file_lst_path() const
{
    return  u_makepath( impex_dir, base_name + "-dbr", "lst" );
}


std::string GEMS3KGenerator::gen_dbr_file_name(int time_point, size_t index) const
{
    char buf[10];
    snprintf( buf, 5, "%4.4ld", index );
    std::string dbr_name = base_name + "-dbr-";
    dbr_name += std::to_string(time_point) + "-" + buf +".";
    dbr_name +=  extension();
    return dbr_name;
}

std::string GEMS3KGenerator::gen_dat_lst_head()
{
    std::stringstream fout;
    fout << mode() << " \"" << gen_dch_file_name() + "\"";
    fout << " \"" << gen_ipm_file_name() << "\" ";
    return fout.str();
}

void GEMS3KGenerator::set_internal_data()
{
    std::string ext;
    u_splitpath( ipmfiles_lst_name, impex_dir, base_name, ext );
    auto pos = base_name.rfind("-");
    if( pos != std::string::npos )
        base_name = base_name.substr(0, pos);
}

void GEMS3KGenerator::load_dat_lst_file()
{
    std::string mode, dbr_name;
    std::fstream f_lst( ipmfiles_lst_name, std::ios::in );
    ErrorIf( !f_lst.good() , ipmfiles_lst_name, " fileopen error");

    f_getline(f_lst, mode, ' ');
    trim(mode);
    get_mode( mode );

    f_getline( f_lst, datach_file_name, ' ');
    trim(datach_file_name, "\"");
    f_getline( f_lst, ipm_file_name, ' ');
    trim(ipm_file_name, "\"");

    while( f_lst.good() )
    {
        f_getline( f_lst, dbr_name, ' ');
        trim( dbr_name );
        trim( dbr_name, "\"");
        if( !dbr_name.empty() )
              databr_file_names.push_back(dbr_name);
    }
    nIV = databr_file_names.size();
}

void GEMS3KGenerator::get_mode( const std::string &str_mode )
{
    io_mode = f_key_value;
    if( str_mode == "-b" )
        io_mode = f_binary;
    else  if( str_mode == "-j" )
        io_mode = f_json;
#ifdef USE_OLD_NLOHMANJSON
    else  if( str_mode == "-n" )
        io_mode = f_nlohmanjson;
#endif
}

std::string GEMS3KGenerator::mode() const
{
    switch( io_mode )
    {
    case f_binary:
        return "-b";
    case f_nlohmanjson:
#ifdef USE_OLD_NLOHMANJSON
        return "-n";
#endif
    case f_json:
        return "-j";
    default:
    case f_key_value:
        break;
    }
    return "-t";
}

size_t GEMS3KGenerator::load_dbr_lst_file( const std::string& dbr_lst_path )
{
    std::string dbr_name, dbr_name_full;
    if( dbr_lst_path.find_first_of("/\\") == std::string::npos )
       dbr_name_full = impex_dir+"/"+dbr_lst_path;
    else
       dbr_name_full = dbr_lst_path;

    std::fstream f_lst( dbr_name_full, std::ios::in );
    ErrorIf( !f_lst.good() , dbr_name_full, " fileopen error");

    databr_file_names.clear();
    while( f_lst.good() )
    {
        f_getline( f_lst, dbr_name, ',' );
        trim( dbr_name );
        trim( dbr_name, "\"");
        if( !dbr_name.empty() )
              databr_file_names.push_back(dbr_name);
    }
    return databr_file_names.size();
}

