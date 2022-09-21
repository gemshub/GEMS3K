#include <fstream>
#include <sstream>
#include "gems3k_impex.h"
#include "v_detail.h"
#include "v_service.h"


GEMS3KGenerator::IOModes GEMS3KGenerator::default_type_f =
        GEMS3KGenerator::f_json;


std::string GEMS3KGenerator::get_dbr_file_lst_path() const
{
    return  u_makepath( impex_dir, base_name + "-dbr", "lst" );
}


std::string GEMS3KGenerator::gen_dbr_name(const std::string the_name, int time_point, size_t index)
{
    char buf[10];
    snprintf( buf, 5, "%4.4zu", index );
    std::string dbr_name = the_name + "-dbr-";
    dbr_name += std::to_string(time_point) + "-" + buf;
    return dbr_name;
}

GEMS3KGenerator::GEMS3KGenerator(const std::string &filepath, long anIV, IOModes file_mode):
    ipmfiles_lst_name(filepath), nIV(anIV), io_mode(file_mode)
{
#ifndef USE_THERMOFUN
    ErrorIf( io_mode>=f_thermofun, ipmfiles_lst_name, " ThermoFun as an option is hidden");
#endif
    set_internal_data();
}

std::string GEMS3KGenerator::gen_dbr_file_name(int time_point, size_t index) const
{
    std::string dbr_name = gen_dbr_name(base_name, time_point, index);
    dbr_name +=  ".";
    dbr_name +=  extension();
    return dbr_name;
}

std::string GEMS3KGenerator::gen_dat_lst_head()
{
    std::stringstream fout;
    fout << mode() << " \"" << gen_dch_file_name() + "\"";
    fout << " \"" << gen_ipm_file_name() << "\" ";
    if( files_mode() >= f_thermofun )
    {
        fout << " \"" << gen_thermofun_file_name() << "\" ";
    }
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
    if( mode[0] == '-' )
    {
        get_mode( mode );
        f_getline( f_lst, datach_file_name, ' ');
        trim(datach_file_name, "\"");
    }
    else
    {
        io_mode = f_key_value;
        datach_file_name = mode;
        trim(datach_file_name, "\"");
    }

    f_getline( f_lst, ipm_file_name, ' ');
    trim(ipm_file_name, "\"");

    if( files_mode() >= f_thermofun )
    {
#ifdef USE_THERMOFUN
        f_getline( f_lst, thermofun_file_name, ' ');
        trim(thermofun_file_name, "\"");
#else
        Error( ipmfiles_lst_name, " ThermoFun as an option is hidden");
#endif
    }

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
    else  if( str_mode == "-f" || str_mode == "-fun" )
        io_mode = f_thermofun;
    else  if( str_mode == "-o" || str_mode == "-fun-kv" )
        io_mode = f_kv_thermofun;
}

std::string GEMS3KGenerator::mode() const
{
    switch( io_mode )
    {
    case f_binary:
        return "-b";
    case f_json:
        return "-j";
    case f_thermofun:
        return "-f";
    case f_kv_thermofun:
        return "-o";
    default:
    case f_key_value:
        break;
    }
    return "-t";
}

size_t GEMS3KGenerator::load_dbr_lst_file( const std::string& dbr_lst_path )
{
    std::string dbr_name, dbr_name_full;
    if( dbr_lst_path.find_first_of("/\\") == std::string::npos && !impex_dir.empty() )
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

