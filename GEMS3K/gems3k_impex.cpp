#include "fstream"
#include "sstream"
#include "gems3k_impex.h"




std::string GEMS3KImpexGenerator::get_dbr_file_lst_path() const
{
    return  u_makepath( impex_dir, base_name + "-dbr", "lst" );
}


std::string GEMS3KImpexGenerator::gen_dbr_file_name(int time_point, size_t index) const
{
    char buf[5];
    snprintf( buf, 5, "%4.4ld", index );
    std::string dbr_name = base_name + "-dbr-";
    dbr_name += std::to_string(time_point) + "-" + buf +".";
    dbr_name +=  extension();
    return dbr_name;
}

std::string GEMS3KImpexGenerator::gen_dat_lst_head()
{
    std::stringstream fout;
    fout << mode() << " \"" << gen_dch_file_name() + "\"";
    fout << " \"" << gen_ipm_file_name() << "\" ";
    return fout.str();
}

void GEMS3KImpexGenerator::set_internal_data()
{
    std::string ext;
    u_splitpath( ipmfiles_lst_name, impex_dir, base_name, ext );
    auto pos = base_name.rfind("-");
    if( pos != std::string::npos )
        base_name = base_name.substr(0, pos);
}

void GEMS3KImpexGenerator::load_dat_lst_file()
{
    std::string mode, dbr_name;
    std::fstream f_lst( ipmfiles_lst_name, std::ios::in );
    ErrorIf( !f_lst.good() , ipmfiles_lst_name, "Fileopen error");
    getline(f_lst, mode, ' ');
    getline(f_lst, datach_file_name, ' ');
    trim(datach_file_name, "\"");
    getline(f_lst, ipm_file_name, ' ');
    trim(ipm_file_name, "\"");

    while( getline(f_lst, dbr_name, ' ') )
    {
        trim( dbr_name );
        trim( dbr_name, "\"");
        if( !dbr_name.empty() )
              databr_file_names.push_back(dbr_name);
    }
    nIV = databr_file_names.size();
    get_mode( mode );
}

void GEMS3KImpexGenerator::get_mode( const std::string &str_mode )
{
    io_mode = f_key_value;
    if( str_mode == "-b" )
        io_mode = f_binary;
    else  if( str_mode == "-j" )
        io_mode = f_json;
}

std::string GEMS3KImpexGenerator::mode() const
{
    switch( io_mode )
    {
    case f_binary:
        return "-b";
    case f_json:
        return "-j";
    default:
    case f_key_value:
        break;
    }
    return "-t";
}
