#ifndef GEMS3K_IMPEX_H
#define GEMS3K_IMPEX_H

#include "sstream"
#include "vector"
#include "v_detail.h"

//#ifndef USE_OLD_KV_IO_FILES
//  const char *dat_ext = "json";
//  const char *dat_filt = "*.json";
//#else
//const char *dat_ext = "dat";
//const char *dat_filt = "*.dat";
//#endif


/// Descripton of data to generate MULTI, DATACH and DATABR files structure prepared from GEMS.
class GEMS3KImpexGenerator
{

public:

    /// These are used io formats
    enum FileIOModes {
        f_binary,
        f_key_value,
        f_json
    };


public:

    /// Constructor
    /// Generate MULTI, DATACH and DATABR files structure prepared from GEMS.
    /// Prints files for separate coupled FMT-GEM programs that use GEMS3K module
    /// \param filepath - IPM work structure file path&name
    /// \param nIV - Number of allocated nodes
    /// \param bin_mode - Write IPM, DCH and DBR files in binary, txt or json mode
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments -Write files with comments for all data entries ( in text mode)
    /// \param addMui - Print internal indices in RMULTS to IPM file for reading into Gems back
    explicit GEMS3KImpexGenerator(  const std::string& filepath, long int anIV, FileIOModes file_mode ):
        ipmfiles_lst_name(filepath), nIV(anIV), io_mode(file_mode)
    {
        std::string ext;
        u_splitpath( ipmfiles_lst_name, impex_dir, base_name, ext );
        auto pos = base_name.rfind("-");
        if( pos != std::string::npos )
            base_name = base_name.substr(0, pos);
    }

    /// Get selected file output mode
    FileIOModes files_mode() const
    {
        return io_mode;
    }

    /// Generate MULTI file name
    std::string get_dbr_file_lst_path() const
    {
        return  u_makepath( impex_dir, base_name + "-dbr", "lst" );
    }

    /// Generate full path
    std::string get_path( const std::string file_name ) const
    {
        return  ( impex_dir.empty() ? file_name : impex_dir+"/"+file_name );
    }

    /// Generate full path for ipm file
    std::string get_ipm_path()
    {
        return  ( impex_dir.empty() ? ipm_file_name : impex_dir+"/"+ipm_file_name );
    }

    /// Generate full path for ipm file
    std::string get_dch_path() const
    {
        return  ( impex_dir.empty() ? datach_file_name : impex_dir+"/"+datach_file_name );
    }

    /// Generate MULTI file name
    std::string gen_ipm_file_name()
    {
        ipm_file_name = base_name + "-ipm." + extension();
        return ipm_file_name;
    }

    /// Generate dataCH file name
    std::string gen_dch_file_name()
    {
        datach_file_name = base_name + "-dch." + extension();
        return datach_file_name;
    }

    /// Generate dataBR file name
    std::string gen_dbr_file_name( int time_point, size_t index ) const
    {
        char buf[5];
        snprintf( buf, 5, "%4.4ld", index );
        std::string dbr_name = base_name + "-dbr-";
        dbr_name += std::to_string(time_point) + "-" + buf +".";
        dbr_name +=  extension();
        return dbr_name;
    }

    /// Generate *-dat.lst file data
    std::string gen_dat_lst_head()
    {
        std::stringstream fout;
        fout << mode() << " \"" << gen_dch_file_name() + "\"";
        fout << " \"" << gen_ipm_file_name() << "\" ";
        return fout.str();
    }

protected:

    /// IPM work structure file path&name
    std::string ipmfiles_lst_name;

    /// Number of allocated nodes
    long int nIV = 1;

    /// Write IPM, DCH and DBR files in binary, txt or json mode)
    FileIOModes io_mode = f_json;

    std::string impex_dir;
    std::string base_name;

    std::string ipm_file_name;
    std::string datach_file_name;
    std::vector<std::string> databr_file_names;

    std::string extension() const
    {
        switch( io_mode )
        {
        case f_binary:
            return "bin";
        case f_json:
            return "json";
        default:
        case f_key_value:
            break;
        }
        return "dat";
    }

    std::string mode() const
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

};



#endif // GEMS3K_IMPEX_H
