
#pragma once

#include <vector>
#include <string>

/// Descripton of data to generate MULTI, DATACH and DATABR files structure prepared from GEMS.
class GEMS3KGenerator
{

public:

    /// These are used io formats
    enum IOModes {
        f_key_value,
        f_binary,
        f_json,
        f_thermofun,     // export thermodynamic data into ThermoFun JSON format file, other in json format
        f_kv_thermofun   // export thermodynamic data into ThermoFun JSON format file, other in key_value format
    };

    /// Default output/exchange file types
    static IOModes default_type_f;

    /// Default output/exchange file extension
    static std::string default_ext()
    {
      return ext( default_type_f );
    }

    static std::string ext( IOModes type )
    {
        switch( type )
        {
        case f_binary:
            return "bin";
        case f_json:
        case f_thermofun:
            return "json";
        default:
        case f_kv_thermofun:
        case f_key_value:
            break;
        }
        return "dat";
    }

    /// Generate MULTI name
    static std::string gen_ipm_name( const std::string the_name )
    {
        return the_name + "-ipm";
    }

    /// Generate dataCH name
    static std::string gen_dch_name( const std::string the_name )
    {
        return the_name + "-dch";
    }

    /// Generate ThermoFun JSON format file name
    static std::string gen_thermofun_name( const std::string the_name )
    {
        return the_name + "-fun";
    }

    /// Generate dataBR name
    static std::string gen_dbr_name( const std::string the_name, int time_point, size_t index );


    /// Constructor
    /// Reads MULTI, DATACH and DATABR files structure prepared from GEMS.
    /// \param filepath - IPM work structure file path&name
    explicit GEMS3KGenerator(  const std::string& filepath ):
        ipmfiles_lst_name( filepath )
    {
        set_internal_data();
        load_dat_lst_file();
    }

    /// Constructor
    /// Generate MULTI, DATACH and DATABR files structure prepared from GEMS.
    /// Prints files for separate coupled FMT-GEM programs that use GEMS3K module
    /// \param filepath - IPM work structure file path&name
    /// \param nIV - Number of allocated nodes
    /// \param bin_mode - Write IPM, DCH and DBR files in binary, txt or json mode
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments -Write files with comments for all data entries ( text mode ) or as "pretty JSON"  ( json mode )
    /// \param addMui - Print internal indices in RMULTS to IPM file for reading into Gems back
    explicit GEMS3KGenerator(  const std::string& filepath, long int anIV, IOModes file_mode );

    /// Get same GEMS3K I/O set name (currently GEM-Selektor asks for the "set" name from .lst file name)
    std::string get_name() const
    {
        return base_name;
    }

    /// Get current extension
    std::string extension() const
    {
       return ext( io_mode );
    }

    /// Get selected file output mode
    IOModes files_mode() const
    {
        return io_mode;
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

    /// Generate full path for dch file
    std::string get_dch_path() const
    {
        return  ( impex_dir.empty() ? datach_file_name : impex_dir+"/"+datach_file_name );
    }

    /// Generate full path for dch file
    std::string get_thermofun_path() const
    {
        return  ( impex_dir.empty() ? thermofun_file_name : impex_dir+"/"+thermofun_file_name );
    }

    /// Generate full path for dbr file
    std::string get_dbr_path( size_t index ) const
    {
        if( index >= databr_file_names.size() )
           return "";
        return  ( impex_dir.empty() ? databr_file_names[index] : impex_dir+"/"+databr_file_names[index] );
    }

    /// Generate MULTI file name
    std::string gen_ipm_file_name()
    {
        ipm_file_name = gen_ipm_name( base_name ) + "." + extension();
        return ipm_file_name;
    }

    /// Generate dataCH file name
    std::string gen_dch_file_name()
    {
        datach_file_name = gen_dch_name( base_name ) + "." + extension();
        return datach_file_name;
    }

    /// Generate ThermoFun JSON format file name
    std::string gen_thermofun_file_name()
    {
        thermofun_file_name = gen_thermofun_name( base_name ) + ".json";
        return thermofun_file_name;
    }

    /// Generate dataBR file name
    std::string gen_dbr_file_name( int time_point, size_t index ) const;

    /// Generate *-dat.lst file data
    std::string gen_dat_lst_head();

    /// Generate dbr lst file path
    std::string get_dbr_file_lst_path() const;

    /// Read dbr lst file
    size_t load_dbr_lst_file( const std::string &dbr_lst_path );

    /// Generate dataBR file name
    std::string get_next_dbr_file( size_t index ) const
    {
       if( index < databr_file_names.size() )
         return ( impex_dir.empty() ? databr_file_names[index] : impex_dir+"/"+databr_file_names[index] );
       else
         return "";
    }

protected:

    /// IPM work structure file path&name
    std::string ipmfiles_lst_name;

    /// Number of allocated nodes
    size_t nIV = 1;

    /// Write IPM, DCH and DBR files in binary, txt or json mode)
    IOModes io_mode = f_json;

    std::string impex_dir;
    std::string base_name;

    std::string ipm_file_name;
    std::string datach_file_name;
    std::string thermofun_file_name;
    std::vector<std::string> databr_file_names;

    std::string mode() const;
    void get_mode( const std::string& str_mode );
    void set_internal_data();
    void load_dat_lst_file();

};

