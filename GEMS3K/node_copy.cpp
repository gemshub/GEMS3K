//--------------------------------------------------------------------
// $Id$
//
/// \file node_copy.cpp
/// Copy constructor for TNode class
/// DATACH and DATABR structures allocations
//
// Copyright (c) 2017-2021 S.Dmytriyeva, D.Kulik
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

#include  <sstream>
#include "node.h"
#include "num_methods.h"
#include "gdatastream.h"
#include "io_keyvalue.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "v_service.h"

// Writes CSD (DATACH structure) to a json string or key-value string
// \param brief_mode - Do not write data items that contain only default values
// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string TNode::datach_to_string( bool with_comments, bool brief_mode ) const
{
    std::stringstream ss;
    write_dch_format_stream( ss, GEMS3KGenerator::default_type_f, with_comments, brief_mode );
    return ss.str();
}

// Reads CSD (DATACH structure) from a json string or key-value string
bool TNode::datach_from_string( const std::string data )
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_dch_format_stream( ss, GEMS3KGenerator::default_type_f );
    return true;
}

// Writes work node (DATABR structure) to a json string or key-value string
// \param brief_mode - Do not write data items that contain only default values
// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string TNode::databr_to_string( bool with_comments, bool brief_mode ) const
{
    std::stringstream ss;
    write_dbr_format_stream( ss, GEMS3KGenerator::default_type_f, with_comments, brief_mode );
    return ss.str();
}

// Reads work node (DATABR structure) from a json string or key-value string
bool TNode::databr_from_string( const std::string data )
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_dbr_format_stream( ss, GEMS3KGenerator::default_type_f );
    return true;
}

void  TNode::read_dbr_format_file( const std::string& dbr_file, GEMS3KGenerator::IOModes  type_f )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
    {
        GemDataStream in_br(dbr_file, std::ios::in|std::ios::binary);
        databr_from_file(in_br);
    }
        break;
    default:
    {
        std::fstream in_br( dbr_file, std::ios::in );
        ErrorIf( !in_br.good() , dbr_file, "DBR_DAT fileopen error");
        read_dbr_format_stream( in_br, type_f );
    }
        break;
    }
}

void  TNode::write_dbr_format_file( const std::string& dbr_file, GEMS3KGenerator::IOModes type_f,
                                    bool with_comments, bool brief_mode )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
    {
        GemDataStream  f_br( dbr_file, std::ios::out|std::ios::binary );
        databr_to_file(f_br);
    }
        break;
    default:
    {
        std::fstream  f_br( dbr_file, std::ios::out);
        write_dbr_format_stream( f_br, type_f, with_comments, brief_mode );
    }
        break;
    }
}

void  TNode::read_dbr_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes  type_f )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonRead in_format( stream, current_input_set_name, "dbr" );
        databr_from_text_file(in_format);
    }
#else
    {
        io_formats::SimdJsonRead in_format( stream, current_input_set_name, "dbr" );
        databr_from_text_file(in_format);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueRead in_format( stream );
        databr_from_text_file(in_format);
    }
        break;
    }
}


void  TNode::write_dbr_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f,
                                    bool with_comments, bool brief_mode ) const
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonWrite out_format( stream, current_output_set_name );
        databr_to_text_file( out_format, with_comments, brief_mode );
    }
#else
    {
        io_formats::SimdJsonWrite out_format( stream, current_output_set_name, with_comments );
        databr_to_text_file( out_format, with_comments, brief_mode );
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueWrite out_format( stream );
        databr_to_text_file( out_format, with_comments, brief_mode );
    }
        break;
    }
}

void  TNode::read_dch_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes  type_f )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonRead in_format( stream, current_input_set_name, "dch" );
        datach_from_text_file(in_format, type_f==GEMS3KGenerator::f_thermofun);
    }
#else
    {
        io_formats::SimdJsonRead in_format( stream, current_input_set_name, "dch" );
        datach_from_text_file(in_format, type_f==GEMS3KGenerator::f_thermofun);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueRead in_format( stream );
        datach_from_text_file(in_format, type_f==GEMS3KGenerator::f_kv_thermofun);
    }
        break;
    }
}


void  TNode::write_dch_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f,
                                    bool with_comments, bool brief_mode ) const
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonWrite out_format( stream, current_output_set_name );
        datach_to_text_file( out_format, type_f==GEMS3KGenerator::f_thermofun, with_comments, brief_mode );
    }
#else
    {
        io_formats::SimdJsonWrite out_format( stream, current_output_set_name, with_comments );
        datach_to_text_file( out_format, type_f==GEMS3KGenerator::f_thermofun, with_comments, brief_mode );
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueWrite out_format( stream );
        datach_to_text_file( out_format, type_f==GEMS3KGenerator::f_kv_thermofun, with_comments, brief_mode );
    }
        break;
    }
}

//-------------------------------------------------------------------
// (1) Initialization of GEM IPM2 data structures in coupled RMT-GEM programs
//  that use GEMS3K module. Also reads in the IPM, DCH and DBR text input files 
//  in key-value, json or binary format. Parameters:
//  ipmfiles_lst_name - name of a text file that contains:
//    " -f | -j | -t |-b <DCH_DAT file name> <IPM_DAT file name> <dataBR file name>
//  dbfiles_lst_name - name of a text file that contains:
//    <dataBR  file name1>, ... , <dataBR file nameN> "
//    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
//    exist in the same directory where the ipmfiles_lst_name file is located.
//    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
//    and must contain valid initial chemical systems (of the same structure
//    as described in the DCH_DAT file) to set up the initial state of the FMT
//    node array. 
//  If -t flag or no flag is specified then all data files must be in key-value text 
//    (ASCII) format (and file names must have .dat extension);
//  If -j and -f flag is specified then all data files must be in JSON format (and file names
//    must have .json extension);
//  if -b flag is specified then all data files are assumed to be binary (little-endian)
//    files.
//-------------------------------------------------------------------
long int  TNode::GEM_init( const char* ipmfiles_lst_name )
{
    clearipmLogError();
    clear_ThermoEngine();

    try
    {
        //  Syntax: -t/-b  "<DCH_DAT file name>"  "<IPM_DAT file name>"
        //       "<DBR_DAT file1 name>" [ ...  "<DBR_DAT fileN name>"]
        GEMS3KGenerator generator( ipmfiles_lst_name );
        current_output_set_name = current_input_set_name = generator.get_name();

        switch( generator.files_mode() )
        {
        case GEMS3KGenerator::f_binary:
        {
            GemDataStream f_ch( generator.get_dch_path(), std::ios::in|std::ios::binary );
            datach_from_file(f_ch);

            GemDataStream f_m( generator.get_ipm_path(), std::ios::in|std::ios::binary );
             multi_ptr()->read_multi(f_m, CSD);
        }
            break;
        default:
        {
            std::fstream f_ch( generator.get_dch_path(), std::ios::in );
            ErrorIf( !f_ch.good() , generator.get_dch_path(), "DCH_DAT fileopen error");
            read_dch_format_stream( f_ch, generator.files_mode() );

            if( generator.files_mode()>=GEMS3KGenerator::f_thermofun )
            {
                load_ThermoEngine(generator.get_thermofun_path());
            }

            std::fstream ff( generator.get_ipm_path(), std::ios::in );
            ErrorIf( !ff.good() , generator.get_ipm_path(), "Fileopen error");
            multi_ptr()->read_ipm_format_stream( ff,generator.files_mode(), CSD, current_input_set_name);

        }
            break;
        }

        // copy intervals for minimization
        pmm->Pai[0] = CSD->Pval[0]/bar_to_Pa;
        pmm->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
        pmm->Pai[2] = getStep( pmm->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
        pmm->Pai[3] = CSD->Ptol/bar_to_Pa;

        pmm->Tai[0] = CSD->TKval[0]-C_to_K;
        pmm->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
        pmm->Tai[2] = getStep( pmm->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
        pmm->Tai[3] = CSD->Ttol;

        pmm->Fdev1[0] = 0.;
        pmm->Fdev1[1] = 1e-6;   // 24/05/2010 must be copied from GEMS3 structure
        pmm->Fdev2[0] = 0.;
        pmm->Fdev2[1] = 1e-6;

        // Reading DBR_DAT file into work DATABR structure from ipmfiles_lst_name
        std::string dbr_file = generator.get_dbr_path( 0 );
        ErrorIf( dbr_file.empty() , ipmfiles_lst_name, " Undefined DBR_DAT file name");
        read_dbr_format_file( dbr_file, generator.files_mode() );
        dbr_file_name = dbr_file;

        // Creating and initializing the TActivity class instance for this TNode instance
        init_into_gems3k();
        node_logger->info("Initialization of system {}", char_array_to_string(pmm->stkey,EQ_RKLEN));
        return 0;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }

    if( ipmfiles_lst_name ) {
        ipmlog_file->error("GEMS3K input : file {}", ipmfiles_lst_name);
    }
    if( !ipmlog_error.empty() ) {
        ipmlog_file->error("GEM_init error: {}",ipmlog_error);
    }
    return 1;
}

void TNode::init_into_gems3k()
{
    //InitReadActivities( mult_in.c_str(),CSD ); // from DCH file in future?
    multi_ptr()->InitalizeGEM_IPM_Data();              // In future, initialize data in TActivity also
    this->InitCopyActivities( CSD, pmm, CNode );
}

//  Parameters:
//  @param dch_json -  DATACH - the Data for CHemistry data structure as a json/key-value string
//  @param ipm_json -  Multi structure as a json/key-value string
//  @param dbr_json -  DATABR - the data bridge structure as a json/key-value string
//  @param fun_json -  ThermoFun data structure as a json string
long int  TNode::GEM_init( std::string dch_json, std::string ipm_json,
                           std::string dbr_json, std::string fun_json)
{
    load_thermodynamic_data = false; // need load thermo
    clearipmLogError();
    clear_ThermoEngine();

    try
    {
       if(fun_json.empty()) {
           GEMS3KGenerator::default_type_f  = GEMS3KGenerator::f_json;
       }
       else {
           GEMS3KGenerator::default_type_f  = GEMS3KGenerator::f_thermofun;
       }
        node_logger->debug("GEMS3KGenerator::default_type_f {}", GEMS3KGenerator::default_type_f);
        node_logger->debug("dch_json {}", dch_json);
        node_logger->debug("dbr_json {}", dbr_json);
        node_logger->debug("ipm_json {}", ipm_json);
        node_logger->debug("fun_json {}", fun_json);

        // This check of data consistency temporarily disabled for perfomance testing
        //if( GEMS3KGenerator::default_type_f  == GEMS3KGenerator::f_json )
        //{
             current_output_set_name = current_input_set_name = extract_string_json( "set", dch_json );
        //     auto ipm_set =  extract_string_json( "set", ipm_json );
        //     auto dbr_set =  extract_string_json( "set", dbr_json );
        //     ErrorIf(  ipm_set!=current_input_set_name,  "GEM_init error", "Multi structure as a json has different set name:  "+ipm_set );
        //     ErrorIf(  dbr_set!=current_input_set_name,  "GEM_init error", "The data bridge structure as a json has different set name:  "+dbr_set );
        //}

        // Reading DCH_DAT data
        datach_from_string(dch_json);

        if( !fun_json.empty() )
        {
            load_ThermoEngine(fun_json);
        }
        // Reading IPM_DAT file into structure MULTI (GEM IPM work structure)
        multi_ptr()->gemipm_from_string( ipm_json, CSD, current_input_set_name );

        // copy intervals for minimization
        pmm->Pai[0] = CSD->Pval[0]/bar_to_Pa;
        pmm->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
        pmm->Pai[2] = getStep( pmm->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
        pmm->Pai[3] = CSD->Ptol/bar_to_Pa;

        pmm->Tai[0] = CSD->TKval[0]-C_to_K;
        pmm->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
        pmm->Tai[2] = getStep( pmm->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
        pmm->Tai[3] = CSD->Ttol;

        pmm->Fdev1[0] = 0.;
        pmm->Fdev1[1] = 1e-6;   // 24/05/2010 must be copied from GEMS3 structure
        pmm->Fdev2[0] = 0.;
        pmm->Fdev2[1] = 1e-6;

        // Reading DBR_DAT file into work DATABR structure from ipmfiles_lst_name
        databr_from_string(dbr_json);

        // Creating and initializing the TActivity class instance for this TNode instance
        init_into_gems3k();
        node_logger->info("Initialization of system {}", std::string(pmm->stkey, 0, EQ_RKLEN));

        return 0;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }

    if( !ipmlog_error.empty() ) {
        ipmlog_file->error("GEM_init error: {}",ipmlog_error);
    }
    return 1;
}

long int  TNode::GEM_write_dbr( std::string& dbr_json )
{
    // Writes work node (DATABR structure) into a json string
    try
    {
        std::stringstream ss;
        write_dbr_format_stream( ss, GEMS3KGenerator::f_json, false, false );
        dbr_json =  ss.str();
        return 0;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }

    ipmlog_file->error("Error Node:{}  time:{}  dt:{}", CNode->NodeHandle, CNode->Tm, CNode->dt);
    ipmlog_file->error("{}", ipmlog_error);
    return 1;
}

void TNode::GEM_write_dbr( const char* fname, GEMS3KGenerator::IOModes type_f, bool with_comments, bool brief_mode )
{
    std::string str_file;
    if( fname == 0)
        str_file = dbr_file_name;//+".out";
    else
        str_file = fname;
    write_dbr_format_file( str_file, type_f, with_comments, brief_mode );
}

// (4) Produces a formatted text file with detailed contents (scalars and arrays) of the GEM IPM work structure.
// This call is useful when GEM_run() returns with a NodeStatusCH value indicating a GEM calculation error
// (see  above).  Another use is for a detailed comparison of a test system calculation after the version upgrade of GEMS3K.
// Parameters: fname   null-terminated (C) string containing a full path to the disk file to be written.
//                     NULL  - the disk file name path stored in the  dbr_file_name  field of the TNode class instance will be used,
//                     extended with ".dump.out".  Usually the dbr_file_name field contains the path to the last input DBR file.
//
void  TNode::GEM_print_ipm( const char* fname )
{
    std::string str_file;
    if( fname == 0)
        str_file = dbr_file_name + ".Dump.out";
    else
        str_file = fname;

     multi_ptr()->to_text_file(str_file.c_str());//profil1->outMultiTxt( str_file.c_str()  );
}


// (5j) Reads another DBR object (with input system composition, T,P etc.) from JSON string \ .
long int TNode::GEM_read_dbr( std::string dbr_json, const bool check_dch_compatibility )
{
    // Reads work node (DATABR structure) from a json string
    try
    {
        if( dbr_json.empty() )
            Error( "GEM_read_dbr" , "Empty input string");

        if( check_dch_compatibility )
        {
            auto dbr_set =  extract_string_json( "set", dbr_json );
            if(  dbr_set!=current_input_set_name )
                return 2;
        }
        std::stringstream ss;
        ss.str(dbr_json);
        read_dbr_format_stream( ss, GEMS3KGenerator::f_json );
        return 0;
     }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }

    ipmlog_file->error("Error Node:{}  time:{}  dt:{}", CNode->NodeHandle, CNode->Tm, CNode->dt);
    ipmlog_file->error("{}", ipmlog_error);
    return 1;

}


// (5) Reads another DBR file (with input system composition, T,P etc.). The DBR file must be compatible with
// the currently loaded IPM and DCH files (see description  of GEM_init() function call).
// Parameters:
//    fname       Null-terminated (C) string containing a full path to the input DBR disk file.
//    binary_f    Flag defining whether the file specified in fname is in text fromat (false or 0, default) or in binary format (true or 1)
// Return values:     0  if successful; 1 if input file(s) has not found been or is corrupt; -1 if internal memory allocation error occurred.
long int  TNode::GEM_read_dbr( const char* fname, GEMS3KGenerator::IOModes type_f )
{
    clearipmLogError();

    try
    {
        read_dbr_format_file( fname, type_f );
        dbr_file_name = fname;
        return 0;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }
    ipmlog_file->error("GEMS3K input : file {}", fname);
    ipmlog_file->error("Error Node:{}  time:{}  dt:{}", CNode->NodeHandle, CNode->Tm, CNode->dt);
    ipmlog_file->error("{}", ipmlog_error);
    return 1;
}

// Copy CSD (DATACH structure) data from other structure.
void TNode::datach_copy( DATACH* otherCSD )
{
    // const data
    copyValues( &CSD->nIC,  &otherCSD->nIC, 14 );
    copyValues( &CSD->Ttol, &otherCSD->Ttol,4 );

    datach_realloc();
    databr_free_internal(CNode);
    databr_realloc();

    //dynamic data
    copyValues( CSD->nDCinPH, otherCSD->nDCinPH, CSD->nPH );
    //   if( CSD->nICb >0 )
     copyValues( CSD->xic, otherCSD->xic, CSD->nICb );
    copyValues( CSD->xdc, otherCSD->xdc, CSD->nDCb );
    copyValues( CSD->xph, otherCSD->xph, CSD->nPHb );

    copyValues( CSD->A, otherCSD->A, CSD->nIC*CSD->nDC );
    copyValues( CSD->ICmm, otherCSD->ICmm, CSD->nIC );
    copyValues( CSD->DCmm, otherCSD->DCmm, CSD->nDC );

    copyValues( CSD->TKval,  otherCSD->TKval,  CSD->nTp );
    copyValues( CSD->Psat,  otherCSD->Psat,CSD->nTp );
    copyValues( CSD->Pval,  otherCSD->Pval,  CSD->nPp );

    copyValues( CSD->ccIC, otherCSD->ccIC, CSD->nIC );
    copyValues( CSD->ccDC, otherCSD->ccDC, CSD->nDC );
    copyValues( CSD->ccPH, otherCSD->ccPH, CSD->nPH );

    if( CSD->ccPH[0] == PH_AQUEL )
       {
         copyValues( CSD->denW,  otherCSD->denW,  5*gridTP() );
         copyValues( CSD->denWg,  otherCSD->denWg,  5*gridTP() );
         copyValues( CSD->epsW, otherCSD->epsW, 5*gridTP() );
         copyValues( CSD->epsWg, otherCSD->epsWg, 5*gridTP() );
       }
    copyValues( CSD->G0,  otherCSD->G0,  CSD->nDC*gridTP() );
    copyValues( CSD->V0,  otherCSD->V0,  CSD->nDC*gridTP() );
    copyValues( CSD->H0,  otherCSD->H0,  CSD->nDC*gridTP() );
    copyValues( CSD->S0, otherCSD->S0, CSD->nDC*gridTP() );
    copyValues( CSD->Cp0, otherCSD->Cp0, CSD->nDC*gridTP() );
    copyValues( CSD->A0, otherCSD->A0, CSD->nDC*gridTP() );
    copyValues( CSD->U0, otherCSD->U0, CSD->nDC*gridTP() );
    if(  CSD->iGrd  )
         copyValues( CSD->DD, otherCSD->DD, CSD->nDCs*gridTP() );

    copyValues( (char *)CSD->ICNL, (char *)otherCSD->ICNL, MaxICN*CSD->nIC );
    copyValues( (char *)CSD->DCNL, (char *)otherCSD->DCNL, MaxDCN*CSD->nDC );
    copyValues( (char *)CSD->PHNL, (char *)otherCSD->PHNL, MaxPHN*CSD->nPH );
}

// Copy node (work DATABR structure) data from other DBR.
void TNode::databr_copy( DATABR* otherCNode )
{
    // const data
    copyValues( &CNode->NodeHandle, &otherCNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
    if( CNode->NodeStatusFMT != No_nodearray )
         copyValues( &CNode->TK, &otherCNode->TK, 32 );
    else
         copyValues( &CNode->TK, &otherCNode->TK, 15 );
#else
    ErrorIf(CNode->NodeStatusFMT != No_nodearray, "databr_copy",
         "Error reading work dataBR structure from binary file (No_nodearray)");
     copyValues( &CNode->TK, &otherCNode->TK, 15 );
#endif
    //dynamic data
    copyValues( CNode->bIC, otherCNode->bIC, CSD->nICb );
    copyValues( CNode->rMB, otherCNode->rMB, CSD->nICb );
    copyValues( CNode->uIC, otherCNode->uIC, CSD->nICb );
    copyValues( CNode->bSP, otherCNode->bSP, CSD->nICb );

    copyValues( CNode->xDC, otherCNode->xDC, CSD->nDCb );
    copyValues( CNode->gam, otherCNode->gam, CSD->nDCb );
    copyValues( CNode->dul, otherCNode->dul, CSD->nDCb );
    copyValues( CNode->dll, otherCNode->dll, CSD->nDCb );

    if( CSD->nAalp >0 )
        copyValues( CNode->aPH, otherCNode->aPH, CSD->nPHb );
    copyValues( CNode->xPH, otherCNode->xPH, CSD->nPHb );
    copyValues( CNode->vPS, otherCNode->vPS, CSD->nPSb );
    copyValues( CNode->mPS, otherCNode->mPS, CSD->nPSb );
    copyValues( CNode->bPS, otherCNode->bPS, CSD->nPSb*CSD->nICb );
    copyValues( CNode->xPA, otherCNode->xPA, CSD->nPSb );
    copyValues( CNode->amru, otherCNode->amru, CSD->nPSb );
    copyValues( CNode->amrl, otherCNode->amrl, CSD->nPSb );
    copyValues( CNode->omPH, otherCNode->omPH, CSD->nPHb );
}


// allocating DataCH structure
void TNode::datach_realloc()
{
  if( CSD->mLook == 1 &&  (CSD->nPp != CSD->nTp) )
     Error( "No-interpolation mode",
           "Different number of points for temperature and pressure ");

 datach_free();
 CSD->nDCinPH = new long int[CSD->nPH];

 if( CSD->nICb >0 )
   CSD->xic = new long int[CSD->nICb];
 else  CSD->xic = 0;
 if( CSD->nDCb >0 )
   CSD->xdc = new long int[CSD->nDCb];
 else  CSD->xdc = 0;
 if( CSD->nPHb >0 )
   CSD->xph = new long int[CSD->nPHb];
 else  CSD->xph = 0;

  CSD->A = new double[CSD->nIC*CSD->nDC];
  CSD->ICmm = new double[CSD->nIC];
  CSD->DCmm = new double[CSD->nDC];
CSD->DCmm[0] = 0.0;   // Added by DK on 03.03.2007

  CSD->TKval = new double[CSD->nTp];
  CSD->Psat = new double[CSD->nTp];
  CSD->Pval = new double[CSD->nPp];

  CSD->denW = new double[ 5*gridTP()];
  CSD->denWg = new double[ 5*gridTP()];
  CSD->epsW = new double[ 5*gridTP()];
  CSD->epsWg = new double[ 5*gridTP()];

  CSD->G0 = new double[CSD->nDC*gridTP()];
  CSD->V0 = new double[CSD->nDC*gridTP()];
  CSD->H0 = new double[CSD->nDC*gridTP()];
  CSD->S0 = new double[CSD->nDC*gridTP()];
  CSD->Cp0 = new double[CSD->nDC*gridTP()];
  CSD->A0 = new double[CSD->nDC*gridTP()];
  CSD->U0 = new double[CSD->nDC*gridTP()];

  if(  CSD->iGrd  )
       CSD->DD = new double[CSD->nDCs*gridTP()];
  else
       CSD->DD = 0;
  CSD->ICNL = new char[CSD->nIC][MaxICN];
  CSD->DCNL = new char[CSD->nDC][MaxDCN];
  CSD->PHNL = new char[CSD->nPH][MaxPHN];

  CSD->ccIC = new char[CSD->nIC];
  CSD->ccDC = new char[CSD->nDC];
  CSD->ccPH = new char[CSD->nPH];
}

// free dynamic memory
void TNode::datach_free()
{
 if( CSD->nDCinPH )
  { delete[] CSD->nDCinPH;
    CSD->nDCinPH = 0;
  }
 if( CSD->xic )
  { delete[] CSD->xic;
    CSD->xic = 0;
  }
 if( CSD->xdc )
  { delete[] CSD->xdc;
    CSD->xdc = 0;
  }
 if( CSD->xph )
  { delete[] CSD->xph;
    CSD->xph = 0;
  }
 if( CSD->A )
  { delete[] CSD->A;
    CSD->A = 0;
  }
 if( CSD->ICmm )
  { delete[] CSD->ICmm;
    CSD->ICmm = 0;
  }
 if( CSD->DCmm )
  { delete[] CSD->DCmm;
    CSD->DCmm = 0;
  }

 if( CSD->TKval )
  { delete[] CSD->TKval;
    CSD->TKval = 0;
  }
 if( CSD->Psat )
  { delete[] CSD->Psat;
    CSD->Psat = 0;
  }
 if( CSD->Pval )
  { delete[] CSD->Pval;
    CSD->Pval = 0;
  }

 if( CSD->denW )
  { delete[] CSD->denW;
    CSD->denW = 0;
  }
 if( CSD->denWg )
  { delete[] CSD->denWg;
    CSD->denWg = 0;
  }
 if( CSD->epsW )
  { delete[] CSD->epsW;
    CSD->epsW = 0;
  }
 if( CSD->epsWg )
  { delete[] CSD->epsWg;
    CSD->epsWg = 0;
  }
 if( CSD->G0 )
  { delete[] CSD->G0;
    CSD->G0 = 0;
  }
 if( CSD->V0 )
  { delete[] CSD->V0;
    CSD->V0 = 0;
  }
 if( CSD->H0 )
  { delete[] CSD->H0;
    CSD->H0 = 0;
  }
 if( CSD->Cp0 )
  { delete[] CSD->Cp0;
    CSD->Cp0 = 0;
  }
  if( CSD->S0 )
  { delete[] CSD->S0;
     CSD->S0 = 0;
  }
  if( CSD->A0 )
  { delete[] CSD->A0;
     CSD->A0 = 0;
  }
  if( CSD->U0 )
  { delete[] CSD->U0;
     CSD->U0 = 0;
  }
  if( CSD->DD )
  { delete[] CSD->DD;
     CSD->DD = 0;
  }

 if( CSD->ICNL )
  { delete[] CSD->ICNL;
    CSD->ICNL = 0;
  }
 if( CSD->DCNL )
  { delete[] CSD->DCNL;
    CSD->DCNL = 0;
  }
 if( CSD->PHNL )
  { delete[] CSD->PHNL;
    CSD->PHNL = 0;
  }

 if( CSD->ccIC )
  { delete[] CSD->ccIC;
    CSD->ccIC = 0;
  }
 if( CSD->ccDC )
  { delete[] CSD->ccDC;
    CSD->ccDC = 0;
  }
 if( CSD->ccPH )
  { delete[] CSD->ccPH;
    CSD->ccPH = 0;
  }
 // delete[] CSD;
}

// Allocates DataBR structure
void TNode::databr_realloc( DATABR * CNode_)
{
  long int j,k;
  CNode_->bIC = new double[CSD->nICb];
  CNode_->rMB = new double[CSD->nICb];
  CNode_->uIC = new double[CSD->nICb];
  CNode_->bSP = new double[CSD->nICb];

  for(  j=0; j<CSD->nICb; j++ )
  {
      CNode_->rMB[j] = 0.;
      CNode_->uIC[j] = 0.;
      CNode_->bSP[j] = 0.;
   }

  CNode_->xDC = new double[CSD->nDCb];
  CNode_->gam = new double[CSD->nDCb];

  for(  j=0; j<CSD->nDCb; j++ )
  {
    CNode_->xDC[j] = 0.;
    CNode_->gam[j] = 1.;
  }

  //  default assignment
 CNode_->dul = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode_->dul[j] = 1.0e6;            // default assignment
 CNode_->dll = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode_->dll[j] = 0.0;              // default assignment

 if( CSD->nAalp >0 )
 {
    CNode_->aPH = new double[CSD->nPHb];
    for(  k=0; k<CSD->nPHb; k++ )
      CNode_->aPH[k] = 0.0;       // default assignment
 }
 else
    CNode_->aPH = 0;

 CNode_->xPH = new double[CSD->nPHb];
 CNode_->omPH = new double[CSD->nPHb];

 for(  k=0; k<CSD->nPHb; k++ )
 {
     CNode_->xPH[k] = 0.0;       // default assignment
     CNode_->omPH[k] = 0.0;
 }

 CNode_->vPS = new double[CSD->nPSb];
 CNode_->mPS = new double[CSD->nPSb];
 CNode_->bPS = new double[CSD->nPSb*CSD->nICb];
 CNode_->xPA = new double[CSD->nPSb];
 CNode_->amru = new double[CSD->nPSb];
 CNode_->amrl = new double[CSD->nPSb];

 for(  k=0; k<CSD->nPSb; k++ )
 {
     CNode_->vPS[k] = 0.0;
     CNode_->mPS[k] = 0.0;
     CNode_->xPA[k] = 0.0;
     for(  j=0; j<CSD->nICb; j++ )
        CNode_->bPS[k*CSD->nICb+j] = 0.0;
     CNode_->amru[k] = 1.0e6;
     CNode_->amrl[k] = 0.0;
 }
}

// free dynamic memory
void TNode::databr_free_internal( DATABR *CNode_ )
{
  if( CNode_ == 0)
    return;

 if( CNode_->bIC )
 { delete[] CNode_->bIC;
   CNode_->bIC = 0;
 }
 if( CNode_->rMB )
 { delete[] CNode_->rMB;
   CNode_->rMB = 0;
 }
 if( CNode_->uIC )
 { delete[] CNode_->uIC;
   CNode_->uIC = 0;
 }

 if( CNode_->xDC )
  { delete[] CNode_->xDC;
    CNode_->xDC = 0;
  }
 if( CNode_->gam )
  { delete[] CNode_->gam;
    CNode_->gam = 0;
  }
 if( CNode_->dul )
   { delete[] CNode_->dul;
     CNode_->dul = 0;
   }
 if( CNode_->dll )
   { delete[] CNode_->dll;
     CNode_->dll = 0;
   }

 if( CNode_->aPH )
 { delete[] CNode_->aPH;
   CNode_->aPH = 0;
 }
 if( CNode_->xPH )
  { delete[] CNode_->xPH;
    CNode_->xPH = 0;
  }
 if( CNode_->omPH )
  { delete[] CNode_->omPH;
    CNode_->omPH = 0;
  }
 if( CNode_->vPS )
  { delete[] CNode_->vPS;
    CNode_->vPS = 0;
  }
 if( CNode_->mPS )
  { delete[] CNode_->mPS;
    CNode_->mPS = 0;
  }
 if( CNode_->bPS )
  { delete[] CNode_->bPS;
    CNode_->bPS = 0;
  }
 if( CNode_->xPA )
  { delete[] CNode_->xPA;
    CNode_->xPA = 0;
  }

 if( CNode_->bSP )
 { delete[] CNode_->bSP;
   CNode_->bSP = 0;
 }
 if( CNode_->amru )
  { delete[] CNode_->amru;
    CNode_->amru = 0;
  }
 if( CNode_->amrl )
  { delete[] CNode_->amrl;
    CNode_->amrl = 0;
  }
}

DATABR * TNode::databr_free( DATABR *CNode_ )
{
    if( CNode_ == 0)
        CNode_ = CNode;
    databr_free_internal(CNode_);
    delete CNode_;
    return NULL;
}

// set default values(zeros) for DATABR structure
void TNode::databr_reset( DATABR *CnNde1, long int level )
{
    //  FMT variables (units or dimensionsless) - to be used for storing them
    //  at the nodearray level = 0.; normally not used in the single-node FMT-GEM coupling
        CnNde1->Tm = 0.;
        CnNde1->dt = 0.;
#ifdef NODEARRAYLEVEL
        CnNde1->Dif = 0.;
        CnNde1->Vt = 0.;
        CnNde1->vp = 0.;
        CnNde1->eps = 0.;
        CnNde1->Km = 0.;
        CnNde1->Kf = 0.;
        CnNde1->S = 0.;
        CnNde1->Tr = 0.;
        CnNde1->h = 0.;
        CnNde1->rho = 0.;
        CnNde1->al = 0.;
        CnNde1->at = 0.;
        CnNde1->av = 0.;
        CnNde1->hDl = 0.;
        CnNde1->hDt = 0.;
        CnNde1->hDv = 0.;
        CnNde1->nto = 0.; //19
#endif
        if(level <1 )
          return;

   CnNde1->NodeHandle = 0;
   CnNde1->NodeTypeHY = normal;
   CnNde1->NodeTypeMT = normal;
   CnNde1->NodeStatusFMT = Initial_RUN;
   CnNde1->NodeStatusCH = NEED_GEM_AIA;
   CnNde1->IterDone = 0;      //6

// Chemical scalar variables
    CnNde1->TK = 0.;
    CnNde1->P = 0.;
    CnNde1->Vs = 0.;
    CnNde1->Vi = 0.;
    CnNde1->Ms = 0.;
    CnNde1->Mi = 0.;
    CnNde1->Gs = 0.;
    CnNde1->Hs = 0.;
    CnNde1->Hi = 0.;
    CnNde1->IC = 0.;
    CnNde1->pH = 0.;
    CnNde1->pe = 0.;
    CnNde1->Eh = 0.; //13

    if( level < 2 )
       return;

// Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
    CnNde1->bIC = 0;
    CnNde1->rMB = 0;
    CnNde1->uIC = 0;
    CnNde1->xDC = 0;
    CnNde1->gam = 0;
   CnNde1->dul = 0;
   CnNde1->dll = 0;
   CnNde1->aPH = 0;
   CnNde1->xPH = 0;
   CnNde1->vPS = 0;
   CnNde1->mPS = 0;
   CnNde1->bPS = 0;
   CnNde1->xPA = 0;
   CnNde1->bSP = 0;
   CnNde1->amru = 0;
   CnNde1->amrl = 0;
   CnNde1->omPH = 0;
}

// set default values(zeros) for DATACH structure
void TNode::datach_reset()
{
    CSD->nIC = 0;
    CSD->nDC = 0;
    CSD->nPH = 0;
    CSD->nPS = 0;
    CSD->nDCs = 0;
    CSD->nTp = 0;
    CSD->nPp = 0;
    CSD->iGrd = 0;
    CSD->nAalp = 0;
    CSD->nICb = 0;
    CSD->nDCb = 0;
    CSD->nPHb = 0;
    CSD->nPSb = 0;
    CSD->mLook = 0;
// Lists = 0; vectors and matrices
    CSD->nDCinPH = 0;
    CSD->xic = 0;
    CSD->xdc = 0;
    CSD->xph = 0;  //18

    CSD->TKval = 0;
    CSD->Psat = 0;
    CSD->Pval = 0;
    CSD->A = 0;
    CSD->Ttol = 0.;
    CSD->Ptol = 0.;
    CSD->dRes1 = 0.;
    CSD->dRes2 = 0.;
    CSD->ICmm = 0;
    CSD->DCmm = 0;
    CSD->DD = 0;
    CSD->denW = 0;
    CSD->epsW = 0;
    CSD->denWg = 0;
    CSD->epsWg = 0;
    CSD->G0 = 0;
    CSD->V0 = 0;
    CSD->S0 = 0;
    CSD->H0 = 0;
    CSD->Cp0 = 0;
    CSD->A0 = 0;
    CSD->U0 = 0;
    CSD->ICNL = 0;
    CSD->DCNL = 0;
    CSD->PHNL = 0;
    CSD->ccIC = 0;
    CSD->ccDC = 0;
    CSD->ccPH = 0;
}

// Test if load thermodynamic data
void TNode::CheckMtparam()
{
    double TK, P, PPa;

    TK = cTK();
    PPa = cP();
    P = PPa/bar_to_Pa;
    //pmp->pTPD = 2;
    node_logger->debug("CheckMtparam T: {} - {}  P: {} - {}", pmm->Tc, TK, pmm->Pc, P);
    if( !load_thermodynamic_data || fabs( pmm->Tc - TK ) > CSD->Ttol
            || fabs( pmm->Pc - P )  > CSD->Ptol/bar_to_Pa  )
    {
        pmm->pTPD = 0;      //T, P is changed
    }
    load_thermodynamic_data = true;
}


void TNode::clear_ThermoEngine()
{
#ifdef USE_THERMOFUN
    // clear previous
    thermo_engine.reset();
    thermo_json_string="";
#endif
}

bool TNode::load_ThermoEngine(const std::string &thermo_file_or_string)
{
    if(thermo_file_or_string.find("\"elements\"")!=std::string::npos) {
        // input string
        thermo_json_string = thermo_file_or_string;
    }
    else {
        // input file
        std::ifstream f_fun(thermo_file_or_string);
        ErrorIf(!f_fun.good(), thermo_file_or_string, "ThermoFun JSON format file fileopen error");
        std::stringstream buffer;
        buffer << f_fun.rdbuf();
        thermo_json_string = buffer.str();
    }

#ifdef USE_THERMOFUN
    thermo_engine.reset(new ThermoFun::ThermoEngine(thermo_file_or_string));
    node_logger->trace("Read ThermoEngine: {}", thermo_file_or_string);
    return true;
#else
    node_logger->warn("Try read ThermoEngine not in USE_THERMOFUN mode {}", thermo_file_or_string);
    return false;
#endif
}


bool TNode::load_all_thermodynamic_from_thermo( double TK, double PPa )
{
#ifdef USE_THERMOFUN
    if( !thermo_engine.get() )
        return false;
    try{
        node_logger->info("Calc ThermoEngine T: {}  P: {}", TK, PPa);
        long int j, jj, k, jb, je=0;
        double G0, P = PPa/bar_to_Pa;

        pmm->T = pmm->Tc = TK;
        pmm->TC = pmm->TCc = TK-C_to_K;
        // new API
        double funT = TK, funP=P*bar_to_Pa;   // T in K, P in Pa

        DATACH  *dCH = pCSD();
        auto water_props = thermo_engine->propertiesSolvent(funT,funP, "H2O@");
        auto water_electro = thermo_engine->electroPropertiesSolvent(funT,funP, "H2O@");

        auto water_vapor = thermo_engine->database().getSubstance("H2O@");
        water_vapor.setMethod_P( ThermoFun::MethodCorrP_Thrift::type::CPM_GAS);

        ThermoFun::Database db2;
        db2.addSubstance(water_vapor);

        ThermoFun::ThermoEngine te2(db2);
        auto water_gas_props = te2.propertiesSolvent(funT,funP, "H2O@");
        auto water_gas_electro = te2.electroPropertiesSolvent(funT,funP, "H2O@");

        pmm->Pc = P;
        if( P < 1e-5 )
        { // Pressure at saturated H2O vapour at given temperature
            P = funP/bar_to_Pa;
            //        long int xT =check_grid_T(TK);
            //        if(xT>= 0)
            //            P = dCH->Psat[xT]/bar_to_Pa;
            //        else
            //            P =  LagranInterp( &PPa, dCH->TKval, dCH->Psat, PPa, TK, dCH->nTp, 1,6 )/bar_to_Pa;
        }

        pmm->P = P;
        pmm->RT = R_CONSTANT * pmm->Tc;
        pmm->FRT = F_CONSTANT/pmm->RT;
        pmm->lnP = log( P );

        pmm->denW[0] = water_props.density.val/1e3;
        pmm->epsW[0] = water_electro.epsilon.val;
        pmm->denW[1] = water_props.densityT.val/1e3;
        pmm->epsW[1] = water_electro.epsilonT.val;
        pmm->denW[2] = water_props.densityTT.val/1e3;
        pmm->epsW[2] = water_electro.epsilonTT.val;
        pmm->denW[3] = water_props.densityP.val*1e2; // /1e3;
        pmm->epsW[3] = water_electro.epsilonP.val;
        pmm->denW[4] = water_props.densityPP.val/1e3;
        pmm->epsW[4] = water_electro.epsilonPP.val;

        pmm->denWg[0] = water_gas_props.density.val/1e3;
        pmm->epsWg[0] = water_gas_electro.epsilon.val;
        pmm->denWg[1] = water_gas_props.densityT.val/1e3;
        pmm->epsWg[1] = water_gas_electro.epsilonT.val;
        pmm->denWg[2] = water_gas_props.densityTT.val/1e3;
        pmm->epsWg[2] = water_gas_electro.epsilonTT.val;
        pmm->denWg[3] = water_gas_props.densityP.val*1e2; // /1e3;
        pmm->epsWg[3] = water_gas_electro.epsilonP.val;
        pmm->denWg[4] = water_gas_props.densityPP.val/1e3;
        pmm->epsWg[4] = water_gas_electro.epsilonPP.val;

#ifdef  USE_THERMO_LOG
        std::fstream f_log("thermodynamic-log.csv", std::ios::out/*|std::ios::app*/ );
        f_log << "\nCalc ThermoEngine;T;" << TK << ";P;" << PPa << "\n";
        f_log << "denW";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pmm->denW[jj]);
        f_log << "\nepsW";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pmm->epsW[jj]);
        f_log << "\ndenWg";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pmm->denWg[jj]);
        f_log << "\nepsWg";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pmm->epsWg[jj]);
#endif

        for( k=0; k<pmm->FI; k++ )
        {
            jb = je;
            je += pmm->L1[k];
            // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
            // depending on the presence of these arrays in DATACH and Multi structures
            for( j=jb; j<je; j++ )
            {
                std::string symbol = std::string(CSD->DCNL[j], 0, MaxDCN);
                auto propAl    = thermo_engine->thermoPropertiesSubstance(funT,funP, symbol);

                G0 = propAl.gibbs_energy.val;
                pmm->Vol[j] = propAl.volume.val*10;
                if( dCH->S0 ) pmm->S0[j] = propAl.entropy.val;
                if( dCH->H0 ) pmm->H0[j] = propAl.enthalpy.val;
                if( dCH->Cp0 ) pmm->Cp0[j] = propAl.heat_capacity_cp.val;
                if( dCH->A0 ) pmm->A0[j] = propAl.helmholtz_energy.val;
                if( dCH->U0 ) pmm->U0[j] = propAl.internal_energy.val;

                pmm->G0[j] = multi_ptr()->ConvertGj_toUniformStandardState(G0, j, k);
#ifdef  USE_THERMO_LOG
                f_log << "\n" << symbol << ";" << floating_point_to_string(G0)
                      << ";" << floating_point_to_string(pmm->G0[j])
                      << ";" << floating_point_to_string(pmm->Vol[j]);
                if( dCH->S0 ) f_log << ";" << floating_point_to_string(pmm->S0[j]);
                if( dCH->H0 ) f_log << ";" << floating_point_to_string(pmm->H0[j]);
                if( dCH->Cp0 ) f_log << ";" << floating_point_to_string(pmm->Cp0[j]);
                if( dCH->A0 ) f_log << ";" << floating_point_to_string(pmm->A0[j]);
                if( dCH->U0 ) f_log << ";" << floating_point_to_string(pmm->U0[j]);
#endif
            }  // j
        } // k
    }
    catch (const std::runtime_error& exception)
    {
        node_logger->error("ThermoEngine error: {}", exception.what());
        Error( "ThermoEngine error:", exception.what());
    }
    return true;
#else
    return false;
#endif

}

