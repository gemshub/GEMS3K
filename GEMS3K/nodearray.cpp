//--------------------------------------------------------------------
// $Id$
//
/// \file nodearray.cpp
/// Implementation of TNodeArray class functionality - advanced
/// interface between GEM IPM and FMT node array
/// working with one DATACH structure and arrays of DATABR structures
//
// Copyright (c) 2004-2012 S.Dmytriyeva, D.Kulik
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
//

#ifndef NO_NODEARRAYLEVEL

#include "nodearray.h"
#include "gdatastream.h"
#include "gems3k_impex.h"
#include "v_service.h"

TNodeArray* TNodeArray::na = nullptr;


TNodeArray::TNodeArray(long int nNod):
    internal_Node(new TNode()), calcNode(internal_Node.get()),
    anNodes(nNod)
{
    sizeN = anNodes;
    sizeM = sizeK =1;
    NodT0 = 0;  // nodes at current time point
    NodT1 = 0;  // nodes at previous time point
    grid  = 0;   // Array of grid point locations, size is anNodes+1
    tcNode = 0;     // Node type codes (see DataBR.h) size anNodes+1
    iaNode = 0;
    na = this;
    allocMemory();
}

TNodeArray::TNodeArray( long int asizeN, long int asizeM, long int asizeK ):
    internal_Node(new TNode()), calcNode(internal_Node.get()),
    sizeN(asizeN), sizeM(asizeM), sizeK(asizeK)
{
    anNodes = asizeN*asizeM*asizeK;
    NodT0 = 0;  // nodes at current time point
    NodT1 = 0;  // nodes at previous time point
    grid  = 0;   // Array of grid point locations, size is anNodes+1
    tcNode = 0;     // Node type codes (see DataBR.h) size anNodes+1
    iaNode = 0;
    na = this;
    allocMemory();
}


TNodeArray::~TNodeArray()
{
    //na = nullptr;
    freeMemory();
}

void TNodeArray::allocMemory()
{
    long int ii;

    // The NodeArray must be allocated here
    /// calcNode = new TNode();

    // alloc memory for data bidge structures
    // did in constructor TNode::allocMemory();

    // alloc memory for all nodes at current time point
    NodT0 = new  DATABRPTR[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        NodT0[ii] = nullptr;

    // alloc memory for all nodes at previous time point
    NodT1 = new  DATABRPTR[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        NodT1[ii] = nullptr;

    // alloc memory for the work array of node types
    tcNode = new char[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        tcNode[ii] = normal;

    // alloc memory for the work array of IA indicators
    iaNode = new bool[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        iaNode[ii] = true;
    // grid ?
}

void TNodeArray::freeMemory()
{
    long int ii;

    if( anNodes )
    {
        if( NodT0 )
            for(  ii=0; ii<anNodes; ii++ )
                if( NodT0[ii] )
                    NodT0[ii] = calcNode->databr_free(NodT0[ii]);
        delete[]  NodT0;
        NodT0 = nullptr;

        if( NodT1 )
            for(  ii=0; ii<anNodes; ii++ )
                if( NodT1[ii] )
                    NodT1[ii] = calcNode->databr_free(NodT1[ii]);
        delete[]  NodT1;
        NodT1 = nullptr;
    }

    if( grid )
        delete[] grid;
    if( tcNode)
        delete[] tcNode;
    if( iaNode )
        delete[] iaNode;

    ///if(calcNode)
    ///   delete calcNode;
}

// To parallelization calculations ========================================================


#ifdef useOMP

#include <omp.h>

//   Here we call a loop on GEM calculations over nodes
//   parallelization should affect this loop only
//   return code   true   Ok
//                 false  Error in GEMipm calculation part
//
bool TNodeArray::CalcIPM_List( const TestModeGEMParam& modeParam, long int start_node, long int end_node, FILE* diffile )
{
    int n;
    long int ii;
    bool iRet = true;
    TNode wrkNode( *calcNode ); // must be copy TNode internal
    DATABRPTR* C0 = pNodT0();
    DATABRPTR* C1 = pNodT1();
    bool* iaN = piaNode();     // indicators for IA in the nodes

    start_node = std::max( start_node, 0L );
    end_node = std::min( end_node, anNodes-1 );


#pragma omp parallel shared( C0, C1, iaN, diffile, iRet ) private( wrkNode, n )
    {
        TNode wrkNode( *calcNode ); // must be copy TNode internal
        n = omp_get_thread_num();
#pragma omp for
        for( ii = start_node; ii<= end_node; ii++) // node iteration
        {
            if( !CalcIPM_Node(  modeParam, &wrkNode, ii, C0, C1, iaN, diffile ) )
            {
#pragma omp atomic write
                iRet = false;
            }

//#pragma omp critical
//            {
//              TNode::node_logger->error(" {} -thread did index: {} ", n, ii);
//            }
        }

    }

    return iRet;
}

//   Here we do a GEM calculation in box ii
//   return code   true   Ok
//                 false  Error in GEMipm calculation part
//
bool TNodeArray::CalcIPM_Node( const TestModeGEMParam& modeParam, TNode* wrkNode,
                               long int ii, DATABRPTR* C0, DATABRPTR* C1, bool* piaN, FILE* diffile )
{
    bool iRet = true;

    long int Mode = SmartMode( modeParam, ii, piaN   );
    bool needGEM = NeedGEMS( wrkNode, modeParam, C0[ii], C1[ii]  );

    if( needGEM )
    {
        long RetCode =  RunGEM( wrkNode, ii, Mode, C1 );
        // checking RetCode from GEM IPM calculation
        if( !(RetCode==OK_GEM_AIA || RetCode == OK_GEM_SIA ))
        {
            std::string err_msg = ErrorGEMsMessage( RetCode,  ii, modeParam.step  );
            iRet = false;

            if( diffile )
            {
#pragma omp critical
                {
                    // write to file here
                    fprintf( diffile, "\nError reported from GEMS3K module\n%s\n",
                             err_msg.c_str() );
                }
            }
        }
    }
    else { // GEM calculation for this node not needed
        C1[ii]->IterDone = 0; // number of GEMIPM iterations is set to 0
    }
    return iRet;
}

#else


//   Here we call a loop on GEM calculations over nodes
//   parallelization should affect this loop only
//   return code   true   Ok
//                 false  Error in GEMipm calculation part
//
bool TNodeArray::CalcIPM_List( const TestModeGEMParam& modeParam, long int start_node, long int end_node, FILE* diffile )
{
    long int ii;
    bool iRet = true;
    DATABRPTR* C0 = pNodT0();
    DATABRPTR* C1 = pNodT1();
    bool* iaN = piaNode();     // indicators for IA in the nodes

    start_node = std::max( start_node, 0L );
    end_node = std::min( end_node, anNodes-1 );


    for( ii = start_node; ii<= end_node; ii++) // node iteration
    {
        if( !CalcIPM_Node(  modeParam, calcNode, ii, C0, C1, iaN, diffile ) )
            iRet = false;
    }

    return iRet;
}

//   Here we do a GEM calculation in box ii
//   return code   true   Ok
//                 false  Error in GEMipm calculation part
//
bool TNodeArray::CalcIPM_Node( const TestModeGEMParam& modeParam, TNode* wrkNode,
                               long int ii, DATABRPTR* C0, DATABRPTR* C1, bool* piaN, FILE* diffile )
{
    bool iRet = true;

    long int Mode = SmartMode( modeParam, ii, piaN   );
    bool needGEM = NeedGEMS( wrkNode, modeParam, C0[ii], C1[ii]  );

    if( needGEM )
    {
        long RetCode =  RunGEM( wrkNode, ii, Mode, C1 );

        // checking RetCode from GEM IPM calculation
        if( !(RetCode==OK_GEM_AIA || RetCode == OK_GEM_SIA ))
        {
            std::string err_msg = ErrorGEMsMessage( RetCode,  ii, modeParam.step  );
            iRet = false;

            if( diffile )
            {
                // write to file here
                fprintf( diffile, "\nError reported from GEMS3K module\n%s\n",
                         err_msg.c_str() );
            }
        }
    }
    else { // GEM calculation for this node not needed
        C1[ii]->IterDone = 0; // number of GEMIPM iterations is set to 0
    }
    return iRet;
}

#endif

void TNodeArray::RunGEM( long int Mode, int nNodes, DATABRPTR* nodeArray, long int* nodeFlags, long int* retCodes )
{
#ifdef useOMP
#pragma omp parallel
#endif
    {
        TNode workNode(*na->getCalcNode());
#ifdef useOMP
#pragma omp for
#endif
        for (long int node=0;node<nNodes;node++)
        {
            if (nodeFlags[node])
                retCodes[node] = na->RunGEM(&workNode, node, Mode, nodeArray);
        }
    }
}

// New init ================================================================

//-------------------------------------------------------------------
// (1) Initialization of GEM IPM2 data structures in coupled FMT-GEM programs
//  that use GEMS3K module. Also reads in the IPM, DCH and DBR text input files.
//  Parameters:
//  ipmfiles_lst_name - name of a text file that contains:
//    " -t/-b <DCH_DAT file name> <IPM_DAT file name> <dataBR file name>
//  dbfiles_lst_name - name of a text file that contains:
//    <dataBR  file name1>, ... , <dataBR file nameN> "
//    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
//    exist in the same directory where the ipmfiles_lst_name file is located.
//    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
//    and must contain valid initial chemical systems (of the same structure
//    as described in the DCH_DAT file) to set up the initial state of the FMT
//    node array. If -t flag or nothing is specified then all data files must
//    be in text (ASCII) format; if -b flag is specified then all data files
//    are  assumed to be binary (little-endian) files.
//  nodeTypes[nNodes] - optional parameter used only on the TNodeArray level,
//    array of node type (fortran) indexes of DBR_DAT files
//    in the ipmfiles_lst_name list. This array (handle for each FMT node),
//    specifies from which DBR_DAT file the initial chemical system should
//    be taken.
//  getNodT1 - optional flag, used only (if true) when reading multiple DBR files
//    after the coupled modeling task interruption in GEM-Selektor
//  This function returns:
//   0: OK; 1: GEM IPM read file error; -1: System error (e.g. memory allocation)
//
//-------------------------------------------------------------------
long int  TNodeArray::GEM_init( const char* ipmfiles_lst_name,
                                const char* dbrfiles_lst_name, long int* nodeTypes, bool getNodT1)
{
    auto ret = calcNode->GEM_init( ipmfiles_lst_name );
    if( ret )
        return ret;

    try
    {
        // Reading DBR_DAT files from dbrfiles_lst_name
        if(dbrfiles_lst_name)
        {
            //  Syntax: -t/-b  "<DCH_DAT file name>"  "<IPM_DAT file name>"
            //       "<DBR_DAT file1 name>" [ ...  "<DBR_DAT fileN name>"]
            GEMS3KGenerator generator( ipmfiles_lst_name );
            InitNodeArray( dbrfiles_lst_name, nodeTypes, getNodT1, generator.files_mode()  );
        }
        else {
            if( nNodes() ==1 )
                setNodeArray( 0 , nullptr  );
            else // undefined TNodeArray
                Error( "GEM_init", "Undefined boundary condition!" );
        }
        return 0;

    }
    catch(TError& err)
    {
        calcNode->ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        calcNode->ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        calcNode->ipmlog_error = "unknown exception";
        return -1;
    }

    if( ipmfiles_lst_name ) {
        TNode::ipmlog_file->error("GEMS3K input : file {}", ipmfiles_lst_name);
    }
    if( !calcNode->ipmLogError().empty() ) {
        TNode::ipmlog_file->error("GEM_init error: {}", calcNode->ipmLogError());
    }
    return 1;
}

long TNodeArray::GEM_init(std::string dch_json, std::string ipm_json, std::vector<std::string> dbr_json, long *nodeTypes)
{

    if(dbr_json.empty())
    {
        calcNode->ipmlog_error = "empty input dbr data";
        return 1;
    }
    auto ret = calcNode->GEM_init(dch_json, ipm_json, dbr_json[0]);
    if( ret )
        return ret;

    try
    {
        for( size_t ii=0; ii<dbr_json.size(); ++ii )
        {
            calcNode->GEM_read_dbr(dbr_json[ii]);
            setNodeArray( ii, nodeTypes  );
        }
        return 0;
    }
    catch(TError& err)
    {
        calcNode->ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        calcNode->ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        calcNode->ipmlog_error = "unknown exception";
        return -1;
    }

    if( !calcNode->ipmLogError().empty() ) {
        TNode::ipmlog_file->error("GEM_init error: {}", calcNode->ipmLogError());
    }
    return 1;
}


// ------------------------------------------------------------------

bool TNodeArray::NeedGEMS( TNode* wrkNode, const TestModeGEMParam& modeParam, DATABR* C0, DATABR* C1  )
{
    bool NeedGEM = false;
    DATACH* CH = wrkNode->pCSD();  // DataCH structure
    double dc;

    if( modeParam.useSIA == S_OFF )
        NeedGEM = true;
    else
    {   C1->IterDone = 0;
        NeedGEM = false;
    }

    // Here we compare this node for current time and for previous time - works for AIA and PIA
    for( long int ic=0; ic < CH->nICb; ic++)    // do we check charge here?
    {
        // It has to be checked on minimal allowed c0 value
        if( C1->bIC[ic] < modeParam.cez )
            C1->bIC[ic] = modeParam.cez; // to prevent loss of Independent Component

        dc = C0->bIC[ic] - C1->bIC[ic];
        if( fabs( dc ) > std::min( modeParam.cdv, (C1->bIC[ic] * 1e-3 ) ))
        {
            NeedGEM = true;  // we still need to recalculate equilibrium
            // in this node because its vector b has changed
        }
    }
    if( CH->ccIC[CH->nICb-1] == IC_CHARGE ) {
        C1->bIC[CH->nICb-1] = 0.;   // zeroing charge off in bulk composition
    }
    return NeedGEM;
}

long int TNodeArray::SmartMode( const TestModeGEMParam& modeParam, long int ii, bool* iaN  )
{
    long int Mode = NEED_GEM_AIA;
    //bool* iaN = piaNode();     // indicators for IA in the nodes

    if( modeParam.useSIA == S_OFF )
        iaN[ii] = true;
    else
        iaN[ii] = false;

    if( modeParam.mode == NEED_GEM_SIA )
    {
        // smart algorithm
        if( iaN[ii] == true )
        {
            Mode = NEED_GEM_AIA;
        }
        else {
            Mode = NEED_GEM_SIA;
            if( modeParam.useSIA == S_ON )   // force loading of primal solution into GEMIPM
                Mode *= -1;            // othervise use internal (old) primal solution
        }
    }

    return Mode;
}

std::string TNodeArray::ErrorGEMsMessage( long int RetCode,  long int ii, long int step  )
{
    std::string err_msg;
    err_msg = " Node= "+std::to_string(ii)+"  Step= "+std::to_string(step)+"\n";

    switch( RetCode )
    {
    case BAD_GEM_AIA:
        err_msg += "Bad GEM result using LPP AIA";
        break;
    case  ERR_GEM_AIA:
        err_msg += "GEM calculation error using LPP AIA";
        break;
    case  BAD_GEM_SIA:
        err_msg += "Bad GEM result using SIA";
        break;
    case  ERR_GEM_SIA:
        err_msg += "GEM calculation error using SIA";
        break;
    case  T_ERROR_GEM:  err_msg +=  "Terminal error in GEMS3K module";
    }

    return  err_msg;
}

//-------------------------------------------------------------------------
// RunGEM()
// GEM IPM calculation of equilibrium state for the iNode node
// from array NodT1. abs(Mode)) - mode of GEMS calculation (NEED_GEM_SIA or NEED_GEM_AIA)
//    if Mode is negative then the loading of primal solution from the node is forced
//    (only in SIA mode)
//  Function returns: NodeStatus code after GEM calculation
//   ( OK_GEM_AIA; OK_GEM_SIA; error codes )
//
//-------------------------------------------------------------------

long int  TNodeArray::RunGEM( TNode* wrkNode, long int  iNode, long int Mode, DATABRPTR* nodeArray )
{
    long int retCode = T_ERROR_GEM;

    // Copy data from the iNode node from array NodT1 to the work DATABR structure
    CopyWorkNodeFromArray( wrkNode, iNode, anNodes, nodeArray );

    // GEM IPM calculation of equilibrium state in MULTI
    wrkNode->pCNode()->NodeStatusCH = std::abs(Mode);
    retCode = CalcNodeServer( wrkNode, iNode, Mode );

    // Copying data for node iNode back from work DATABR structure into the node array
    //   if( retCode == OK_GEM_AIA ||
    //       retCode == OK_GEM_PIA  )
    MoveWorkNodeToArray( wrkNode, iNode, anNodes, nodeArray );

    return retCode;
}

long TNodeArray::CalcNodeServer(TNode* wrkNode, long ,  long int Mode)
{
    bool uPrimalSol = false;
    if( Mode < 0 || std::abs(Mode) == NEED_GEM_SIA )
        uPrimalSol = true;
    return  wrkNode->GEM_run( uPrimalSol );
}

//-------------------------------------------------------------------
// Initialization of TNodeArray data structures. Reads in the DBR text input files and
// copying data from work DATABR structure into the node array
// (as specified in nodeTypes array, ndx index of dataBR files in
//    the dbrfiles_lst_name list).
//   type_f    defines if the file is in binary format (1), in text format (0) or in json format (2)
//-------------------------------------------------------------------
void  TNodeArray::InitNodeArray( const char *dbrfiles_lst_name,
                                 long int *nodeTypes, bool getNodT1, GEMS3KGenerator::IOModes type_f  )
{
    int i;
    std::string datachbr_fn;
    std::string lst_in = dbrfiles_lst_name;
    std::string dir_path = u_getpath( lst_in );
    if( !dir_path.empty() )
        dir_path += "/";

    //  open file stream for the file names list file
    std::fstream f_lst( lst_in, std::ios::in );
    ErrorIf( !f_lst.good() , lst_in, "Fileopen error");

    // Prepare for reading DBR_DAT files
    i = 0;
    while( !f_lst.eof() )  // For all DBR_DAT files listed
    {
        pVisor_Message( false, i, nNodes() );

        // Reading DBR_DAT file into work DATABR structure
        getline( f_lst, datachbr_fn, ',');
        trim( datachbr_fn );
        trim( datachbr_fn, "\"" );

        std::string dbr_file = dir_path + datachbr_fn;
        calcNode->read_dbr_format_file( dbr_file,  type_f );

        // Unpacking work DATABR structure into MULTI (GEM IPM work structure): uses DATACH
        //    unpackDataBr();

        if( getNodT1 )  // optional parameter used only when reading multiple
            // DBR files after coupled modeling task interruption in GEM-Selektor
        {
            setNodeArray( dbr_file, i, type_f );
        }
        else
        {
            // Copying data from work DATABR structure into the node array
            // (as specified in nodeTypes array)
            setNodeArray( i, nodeTypes  );
        }
        i++;
    }  // end while()

    pVisor_Message( true );

    ErrorIf( i==0, datachbr_fn.c_str(), "No DBR_DAT files read!" );
    checkNodeArray( i, nodeTypes, datachbr_fn.c_str()  );
}

void  TNodeArray::checkNodeArray(
        long int i, long int* nodeTypes, const char*  datachbr_file )
{
    if(nodeTypes)
        for( long int ii=0; ii<anNodes; ii++)
            if(   nodeTypes[ii]<0 || nodeTypes[ii] >= i )
            {
                TNode::node_logger->warn("{} {} i= {}", anNodes, nodeTypes[ii], i);
                Error( datachbr_file, "Undefined boundary condition!" );
            }
}

//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array NodT0
// and read DATABR structure into the node array NodT1 from file
// dbr_file
//
//-------------------------------------------------------------------

void  TNodeArray::setNodeArray( std::string& dbr_file, long int ndx, GEMS3KGenerator::IOModes type_f )
{
    replace( dbr_file, "dbr-0-","dbr-1-" );
    calcNode->read_dbr_format_file( dbr_file,  type_f );

    NodT0[ndx] = allocNewDBR( calcNode);
    NodT1[ndx] = allocNewDBR( calcNode);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT0);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT1);
}


// Writing dataCH, dataBR structure to binary/text files
// and other necessary GEM2MT files
std::string TNodeArray::genGEMS3KInputFiles(  const std::string& filepath, ProcessProgressFunction message,
                                              long int nIV, GEMS3KGenerator::IOModes type_f, bool brief_mode,
                                              bool with_comments,   bool putNodT1, bool addMui )
{
    std::fstream fout_dat_lst;
    std::fstream fout_dbr_lst;
    GEMS3KGenerator generator( filepath, nIV, type_f );
    calcNode->current_output_set_name = generator.get_name();

    // open *-dat.lst
    fout_dat_lst.open( filepath, std::ios::out );
    fout_dat_lst << generator.gen_dat_lst_head();

    // open *-dbr.lst
    std::string dbr_lst_file_path = generator.get_dbr_file_lst_path();
    fout_dbr_lst.open( dbr_lst_file_path, std::ios::out);

    switch( generator.files_mode() )
    {
    case GEMS3KGenerator::f_binary:
    {
        GemDataStream  ff( generator.get_ipm_path(), std::ios::out|std::ios::binary );
        calcNode->multi_ptr()->out_multi( ff  );

        GemDataStream  f_ch( generator.get_dch_path(), std::ios::out|std::ios::binary);
        dbr_dch_api::datach_to_file(calcNode->CSD, f_ch);
    }
        break;
    default:
    {
        std::fstream ff( generator.get_ipm_path(), std::ios::out );
        ErrorIf( !ff.good(), generator.get_ipm_path(), "Fileopen error");
        calcNode->multi_ptr()->write_ipm_format_stream( ff, generator.files_mode(), addMui, with_comments, brief_mode, calcNode->output_set_name() );

        std::fstream  f_ch( generator.get_dch_path(), std::ios::out);
        dbr_dch_api::write_dch_format_stream(calcNode->current_output_set_name, calcNode->CSD, f_ch, generator.files_mode(), with_comments, brief_mode );

        if(generator.files_mode()>=GEMS3KGenerator::f_thermofun )
        {
            std::fstream  f_fun( generator.get_thermofun_path(), std::ios::out);
            calcNode->write_ThermoFun_format_stream(f_fun, false);
        }
    }
        break;
    }

    nIV = std::min( nIV, nNodes() );
    bool first = true;
    for( long int ii = 0; ii < nIV; ii++ )
    {
        if( !NodT0[ii] )
            continue;

        message( "Writing to disk a set of node array files from interrupted RMT task. \nPlease, wait...", ii );
        CopyWorkNodeFromArray( calcNode, ii, anNodes, NodT0 );
        auto dbr_file_name = generator.gen_dbr_file_name( 0, ii );
        calcNode->write_dbr_format_file( generator.get_path(dbr_file_name), generator.files_mode(), with_comments, brief_mode );

        fout_dat_lst << " \"" << dbr_file_name << "\"";
        if( !first )
            fout_dbr_lst << ",";
        fout_dbr_lst << " \"" << dbr_file_name << "\"";
        first = false;

        if( putNodT1 && NodT1[ii]) // put NodT1[ii] data
        {

            CopyWorkNodeFromArray( calcNode, ii, anNodes, NodT1 );
            auto dbr2_file_name = generator.gen_dbr_file_name( 1, ii );
            calcNode->write_dbr_format_file( generator.get_path(dbr2_file_name), generator.files_mode(), with_comments, brief_mode );
        }
    } // ii

    // Add full multy to test
    auto multi_txt = filepath+".txt";
    calcNode->multi_ptr()->to_text_file(multi_txt.c_str(), false);
    return dbr_lst_file_path;
}


void  TNodeArray::GEMS3k_write_dbr( const char* fname, GEMS3KGenerator::IOModes  type_f,
                                    bool with_comments, bool brief_mode )
{
    calcNode->packDataBr();
    calcNode->GEM_write_dbr( fname,  type_f, with_comments, brief_mode );
}

void  TNodeArray::GEMS3k_write_dbr( long int ndx, const char* fname, GEMS3KGenerator::IOModes  type_f,
                                    bool with_comments, bool brief_mode )
{
    CopyWorkNodeFromArray( calcNode, ndx, anNodes, pNodT1() );
    calcNode->GEM_write_dbr( fname,  type_f, with_comments, brief_mode );
}

long int  TNodeArray::GEMS3k_write_dbr( long int ndx, std::string& dbr_json )
{
    CopyWorkNodeFromArray( calcNode, ndx, anNodes, pNodT1() );
    return calcNode->GEM_write_dbr( dbr_json );
}

long int   TNodeArray::GEMS3k_read_dbr( long int ndx, std::string& dbr_file, GEMS3KGenerator::IOModes  type_f )
{
    auto ret = calcNode->GEM_read_dbr( dbr_file,  type_f );
    if( ret )
        return ret;
    if( !NodT0[ndx] )
        NodT0[ndx] = allocNewDBR( calcNode);
    if( !NodT1[ndx] )
        NodT1[ndx] = allocNewDBR( calcNode);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT0);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT1);
    return 0;
}

long int   TNodeArray::GEMS3k_read_dbr( long int ndx, std::string dbr_json, const bool check_dch_compatibility )
{
    auto ret = calcNode->GEM_read_dbr( dbr_json, check_dch_compatibility );
    if( ret )
        return ret;
    if( !NodT0[ndx] )
        NodT0[ndx] = allocNewDBR( calcNode);
    if( !NodT1[ndx] )
        NodT1[ndx] = allocNewDBR( calcNode);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT0);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT1);
    return 0;
}

//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array
// (as specified in nodeTypes array, ndx index of dataBR files in
//    the ipmfiles_lst_name list).
//
//-------------------------------------------------------------------
void  TNodeArray::setNodeArray( long int ndx, long int* nodeTypes  )
{

    for( long int ii=0; ii<anNodes; ii++)
        if(  (!nodeTypes && ndx==0) ||
             ( nodeTypes && (nodeTypes[ii] == ndx/*i+1*/ )) )
        {
            calcNode->pCNode()->NodeHandle = ndx/*(i+1)*/;
            NodT0[ii] = allocNewDBR( calcNode);
            NodT1[ii] = allocNewDBR( calcNode);

            MoveWorkNodeToArray( calcNode, ii, anNodes, NodT0);
            //  CopyWorkNodeFromArray( calcNode, ii, anNodes,NodT0);
            MoveWorkNodeToArray( calcNode, ii, anNodes, NodT1);
            //  CopyWorkNodeFromArray( calcNode, ii, anNodes,NodT1);
        }
}

#endif

//-----------------------End of nodearray.cpp--------------------------

