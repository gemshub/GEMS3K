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

#include <cmath>
#ifdef _MSC_VER
#include <stdint.h>
#else
#include <unistd.h>
#endif
#include <algorithm>
#include "v_detail.h"

#ifdef NODEARRAYLEVEL
#include "nodearray.h"
#include "gdatastream.h"
#include "io_keyvalue.h"
#include "io_nlohmann.h"

TNodeArray::TNodeArray( long int nNod  ):
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
    allocMemory();
    na = this;
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
    allocMemory();
    na = this;
}


TNodeArray::~TNodeArray()
{
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
    { if( NodT0 )
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
    TNode wrkNode;//( calcNode ); // must be copy TNode internal
    DATABRPTR* C0 = pNodT0();
    DATABRPTR* C1 = pNodT1();
    bool* iaN = piaNode();     // indicators for IA in the nodes

    start_node = max( start_node, 0L );
    end_node = min( end_node, anNodes );


#pragma omp parallel shared( C0, C1, iaN, diffile, iRet ) private( wrkNode, n )
    {
        TNode wrkNode( calcNode ); // must be copy TNode internal
        n = omp_get_thread_num();
#pragma omp for
        for( ii = start_node; ii<= end_node; ii++) // node iteration
        {
            if( !CalcIPM_Node(  modeParam, wrkNode, ii, C0, C1, iaN, diffile ) )
            {
#pragma omp atomic write
                iRet = false;
            }

            //#pragma omp critical
            //{
            //   cout << n << "-thread did index: " << ii << endl;
            //}
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

    start_node =std:: max( start_node, 0L );
    end_node = std::min( end_node, anNodes );


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

    if( needGEM ) {
        long RetCode =  RunGEM( wrkNode, ii, Mode, C1 );

        // checking RetCode from GEM IPM calculation
        if( !(RetCode==OK_GEM_AIA || RetCode == OK_GEM_SIA ))
        {
            std::string err_msg = ErrorGEMsMessage( RetCode,  ii, modeParam.step  );
            iRet = false;
            if( diffile )
            {
                // write to file here
                fprintf( diffile, "\nError reported from GEMS3K module\n%s\n",  err_msg.c_str() );
            }
        }
    }
    else { // GEM calculation for this node not needed
        C1[ii]->IterDone = 0; // number of GEMIPM iterations is set to 0
    }
    return iRet;
}

#endif


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

    // cout << ipmfiles_lst_name << "  " << dbrfiles_lst_name << endl;
    std::string curPath = ""; //current reading file path
    std::fstream f_log(TNode::ipmLogFile.c_str(), std::ios::out|std::ios::app );
    try
    {
        bool binary_f = false;

        //  open file stream for the file names list file
        std::fstream f_lst( ipmfiles_lst_name, std::ios::in );
        ErrorIf( !f_lst.good() , ipmfiles_lst_name, "Fileopen error");

        std::string datachbr_fn;
        f_getline( f_lst, datachbr_fn, ' ');

        //Testing flag "-t" or "-b" (by default "-t")   // use binary or text files
        size_t pos = datachbr_fn.find( '-');
        if( pos != std::string::npos )
        {
            if( datachbr_fn[pos+1] == 'b' )
                binary_f = true;
            f_getline( f_lst, datachbr_fn, ' ');
        }
        // Reading DBR_DAT files from dbrfiles_lst_name
        if(  dbrfiles_lst_name )
            InitNodeArray( dbrfiles_lst_name, nodeTypes, getNodT1, binary_f  );
        else
            if( nNodes() ==1 )
                setNodeArray( 0 , nullptr  );
            else // undefined TNodeArray
                Error( "GEM_init", "GEM_init() error: Undefined boundary condition!" );
        return 0;

    }
    catch(TError& err) {
        if( !curPath.empty() )
            f_log << "GEMS3K input : file " << curPath.c_str() << std::endl;
        f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
    }
    catch(...) {
        return -1;
    }
    return 1;
}


//-------------------------------------------------------------------
// Initialization of TNodeArray data structures. Reads in the DBR text input files and
// copying data from work DATABR structure into the node array
// (as specified in nodeTypes array, ndx index of dataBR files in
//    the dbrfiles_lst_name list).
//
//-------------------------------------------------------------------
void  TNodeArray::InitNodeArray( const char *dbrfiles_lst_name,
                                 long int *nodeTypes, bool getNodT1, bool binary_f  )
{
    int i;
    std::string datachbr_fn;

    std::string curPath = ""; //current reading file path
    std::string lst_in = dbrfiles_lst_name;
    std::string Path = "";

    // Get path
    size_t pos = lst_in.rfind("\\");
    if( pos == std::string::npos )
        pos = lst_in.rfind("/");
    else
        pos = std::max(pos, lst_in.rfind("/") );

    if( pos < std::string::npos )
        Path = lst_in.substr(0, pos+1);

    //  open file stream for the file names list file
    std::fstream f_lst( lst_in.c_str(), std::ios::in );
    ErrorIf( !f_lst.good() , lst_in.c_str(), "Fileopen error");


    // Prepare for reading DBR_DAT files
    i = 0;
    while( !f_lst.eof() )  // For all DBR_DAT files listed
    {

        pVisor_Message(false, i, nNodes() );

        // Reading DBR_DAT file into work DATABR structure
        if( i )  // Comma only after the first DBR_DAT file!
            f_getline( f_lst, datachbr_fn, ',');
        else
            f_getline( f_lst, datachbr_fn, ' ');

        std::string dbr_file = Path + datachbr_fn;
        curPath = dbr_file;
        if( binary_f )
        {
            GemDataStream in_br(dbr_file, std::ios::in|std::ios::binary);
            calcNode->databr_from_file(in_br);
        }
        else
        {   std::fstream in_br(dbr_file.c_str(), std::ios::in );
            ErrorIf( !in_br.good() , datachbr_fn.c_str(),
                     "DBR_DAT fileopen error");
            calcNode->databr_from_text_file(in_br);
        }
        curPath = "";
        // Unpacking work DATABR structure into MULTI (GEM IPM work structure): uses DATACH
        //    unpackDataBr();

        if( getNodT1 )  // optional parameter used only when reading multiple
            // DBR files after coupled modeling task interruption in GEM-Selektor
        {
            setNodeArray( dbr_file, i, binary_f );
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

    ErrorIf( i==0, datachbr_fn.c_str(), "GEM_init() error: No DBR_DAT files read!" );
    checkNodeArray( i, nodeTypes, datachbr_fn.c_str()  );
}

void  TNodeArray::checkNodeArray(
        long int i, long int* nodeTypes, const char*  datachbr_file )
{
    if(nodeTypes)
        for( long int ii=0; ii<anNodes; ii++)
            if(   nodeTypes[ii]<0 || nodeTypes[ii] >= i )
            {
                std::cout << anNodes << " " << nodeTypes[ii] << " i = " << i<< std::endl;
                Error( datachbr_file,
                       "GEM_init() error: Undefined boundary condition!" );
            }
}

//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array NodT0
// and read DATABR structure into the node array NodT1 from file
// dbr_file
//
//-------------------------------------------------------------------

void  TNodeArray::setNodeArray( std::string& dbr_file, long int ndx, bool binary_f )
{
    replace(dbr_file, "dbr-0-","dbr-1-");
    if( binary_f )
    {
        GemDataStream in_br(dbr_file, std::ios::in|std::ios::binary);
        calcNode->databr_from_file(in_br);
    }
    else
    {   std::fstream in_br(dbr_file.c_str(), std::ios::in );
        ErrorIf( !in_br.good() , dbr_file.c_str(),
                 "DataBR Fileopen error");
        calcNode->databr_from_text_file(in_br);
    }

    NodT0[ndx] = allocNewDBR( calcNode);
    NodT1[ndx] = allocNewDBR( calcNode);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT0);
    MoveWorkNodeToArray(calcNode, ndx, anNodes, NodT1);
}

//// Writing dataCH, dataBR structure to binary/text files
//// and other necessary GEM2MT files
//std::string TNodeArray::PutGEM2MTFiles(  QWidget* par, long int nIV,
//                                     bool bin_mode, bool brief_mode, bool with_comments,
//                                     bool putNodT1, bool addMui )
//{
//    // Get name of filenames structure
//    std::string path = std::string( rt[RT_SYSEQ].FldKey(2), 0, rt[RT_SYSEQ].FldLen(2));;
//    strip(path);
//    if( bin_mode )
//        path += "-bin.lst";
//    else
//        path += "-dat.lst";

//AGAIN:
//    // open file to output
//    if( vfChooseFileSave(par, path, "Please, enter IPM work structure file name", "*.lst" ) == false )
//        return "";

//    if( !access(path.c_str(), 0 ) ) //file exists
//        switch( vfQuestion3( par, path.c_str(), "This set of files exists!",
//                             "&Overwrite", "&Rename", "&Cancel") )
//        {
//        case VF3_2:
//            goto AGAIN;
//        case VF3_1:
//            break;
//        case VF3_3:
//            return path;
//        }

//    ProcessProgressFunction messageF = [nIV, par](const std::string& message, long point){
//        return  pVisor->Message( par, "GEM2MT node array",  message.c_str() , point, nIV );
//    };
//    genGEMS3KInputFiles(  path, messageF, nIV, bin_mode, brief_mode, with_comments,
//                          putNodT1, addMui );

//    pVisor->CloseMessage();
//    return path;
//}

// Writing dataCH, dataBR structure to binary/text files
// and other necessary GEM2MT files
std::string TNodeArray::genGEMS3KInputFiles(  const std::string& filepath, ProcessProgressFunction message,
                                              long int nIV, bool bin_mode, bool brief_mode, bool with_comments,
                                              bool putNodT1, bool addMui )
{
    std::fstream fout;
    std::fstream fout2;
    std::string Path_;
    std::string dir;
    std::string name;
    std::string newname;
    std::string path;
    char buf[20];

    path = filepath;
    u_splitpath( path, dir, name, newname );

    // get name
    auto pos = name.rfind("-");
    if( pos != std::string::npos )
        name = name.substr(0, pos);

    // making special files
    // put data to pmfiles-bin.lst file
    if( bin_mode )
    {
        fout.open(path.c_str(), std::ios::out);
        fout << "-b \"" << name.c_str() << "-dch.bin\"";
        fout << " \"" << name.c_str() << ".ipm\" ";
    }
    // put data to pmfiles-dat.lst file
    else
    {   fout.open(path.c_str(), std::ios::out);
        fout << "-t \"" << name.c_str() << "-dch."<<dat_ext << "\"";
        fout << " \"" << name.c_str() << "-ipm."<<dat_ext << "\" ";
    }

    std::string path2 = name;
    path2 += "-dbr";
    path2 = u_makepath( dir, path2, "lst" );
    fout2.open(path2.c_str(), std::ios::out);

    if( bin_mode )
    {
        //  putting MULTI to binary file
        Path_ = u_makepath( dir, name, "ipm" );
        GemDataStream  ff(Path_, std::ios::out|std::ios::binary);
        calcNode->multi->out_multi( ff, Path_  );
    }
    else
    {
        // output MULTI to txt file
        newname = name+"-ipm";
        Path_ = u_makepath( dir, newname, dat_ext );
        std::fstream ff( Path_.c_str(), std::ios::out );
        ErrorIf( !ff.good() , Path_, "Fileopen error");
        calcNode->multi->to_text_file_gemipm( ff, addMui, with_comments, brief_mode );
    }

    // out dataCH to binary file
    newname = name+"-dch";
    if( bin_mode )
    {  Path_ = u_makepath( dir, newname, "bin" );
        GemDataStream  f_ch1(Path_, std::ios::out|std::ios::binary);
        calcNode->datach_to_file(f_ch1);
        f_ch1.close();
    }
    // out dataCH to text file
    else
    {  //newname = name+"-dch";
        Path_ = u_makepath( dir, newname, dat_ext );
        std::fstream  f_ch2(Path_.c_str(), std::ios::out);
        calcNode->datach_to_text_file(f_ch2, with_comments, brief_mode, Path_.c_str() );
        f_ch2.close();
    }

    nIV = std::min( nIV, nNodes() );
    bool first = true;
    for( long int ii = 0; ii < nIV; ii++ )
    {
        if( !NodT0[ii] )
            continue;

        message( "Writing to disk a set of node array files from interrupted RMT task. \nPlease, wait...", ii );
        // Save databr
        CopyWorkNodeFromArray( calcNode, ii, anNodes, NodT0 );

        sprintf( buf, "%4.4ld", ii );
        // dataBR files - binary
        if( bin_mode )
        {
            newname =  name + + "-dbr-0-"  + buf;
            Path_ = u_makepath( dir, newname, "bin" );
            GemDataStream  f_br1(Path_, std::ios::out|std::ios::binary);
            calcNode->databr_to_file(f_br1);
            f_br1.close();
            if( first )
                fout << " \"" << newname.c_str() << ".bin\"";
            if( !first )
                fout2 << ",";
            fout2 << " \"" << newname.c_str() << ".bin\"";
        }
        else
        {
            newname = name + "-dbr-0-" + buf;
            Path_ = u_makepath( dir, newname, dat_ext );
            std::fstream  f_br2(Path_.c_str(), std::ios::out);
            calcNode->databr_to_text_file(f_br2, with_comments, brief_mode, Path_.c_str() );
            f_br2.close();
            if( first )
                fout << " \"" << newname.c_str() << "."<<dat_ext << "\"";
            if( !first )
                fout2 << ",";
            fout2 << " \"" << newname.c_str() << "."<<dat_ext << "\"";
        }
        first = false;

        if( putNodT1 && NodT1[ii]) // put NodT1[ii] data
        {

            // Save databr
            CopyWorkNodeFromArray( calcNode, ii, anNodes, NodT1 );

            // dataBR files - binary
            if( bin_mode )
            {
                newname =  name +  "-dbr-1-"  + buf;
                Path_ = u_makepath( dir, newname, "bin" );
                GemDataStream  f_br1(Path_, std::ios::out|std::ios::binary);
                calcNode->databr_to_file(f_br1);
                f_br1.close();
                //         fout << ", \"" << newname.c_str() << ".bin\"";
            }
            else
            {
                newname = name + "-dbr-1-" + buf;
                Path_ = u_makepath( dir, newname, dat_ext );
                std::fstream  f_br2(Path_.c_str(), std::ios::out);
                calcNode->databr_to_text_file(f_br2, with_comments, brief_mode, Path_.c_str() );
                f_br2.close();
                //         fout << ", \"" << newname.c_str() << "."<< dat_ext << "\"";
            }
        }
    } // ii

    return path2; // dbr list file
}



void  TNodeArray::GEMS3k_write_dbr( const char* fname,  bool binary_f,
                                    bool with_comments, bool brief_mode )
{
    calcNode->packDataBr();
    calcNode->GEM_write_dbr( fname,  binary_f, with_comments, brief_mode );
}

#endif
//-----------------------End of nodearray_new.cpp--------------------------
