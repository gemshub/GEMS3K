//--------------------------------------------------------------------
// $Id: nodearray.cpp 684 2005-11-23 13:17:15Z gems $
//
// C/C++ interface between GEM IPM and FMT node array
// Working whith DATACH and DATABR structures
//
// Copyright (C) 2004,2007 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef NOPARTICLEARRAY
#include "particlearray.h"
#endif

#include "nodearray.h"
#include "gdatastream.h"
#include "num_methods.h"
#include <math.h>

#ifndef __unix
#include <io.h>
#endif

#ifndef IPMGEMPLUGIN
  #include "service.h"
  #include "visor.h"
#else
  istream& f_getline(istream& is, gstring& str, char delim);
#endif

TNodeArray* TNodeArray::na;

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

int  TNodeArray::RunGEM( int  iNode, int Mode )
{

bool uPrimalSol = false;  
  if( Mode < 0 || (short)abs(Mode) == NEED_GEM_SIA )
	  uPrimalSol = true;
	  
// Copy data from the iNode node from array NodT1 to the work DATABR structure
   CopyWorkNodeFromArray( iNode, anNodes, NodT1 );

// GEM IPM calculation of equilibrium state in MULTI
  pCNode()->NodeStatusCH = (short)abs(Mode);
  int retCode = GEM_run( uPrimalSol );

// Copying data for node iNode back from work DATABR structure into the node array
//   if( retCode == OK_GEM_AIA ||
//       retCode == OK_GEM_PIA  )
   MoveWorkNodeToArray( iNode, anNodes, NodT1 );

   return retCode;
}

void  TNodeArray::checkNodeArray(
    int i, int* nodeTypes, const char*  datachbr_file )
{
 if(nodeTypes)
   for( int ii=0; ii<anNodes; ii++)
     if(   nodeTypes[ii]<0 || nodeTypes[ii] >= i )
        {
	  cout << anNodes << " " << nodeTypes[ii] << " i = " << i<< endl;
	Error( datachbr_file,
          "GEM_init() error: Undefined boundary condition!" );
	}
}

//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array
// (as specified in nodeTypes array, ndx index of dataBR files in
//    the ipmfiles_lst_name list).
//
//-------------------------------------------------------------------
void  TNodeArray::setNodeArray( int ndx, int* nodeTypes  )
{
   for( int ii=0; ii<anNodes; ii++)
            if(  (!nodeTypes && ndx==0) ||
              ( nodeTypes && (nodeTypes[ii] == ndx/*i+1*/ )) )
                  {    pCNode()->NodeHandle = (short)ndx/*(i+1)*/;
                       MoveWorkNodeToArray(ii, anNodes, NodT0);
                       CopyWorkNodeFromArray(ii, anNodes,NodT0);
                       MoveWorkNodeToArray(ii, anNodes, NodT1);
                       CopyWorkNodeFromArray(ii, anNodes,NodT1);
                   }
}

#ifndef IPMGEMPLUGIN

//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array NodT0
// and read DATABR structure into the node array NodT1 from file
// dbr_file
//
//-------------------------------------------------------------------

void  TNodeArray::setNodeArray( gstring& dbr_file, int ndx, bool binary_f )
{
   MoveWorkNodeToArray(ndx, anNodes, NodT0);
   dbr_file = dbr_file.replace("dbr-0-","dbr-1-");
   if( binary_f )
   {
      GemDataStream in_br(dbr_file, ios::in|ios::binary);
      databr_from_file(in_br);
   }
   else
   {   fstream in_br(dbr_file.c_str(), ios::in );
         ErrorIf( !in_br.good() , dbr_file.c_str(),
                    "DataBR Fileopen error");
       databr_from_text_file(in_br);
    }
   MoveWorkNodeToArray(ndx, anNodes, NodT1);
}

// Writing dataCH, dataBR structure to binary/text files
// and other necessary GEM2MT files
gstring TNodeArray::PutGEM2MTFiles( QWidget* par, int nIV,
      bool addMui, bool bin_mode, bool putNodT1  )
{
  fstream fout;
  gstring Path_;
  gstring dir;
  gstring name;
  gstring newname;
  gstring path;
  char buf[20];

  // Get name of filenames structure
   path = gstring( rt[RT_SYSEQ].FldKey(2), 0, rt[RT_SYSEQ].FldLen(2));;
   path.strip();
   if( bin_mode )
     path += "-bin.lst";
   else
     path += "-dat.lst";

AGAIN:
      // open file to output
   if( vfChooseFileSave(par, path,
          "Please, enter IPM work structure file name", "*.lst" ) == false )
               return "";
   u_splitpath( path, dir, name, newname );
   if( !access(path.c_str(), 0 ) ) //file exists
     switch( vfQuestion3( par, name.c_str(),
        "This set of files exists!",
             "&Overwrite", "&Rename", "&Cancel") )
            {
            case VF3_2:
                goto AGAIN;
            case VF3_1:
                break;
            case VF3_3:
                return path;
            }

// get name
   size_t pos = name.rfind("-");
   if( pos != gstring::npos )
      name = name.substr(0, pos);

// making special files
// put data to pmfiles-bin.lst file
   if( bin_mode )
   {
       fout.open(path.c_str(), ios::out);
       fout << "-b \"" << name.c_str() << "-dch.bin\"";
       fout << " \"" << name.c_str() << ".ipm\" ";
   }
// put data to pmfiles-dat.lst file
   else
   {   fout.open(path.c_str(), ios::out);
       fout << "-t \"" << name.c_str() << "-dch.dat\"";
       fout << " \"" << name.c_str() << "-ipm.dat\" ";
   }

  if( bin_mode )
  {
//  putting MULTI to binary file
    Path_ = u_makepath( dir, name, "ipm" );
    GemDataStream  ff(Path_, ios::out|ios::binary);
    TProfil::pm->outMulti( ff, Path_  );
  }
 else
  {
// output MULTI to txt file
    newname = name+"-ipm";
    Path_ = u_makepath( dir, newname, "dat" );
    TProfil::pm->outMulti( Path_, addMui  );
  }

// out dataCH to binary file
   newname = name+"-dch";
   if( bin_mode )
   {  Path_ = u_makepath( dir, newname, "bin" );
      GemDataStream  f_ch1(Path_, ios::out|ios::binary);
      datach_to_file(f_ch1);
      f_ch1.close();
   }
// out dataCH to text file
   else
   {  //newname = name+"-dch";
      Path_ = u_makepath( dir, newname, "dat" );
      fstream  f_ch2(Path_.c_str(), ios::out);
      datach_to_text_file(f_ch2);
      f_ch2.close();
   }

 nIV = min( nIV, nNodes() );
 bool first = true;
 for( int ii = 0; ii < nIV; ii++ )
 {
   if( !NodT0[ii] )
      continue;

   pVisor->Message( par, "GEM2MT node array",
      "Writing to disk a set of node array files from interrupted RMT task.\n"
           "Please, wait...", ii, nIV );
   // Save databr
   CopyWorkNodeFromArray( ii, anNodes, NodT0 );

   sprintf( buf, "%4.4d", ii );
   // dataBR files - binary
   if( bin_mode )
    {
       newname =  name + + "-dbr-0-"  + buf;
       Path_ = u_makepath( dir, newname, "bin" );
       GemDataStream  f_br1(Path_, ios::out|ios::binary);
       databr_to_file(f_br1);
       f_br1.close();
       if( !first )
          fout << ",";
       fout << " \"" << newname.c_str() << ".bin\"";
     }
     else
     {
        newname = name + "-dbr-0-" + buf;
        Path_ = u_makepath( dir, newname, "dat" );
        fstream  f_br2(Path_.c_str(), ios::out);
        databr_to_text_file(f_br2);
        f_br2.close();
        if( !first )
           fout << ",";
        fout << " \"" << newname.c_str() << ".dat\"";
     }
     first = false;

   if( putNodT1 && NodT1[ii]) // put NodT1[ii] data
   {

       // Save databr
      CopyWorkNodeFromArray( ii, anNodes, NodT1 );

      // dataBR files - binary
      if( bin_mode )
      {
         newname =  name + + "-dbr-1-"  + buf;
         Path_ = u_makepath( dir, newname, "bin" );
         GemDataStream  f_br1(Path_, ios::out|ios::binary);
         databr_to_file(f_br1);
         f_br1.close();
//         fout << ", \"" << newname.c_str() << ".bin\"";
      }
     else
      {
         newname = name + "-dbr-1-" + buf;
         Path_ = u_makepath( dir, newname, "dat" );
         fstream  f_br2(Path_.c_str(), ios::out);
         databr_to_text_file(f_br2);
         f_br2.close();
//         fout << ", \"" << newname.c_str() << ".dat\"";
      }
   }
 } // ii
 pVisor->CloseMessage();
 return path;
}

#endif

//---------------------------------------------------------//

void TNodeArray::allocMemory()
{
  int ii;

#ifndef NOPARTICLEARRAY
TParticleArray::pa = 0;
#endif

// alloc memory for data bidge structures
// did in constructor TNode::allocMemory();

// alloc memory for all nodes at current time point
    NodT0 = new  DATABRPTR[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        NodT0[ii] = 0;

// alloc memory for all nodes at previous time point
    NodT1 = new  DATABRPTR[anNodes];
    for(  ii=0; ii<anNodes; ii++ )
        NodT1[ii] = 0;
    
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
   int ii;

   if( anNodes )
   { if( NodT0 )
       for(  ii=0; ii<anNodes; ii++ )
        if( NodT0[ii] )
           NodT0[ii] = databr_free(NodT0[ii]);
     delete[]  NodT0;
     NodT0 = 0;
     if( NodT1 )
       for(  ii=0; ii<anNodes; ii++ )
        if( NodT1[ii] )
          NodT1[ii] = databr_free(NodT1[ii]);
     delete[]  NodT1;
     NodT1 = 0;
   }

  if( grid )
     delete[] grid;
  if( tcNode)
     delete[] tcNode;
  if( iaNode )
	 delete[] iaNode; 
}

#ifndef IPMGEMPLUGIN

TNodeArray::TNodeArray( int nNod, MULTI *apm  ):TNode( apm )
{
    anNodes = nNod;
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

TNodeArray::TNodeArray( int asizeN, int asizeM, int asizeK, MULTI *apm  ):
TNode( apm ), sizeN(asizeN), sizeM(asizeM), sizeK(asizeK)
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


#else

TNodeArray::TNodeArray( int nNod  ):
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

TNodeArray::TNodeArray( int asizeN, int asizeM, int asizeK ):
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

#endif


TNodeArray::~TNodeArray()
{
   freeMemory();
}


// Copying data for node ii from node array into work DATABR structure
void TNodeArray::CopyWorkNodeFromArray( int ii, int nNodes, DATABRPTR* arr_BR )
{
  // from arr_BR[ii] to pCNode() structure
  if( ii < 0 || ii>= nNodes )
    return;
  // memory must be allocated before

  memcpy( &pCNode()->NodeHandle, &arr_BR[ii]->NodeHandle, 6*sizeof(short));
  memcpy( &pCNode()->TC, &arr_BR[ii]->TC, 32*sizeof(double));
// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
// exchange of values occurs through lists of indices, e.g. xDC, xPH
  memcpy( pCNode()->xDC, arr_BR[ii]->xDC, pCSD()->nDCb*sizeof(double) );
  memcpy( pCNode()->gam, arr_BR[ii]->gam, pCSD()->nDCb*sizeof(double) );
  if( CSD->nAalp >0 )
    memcpy( pCNode()->aPH, arr_BR[ii]->aPH, pCSD()->nPHb*sizeof(double) );
  else  pCNode()->aPH = 0;
  memcpy( pCNode()->xPH, arr_BR[ii]->xPH, pCSD()->nPHb*sizeof(double) );
  memcpy( pCNode()->vPS, arr_BR[ii]->vPS, pCSD()->nPSb*sizeof(double) );
  memcpy( pCNode()->mPS, arr_BR[ii]->mPS, pCSD()->nPSb*sizeof(double) );

  memcpy( pCNode()->bPS, arr_BR[ii]->bPS,
                          pCSD()->nPSb*pCSD()->nICb*sizeof(double) );
  memcpy( pCNode()->xPA, arr_BR[ii]->xPA, pCSD()->nPSb*sizeof(double) );
  memcpy( pCNode()->dul, arr_BR[ii]->dul, pCSD()->nDCb*sizeof(double) );
  memcpy( pCNode()->dll, arr_BR[ii]->dll, pCSD()->nDCb*sizeof(double) );
  memcpy( pCNode()->bIC, arr_BR[ii]->bIC, pCSD()->nICb*sizeof(double) );
  memcpy( pCNode()->rMB, arr_BR[ii]->rMB, pCSD()->nICb*sizeof(double) );
  memcpy( pCNode()->uIC, arr_BR[ii]->uIC, pCSD()->nICb*sizeof(double) );
  pCNode()->dRes1 = 0;
}

// Copying data for node iNode back from work DATABR structure into the node array
void TNodeArray::MoveWorkNodeToArray( int ii, int nNodes, DATABRPTR* arr_BR )
{
  if( ii < 0 || ii>= nNodes )
    return;
  if( arr_BR[ii] )
  {
       arr_BR[ii] = databr_free( arr_BR[ii]);
       // delete[] arr_BR[ii];
  }
  arr_BR[ii] = pCNode();
// alloc new memory
  TNode::CNode = new DATABR;
  databr_realloc();
  memset( &pCNode()->TC, 0, 32*sizeof(double));
}

void TNodeArray::CopyNodeFromTo( int ndx, int nNod,
                       DATABRPTR* arr_From, DATABRPTR* arr_To )
{
  if( !arr_From || !arr_To )
      return;
  CopyWorkNodeFromArray( ndx, nNod, arr_From );
  MoveWorkNodeToArray( ndx,  nNod, arr_To );
}

//---------------------------------------------------------
// Methods for working with node arrays (access to data from DBR)

// Calculate phase (carrier) mass, g  of single component phase
double TNodeArray::get_mPH( int ia, int nodex, int PHx )
{
  int DCx = Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
  double val=0.;

  if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
  {
    val = pCSD()->DCmm[DCx];
    if( ia == 0)
     val *= pNodT0()[nodex]->xDC[DC_xCH_to_xDB(DCx)];
    else
     val *= pNodT1()[nodex]->xDC[DC_xCH_to_xDB(DCx)];
  }

  return val;
}

// Calculate phase volume, cm3/mol  of single component phase
double TNodeArray::get_vPH( int ia, int nodex, int PHx )
{
  int DCx = Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
  double val=0.;

  if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
  {
     double T, P;
     if( ia == 0 )
     {
      T = pNodT0()[(nodex)]->TC;
      P = pNodT0()[(nodex)]->P;
      val = pNodT0()[nodex]->xDC[DC_xCH_to_xDB(DCx)]; // number of moles
     }
     else
     {
      T = pNodT1()[(nodex)]->TC;
      P = pNodT1()[(nodex)]->P;
      val = pNodT1()[nodex]->xDC[DC_xCH_to_xDB(DCx)];
     }
     val *= DC_V0_TP( DCx, T, P );
  }
  return val;
}


// Calculate bulk compositions  of single component phase
double TNodeArray::get_bPH( int ia, int nodex, int PHx, int ICx )
{
  int DCx = Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
  double val=0.;

  if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
  {
    val = pCSD()->A[ pCSD()->xIC[ICx] + DCx * pCSD()->nIC];
    if( ia == 0)
     val *= pNodT0()[nodex]->xDC[DC_xCH_to_xDB(DCx)];
    else
     val *= pNodT1()[nodex]->xDC[DC_xCH_to_xDB(DCx)];
  }

  return val;
}

//---------------------------------------------------------
// working with grid

// Set grid coordinate array use predefined array aGrid
// or set up regular scale
void TNodeArray::SetGrid( float aSize[3], float (*aGrid)[3] )
{
  int i, j, k, ndx;
  size.x =  aSize[0];
  size.y =  aSize[1];
  size.z =  aSize[2];

  LOCATION delta( size.x/sizeN, size.y/sizeM, size.z/sizeK );
  if( !grid )
      grid = new LOCATION[ anNodes ];

  for( i = 0; i < sizeN; i++ )
    for( j = 0; j < sizeM; j++ )
      for( k = 0; k < sizeK; k++ )
      {
         ndx = iNode( i, j, k );
         if( aGrid )
         {
           grid[ndx].x = aGrid[ndx][0];
           grid[ndx].y = aGrid[ndx][1];
           grid[ndx].z = aGrid[ndx][2];
         }
         else
         {
           grid[ndx].x =delta.x*i;
           grid[ndx].y =delta.y*j;
           grid[ndx].z =delta.z*k;
         }
      }
   // outside limits settted in size
}

// test location in node
bool TNodeArray::isLocationInNode( int iNode, LOCATION cxyz ) const
{
  int i1, j1, k1;
  i1 = indN( iNode );
  j1 = indM( iNode );
  k1 = indK( iNode );
  return isLocationInNode( i1,j1, k1, cxyz ) ;
}

// test location in node
bool TNodeArray::isLocationInNode( int ii, int jj, int kk, LOCATION cxyz ) const
{
  LOCATION maxl;

  if( ii<0 || ii>= sizeN || jj<0 ||
      jj >= sizeM || kk <0 || kk >= sizeK )
    return false;

  int ndx = iNode( ii, jj, kk );
  maxl = getGrid( ii+1, jj+1, kk+1 ); // only for rectangular
  // must be changed
  //  x = const, find new y,z srez i
  // analiz pryamougol`nika pri y1 == const, poisk z21 i z22
  // analiz otrezka po z2
  if( grid[ndx].x <= cxyz.x &&
     ( cxyz.x < maxl.x || ( cxyz.x <= maxl.x && size.x == maxl.x ) )&&
      grid[ndx].y <= cxyz.y &&  cxyz.y <= maxl.y &&
      grid[ndx].z <= cxyz.z &&  cxyz.z <= maxl.z )
           return true;

   return false; // location behind the node
}

// Finds a node absolute index for the current
// point location (uses grid coordinate array grid[])
// performance-important functions to be used e.g. in particle tracking methods
int TNodeArray::FindNodeFromLocation( LOCATION cxyz, int old_node ) const
{
//  LOCATION maxl;
  int i, j, k/*, ndx*/;

 if( old_node == -1 )
 { // check all nodes
   for( i = 0; i < sizeN; i++ )
    for( j = 0; j < sizeM; j++ )
      for( k = 0; k < sizeK; k++ )
      {
        if(  isLocationInNode( i, j, k, cxyz ) )
         return iNode( i, j, k );
      }
 }
 else // check only nearest nodes
 {
   int i1, j1, k1;
   i1 = indN( old_node );
   j1 = indM( old_node );
   k1 = indK( old_node );

   for( i = i1-1; i <= i1+1; i++ )
     for( j = j1-1; j <= j1+1; j++ )
       for( k = k1-1; k <= k1+1; k++ )
       {
          if(  isLocationInNode( i, j, k, cxyz ) )
              return iNode( i, j, k );
       }
 }
 return -1; // behind region
}

// get current node location
// if iN, jN or kN more then corresponding sizeN, sizeM, sizeK
// return size of system
LOCATION TNodeArray::getGrid( int iN, int jN, int kN ) const
{
  LOCATION loc;
  int i1, j1, k1;

// only for test
  if( iN < 0 || iN > sizeN ||
      jN < 0 || jN > sizeM ||
      kN < 0 || kN > sizeK  )
   Error( "", "getGrid - programm error");

  if( iN == sizeN )
    i1 = iN-1;
  else i1 = iN;
  if( jN == sizeM )
     j1 = jN-1;
  else j1 = jN;
  if( kN == sizeK )
   k1 = kN-1;
   else k1 = kN;

  loc = grid[ iNode( i1, j1, k1)];
  if( i1 != iN )
    loc.x = size.x;
  if( j1 != jN ) loc.y = size.y;
  if( k1 != kN ) loc.z = size.z;

  return loc;
}

// get 3D sizes for node (  from cxyz[0] - to cxyz[1] )
// only for rectangular -  must be changed
// for any must be LOCATION cxyz[8]
void TNodeArray::GetNodeSizes( int ndx, LOCATION cxyz[2] )
{
  LOCATION maxl;
  int i, j, k;

  i = indN( ndx );
  j = indM( ndx );
  k = indK( ndx );

  cxyz[0] = grid[ndx];
  cxyz[1] = getGrid( i+1, j+1, k+1 );
  /*
  cxyz[1] = getGrid( i, j, k+1 );
  cxyz[2] = getGrid( i, j+1, k );
  cxyz[3] = getGrid( i, j+1, k+1 );
  cxyz[4] = getGrid( i+1, j, k );
  cxyz[5] = getGrid( i+1, j, k+1 );
  cxyz[6] = getGrid( i+1, j+1, k );
  cxyz[7] = getGrid( i+1, j+1, k+1 );
  */
}

// get full mass of all particles of the type ptype in node ndx
// ndx    -  (absolute) index of the node
// ptype  -  particle type index ( 1 to 255 )
// tcode  -  particle transport mechanism code (see enum PTCODE)
// ips   - DataBr index of phase or species to which this particle is connected
double TNodeArray::GetNodeMass( int ndx,
       char /*ptype*/, char tcode, unsigned char iips )
{
   double mass = 0.;
//   DATABR* dbr = NodT1[ndx];
   DATACH* dch = CSD;
   int xWatCH, ips = (int)iips;

     switch( tcode )
     {
        case DISSOLVED: // mass of dissolved matter in aqueous solution
                        xWatCH = dch->nDCinPH[dch->xPH[0]]-1; // CH index of water
//                       mass = node1_mPS(ndx,ips); // - node1_xPA(ndx,ips)*dch->DCmm[xWatCH];
                        mass = node1_xPA(ndx,ips)*dch->DCmm[xWatCH]; // Mass of aq-solvent
                        break;
        case ADVECTIVE: // mass of aq solution
                        mass = node1_mPS( ndx, ips );
                        break;
        case COLLOID:  // mass of phase - solid solution, sorption or pure solid
                       if( ips < dch->nPSb )
                          mass = node1_mPS( ndx, ips );
                        else
                          mass = node1_mPH( ndx, ips );
                       break;
        case DIFFUSIVE: // mass of the diffusing species
                        mass = nodeCH_DCmm( ips ) * node1_xDC( ndx, ips );
                        break;
        default: break;
     }
   return mass;
}

// move mass m_v from node ndx_from to node ind_to, one particle move
// ndx_from    -  (absolute) index of the old node
// ndx_to     -  (absolute) index of the new  node
// type  -  particle type index ( 1 to 255 )
// COmpMode: true: transport of DCs; false - transport of ICs 
// tcode  -  particle transport mechanism code (see enum PTCODE)
// iips   - DataBr index of phase or species to which this particle is connected
// m_v -  mass or volume of the particle (depending on ptype and mmode)
void TNodeArray::MoveParticleMass( int ndx_from, int ndx_to,
       char /*type*/, char CompMode, char tcode, unsigned char iips, double m_v )
{
   double mass = 0., coeff, mol, mWat, fmolal, aji;
   DATABR* dbr = NodT1[ndx_from];
   DATACH* dch = CSD;
   int xWatCH=0, ic, ips = (int)iips;
   if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
   {
	   xWatCH = CSD->nDCinPH[CSD->xPH[0]]-1; // CH index of water
//	   mWat = node1_xDC( ndx_from, xWatCH )* CSD->DCmm[xWatCH]; 
	   mWat = node1_xPA(ndx_from, ips) * CSD->DCmm[xWatCH];  // Mass of water-solvent
	   fmolal = 1.0; // 1000./mWat;              // molality conversion factor
   }

   switch( tcode )
   {
    case DISSOLVED: // mass of dissolved matter in aqueous solution
    				mass = mWat;  // trying normalization over mass of water-solvent
//    				mass = 1000/fmolal; // grams of water in the node
//    				fmolal = 1.0;
    				//    	mass = node1_mPS(ndx_from,ips) - mWat;
                   break;
    case ADVECTIVE: // mass of liquid phase for full advection
                   mass = node1_mPS( ndx_from, ips );
                    //  - node1_xPA(ndx_from,ips) * dch->DCmm[xWatCH];
                   break;
    case COLLOID:  // colloid particles
                   if( ips < dch->nPSb )
                       mass = node1_mPS( ndx_from, ips );
                   else
                       mass = node1_mPH( ndx_from, ips );
                   break;
    case DIFFUSIVE: // gets the mass of diffusing species in the phase
                   mass = nodeCH_DCmm( ips ) * node1_xDC( ndx_from, ips );
                   break;
    }
   coeff = m_v/mass; // mass of particle/mass of phase (solvent). Is this reasonable? 

  if( CompMode == true )
  { // Moving dependent components 
	for(short jc=0; jc < CSD->nDC; jc++ )
	{
	  mol = 0.; // moles of DC transported in the particle
	  switch( tcode )
	  {
        case DISSOLVED: // moving only dissolved DC (-H2O)
//             if( jc == xWatCH )
//            	 continue;  // H2O is ignored - not moved with the particle 
	    case ADVECTIVE: // moving DC of the whole aq phase
             if( jc > xWatCH )
            	 continue;     // ignoring non-aqueous species 
	    	 mol = node1_xDC( ndx_from, jc ) * coeff * fmolal;
             break;                	 
//             ( node1_bPS( ndx_from, ips, ie )
//                   - nodeCH_A( xWatCH, ie)
//                   * node1_xPA(ndx_from,ips)) * coeff;
	     case COLLOID:  // moving DC of solid particle - to be completed!     
// 	    	 if( ips < dch->nPSb )
//                  mol = node1_bPS( ndx_from, ips, ie ) * coeff;
//             else
//                  mol = node1_bPH( ndx_from, ips, ie ) * coeff;
	    	 break;
	     case DIFFUSIVE: // moving DC - a diffusing species
             if( jc != ips )
            	 continue;     // ignoring other diffusive species 
	    	 mol = node1_xDC( ndx_from, jc ) * coeff * fmolal;
             break;                	 
	  }
	  if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
	      mol /= fmolal;       // back from molality to moles
      
      if( fabs(mol) > 1e-20 ) // mtp->cdv ) // Threshold for DC change carried over in the particle  
      { 	   
    	  if( NodT1[ndx_from]->NodeTypeHY != NBC3source )
    	  {
    		  node1_xDC( ndx_from, jc ) -= mol;   // Correcting species amount in source node at T1 
    		  for( ic=0; ic<CSD->nICb; ic++)  // incrementing independent components
    		  {
    			  aji = nodeCH_A( jc, ic );
    			  if( aji )
    				  node1_bIC(ndx_from, ic) -= aji * mol;
    		  }
    	  }
 	  	  if( ndx_to >= 0 && ndx_to < anNodes )
 	  	  {
 	  		  if( NodT1[ndx_to]->NodeTypeHY != NBC3source )
 	  		  {
 	  			  node1_xDC( ndx_to, jc ) += mol;
 	  			  for( ic=0; ic<CSD->nICb; ic++)  // incrementing independent components
 	  			  {
 	  				  aji = nodeCH_A( jc, ic );
 	  				  if( aji )
 	  					  node1_bIC(ndx_to, ic) += aji * mol;
 	  			  }	
 	  		  }	 
 	  	  }
 	  	  else
 	  		  if(dbr->NodeTypeHY != NBC3sink  && dbr->NodeTypeHY != NBC3source)
 	  			cout << "W002MTRW " << "Warning: Particle jumped outside the domain" << endl; 
// 	  			  Error( "W002MTRW", "Warning: Particle jumped outside the domain" );
//	   } 
     } 	  
	} // loop jc  	  
  }
  else {  
	      // Transport of independent components 
   for(short ie=0; ie < CSD->nICb; ie++ )
   {
     mol = 0.; // moles of IC in the particle
     switch( tcode )
     {
        case DISSOLVED: // moving only dissolved IC (-H2O)
                        mol = ( node1_bPS( ndx_from, ips, ie )
//                              - nodeCH_A( xWatCH, ie)* node1_xPA(ndx_from,ips)
                              ) * coeff * fmolal;
                        break;
        case ADVECTIVE: // moving IC of the whole aq phase
                        mol = ( node1_bPS( ndx_from, ips, ie )
                        // - nodeCH_A( xWatCH, ie) * node1_xPA(ndx_from,ips)
                        ) * coeff * fmolal;
                        break;
        case COLLOID:   // moving IC of solid particle
                        if( ips < dch->nPSb )
                             mol = node1_bPS( ndx_from, ips, ie ) * coeff;
                        else
                             mol = node1_bPH( ndx_from, ips, ie ) * coeff;
                        break;
        case DIFFUSIVE: // moving IC of diffusing species
                        mol = node1_xDC( ndx_from, ips ) * nodeCH_A( ips, ie )
                                * coeff * fmolal;
                        break;
     }
     if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
    	 mol /= fmolal; 
     if( dbr->NodeTypeHY != NBC3source )
       dbr->bIC[ie] -= mol;
     if( ndx_to >= 0 && ndx_to < anNodes )
     {
         if( NodT1[ndx_to]->NodeTypeHY != NBC3source )
           NodT1[ndx_to]->bIC[ie] += mol;
     }
     else
         if(dbr->NodeTypeHY != NBC3sink  && dbr->NodeTypeHY != NBC3source)
        	 cout << "W002MTRW " << "Warning: Particle jumped outside the domain" << endl;        	                     
//        	 Error( "W002MTRW", "Warning: Particle jumped outside the domain" );

   } // loop ie 
  } // else
   // End of function
}

//==========================================================================
// Data collection for monitoring differences

// Prints difference increments in a all nodes (cells) for time point t / at
void TNodeArray::logDiffsIC( FILE* diffile, int t, double at, int nx, int every_t )
{
  double dc;
  int i, ie;

  if( t % every_t )
    return;

  fprintf( diffile, "\nStep= %-8d  Time= %-12.4g\nNode#   ", t, at );
  for( ie=0; ie < int(pCSD()->nICb); ie++ )
    fprintf( diffile, "%-12.4s ", pCSD()->ICNL[ie] );
  for (i=0; i<nx; i++)    // node iteration
  {
     fprintf( diffile, "\n%5d   ", i );
     for( ie=0; ie < int(pCSD()->nICb); ie++ )
     {
        dc = NodT1[i]->bIC[ie] - NodT0[i]->bIC[ie];
        fprintf( diffile, "%-12.4g ", dc );
     }
  }
  fprintf( diffile, "\n" );
}

// Data collection for monitoring 1D profiles in debugging FMT models
// Prints dissolved elemental molarities in all cells for time point t / at
void TNodeArray::logProfileAqIC( FILE* logfile, int t, double at, int nx, int every_t )
{
  double pm;
  int i, ie;
  if( t % every_t )
    return;
  fprintf( logfile, "\nStep= %-8d  Time= %-12.4g     Dissolved IC total concentrations, M\n", t, at/(365*86400) );
  fprintf(logfile, "%s","Node#   ");
  for( ie=0; ie < int(pCSD()->nICb); ie++ )
    fprintf( logfile, "%-12.4s ", pCSD()->ICNL[ie] );
  for (i=0; i<nx; i++)    // node iteration
  {
     fprintf( logfile, "\n%5d   ", i );
     for( ie=0; ie < int(pCSD()->nICb); ie++ )
     {
       pm = NodT1[i]->bPS[ie]/NodT1[i]->vPS[0]*1000.;  // Assumes there is aq phase!
                 // total dissolved element molarity
       fprintf( logfile, "%-12.4g ", pm );
     }
  }
  fprintf( logfile, "\n" );
}

// Data collection for monitoring 1D profiles
// Prints total elemental amounts in all cells for time point t / at
void TNodeArray::logProfileTotIC( FILE* logfile, int t, double at, int nx, int every_t )
{
  double pm;
  int i, ie;
  if( t % every_t )
    return;
  fprintf( logfile, "\nStep= %-8d  Time= %-12.4g     Bulk IC amounts, moles\n", t, at/(365*86400) );
  fprintf(logfile, "%s","Node#   ");
  for( ie=0; ie < int(pCSD()->nICb); ie++ )
    fprintf( logfile, "%-12.4s ", pCSD()->ICNL[ie] );
  for (i=0; i<nx; i++)    // node iteration
  {
     fprintf( logfile, "\n%5d   ", i );
     for( ie=0; ie < int(pCSD()->nICb); ie++ )
     {
       pm = NodT1[i]->bIC[ie];
       fprintf( logfile, "%-12.4g ", pm );
     }
  }
  fprintf( logfile, "\n" );
}

// Prints amounts of reactive phases in all cells for time point t / at
void TNodeArray::logProfilePhMol( FILE* logfile, int t, double at, int nx, int every_t )
{
  double pm;
  int i, ip;
  if( t % every_t )
    return;
  fprintf( logfile, "\nStep= %-8d  Time= %-12.4g     Amounts of reactive phases, moles\n", t, at/(365*86400) );
  fprintf(logfile, "%s","Node#   ");
  for( ip=0; ip < int(pCSD()->nPHb); ip++ )
    fprintf( logfile, "%-12.12s ", pCSD()->PHNL[ip] );
  for (i=0; i<nx; i++)    // node iteration
  {
     fprintf( logfile, "\n%5d   ", i );
     for( ip=0; ip < int(pCSD()->nPHb); ip++ )
     {
       pm = NodT1[i]->xPH[ip];
       fprintf( logfile, "%-12.4g ", pm );
     }
#ifndef NOPARTICLEARRAY
    if( TParticleArray::pa )
     {       // printing number of particles in nodes
       TParticleArray * ppa = TParticleArray::pa;
       int npa;
       for( int jp=0; jp < ppa->nPTypes(); jp++ )
       {
            npa = ppa->getNPnum( i, jp);   // number of particles in the node
            fprintf( logfile, "%-8d ", npa );
       }
     }
#endif
  }
  fprintf( logfile, "\n" );
}

// Prints dissolved species molarities in all cells for time point t / at
void TNodeArray::logProfileAqDC( FILE* logfile, int t, double at, int nx, int every_t )
{
  double pm;
  int i, is;
  if( t % every_t )
    return;
  fprintf( logfile, "\nStep= %-8d  Time= %-12.4g     Dissolved species concentrations, M\n",
        t, at/(365*86400) );
  fprintf(logfile, "%s","Node#   ");
  for( is=0; is < int(pCSD()->nDCb); is++ )
    fprintf( logfile, "%-12.4s ", pCSD()->DCNL[is] );
  for (i=0; i<nx; i++)    // node iteration
  {
     fprintf( logfile, "\n%5d   ", i );
     for( is=0; is < int(pCSD()->nDCinPH[0]); is++ )
     {
       pm = NodT1[i]->xDC[is]/NodT1[i]->vPS[0]*1000.;  // Assumes there is aq phase!
                 // dissolved species molarity
       fprintf( logfile, "%-12.4g ", pm );
     }
  }
  fprintf( logfile, "\n" );
}

//-----------------------End of nodearray.cpp--------------------------

