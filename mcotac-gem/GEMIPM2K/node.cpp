//-------------------------------------------------------------------
// $Id: node.cpp 684 2005-11-23 13:17:15Z gems $
//
// Implementation of Tnode class including initialization and
// execution of GEMIPM2 kernel

// Works whith DATACH and DATABR structures
//
// Copyright (C) 2005-2007 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#include "node.h"
#include "gdatastream.h"
#include "num_methods.h"
#include <math.h>
#include <algorithm>

#ifndef __unix
#include <io.h>
#endif

#ifndef IPMGEMPLUGIN
  #include "service.h"
  #include "visor.h"
#else
  istream& f_getline(istream& is, gstring& str, char delim);
#endif

TNode* TNode::na;

// This function proves whether the given Tc and P values fit within
// the interpolation regions.
// Result is returned in the ok field in this
// TNode class instance (true if T and P fit, otherwise false)
bool  TNode::check_TP( double& Tc, double& P )
{
   bool ok = true;
   double T_ = Tc, P_ = P;

   if( Tc < CSD->TCval[0] )
   { ok = false;
     Tc = CSD->TCval[0];
   }
   if( Tc > CSD->TCval[CSD->nTp-1] )
   { ok = false;
     Tc = CSD->TCval[CSD->nTp-1];
   }

  if( !ok )
  {
    fstream f_log("ipmlog.txt", ios::out|ios::app );
    f_log << "In node "<< CNode->NodeHandle << "  Given T= "<<  T_ <<
             "  is beyond the range for thermodynamic data;" <<
             " set to T= " << Tc << endl;
  }

  ok = true;
  if( P < CSD->Pval[0] )
  { ok = false;
    P = CSD->Pval[0];
  }
  if( P > CSD->Pval[CSD->nPp-1] )
  { ok = false;
    P = CSD->Pval[CSD->nPp-1];
  }

  if( !ok )
  {
    fstream f_log("ipmlog.txt", ios::out|ios::app );
    f_log << "In node "<< CNode->NodeHandle << "  Given P= "<<  P_ <<
           "  is beyond the range for thermodynamic data;" <<
           " set to P= " << P << endl;
  }
  return ok;
}

//-------------------------------------------------------------------------
// GEM_run()
// GEM IPM calculation of equilibrium state for the work node.
//   mode of GEM calculation is taken from the DATABR work structure
//
//  if uPrimalSol == true then the primal solution (vectors x, gamma, IC etc.) 
//    will be unpacked before GEMIPM calculation - as an option for the PIA mode 
//    with previous GEM solution taken from the same node. 
//  If uPrimalSol == false then the primal solution data will not be unpacked 
//    into the MULTI structure (AIA mode or PIA mode with primal solution retained 
//    in the MULTI structure from any previous IPM calculation)
//
//   Function returns: NodeStatusCH code from DATABR structure
//   ( OK; GEMIPM2K calculation error; system error )
//
int TNode::GEM_run( bool uPrimalSol )
{
//  fstream f_log("ipmlog.txt", ios::out|ios::app );
  CalcTime = 0.0;
  PrecLoops = 0; NumIterFIA = 0; NumIterIPM = 0;
//
  try
  {
// f_log << " GEM_run() begin Mode= " << p_NodeStatusCH endl;
//---------------------------------------------
// Checking T and P  for interpolation intervals
   check_TP( CNode->TC, CNode->P);
// Unpacking work DATABR structure into MULTI (GEM IPM structure): uses DATACH
// setting up up PIA or AIA mode
   if( CNode->NodeStatusCH == NEED_GEM_PIA )
   {
	   pmm->pNP = 1;
	   unpackDataBr( uPrimalSol );
   }
   else { 
	   pmm->pNP = 0; // As default setting AIA mode
	   unpackDataBr( false );
   }
   // GEM IPM calculation of equilibrium state in MULTI
   CalcTime = TProfil::pm->calcMulti( PrecLoops, NumIterFIA, NumIterIPM );
// Extracting and packing GEM IPM results into work DATABR structure
    packDataBr();
    CNode->IterDone = NumIterFIA+NumIterIPM;
//**************************************************************
// only for testing output results for files
//    GEM_write_dbr( "calculated_dbr.dat",  false );
//    GEM_printf( "calc_multi.ipm", "calculated_dbr.dat", "calculated.dbr" );
// *********************************************************

   if( CNode->NodeStatusCH  == NEED_GEM_AIA )
         CNode->NodeStatusCH = OK_GEM_AIA;
   else
         CNode->NodeStatusCH = OK_GEM_PIA;

   }
   catch(TError& err)
    {
     fstream f_log("ipmlog.txt", ios::out|ios::app );
     f_log << err.title.c_str() << ": " << err.mess.c_str() << endl;
     if( CNode->NodeStatusCH  == NEED_GEM_AIA )
       CNode->NodeStatusCH = BAD_GEM_AIA;
     else
       CNode->NodeStatusCH = BAD_GEM_PIA;

    }
    catch(...)
    {
     fstream f_log("ipmlog.txt", ios::out|ios::app );
     f_log << "gems2: Unknown exception: GEM calculation aborted" << endl;
       if( CNode->NodeStatusCH  == NEED_GEM_AIA )
         CNode->NodeStatusCH = ERR_GEM_AIA;
       else
         CNode->NodeStatusCH = ERR_GEM_PIA;
    }
   return CNode->NodeStatusCH;
}

//-----------------------------------------------------------------------
// Returns calculation time after the last GEM_run() call
//
double TNode::GEM_CalcTime()
{
  return CalcTime;
}

// Returns total number of FIA + IPM iterations after the last call to GEM_run()
// More detailed info is returned via parameters by reference:
//    PrecLoops:  Number of performed IPM-2 precision refinement loops
//    NumIterFIA: Total Number of performed FIA entry iterations
//    NumIterIPM: Total Number of performed IPM main iterations
int TNode::GEM_Iterations( int& PrecLoops_, int& NumIterFIA_, int& NumIterIPM_ )
{
	PrecLoops_ = PrecLoops;
	NumIterFIA_ = NumIterFIA;
	NumIterIPM_ = NumIterIPM;
	return NumIterFIA+NumIterIPM;
}

// ----------------------------------------------------------------------
// reads work node (DATABR structure) from a  file
int  TNode::GEM_read_dbr( const char* fname, bool binary_f )
{
  try
  {
    if( binary_f )
	{
       gstring str_file = fname;
	   GemDataStream in_br(str_file, ios::in|ios::binary);
       databr_from_file(in_br);
	}
   else
   {   fstream in_br(fname, ios::in );
       ErrorIf( !in_br.good() , fname, "DataBR Fileopen error");
       databr_from_text_file(in_br);
   }

    dbr_file_name = fname;

  } catch(TError& /*err*/)
    {
      return 1;
    }
    catch(...)
    {
      return -1;
    }
  return 0;
}

//-------------------------------------------------------------------
// GEM_init()
// reads in the data from IPM_DAT, DCH_DAT, DBR_DAT files (prepared by hand or
// using the GEMS-PSI GEM2MT module).
//  Parameters:
//  ipmfiles_lst_name - name of a text file that contains:
//    " -t/-b <DCH_DAT file name> <IPM_DAT file name> <dataBR file name1>,
//      ... , <dataBR file nameN> "
//    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
//    exist in the same directory where the ipmfiles_lst_name file is located.
//    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
//    and must contain valid initial chemical systems (of the same structure
//    as described in the DCH_DAT file) to set up the initial state of the FMT
//    node array. If -t flag or nothing is specified then all data files must
//    be in text (ASCII) format; if -b flag is specified then all data files
//    are  assumed to be binary (little-endian) files.
//  nodeTypes[nNodes] - array of node type (fortran) indexes of DBR_DAT files
//    in the ipmfiles_lst_name list. This array (handle for each FMT node),
//    specifies from which DBR_DAT file the initial chemical system should
//    be taken.
//  This function returns:
//   0: OK; 1: GEM IPM read file error; -1: System error (e.g. memory allocation)
//
//-------------------------------------------------------------------

int  TNode::GEM_init( const char* ipmfiles_lst_name,
                          int* nodeTypes, bool getNodT1)
{
  int i;
#ifdef IPMGEMPLUGIN
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
    {
#else
      size_t npos = gstring::npos;
#endif
//     bool binary_mult = true;
     bool binary_f = false;
     gstring lst_in = ipmfiles_lst_name;

// Get path
      size_t pos = lst_in.rfind("/");
      gstring Path = "";
      if( pos < npos )
       Path = lst_in.substr(0, pos+1);

//  open file stream for the file names list file
      fstream f_lst( lst_in.c_str(), ios::in );
      ErrorIf( !f_lst.good() , lst_in.c_str(), "Fileopen error");

      gstring datachbr_fn;
      f_getline( f_lst, datachbr_fn, ' ');

//  Syntax: -t/-b  "<DCH_DAT file name>"  "<IPM_DAT file name>"
//       "<DBR_DAT file1 name>" [, ... , "<DBR_DAT fileN name>"]
// Rearranged in logical shape by KD on 12.01.2007

//Testing flag "-t" or "-b" (by default "-t")   // use binary or text files
      pos = datachbr_fn.find( '-');
      if( pos != /*gstring::*/npos )
      {
         if( datachbr_fn[pos+1] == 'b' )
            binary_f = true;
         f_getline( f_lst, datachbr_fn, ' ');
      }

      // Reading name of DCH_DAT file
      gstring dat_ch = Path + datachbr_fn;

      // Reading name of IPM_DAT file for structure MULTI (GEM IPM work structure)
      f_getline( f_lst, datachbr_fn, ' ');
      gstring mult_in = Path + datachbr_fn;

// Reading DCH_DAT file in binary or text format
      if( binary_f )
      {  GemDataStream f_ch(dat_ch, ios::in|ios::binary);
         datach_from_file(f_ch);
       }
      else
      { fstream f_ch(dat_ch.c_str(), ios::in );
         ErrorIf( !f_ch.good() , dat_ch.c_str(), "DCH_DAT fileopen error");
         datach_from_text_file(f_ch);
      }

// Reading IPM_DAT file into structure MULTI (GEM IPM work structure)
if( binary_f )
 {
   GemDataStream f_m(mult_in, ios::in|ios::binary);
#ifdef IPMGEMPLUGIN
    profil->readMulti(f_m);
#else
    TProfil::pm->readMulti(f_m);
#endif
  }
  else
  {
#ifdef IPMGEMPLUGIN
        profil->readMulti(mult_in.c_str());
#endif
  }

// Prepare for reading DBR_DAT files
     i = 0;
     while( !f_lst.eof() )  // For all DBR_DAT files listed
     {

#ifndef IPMGEMPLUGIN
   pVisor->Message( 0, "GEM2MT node array",
      "Reading from disk a set of node array files to resume an interrupted RMT task.\n"
           "Please, wait...", i, nNodes() );
#endif

// Reading DBR_DAT file into work DATABR structure
         if( i )  // Comma only after the first DBR_DAT file!
            f_getline( f_lst, datachbr_fn, ',');
         else
            f_getline( f_lst, datachbr_fn, ' ');

         gstring dbr_file = Path + datachbr_fn;
         if( binary_f )
         {
             GemDataStream in_br(dbr_file, ios::in|ios::binary);
             databr_from_file(in_br);
          }
         else
          {   fstream in_br(dbr_file.c_str(), ios::in );
                 ErrorIf( !in_br.good() , datachbr_fn.c_str(),
                    "DBR_DAT fileopen error");
               databr_from_text_file(in_br);
          }
          if(!i) 
        	  dbr_file_name = dbr_file;
// Unpacking work DATABR structure into MULTI (GEM IPM work structure): uses DATACH
//    unpackDataBr();

#ifndef IPMGEMPLUGIN
        if( getNodT1 )
        {
           setNodeArray( dbr_file, i, binary_f );
        }
        else
#endif
        {
// Copying data from work DATABR structure into the node array
// (as specified in nodeTypes array)
           setNodeArray( i, nodeTypes  );
         }
          i++;
     }  // end while()
#ifndef IPMGEMPLUGIN
   pVisor->CloseMessage();
#endif

    ErrorIf( i==0, datachbr_fn.c_str(), "GEM_init() error: No DBR_DAT files read!" );
    checkNodeArray( i, nodeTypes, datachbr_fn.c_str()  );

   return 0;

#ifdef IPMGEMPLUGIN
    }
    catch(TError& err)
    {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
    }
    catch(...)
    {
        return -1;
    }
    return 1;
#endif
}

//-----------------------------------------------------------------
// work with lists

// Return ICH index of IC by Name or -1 if name not found
int TNode::IC_name_to_x( const char *Name )
{
  uint len = strlen( Name );
  len =  min(len,MaxICN);

  for(int ii = 0; ii<CSD->nIC; ii++ )
       if(!memcmp(Name, CSD->ICNL[ii], len ))
        if( len == MaxICN || CSD->ICNL[ii][len] == ' ' || CSD->ICNL[ii][len] == '\0' )
         return ii;
  return -1;
}

// Return DCH index of DC by Name or -1 if name not found
int TNode::DC_name_to_x( const char *Name )
{
  uint len = strlen( Name );
  len =  min(len,MaxDCN);

  for(int ii = 0; ii<CSD->nDC; ii++ )
       if(!memcmp(Name, CSD->DCNL[ii], min(len,MaxDCN)))
        if( len == MaxDCN || CSD->DCNL[ii][len] == ' ' || CSD->DCNL[ii][len] == '\0' )
         return ii;
  return -1;
}

// Return DCH index of Ph by Name or -1 if name not found
int TNode::Ph_name_to_x( const char *Name )
{
  uint len = strlen( Name );
  len =  min(len,MaxPHN);

  for(int ii = 0; ii<CSD->nPH; ii++ )
       if(!memcmp(Name, CSD->PHNL[ii], min(len,MaxPHN)))
        if( len == MaxPHN || CSD->PHNL[ii][len] == ' ' || CSD->PHNL[ii][len] == '\0' )
         return ii;
  return -1;
}

// Return for IComp DBR index from DCH index or -1 if IComp is not used in the data bridge
int TNode::IC_xCH_to_xDB( const int xCH )
{
  for(int ii = 0; ii<CSD->nICb; ii++ )
       if( CSD->xIC[ii] == xCH )
         return ii;
  return -1;
}

// Return for DComp DBR index from DCH index or -1 if not used in the data bridge
int TNode::DC_xCH_to_xDB( const int xCH )
{
  for(int ii = 0; ii<CSD->nDCb; ii++ )
       if( CSD->xDC[ii] == xCH )
         return ii;
  return -1;
}

// Return for Ph DBR index from DCH index or -1 if not used in the data bridge
int TNode::Ph_xCH_to_xDB( const int xCH )
{
  for(int ii = 0; ii<CSD->nPHb; ii++ )
       if( CSD->xPH[ii] == xCH )
         return ii;
  return -1;
}

// Converts the Phase DCH index into the DC DCH index
 int  TNode::Phx_to_DCx( const int Phx )
 {
   int k, DCx = 0;
   for( k=0; k<CSD->nPHb; k++ )
   {
     if( k == Phx )
      break;
     DCx += CSD->nDCinPH[ k];
   }
   return DCx;
 }

// Converts the Phase DCH index into the DC DCH index (1-st)
// returns into nDCinPh number of DC included into Phx phase
 int  TNode::PhtoDC_DCH( const int Phx, int& nDCinPh )
 {
   int k, DCx = 0;
   for( k=0; k<CSD->nPHb; k++ )
   {
     if( k == Phx )
      break;
     DCx += CSD->nDCinPH[ k];
   }
   nDCinPh = CSD->nDCinPH[k];
   return DCx;
 }

// Converts the Phase DBR index into the DC DBR index (1-st selected )
// returns into nDCinPh number of DC selected into Phx phase
 int  TNode::PhtoDC_DBR( const int Phx, int& nDCinPh )
 {
   int DCx, DCxCH, PhxCH, nDCinPhCH;

   PhxCH = Ph_xBR_to_xCH( Phx );
   DCxCH = PhtoDC_DCH( PhxCH, nDCinPhCH );

   DCx = -1;
   nDCinPh = 0;
   for(int ii = 0; ii<CSD->nDCb; ii++ )
   {
      if( CSD->xDC[ii] >= DCxCH )
      {
        if( CSD->xDC[ii] >= DCxCH+nDCinPhCH  )
          break;
        nDCinPh++;
        if( DCx == -1)
          DCx = ii;
      }
   }
   return DCx;
 }

 // Test Tc as grid point for the interpolation of thermodynamic data
 // Return index in grid array or -1
 int  TNode::check_grid_T( double& Tc )
 {
   int jj;
   for( jj=0; jj<CSD->nTp; jj++)
    if( fabs( Tc - CSD->TCval[jj] ) < CSD->Ttol )
    {
       Tc = CSD->TCval[jj];
       return jj;
    }
   return -1;
 }

 // Test P as grid point for the interpolation of thermodynamic data
 // Return index in grid array or -1
 int  TNode::check_grid_P( double& P )
 {
   int jj;
   for( jj=0; jj<CSD->nPp; jj++)
    if( fabs( P - CSD->Pval[jj] ) < CSD->Ptol )
    {
      P = CSD->Pval[jj];
      return jj;
    }
   return -1;
 }

  // Test Tc and P as grid point for the interpolation of thermodynamic data
 // Return index in grid matrix or -1
  int  TNode::check_grid_TP(  double& Tc, double& P )
  {
    int xT, xP, ndx=-1;

    xT = check_grid_T( Tc );
    xP = check_grid_P( P );
    if( xT >=0 && xP>= 0 )
     ndx =  xP * CSD->nTp + xT;
    return ndx;
  }

 // Access to interpolated G0 from DCH structure ( xCH the DC DCH index)
  double  TNode::DC_G0_TP( const int xCH, double& Tc, double& P )
  {
    int xTP, jj;
    double G0;

    check_TP( Tc, P );
    xTP = check_grid_TP( Tc, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       G0 = CSD->G0[ jj + xTP ];
    else
       G0 = LagranInterp( CSD->Pval, CSD->TCval, CSD->G0+jj,
                        P, Tc, CSD->nTp, CSD->nPp, 1 );
    return G0;
 }

  // Access to interpolated V0 from DCH structure ( xCH the DC DCH index)
  double  TNode::DC_V0_TP( const int xCH, double& Tc, double& P )
  {
    int xTP, jj;
    double V0;

    check_TP( Tc, P );
    xTP = check_grid_TP( Tc, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       V0 = CSD->V0[ jj + xTP ];
    else
       V0 = LagranInterp( CSD->Pval, CSD->TCval, CSD->V0+jj,
                        P, Tc, CSD->nTp, CSD->nPp, 1 );
    return V0;
}

 // Retrieval of Phase Volume ( xBR the Ph DBR index)
 double  TNode::Ph_Volume( const int xBR )
 {
   double vol;
   if( xBR < CSD->nPSb )
    vol = CNode->vPS[xBR];
   else
   {
     int xDC = Phx_to_DCx( Ph_xBR_to_xCH( xBR ));
     vol = DC_V0_TP( xDC, CNode->TC, CNode->P );
     vol *= CNode->xDC[DC_xCH_to_xDB(xDC)] *10.;
   }
   return vol;
 }

  // Retrieval of Phase mass ( xBR the Ph DBR index)
  double  TNode::Ph_Mass( const int xBR )
  {
     double mass;
     if( xBR < CSD->nPSb )
        mass = CNode->mPS[xBR];
     else
     {
        int xDC = Phx_to_DCx( Ph_xBR_to_xCH( xBR ));
        mass = CNode->xDC[ DC_xCH_to_xDB(xDC) ] * CSD->DCmm[xDC];
     }
    return mass;
  }

  // Retrieval of Phase composition ( xBR the Ph DBR index)
  void  TNode::Ph_BC( const int xBR, double* ARout )
  {
    int ii;
    if( !ARout )
      ARout = new double[ CSD->nICb ];

    if( xBR < CSD->nPSb )
       for( ii=0; ii<pCSD()->nICb; ii++ )
         ARout[ii] = CNode->bPS[ xBR * CSD->nICb + ii ];
    else
    {
      int DCx = Phx_to_DCx( Ph_xBR_to_xCH(xBR) );
      for( ii=0; ii<pCSD()->nICb; ii++ )
      {
         ARout[ii] = CSD->A[ IC_xBR_to_xCH(ii) + DCx * CSD->nIC];
         ARout[ii] *= CNode->xDC[ DC_xCH_to_xDB(DCx) ];
      }
    }
  }

//---------------------------------------------------------//

void TNode::allocMemory()
{
// memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;

    memset( CSD, 0, sizeof(DATACH) );
    memset( CNode, 0, sizeof(DATABR) );

#ifdef IPMGEMPLUGIN
// internal structures
    multi = new TMulti();
    pmm = multi->GetPM();
    profil = new TProfil( multi );
    TProfil::pm = profil;
#endif
}

void TNode::freeMemory()
{
   datach_free();
   CSD = 0;
   CNode = databr_free( CNode );

#ifdef IPMGEMPLUGIN
  delete multi;
  delete profil;
#endif
}

#ifndef IPMGEMPLUGIN

// Makes start DATACH and DATABR data from GEMS internal data (MULTI and other)
void TNode::MakeNodeStructures(
    short anICb,       // number of stoichiometry units (<= nIC) used in the data bridge
    short anDCb,      	// number of DC (chemical species, <= nDC) used in the data bridge
    short anPHb,     	// number of phases (<= nPH) used in the data bridge
    short* axIC,   // ICNL indices in DATABR IC vectors [nICb]
    short* axDC,   // DCNL indices in DATABR DC list [nDCb]
    short* axPH,   // PHNL indices in DATABR phase vectors [nPHb]
    float* Tai, float* Pai,
    short nTp_, short nPp_, float Ttol_, float Ptol_  )
{
  short ii;
  TCIntArray aSelIC;
  TCIntArray aSelDC;
  TCIntArray aSelPH;

// make lists
  for( ii=0; ii<anICb; ii++)
     aSelIC.Add( axIC[ii] );
  for( ii=0; ii<anDCb; ii++)
     aSelDC.Add( axDC[ii] );
  for( ii=0; ii<anPHb; ii++)
     aSelPH.Add( axPH[ii] );

// set default data and realloc arrays
   makeStartDataChBR( aSelIC, aSelDC, aSelPH,
                      nTp_, nPp_, Ttol_, Ptol_, Tai, Pai );
}

// Make start DATACH and DATABR data from GEMS internal data (MULTI and other)
void TNode::MakeNodeStructures( QWidget* par, bool select_all,
    float *Tai, float *Pai,
    short nTp_, short nPp_, float Ttol_, float Ptol_  )
{

  //MULTI *mult = TProfil::pm->pmp;
  TCStringArray aList;
  TCIntArray aSelIC;
  TCIntArray aSelDC;
  TCIntArray aSelPH;

// select lists
    aList.Clear();
    for(int ii=0; ii< pmm->N; ii++ )
    {  if( select_all )
         aSelIC.Add( ii );
       else
         aList.Add( gstring( pmm->SB[ii], 0, MAXICNAME+MAXSYMB));
    }
    if( !select_all  )
      aSelIC = vfMultiChoice(par, aList,
          "Please, mark independent components for selection into DataBridge");

    aList.Clear();
    for(int ii=0; ii< pmm->L; ii++ )
    {  if( select_all )
         aSelDC.Add( ii );
       else
       aList.Add( gstring( pmm->SM[ii], 0, MAXDCNAME));
    }
    if( !select_all  )
       aSelDC = vfMultiChoice(par, aList,
         "Please, mark dependent components for selection into DataBridge");

    aList.Clear();
    for(int ii=0; ii< pmm->FI; ii++ )
    {  if( select_all )
         aSelPH.Add( ii );
       else
       aList.Add( gstring( pmm->SF[ii], 0, MAXPHNAME+MAXSYMB));
    }
    if( !select_all  )
       aSelPH = vfMultiChoice(par, aList,
         "Please, mark phases for selection into DataBridge");


// set default data and realloc arrays
   makeStartDataChBR( aSelIC, aSelDC, aSelPH,
                      nTp_, nPp_, Ttol_, Ptol_, Tai, Pai );
}


// Writing dataCH structure to binary file
void TNode::makeStartDataChBR(
  TCIntArray& selIC, TCIntArray& selDC, TCIntArray& selPH,
  short nTp_, short nPp_, float Ttol_, float Ptol_,
  float *Tai, float *Pai )
{
// set sizes for DataCh
  uint ii;
  short i1;
// reallocates memory for     DATACH  *CSD;  and  DATABR  *CNode;
  if( !CSD )
     CSD = new DATACH;
  if( !CNode )
     CNode = new DATABR;

  CSD->nIC = pmm->N;
  CSD->nDC = pmm->L;
  CSD->nDCs = pmm->Ls;
  CSD->nPH = pmm->FI;
  CSD->nPS = pmm->FIs;
  CSD->nTp = nTp_;
  CSD->nPp = nPp_;
  if( pmm->Aalp )
    CSD->nAalp = 1;
  else
    CSD->nAalp = 0;
  CSD->iGrd = 0;
  if ( pmm->H0 )
    CSD->iGrd = 1;
  if ( pmm->S0 )
    CSD->iGrd = 2;
  if ( pmm->Cp0 )
    CSD->iGrd = 3;

// These dimensionalities define sizes of dynamic data in DATABR structure!!!

  CSD->nICb = (short)selIC.GetCount();
  CSD->nDCb = (short)selDC.GetCount();
  CSD->nPHb = (short)selPH.GetCount();
  CSD->nPSb = 0;
  for( ii=0; ii< selPH.GetCount(); ii++, CSD->nPSb++ )
   if( selPH[ii] >= pmm->FIs )
       break;
  CSD->uRes1 = 0;
  CSD->dRes1 = 0.;
  CSD->dRes2 = 0.;

  CSD->Ttol = Ttol_;
  CSD->Ptol = Ptol_;

// realloc structures DataCh&DataBr

  datach_realloc();
  databr_realloc();

// set dynamic data to DataCH

  for( ii=0; ii< selIC.GetCount(); ii++ )
    CSD->xIC[ii] = (short)selIC[ii];
  for( ii=0; ii< selDC.GetCount(); ii++ )
    CSD->xDC[ii] = (short)selDC[ii];
  for( ii=0; ii< selPH.GetCount(); ii++ )
    CSD->xPH[ii] = (short)selPH[ii];

  for( i1=0; i1< CSD->nIC*CSD->nDC; i1++ )
    CSD->A[i1] = pmm->A[i1];

  for( i1=0; i1< CSD->nPH; i1++ )
  {
    CSD->nDCinPH[i1] = pmm->L1[i1];
    CSD->ccPH[i1] = pmm->PHC[i1];
    memcpy( CSD->PHNL[i1], pmm->SF[i1]+4 , MaxPHN*sizeof(char));
  }
  for( i1=0; i1< CSD->nIC; i1++ )
  {
    CSD->ICmm[i1] = pmm->Awt[i1];
    CSD->ccIC[i1] = pmm->ICC[i1];
    memcpy( CSD->ICNL[i1], pmm->SB[i1] , MaxICN*sizeof(char));
  }
  for( i1=0; i1< CSD->nDC; i1++ )
  {
    CSD->DCmm[i1] = pmm->MM[i1];
    CSD->ccDC[i1] = pmm->DCC[i1];
    memcpy( CSD->DCNL[i1], pmm->SM[i1] , MaxDCN*sizeof(char));
  }
//  memset( CSD->DD, 0, CSD->nDCs*sizeof(double));

// set default data to DataBr

   CNode->NodeHandle = 0;
   CNode->NodeTypeHY = normal;
   CNode->NodeTypeMT = normal;
   CNode->NodeStatusFMT = Initial_RUN;
//   CNode->NodeStatusCH = NEED_GEM_AIA;
   if( pmm->pNP == 0 )
    CNode->NodeStatusCH = NEED_GEM_AIA;
  else
     CNode->NodeStatusCH = NEED_GEM_PIA;

   CNode->IterDone = 0;

   memset( &CNode->TC, 0, 32*sizeof(double));
   CNode->TC = pmm->TCc; //25
   CNode->P = pmm->Pc; //1
   CNode->Ms = pmm->MBX; // in kg

// arrays
   for( i1=0; i1<CSD->nICb; i1++ )
    CNode->bIC[i1] = pmm->B[ CSD->xIC[i1] ];

   for( i1=0; i1<CSD->nDCb; i1++ )
   {
     CNode->dul[i1] = pmm->DUL[ CSD->xDC[i1] ];
     CNode->dll[i1] = pmm->DLL[ CSD->xDC[i1] ];
    }

   if( CSD->nAalp >0 )
      for( i1=0; i1< CSD->nPHb; i1++ )
        CNode->aPH[i1] = pmm->Aalp[CSD->xPH[i1]];

// puts calculated & dynamic data to DataBR
   packDataBr();

// must be changed to matrix structure  ???????
// setted CSD->nPp*CSD->nTp = 1
   for( i1=0; i1<CSD->nTp; i1++ )
    CSD->TCval[i1] = Tai[i1];
   for( i1=0; i1<CSD->nPp; i1++ )
    CSD->Pval[i1] = Pai[i1];

   G0_V0_H0_Cp0_DD_arrays();

   if(  CSD->iGrd > 3 )
     for( i1=0; i1< CSD->nDCs*CSD->nPp*CSD->nTp; i1++ )
       CSD->DD[i1] = 0.;
}

void TNode::G0_V0_H0_Cp0_DD_arrays()
{

  double cT, cP;
  double *G0, *V0, *H0, *Cp0, *S0, roW, epsW;

  G0 =  new double[TProfil::pm->mup->L];
  V0 =  new double[TProfil::pm->mup->L];
  if ( pmm->H0 )
    H0 =  new double[TProfil::pm->mup->L];
  else
    H0 = 0;
  if ( pmm->Cp0 )
    Cp0 = new double[TProfil::pm->mup->L];
  else
    Cp0 = 0;
  if ( pmm->S0 )
      S0 = new double[TProfil::pm->mup->L];
  else
      S0 = 0;

  for( int ii=0; ii<CSD->nTp; ii++)
  {
    cT = CSD->TCval[ii];
    for( int jj=0; jj<CSD->nPp; jj++)
    {
      cP = CSD->Pval[jj];
     // calculates new G0, V0, H0, Cp0, S0
     TProfil::pm->LoadFromMtparm( cT, cP, G0, V0, H0, S0, Cp0, roW, epsW );
     CSD->roW[ jj * CSD->nTp + ii] = roW;
     CSD->epsW[ jj * CSD->nTp + ii] = epsW;
     // copy to arrays
     for(int kk=0; kk<CSD->nDC; kk++)
      {
         int ll = ( kk * CSD->nPp + jj) * CSD->nTp + ii;
         CSD->G0[ll] =  G0[pmm->muj[kk]]; //
         CSD->V0[ll] =  V0[pmm->muj[kk]];
         if(  CSD->iGrd > 0 )
         {  if ( H0 )
              CSD->H0[ll] = H0[pmm->muj[kk]];
            else
              CSD->H0[ll] = 0.;
         }
         if(  CSD->iGrd > 2 )
         { if ( Cp0 )
             CSD->Cp0[ll] = Cp0[pmm->muj[kk]];
           else
             CSD->Cp0[ll] = 0.;
         }
         if(  CSD->iGrd > 1 )
         {
            if ( S0 )
               CSD->S0[ll] = S0[pmm->muj[kk]];
            else
             CSD->S0[ll] = 0.;
         }
     }
    }
  }
  // free memory
  delete[] G0;
  delete[] V0;
  if( H0 )
   delete[] H0;
  if( Cp0 )
   delete[] Cp0;
  if( S0 )
   delete[] S0;
}

TNode::TNode( MULTI *apm  )
{
    pmm = apm;
    CSD = 0;
    CNode = 0;
    allocMemory();
    na = this;
    dbr_file_name = "dbr_file_name";
}

#else

TNode::TNode()
{
  CSD = 0;
  CNode = 0;
  allocMemory();
  na = this;
  dbr_file_name = "dbr_file_name";
}

#endif


TNode::~TNode()
{
   freeMemory();
   na = 0;
}

// Extracting and packing GEM IPM results into work DATABR structure
void TNode::packDataBr()
{
 short ii;

// set default data to DataBr
#ifndef IPMGEMPLUGIN
   CNode->NodeHandle = 0;
//   CNode->NodeTypeHY = normal;
   CNode->NodeTypeMT = normal;
   CNode->NodeStatusFMT = Initial_RUN;
   //   CNode->NodeStatusCH = NEED_GEM_AIA;
   if( pmm->pNP == 0 )
    CNode->NodeStatusCH = NEED_GEM_AIA;
  else
     CNode->NodeStatusCH = NEED_GEM_PIA;
//#else
//
 // numbers
//  if( pmm->pNP == 0 )
//    CNode->NodeStatusCH = OK_GEM_AIA;
//  else
//    CNode->NodeStatusCH = OK_GEM_PIA;
//
#endif

   CNode->TC = pmm->TCc; //25
   CNode->P = pmm->Pc; //1
//   CNode->IterDone = pmm->IT;
   CNode->IterDone = pmm->ITF+pmm->ITG;   // Now complete number of FIA and IPM iterations
// values
  CNode->Vs = pmm->VXc*1.e-6; // from cm3 to m3
  CNode->Gs = pmm->FX;
  CNode->Hs = pmm->HXc;
  CNode->IC = pmm->IC;
  CNode->pH = pmm->pH;
  CNode->pe = pmm->pe;
//  CNode->Eh = pmm->FitVar[3];  Bugfix 19.12.2006  KD
  CNode->Eh = pmm->Eh;
  CNode->Ms = pmm->MBX;

  // arrays
   for( ii=0; ii<CSD->nPHb; ii++ )
   {  CNode->xPH[ii] = pmm->XF[ CSD->xPH[ii] ];
      if( CSD->nAalp >0 )
       CNode->aPH[ii] = pmm->Aalp[ CSD->xPH[ii] ];//??? only insert
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
   {   CNode->vPS[ii] = pmm->FVOL[ CSD->xPH[ii] ];
       CNode->mPS[ii] = pmm->FWGT[ CSD->xPH[ii] ];
       CNode->xPA[ii] = pmm->XFA[ CSD->xPH[ii] ];
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
   for(short jj=0; jj<CSD->nICb; jj++ )
   { int   new_ndx= (ii*CSD->nICb)+jj,
           mul_ndx = ( CSD->xPH[ii]*CSD->nIC )+ CSD->xIC[jj];
     CNode->bPS[new_ndx] = pmm->BF[ mul_ndx ];
   }
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
      CNode->xDC[ii] = pmm->X[ CSD->xDC[ii] ];
      CNode->gam[ii] = pmm->Gamma[ CSD->xDC[ii] ];
      CNode->dul[ii] = pmm->DUL[ CSD->xDC[ii] ];//??? only insert
      CNode->dll[ii] = pmm->DLL[ CSD->xDC[ii] ];//??? only insert
   }
   for( ii=0; ii<CSD->nICb; ii++ )
   {  CNode->bIC[ii] = pmm->B[ CSD->xIC[ii] ];//??? only insert
      CNode->rMB[ii] = pmm->C[ CSD->xIC[ii] ];
      CNode->uIC[ii] = pmm->U[ CSD->xIC[ii] ];
   }
}

// Unpacking work DATABR structure into MULTI
//(GEM IPM work structure): uses DATACH
//  if uPrimalSol is true then the primal solution (vectors x, gamma, IC etc.) 
//  will be unpacked - as an option for PIA mode with previous GEM solution from 
//  the same node. 
//  If uPrimalSol = false then the primal solution data will not be unpacked 
//  into the MULTI structure (AIA mode or PIA mode with primal solution retained 
//    in the MULTI structure from previous IPM calculation)
void TNode::unpackDataBr( bool uPrimalSol )
{
 short ii;
 double Gamm;
// numbers

//  if( CNode->NodeStatusCH >= NEED_GEM_PIA )
//   pmm->pNP = 1;
//  else
//   pmm->pNP = 0; //  NEED_GEM_AIA;
  CNode->IterDone = 0;
  pmm->TCc = CNode->TC;
  pmm->Tc = CNode->TC+C_to_K;
  pmm->Pc  = CNode->P;
  // Obligatory arrays - always unpacked!
  for( ii=0; ii<CSD->nDCb; ii++ )
  {
    pmm->DUL[ CSD->xDC[ii] ] = CNode->dul[ii];
    pmm->DLL[ CSD->xDC[ii] ] = CNode->dll[ii];
  }
  for( ii=0; ii<CSD->nICb; ii++ )
    pmm->B[ CSD->xIC[ii] ] = CNode->bIC[ii];
  for( ii=0; ii<CSD->nPHb; ii++ )
  {
    if( CSD->nAalp >0 )
        pmm->Aalp[ CSD->xPH[ii] ] = CNode->aPH[ii];
  }
 
 if( !uPrimalSol )
 {    //  Using primal solution retained in the MULTI structure instead - 
     // the primal solution data from the DATABR structure are not unpacked
   pmm->IT = 0;	
 }
 else {   // Unpacking primal solution provided in the node DATABR structure 	
  pmm->MBX = CNode->Ms;
  pmm->IC = CNode->IC;
//  pmm->FitVar[3] = CNode->Eh;  Bugfix 19.12.2006  KD
  pmm->Eh = CNode->Eh;
  for( ii=0; ii<CSD->nDCb; ii++ )
  /*    pmm->X[ CSD->xDC[ii] ] = */
        pmm->Y[ CSD->xDC[ii] ] = CNode->xDC[ii];
  for( ii=0; ii<CSD->nDCb; ii++ )
  {
     pmm->lnGam[ CSD->xDC[ii] ] = log( CNode->gam[ii] );
  //       Gamm = CNode->gam[ii];
  //      pmm->Gamma[ CSD->xDC[ii] ] = Gamm;
  //      pmm->lnGmo[ CSD->xDC[ii] ] = pmm->lnGam[ CSD->xDC[ii] ] = log(Gamm);
  }
  for( ii=0; ii<CSD->nPSb; ii++ )
   pmm->FVOL[ CSD->xPH[ii] ] = CNode->vPS[ii];
  for( ii=0; ii<CSD->nPSb; ii++ )
   pmm->FWGT[ CSD->xPH[ii] ] = CNode->mPS[ii];

  for( ii=0; ii<CSD->nPHb; ii++ )
  {
    pmm->XF[ CSD->xPH[ii] ] =
    pmm->YF[ CSD->xPH[ii] ] = CNode->xPH[ii];
  } 
  
  for( int k=0; k<CSD->nPSb; k++ )
  for(int i=0; i<CSD->nICb; i++ )
  { int dbr_ndx= (k*CSD->nICb)+i,
          mul_ndx = ( CSD->xPH[k]*CSD->nIC )+ CSD->xIC[i];
    pmm->BF[ mul_ndx ] = CNode->bPS[dbr_ndx];
  }

  for( ii=0; ii<CSD->nPSb; ii++ )
   pmm->XFA[ CSD->xPH[ii] ] = pmm->YFA[ CSD->xPH[ii] ] = CNode->xPA[ii];

  for( ii=0; ii<CSD->nICb; ii++ )
   pmm->C[ CSD->xIC[ii] ] = CNode->rMB[ii];
  for( ii=0; ii<CSD->nICb; ii++ )
   pmm->U[ CSD->xIC[ii] ] = CNode->uIC[ii];
 }
//  End
}

// (5) For interruption/debugging
// Writes work node (DATABR structure) into a file path name fname
// Parameter binary_f defines if the file is to be written in binary
// format (true or 1, good for interruption of coupled modeling task
// if called in loop for each node), or in text format
// (false or 0, default)
//
   void  TNode::GEM_write_dbr( const char* fname, bool binary_f, bool with_comments )
   {
       gstring str_file;
       if( fname == 0)	   
    	   str_file = dbr_file_name+".out";
       else
           str_file = fname;
       
	   if( binary_f )
           {
            // gstring str_file = fname;
              GemDataStream out_br(str_file, ios::out|ios::binary);
              databr_to_file(out_br);
           }
      else
      {  fstream out_br(str_file.c_str(), ios::out );
         ErrorIf( !out_br.good() , str_file.c_str(), "DataBR text make error");
         databr_to_text_file(out_br, with_comments );
      }
   }

// (5a) For detailed examination of GEM work data structure:
// writes GEMIPM internal MULTI data structure into text file
// path name fname in debugging format (different from MULTI input format).
// This file cannot be read back with GEM_init()!
//
   void  TNode::GEM_print_ipm( const char* fname )
   {
     gstring str_file;  
     if( fname == 0)	   
    	   str_file = dbr_file_name + ".Dump.out";
     else
           str_file = fname;
      
	   TProfil::pm->outMultiTxt( str_file.c_str()  );
   }

#ifdef IPMGEMPLUGIN

// calculation mode: passing input GEM data changed on previous FMT iteration
//                   into the work DATABR structure
void TNode::GEM_from_MT(
   short  p_NodeHandle,   // Node identification handle
   short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double p_TC,      // Temperature T, K                         +       -      -
   double p_P,      // Pressure P, bar                          +       -      -
   double p_Vs,     // Volume V of reactive subsystem, cm3      -       -      +
   double p_Ms,     // Mass of reactive subsystem, kg           -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
   double *p_dul,   // upper kinetic restrictions [nDCb]        +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]        +       -      -
   double *p_aPH  // Specific surface areas of phases (m2/g)    +       -      -
)
{
  int ii;
  bool useSimplex = false;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->TC = p_TC;
  CNode->P = p_P;
  CNode->Vs = p_Vs;
  CNode->Ms = p_Ms;
// Checking if no-simplex IA is Ok
   for( ii=0; ii<CSD->nICb; ii++ )
   {  //  SD 11/02/05 for test
      //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
       //     useSimplex = true;
     CNode->bIC[ii] = p_bIC[ii];
   }
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
     CNode->dul[ii] = p_dul[ii];
     CNode->dll[ii] = p_dll[ii];
   }
    if( CSD->nAalp >0 )
     for( ii=0; ii<CSD->nPHb; ii++ )
         CNode->aPH[ii] = p_aPH[ii];
   if( useSimplex && CNode->NodeStatusCH == NEED_GEM_PIA )
     CNode->NodeStatusCH = NEED_GEM_AIA;
   // Switch only if PIA is ordered, leave if simplex is ordered (KD)
}

// readonly mode: passing input GEM data to FMT
void TNode::GEM_restore_MT(
   short  &p_NodeHandle,   // Node identification handle
   short  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double &p_TC,     // Temperature T, K                         +       -      -
   double &p_P,      // Pressure P, bar                          +       -      -
   double &p_Vs,     // Volume V of reactive subsystem, cm3      -       -      +
   double &p_Ms,     // Mass of reactive subsystem, kg           -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]           +       -      -
   double *p_dul,    // upper kinetic restrictions [nDCb]        +       -      -
   double *p_dll,    // lower kinetic restrictions [nDCb]        +       -      -
   double *p_aPH     // Specific surface areas of phases (m2/g)  +       -      -
   )
{
 int ii;
  p_NodeHandle = CNode->NodeHandle;
  p_NodeStatusCH = CNode->NodeStatusCH;
  p_TC = CNode->TC;
  p_P = CNode->P;
  p_Vs = CNode->Vs;
  p_Ms = CNode->Ms;
// Checking if no-simplex IA is Ok
   for( ii=0; ii<CSD->nICb; ii++ )
     p_bIC[ii] = CNode->bIC[ii];
   for( ii=0; ii<CSD->nDCb; ii++ )
   {  p_dul[ii] = CNode->dul[ii];
      p_dll[ii] = CNode->dll[ii];
   }
   if( CSD->nAalp >0 )
     for( ii=0; ii<CSD->nPHb; ii++ )
        p_aPH[ii] = CNode->aPH[ii];
}

// Copying results that must be returned into the FMT part into MAIF_CALC parameters
void TNode::GEM_to_MT(
       short &p_NodeHandle,    // Node identification handle
       short &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
       short &p_IterDone,      // Number of iterations performed by GEM IPM
                         //                                     GEM input output  FMT control
       // Chemical scalar variables
       double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
       double &p_Ms,    // Mass of reactive subsystem, kg          -      -      +     +
       double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
       double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
       double &p_IC,    // Effective molal aq ionic strength       -      -      +     +
       double &p_pH,    // pH of aqueous solution                  -      -      +     +
       double &p_pe,    // pe of aqueous solution                  -      -      +     +
       double &p_Eh,    // Eh of aqueous solution, V               -      -      +     +
       // Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
       // exchange of values occurs through lists of indices, e.g. xDC, xPH
       double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
       double  *p_uIC,  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
       double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
       double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
       double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
       double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
       double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
       double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
       double  *p_xPA  // amount of carrier in phases  [nPSb] ??       -      -      +     +
)
{
   int ii;
   p_NodeHandle = CNode->NodeHandle;
   p_NodeStatusCH = CNode->NodeStatusCH;
   p_IterDone = CNode->IterDone;

   p_Vs = CNode->Vs;
   p_Ms = CNode->Ms;
   p_Gs = CNode->Gs;
   p_Hs = CNode->Hs;
   p_IC = CNode->IC;
   p_pH = CNode->pH;
   p_pe = CNode->pe;
   p_Eh = CNode->Eh;

  for( ii=0; ii<CSD->nICb; ii++ )
  {
    p_rMB[ii] = CNode->rMB[ii];
    p_uIC[ii] = CNode->uIC[ii];
  }
  for( ii=0; ii<CSD->nDCb; ii++ )
  {
    p_xDC[ii] = CNode->xDC[ii];
    p_gam[ii] = CNode->gam[ii];
  }
  for( ii=0; ii<CSD->nPHb; ii++ )
    p_xPH[ii] = CNode->xPH[ii];
  for( ii=0; ii<CSD->nPSb; ii++ )
  {
    p_vPS[ii] = CNode->vPS[ii];
    p_mPS[ii] = CNode->mPS[ii];
    p_xPA[ii] = CNode->xPA[ii];
  }
  for( ii=0; ii<CSD->nPSb*CSD->nICb; ii++ )
    p_bPS[ii] = CNode->bPS[ii];
}

// Overloaded variant - takes input to bIC vector also from the speciation
//     vector xDC.     Added by DK on 09.07.2007
// calculation mode: passing input GEM data changed on previous FMT iteration
//                   into the work DATABR structure
void TNode::GEM_from_MT(
   short  p_NodeHandle,   // Node identification handle
   short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double p_TC,      // Temperature T, K                         +       -      -
   double p_P,      // Pressure P, bar                          +       -      -
   double p_Vs,     // Volume V of reactive subsystem, cm3      -       -      +
   double p_Ms,     // Mass of reactive subsystem, kg           -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
   double *p_dul,   // upper kinetic restrictions [nDCb]        +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]        +       -      -
   double *p_aPH,  // Specific surface areas of phases (m2/g)    +       -      -
   double *p_xDC  // Optional: mole amounts of DCs [nDCb] - will be convoluted
                    // and added to the bIC GEM input vector
)
{
  int ii;
  bool useSimplex = false;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->TC = p_TC;
  CNode->P = p_P;
  CNode->Vs = p_Vs;
  CNode->Ms = p_Ms;
// Checking if no-simplex IA is Ok
   for( ii=0; ii<CSD->nICb; ii++ )
   {
     CNode->bIC[ii] = p_bIC[ii];
   }
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
     CNode->dul[ii] = p_dul[ii];
     CNode->dll[ii] = p_dll[ii];
   }
    if( CSD->nAalp >0 )
     for( ii=0; ii<CSD->nPHb; ii++ )
         CNode->aPH[ii] = p_aPH[ii];
   if( useSimplex && CNode->NodeStatusCH == NEED_GEM_PIA )
     CNode->NodeStatusCH = NEED_GEM_AIA;
   // Switch only if PIA is ordered, leave if simplex is ordered (KD)

   // Optional part - convolution of xDC vector into bIC vector
   if( p_xDC )
   {  int jj;
      // Correction of bIC vector by convoluting the amounts of DCs
      for( jj=0; jj<CSD->nDCb; jj++ )
        if( p_xDC[jj] )
          for( ii=0; ii<CSD->nICb; ii++ )
            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
   }
}

// Overloaded variant - uses xDC and gam vectors as old primal solution
// for the node in GEM IPM2 input when NEED_GEM_PIA flag is set for calculation
// Important! This variant works only when DATACH contains a full list of DCs
// with passed through the DATABR structure.
// added by DK on 17.09.2007
// calculation mode: passing input GEM data changed on previous FMT iteration
//                   into the work DATABR structure
void TNode::GEM_from_MT(
   short  p_NodeHandle,   // Node identification handle
   short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double p_TC,      // Temperature T, K                         +       -      -
   double p_P,      // Pressure P, bar                          +       -      -
   double p_Vs,     // Volume V of reactive subsystem, cm3      -       -      +
   double p_Ms,     // Mass of reactive subsystem, kg           -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
   double *p_dul,   // upper kinetic restrictions [nDCb]        +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]        +       -      -
   double *p_aPH,  // Specific surface areas of phases (m2/g)    +       -      -
   double *p_xDC,  // Amounts of DCs [nDCb] - old primal soln.  +      -      -
   double *p_gam   // DC activity coeffs [nDCb] - old primal s. +      -      -
)
{
  int ii;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->TC = p_TC;
  CNode->P = p_P;
  CNode->Vs = p_Vs;
  CNode->Ms = p_Ms;
// Checking if no-simplex IA is Ok
   for( ii=0; ii<CSD->nICb; ii++ )
   {
     CNode->bIC[ii] = p_bIC[ii];
   }
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
     CNode->dul[ii] = p_dul[ii];
     CNode->dll[ii] = p_dll[ii];
   }
    if( CSD->nAalp >0 )
     for( ii=0; ii<CSD->nPHb; ii++ )
         CNode->aPH[ii] = p_aPH[ii];
         
   // Optional part - copying old primal solution from p_xDC and p_gam vectors
   if( p_xDC && p_gam )
   {
      for( ii=0; ii<CSD->nDCb; ii++ )
      {
        CNode->xDC[ii] = p_xDC[ii];
        CNode->gam[ii] = p_gam[ii];
      }
   }
   else if( CNode->NodeStatusCH == NEED_GEM_PIA )
            CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal
                                                  // provided!
}

#endif
//-----------------------End of node.cpp--------------------------



