//-------------------------------------------------------------------
// $Id: node.cpp 684 2005-11-23 13:17:15Z gems $
//
// Implementation of TNode class including initialization and
// execution of GEMIPM2 kernel

// Works with DATACH and DATABR structures
//
// Copyright (C) 2005,2009 S.Dmytriyeva, D.Kulik, G.Kosakowski
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

  istream& f_getline(istream& is, gstring& str, char delim);

// TNode* TNode::na;

// Conversion factors
const double bar_to_Pa = 1e5,
               m3_to_cm3 = 1e6,
               kg_to_g = 1e3;


// Checks if given temperature TK and pressure P fit within the interpolation
//intervals of the DATACH lookup arrays (returns true) or not (returns false)
bool  TNode::check_TP( double TK, double P )
{
   bool okT = true, okP = true;
   double T_=TK, P_=P;

   if( TK <= CSD->TKval[0] - CSD->Ttol )
   { 				// Lower boundary of T interpolation interval
	 okT = false;
     T_ = CSD->TKval[0] - CSD->Ttol;
   }
   if( TK >= CSD->TKval[CSD->nTp-1] + CSD->Ttol )
   {
	 okT = false;
     T_ = CSD->TKval[CSD->nTp-1] + CSD->Ttol;
   }
   if( okT == false )
   {
     fstream f_log("ipmlog.txt", ios::out|ios::app );
     f_log << "In node "<< CNode->NodeHandle << ",  Given TK= "<<  TK <<
             "  is beyond the interpolation range for thermodynamic data near boundary T_= "
     		<< T_ << endl;
   }

   if( P <= CSD->Pval[0] - CSD->Ptol )
   {
	  okP = false;
      P_ = CSD->Pval[0] - CSD->Ptol;
   }
   if( P >= CSD->Pval[CSD->nPp-1] + CSD->Ptol )
   {
	  okP = false;
      P_ = CSD->Pval[CSD->nPp-1] + CSD->Ptol;
   }
   if( !okP )
   {
     fstream f_log("ipmlog.txt", ios::out|ios::app );
     f_log << "In node "<< CNode->NodeHandle << ", Given P= "<<  P <<
           "  is beyond the interpolation range for thermodynamic data near boundary P_= "
           << P_ << endl;
   }
   return ( okT && okP );
}

//-------------------------------------------------------------------------------------------------------------------------------
// (2) Main call for GEM IPM calculations using the input bulk composition, temperature, pressure
//   and metastability constraints provided in the work instance of DATABR structure.
//   Actual calculation will be performed only when dBR->NodeStatusCH == NEED_GEM_SIA (5) or dBR->NodeStatusCH = NEED_GEM_AIA (1).
//   By other values of NodeStatusCH, no calculation will be performed and the status will remain unchanged.
//  In "smart initial approximation" (SIA) mode, the program can automatically switch into the "automatic initial
//  approximation" (AIA) mode and return  OK_GEM_AIA instead of OK_GEM_SIA.
//  Parameter:
//   uPrimalSol  flag to define the mode of GEM smart initial approximation
//               (only if dBR->NodeStatusCH = NEED_GEM_SIA has been set before GEM_run() call).
//               false  (0) -  use speciation and activity coefficients from previous GEM_run() calculation
//               true  (1)  -  use speciation provided in the DATABR memory structure (e.g. after reading the DBR file)
//  Return values:    NodeStatusCH  (the same as set in dBR->NodeStatusCH). Possible values (see "databr.h" file for the full list)
long int TNode::GEM_run( bool uPrimalSol )
{
  CalcTime = 0.0;
  PrecLoops = 0; NumIterFIA = 0; NumIterIPM = 0;
//
  try
  {
// f_log << " GEM_run() begin Mode= " << p_NodeStatusCH endl;
//---------------------------------------------
// Checking T and P  for interpolation intervals
   check_TP( CNode->TK, CNode->P);
// Unpacking work DATABR structure into MULTI (GEM IPM structure): uses DATACH
// setting up up PIA or AIA mode
   if( CNode->NodeStatusCH == NEED_GEM_SIA )
   {
	   pm.pNP = 1;
	   unpackDataBr( uPrimalSol );
   }
   else if( CNode->NodeStatusCH == NEED_GEM_AIA )
	     {  pm.pNP = 0; // As default setting AIA mode
	        unpackDataBr( false );
         }
        else
	       return CNode->NodeStatusCH;
   // GEM IPM calculation of equilibrium state
   CalcTime = ComputeEquilibriumState( PrecLoops, NumIterFIA, NumIterIPM , this);
// Extracting and packing GEM IPM results into work DATABR structure
    packDataBr();
    CNode->IterDone = NumIterFIA+NumIterIPM;
//**************************************************************
// only for testing output results for files
//    GEM_write_dbr( "calculated_dbr.dat",  false );
//    GEM_print_ipm( "calc_multi.ipm" );
// *********************************************************

    // test error result GEM IPM calculation of equilibrium state in MULTI
    long int erCode = testMulti();

    if( erCode )
    {
        if( CNode->NodeStatusCH  == NEED_GEM_AIA )
          CNode->NodeStatusCH = BAD_GEM_AIA;
        else
          CNode->NodeStatusCH = BAD_GEM_SIA;
    }
    else
    {
      if( CNode->NodeStatusCH  == NEED_GEM_AIA )
          CNode->NodeStatusCH = OK_GEM_AIA;
      else
         CNode->NodeStatusCH = OK_GEM_SIA;
    }

   }
   catch(TError& err)
   {
        if( pa.p.PSM  )
	{
		fstream f_log("ipmlog.txt", ios::out|ios::app );
        f_log << "Error Node:" << CNode->NodeHandle << ":time:" << CNode->Tm << ":dt:" << CNode->dt<< ": " <<
          err.title.c_str() << ":" << endl;
       if( pa.p.PSM >= 2  )
          f_log  << err.mess.c_str() << endl;
	}
    if( CNode->NodeStatusCH  == NEED_GEM_AIA )
      CNode->NodeStatusCH = ERR_GEM_AIA;
    else
      CNode->NodeStatusCH = ERR_GEM_SIA;

   }
   catch(...)
   {
    fstream f_log("ipmlog.txt", ios::out|ios::app );
    f_log << "Node:" << CNode->NodeHandle << ":time:" << CNode->Tm << ":dt:" << CNode->dt<< ": "
                << "gems3: Unknown exception: GEM calculation aborted" << endl;
    CNode->NodeStatusCH = T_ERROR_GEM;
    }
   return CNode->NodeStatusCH;
}



// Returns GEMIPM2 calculation time in seconds elapsed during the last call of GEM_run() -
// can be used for monitoring the performance of calculations.
// Return value:  double number, may contain 0.0 if the calculation time is less than the
//                internal time resolution of C/C++ function
double TNode::GEM_CalcTime()
{
  return CalcTime;
}

// To obtain the number of GEM IPM2 iterations performed during the last call of GEM_run() e.g. for monitoring the
// performance of GEMIPM2K in AIA or SIA modes, or for problem diagnostics.
// Parameters:  long int variables per reference (must be allocated before calling GEM_Iterations(), previous values will be lost. See Return values.
// Return values:
//   Function         Total number of EFD + IPM iterations from the last call to GEM_run()
//   PrecLoops        Number of performed IPM-2 precision refinement loops
//   NumIterFIA       Total number of performed MBR() iterations to obtain a feasible initial approximation for the IPM algorithm.
//   NumIterIPM       Total number of performed IPM main descent algorithm iterations.
long int TNode::GEM_Iterations( long int& PrecLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
	PrecLoops_ = PrecLoops;
	NumIterFIA_ = NumIterFIA;
	NumIterIPM_ = NumIterIPM;
	return NumIterFIA+NumIterIPM;
}

// (5) Reads another DBR file (with input system composition, T,P etc.). The DBR file must be compatible with
// the currently loaded IPM and DCH files (see description  of GEM_init() function call).
// Parameters:
//    fname       Null-terminated (C) string containing a full path to the input DBR disk file.
//    binary_f    Flag defining whether the file specified in fname is in text fromat (false or 0, default) or in binary format (true or 1)
// Return values:     0  if successful; 1 if input file(s) has not found been or is corrupt; -1 if internal memory allocation error occurred.
long int  TNode::GEM_read_dbr( const char* fname, bool binary_f )
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

  } catch(TError& err)
    {
	  fstream f_log("ipmlog.txt", ios::out|ios::app );
      f_log << "Error Node:" << CNode->NodeHandle << ":time:" << CNode->Tm << ":dt:" << CNode->dt<< ": " <<
        err.title.c_str() << ":" <<  err.mess.c_str() << endl;
      return 1;
    }
    catch(...)
    {
      return -1;
    }
  return 0;
}

//-------------------------------------------------------------------
// (1) Initialization of GEM IPM2 data structures in coupled FMT-GEM programs
//  that use GEMIPM2K module. Also reads in the IPM, DCH and DBR text input files.
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
long int  TNode::GEM_init( const char* ipmfiles_lst_name,
                          long int* nodeTypes, bool /*getNodT1*/)
{
  long int i;
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
    {
     bool binary_f = false;
     gstring lst_in = ipmfiles_lst_name;
     gstring Path = "";         // was " "   fixed 10.12.2009 by DK
// Get path
#ifdef _WIN32
      size_t pos = lst_in.rfind("\\");// HS keep this on windows
#else      
      size_t pos = lst_in.rfind("/"); // HS keep this on linux
#endif
      if( pos < npos )
      Path = lst_in.substr(0, pos+1);

//  open file stream for the file names list file
      fstream f_lst( lst_in.c_str(), ios::in );
      ErrorIf( !f_lst.good() , lst_in.c_str(), "Fileopen error");

      gstring datachbr_fn;
      f_getline( f_lst, datachbr_fn, ' ');

//  Syntax: -t/-b  "<DCH_DAT file name>"  "<IPM_DAT file name>"
//       "<DBR_DAT file1 name>" [, ... , "<DBR_DAT fileN name>"]

//Testing flag "-t" or "-b" (by default "-t")   // use binary or text files
      pos = datachbr_fn.find( '-');
      if( pos != /*gstring::*/npos )
      {
         if( datachbr_fn[pos+1] == 'b' )
            binary_f = true;
         f_getline( f_lst, datachbr_fn, ' ');
      }
 //     f_getline( f_lst, datachbr_fn, ' ');
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
    readMulti( f_m, this);
  }
  else
  {
    readMulti( this, mult_in.c_str());
  }

// Prepare for reading DBR_DAT files
     i = 0;
     while( !f_lst.eof() )  // For all DBR_DAT files listed
     {


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

        {
// Copying data from work DATABR structure into the node array
// (as specified in nodeTypes array)
           setNodeArray( i, nodeTypes  );
         }
          i++;
     }  // end while()

    ErrorIf( i==0, datachbr_fn.c_str(), "GEM_init() error: No DBR_DAT files read!" );
    checkNodeArray( i, nodeTypes, datachbr_fn.c_str()  );

   return 0;

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
}

//-----------------------------------------------------------------
// work with lists

//Returns DCH index of IC given the IC Name string (null-terminated)
// or -1 if no such name was found in the DATACH IC name list
long int TNode::IC_name_to_xCH( const char *Name )
{
  long int ii, len = strlen( Name );
  len =  min(len,MaxICN);

  for(ii = 0; ii<CSD->nIC; ii++ )
       if(!memcmp(Name, CSD->ICNL[ii], len ))
        if( len == MaxICN || CSD->ICNL[ii][len] == ' ' || CSD->ICNL[ii][len] == '\0' )
         return ii;
  return -1;
}

// Returns DCH index of DC given the DC Name string
// or -1 if no such name was found in the DATACH DC name list
 long int TNode::DC_name_to_xCH( const char *Name )
 {
  long int ii, len = strlen( Name );
  len =  min(len,MaxDCN);

  for( ii = 0; ii<CSD->nDC; ii++ )
       if(!memcmp(Name, CSD->DCNL[ii], min(len,MaxDCN)))
        if( len == MaxDCN || CSD->DCNL[ii][len] == ' ' || CSD->DCNL[ii][len] == '\0' )
         return ii;
  return -1;
 }

// Returns DCH index of Phase given the Phase Name string
// or -1 if no such name was found in the DATACH Phase name list
long int TNode::Ph_name_to_xCH( const char *Name )
{
  long int ii, len = strlen( Name );
  len =  min(len,MaxPHN);

  for( ii = 0; ii<CSD->nPH; ii++ )
       if(!memcmp(Name, CSD->PHNL[ii], min(len,MaxPHN)))
        if( len == MaxPHN || CSD->PHNL[ii][len] == ' ' || CSD->PHNL[ii][len] == '\0' )
         return ii;
  return -1;
}

// Converts the IC DCH index into the IC DBR index
// or returns -1 if this IC is not used in the data bridge
long int TNode::IC_xCH_to_xDB( const long int xCH )
{
  for(long int ii = 0; ii<CSD->nICb; ii++ )
       if( CSD->xic[ii] == xCH )
         return ii;
  return -1;
}

// Converts the DC DCH index into the DC DBR index
// or returns -1 if this DC is not used in the data bridge
long int TNode::DC_xCH_to_xDB( const long int xCH )
{
  for(long int ii = 0; ii<CSD->nDCb; ii++ )
       if( CSD->xdc[ii] == xCH )
         return ii;
  return -1;
}

// Converts the Phase DCH index into the Phase DBR index
// or returns -1 if this Phase is not used in the data bridge
long int TNode::Ph_xCH_to_xDB( const long int xCH )
{
  for(long int ii = 0; ii<CSD->nPHb; ii++ )
       if( CSD->xph[ii] == xCH )
         return ii;
  return -1;
}

// Returns the DCH index of the first DC belonging to the phase with DCH index Phx
 long int  TNode::Phx_to_DCx( const long int Phx )
 {
   long int k, DCx = 0;
   for( k=0; k<CSD->nPHb; k++ )
   {
     if( k == Phx )
      break;
     DCx += CSD->nDCinPH[ k];
   }
   return DCx;
 }

 // Returns the DCH index of the Phase to which the Dependent Component with index xCH belongs
  long int  TNode::DCtoPh_DCH( const long int xdc )
  {
    long int k, DCx = 0;
    for( k=0; k<CSD->nPHb; k++ )
    {
    	DCx += CSD->nDCinPH[ k];
        if( xdc < DCx )
          break;
    }
    return k;
  }


 // Returns the DCH index of the first DC belonging to the phase with DCH index Phx,
 // plus returns through the nDCinPh (reference) parameter the number of DCs included into this phase
 long int  TNode::PhtoDC_DCH( const long int Phx, long int& nDCinPh )
 {
   long int k, DCx = 0;
   for( k=0; k<CSD->nPHb; k++ )
   {
     if( k == Phx )
      break;
     DCx += CSD->nDCinPH[ k];
   }
   nDCinPh = CSD->nDCinPH[k];
   return DCx;
 }

 // Returns the DBR index of the Phase to which the  Dependent Component with index xBR belongs
  long int  TNode::DCtoPh_DBR( const long int xBR )
  {
    long int DCxCH = DC_xDB_to_xCH( xBR );
    long int PhxCH = DCtoPh_DCH( DCxCH );
    return Ph_xCH_to_xDB(PhxCH);
  }

// Returns the DBR index of the first DC belonging to the phase with DBR index Phx,
//plus returns through the nDCinPh (reference) parameter the number of DCs included into DBR for this phase
 long int  TNode::PhtoDC_DBR( const long int Phx, long int& nDCinPh )
 {
   long int ii, DCx, DCxCH, PhxCH, nDCinPhCH;

   PhxCH = Ph_xDB_to_xCH( Phx );
   DCxCH = PhtoDC_DCH( PhxCH, nDCinPhCH );

   DCx = -1;
   nDCinPh = 0;
   for( ii = 0; ii<CSD->nDCb; ii++ )
   {
      if( CSD->xdc[ii] >= DCxCH )
      {
        if( CSD->xdc[ii] >= DCxCH+nDCinPhCH  )
          break;
        nDCinPh++;
        if( DCx == -1)
          DCx = ii;
      }
   }
   return DCx;
 }

 // Test TK as lying in the vicinity of a grid point for the interpolation of thermodynamic data
 // Return index of the node in lookup array or -1
 long int  TNode::check_grid_T( double TK )
 {
   long int jj;
   for( jj=0; jj<CSD->nTp; jj++)
     if( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol )
        return jj;
   return -1;
 }

 // Test P as lying in the vicinity of a grid point for the interpolation of thermodynamic data
 // Return index of the node in lookup array or -1
 long int  TNode::check_grid_P( double P )
 {
   long int jj;
   for( jj=0; jj<CSD->nPp; jj++)
     if( fabs( P - CSD->Pval[jj] ) < CSD->Ptol )
       return jj;
   return -1;
 }

 // Tests TK and P as a grid point for the interpolation of thermodynamic data using DATACH
 // lookup arrays. Returns -1L if interpolation is needed, or 1D index of the lookup array element
 // if TK and P fit within the respective tolerances.
 // For producing lookup arrays (in GEMS), we recommend using step for temperature less or equal to 10 degrees
 // in order to assure good accuracy of interpolation especially for S0 and Cp0 of aqueous species.
  long int  TNode::check_grid_TP(  double TK, double P )
  {
    long int xT, xP, ndx=-1;

    xT = check_grid_T( TK );
    xP = check_grid_P( P );
    if( xT >=0 && xP>= 0 )
     ndx =  xP * CSD->nTp + xT;
    return ndx;
  }

  //Retrieves (interpolated) molar Gibbs energy G0(P,TK) value for Dependent Component
  //from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
  // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
  // Parameter norm defines in wnich units the value is returned: false - in J/mol; true (default) - in mol/mol
   double TNode::DC_G0(const long int xCH, const double P, const double TK,  bool norm )
   {
    long int xTP, jj;
    double G0;

    if( check_TP( TK, P ) == false )
    	return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       G0 = CSD->G0[ jj + xTP ];
    else
       G0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->G0+jj,
               P, TK, CSD->nTp, CSD->nPp, 6 );

    if( norm )
      return G0/(R_CONSTANT * (TK));
    else
      return G0;
   }

   // Retrieves (interpolated, if necessary) molar volume V0(P,TK) value for Dependent Component (in J/Pa)
   // from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_V0(const long int xCH, const double P, const double TK)
   {
    long int xTP, jj;
    double V0;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       V0 = CSD->V0[ jj + xTP ];
    else
       V0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->V0+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return V0;
   }


   // Retrieves (interpolated) molar enthalpy H0(P,TK) value for Dependent Component (in J/mol)
   // from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_H0(const long int xCH, const double P, const double TK)
   {
    long int xTP, jj;
    double H0;

    if( check_TP( TK, P ) == false )
    	return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       H0 = CSD->H0[ jj + xTP ];
    else
       H0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->H0+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return H0;
   }

   // Retrieves (interpolated) absolute molar enropy S0(P,TK) value for Dependent Component (in J/K/mol)
   // from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_S0(const long int xCH, const double P, const double TK)
   {
    long int xTP, jj;
    double s0;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       s0 = CSD->S0[ jj + xTP ];
    else
       s0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->S0+jj,
                P, TK, CSD->nTp, CSD->nPp, 4 );
    return s0;
   }

   // Retrieves (interpolated) constant-pressure heat capacity Cp0(P,TK) value for Dependent Component (in J/K/mol)
   // from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_Cp0(const long int xCH, const double P, const double TK)
   {
    long int xTP, jj;
    double cp0;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       cp0 = CSD->Cp0[ jj + xTP ];
    else
       cp0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->Cp0+jj,
                P, TK, CSD->nTp, CSD->nPp, 3 );
    return cp0;
   }

   // Retrieves (interpolated) Helmholtz energy  of Dependent Component (in J/mol)
   // from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_A0(const long int xCH, const double P, const double TK)
   {
    long int xTP, jj;
    double a0;

    if( check_TP( TK, P ) == false )
    	return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
       a0 = CSD->A0[ jj + xTP ];
    else
       a0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->A0+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return a0;
   }

   // Retrieves (interpolated) Internal energy of  Dependent Component (in J/mol)
   // from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
   // or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DC_U0(const long int xCH, const double P, const double TK)
   {
       long int xTP, jj;
       double u0;

       if( check_TP( TK, P ) == false )
       	return 7777777.;

       xTP = check_grid_TP( TK, P );
       jj =  xCH * CSD->nPp * CSD->nTp;

       if( xTP >= 0 )
          u0 = CSD->U0[ jj + xTP ];
       else
          u0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->U0+jj,
                   P, TK, CSD->nTp, CSD->nPp, 5 );
       return u0;
  }


   // Retrieves (interpolated) dielectric constant of liquid water at (P,TK) from the DATACH structure or 0.0,
   // if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::EpsH2Ow(const double P, const double TK)
   {
    long int xTP, jj;
    double epsW;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
    	epsW = CSD->epsW[ jj + xTP ];
    else
    	epsW = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return epsW;
   }

   // Retrieves (interpolated) density of liquid water (in kg/m3) at (P,TK) from the DATACH structure or 0.0,
   // if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DenH2Ow(const double P, const double TK)
   {
    long int xTP, jj;
    double denW;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
    	denW = CSD->denW[ jj + xTP ];
    else
    	denW = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return denW;
   }

   // Retrieves (interpolated) dielectric constant of H2O vapor at (P,TK) from the DATACH structure or 0.0,
   // if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::EpsH2Og(const double P, const double TK)
   {
     long int xTP, jj;
     double epsWg;

     if( check_TP( TK, P ) == false )
     	return 0.;

     xTP = check_grid_TP( TK, P );
     jj = 0; // 0 * CSD->nPp * CSD->nTp;

     if( xTP >= 0 )
     	epsWg = CSD->epsWg[ jj + xTP ];
     else
     	epsWg = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsWg+jj,
                 P, TK, CSD->nTp, CSD->nPp, 5 );
     return epsWg;
    }

   // Retrieves (interpolated) density of H2O vapor (in kg/m3) at (P,TK) from the DATACH structure or 0.0,
   // if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
   double TNode::DenH2Og(const double P, const double TK)
   {
    long int xTP, jj;
    double denWg;

    if( check_TP( TK, P ) == false )
    	return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * CSD->nPp * CSD->nTp;

    if( xTP >= 0 )
    	denWg = CSD->denWg[ jj + xTP ];
    else
    	denWg = LagranInterp( CSD->Pval, CSD->TKval, CSD->denWg+jj,
                P, TK, CSD->nTp, CSD->nPp, 5 );
    return denWg;
   }

 //Retrieves the current phase volume in m3 ( xph is DBR phase index) in the reactive sub-system.
 // Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
 double  TNode::Ph_Volume( const long int xBR )
 {
   double vol;
   if( xBR < CSD->nPSb )
    vol = CNode->vPS[xBR];
   else
   {
     long int xDC = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
     vol = DC_V0( xDC, CNode->P, CNode->TK );
     vol *= CNode->xDC[DC_xCH_to_xDB(xDC)];
   }
   return vol;
 }

  // Retrieves the phase mass in kg ( xph is DBR phase index).
  // Works for multicomponent and for single-component phases. Returns 0.0 if phase amount is zero.
  double  TNode::Ph_Mass( const long int xBR )
  {
     double mass;
     if( xBR < CSD->nPSb )
        mass = CNode->mPS[xBR];
     else
     {
        long int xDC = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
        mass = CNode->xDC[ DC_xCH_to_xDB(xDC) ] * CSD->DCmm[xDC];
     }
    return mass;
  }

  // Retrieves the phase saturation index ( xph is DBR phase index). Works for multicomponent and for
  // single-component phases. Returns 0.0 if phase amount is zero.
  double TNode::Ph_SatInd(const long int xph )
  {
    double SatInd=0.;
    long int jj, dcx1, Ndc;
    dcx1 = PhtoDC_DBR( xph, Ndc );

    if( xph < CSD->nPSb )
	{
        for( jj=dcx1; jj<Ndc+dcx1; jj++)
        	SatInd +=  Get_aDC( jj )/Get_gDC(jj);
	}
	else
	{
	  SatInd = Get_aDC( dcx1 );
	}
    return SatInd;
  }

  // Retrieval of the phase bulk composition ( xph is DBR phase index) into memory indicated by
  // ARout (array of at least [dCH->nICb elements]). Returns pointer to ARout which may also be
  // allocated inside of Ph_BC() in the case if parameter ARout = NULL is specified;
  // to avoid a memory leak, you will have to free this memory wherever appropriate.
  // This function works for multicomponent and for single-component phases
  double *TNode::Ph_BC( const long int xBR, double* ARout )
  {
    long int ii;
    if( !ARout )
      ARout = new double[ CSD->nICb ];   // Potential memory leak ! ! ! ! ! ! ! !

    if( xBR < CSD->nPSb )
       for( ii=0; ii<CSD->nICb; ii++ )
         ARout[ii] = CNode->bPS[ xBR * CSD->nICb + ii ];
    else
    {
      long int DCx = Phx_to_DCx( Ph_xDB_to_xCH(xBR) );
      for( ii=0; ii<CSD->nICb; ii++ )
      {
         ARout[ii] = CSD->A[ IC_xDB_to_xCH(ii) + DCx * CSD->nIC];
         ARout[ii] *= CNode->xDC[ DC_xCH_to_xDB(DCx) ];
      }
    }
    return ARout;
  }

  // Retrieval of (dual-thermodynamic) chemical potential of the DC (xdc is the DC DBR index).
  // Parameter norm defines the scale: if true (1) then in mol/mol, otherwise in J/mol
  double TNode::Get_muDC( const long int xdc, bool norm )
  {	long int xCH, ii;
	double muDC = 0;

	xCH = DC_xDB_to_xCH(xdc);
    for( ii=0; ii<CSD->nICb; ii++ )
           muDC += CSD->A[  xCH * CSD->nIC + IC_xDB_to_xCH(ii) ] * CNode->uIC[ ii ];

    if( norm )
      return muDC;
    else
      return muDC*(R_CONSTANT * (CNode->TK));
  }

  //Retrieval of (dual-thermodynamic) activity of the DC (xdc is the DC DBR index)
  double TNode::Get_aDC( const long int xdc )
   {
	 double Mj  = Get_muDC( xdc, true );
	 double Mj0 = DC_G0( DC_xDB_to_xCH(xdc), CNode->P, CNode->TK,  true );
     return exp(Mj-Mj0);
	 // return 	pow(10.0,pm.Y_la[xCH]);
  }

  // Retrieves concentration of Dependent Component (xdc is the DC DBR index) in its phase
  // in the respective concentration scale. For aqueous species, molality is returned;
  // for gas species, partial pressure; for surface complexes - density in mol/m2; 
  // for species in other phases - mole fraction. If DC has zero amount, the function returns 0.0.
  double TNode::Get_cDC( const long int xdc )
  {
    long int xph = DCtoPh_DBR( xdc);
    long int DCxCH = DC_xDB_to_xCH(xdc);
	double DCcon = 0.;

    switch( CSD->ccDC[DCxCH] )
	    {
	     case DC_SCP_CONDEN:

	     case DC_AQ_SOLVENT:
	     case DC_AQ_SOLVCOM:

         case DC_SOL_IDEAL:
	     case DC_SOL_MINOR:
	     case DC_SOL_MAJOR:

	     case DC_PEL_CARRIER:
	     case DC_SUR_MINAL:
	     case DC_SUR_CARRIER:
	    	                  DCcon =  CNode->xDC[xdc]/CNode->xPH[xph];  //pm.Wx[xCH];
	                          break;
	      case DC_GAS_COMP:
	      case DC_GAS_H2O:
	      case DC_GAS_CO2:
	      case DC_GAS_H2:
	      case DC_GAS_N2:    DCcon =  CNode->xDC[xdc]/CNode->xPH[xph]*CNode->P;  //pm.Wx[xCH]*(pm.P*bar_to_Pa);
	                          break;
	     case DC_AQ_PROTON:
	     case DC_AQ_SPECIES:
	     case DC_AQ_SURCOMP:

	      case DC_SUR_GROUP:
	      case DC_SSC_A0:
	      case DC_SSC_A1:
	      case DC_SSC_A2:
	      case DC_SSC_A3:
	      case DC_SSC_A4:
	      case DC_WSC_A0:
	      case DC_WSC_A1:
	      case DC_WSC_A2:
	      case DC_WSC_A3:
	      case DC_WSC_A4:
	      case DC_SUR_COMPLEX:
	      case DC_SUR_IPAIR:
	      case DC_IESC_A:
	      case DC_IEWC_B:
	          {
	            double MMC = 0., Factor;
	            long int PhxCH_aquel = Ph_xDB_to_xCH(0);
	            long int H2Obr, jj;

	            if(CSD->ccPH[PhxCH_aquel] != PH_AQUEL )
	            	break;                                  // no aquel phase in DATABR

	            if( CNode->xPA[0] > 1e-19 )
	            {
	                for(jj=0; jj<CSD->nDCinPH[PhxCH_aquel]; jj++)
	                    if( CSD->ccDC[jj] == DC_AQ_SOLVENT || CSD->ccDC[jj] == DC_AQ_SOLVCOM )
	                    {  H2Obr = DC_xCH_to_xDB(jj);
	                       if( H2Obr >= 0 )
	                      	  MMC += CSD->DCmm[jj]*CNode->xDC[H2Obr]/CNode->xPA[0];
	                    }
	            }
	            else MMC=18.01528; // Assuming water-solvent
	            if( (CNode->xPA[0] > 1e-19) && (MMC > 1e-19) )
	                Factor = 1./MMC/CNode->xPA[0]; // molality
	            else Factor = 0.0;

	    	    DCcon =  CNode->xDC[xdc]*Factor;  //pm.Y_m[xCH];
	          }
	          break;
	      default:
	          break; // error in DC class code
	      }
	   return DCcon;
  }

  // Access to equilibrium properties of phases and components using DATACH indexation

  // Retrieves the current (dual-thermodynamic) activity of DC (xCH is DC DCH index)
  // directly from GEM IPM work structure. Also activity of a DC not included into DATABR list
  // can be retrieved. If DC has zero amount, its dual-thermodynamic activity is returned anyway.
  // For single condensed phase component, this value has a meaning of the saturation index,
  // also in the presence of metastability constraint(s).
  double TNode::DC_a(const long int xCH)
  {
	 //double Mj  = DC_mu( xCH, true );
	 //double Mj0 = DC_G0( xCH, CNode->P, CNode->TK,  true );
	 //return (Mj-Mj0)/2.302585093;
	 return 	pow(10.0,pm.Y_la[xCH]);
  }

    
  // Functions needed by GEM_FIT. Setting parameters for activity coefficient models.
    //     aIPx     = pm.IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
    //                               // -> NPar x MaxOrd   added 07.12.2006   KD
    //     aIPc     = pm.PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    //     aDCc     = pm.DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    //     NComp    = pm.L1[k];          // Number of components in the phase
    //     NPar     = pm.LsMod[k*3];      // Number of interaction parameters
    //     NPcoef   = pm.LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    //     MaxOrd   = pm.LsMod[k*3+1];  // max. parameter order (cols in IPx)
    //     NP_DC    = pm.LsMdc[k]; // Number of non-ideality coeffs per one DC in multicomponent phase
  void TNode::Get_IPc_IPx_DCc_indices( long* index_phase_aIPx, long* index_phase_aIPc, long* index_phase_aDCc, long* index_phase)
  {
    long ip_IPx=0; long ip_IPc=0; long ip_DCc=0;
    for(int k=0;k<((*index_phase)+1);k++)
    {
        *index_phase_aIPx = ip_IPx;
        ip_IPx          += pm.LsMod[k*3] * pm.LsMod[k*3+1];
        *index_phase_aIPc = ip_IPc;
        ip_IPc          += pm.LsMod[k*3] * pm.LsMod[k*3+2];
        *index_phase_aDCc = ip_DCc;
        ip_DCc          += pm.LsMdc[k] * pm.L1[k];
    }
  }

  void TNode::Get_NPar_NPcoef_MaxOrd_NComp_NP_DC ( long* NPar, long* NPcoef, long* MaxOrd, long* NComp, long* NP_DC, long* index_phase )
  {
    *NPar   = pm.LsMod[(*index_phase)*3];
    *NPcoef = pm.LsMod[(*index_phase)*3+2];
    *MaxOrd = pm.LsMod[(*index_phase)*3+1];
    *NComp  = pm.L1[(*index_phase)];
    *NP_DC  = pm.LsMdc[(*index_phase)];
    cout<<"*NPcoef in node.cpp = "<<*NPcoef<<endl;
  }

  void TNode::Set_aIPc ( double* aIPc, long* index_phase_aIPc, long *index_phase )
  { 
    int rc, NPar, NPcoef;
    NPar = pm.LsMod[(*index_phase)*3];
    NPcoef =  pm.LsMod[(*index_phase)*3+2];
    for ( rc=0;rc<(NPar*NPcoef);rc++ )
    {
        (pm.PMc[(*index_phase_aIPc)+rc]) = aIPc[rc];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
    }
  }

  void TNode::Get_aIPc ( double *aIPc, long* index_phase_aIPc, long* index_phase )
  { 
    int i; long NPar, NPcoef;
    NPar   = pm.LsMod[(*index_phase)*3];
    NPcoef = pm.LsMod[(*index_phase)*3+2];
    i = 0;
    while (i<(NPar*NPcoef))
    {
      *(aIPc + i)   = pm.PMc[(*index_phase_aIPc) + i];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
      i++;
    }
  }

  void TNode::Get_aIPx ( long* aIPx, long* index_phase_aIPx, long* index_phase )
  {
    int i; long NPar, MaxOrd;
    NPar   = pm.LsMod[(*index_phase)*3];
    MaxOrd = pm.LsMod[(*index_phase)*3+1];
    i = 0;
    while (i<(NPar*MaxOrd))
    {
      *(aIPx + i)   = pm.IPx[(*index_phase_aIPx) + i];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
      i++;
    }
  }

  void TNode::Set_aDCc ( const double* aDCc, long* index_phase_aDCc, long* index_phase )
  { 
    int rc, NComp, NP_DC;
    NComp = pm.L1[(*index_phase)];
    NP_DC = pm.LsMdc[(*index_phase)];
    for ( rc=0;rc<NComp*NP_DC;rc++ )
    {
        (pm.DMc[(*index_phase_aDCc)+rc]) = aDCc[rc];		// end-member param coeffs, NComp * NP_DC
    }
  }

  void TNode::Get_aDCc ( double* aDCc, long* index_phase_aDCc, long* index_phase )
  {
    int i; long NComp, NP_DC;
    NComp = pm.L1[(*index_phase)];
    NP_DC = pm.LsMdc[(*index_phase)];
    i = 0;
    while (i<(NComp*NP_DC))
    {
      *(aDCc + i)   = pm.DMc[(*index_phase_aDCc)+i];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
      i++;
    }
  }

  void TNode::Set_Tk   ( double* T_k)
  {
      CNode->TK = *T_k;
  }

  void TNode::Set_Pb   ( double* P_b)
  {
      CNode->P = *P_b;
  }


  // Retrieves the current concentration of Dependent Component (xCH is DC DCH index) in its
  // phase directly from GEM IPM work structure.Also activity of a DC not included into 
  // DATABR list can be retrieved. For aqueous species, molality is returned; for gas species, 
  // partial pressure; for surface complexes - density in mol/m2; for species in other phases - 
  // mole fraction. If DC has zero amount, the function returns 0.0. 
  double TNode::DC_c(const long int xCH)
  {
    double DCcon = 0.;
	switch( pm.DCC[xCH] )
    {
     case DC_SCP_CONDEN: DCcon =  pm.Wx[xCH];
                          break;
     case DC_AQ_PROTON:
     case DC_AQ_SPECIES:
     case DC_AQ_SURCOMP: DCcon =  pm.Y_m[xCH];
                          break;
     case DC_AQ_SOLVENT:
     case DC_AQ_SOLVCOM: DCcon =  pm.Wx[xCH];
                          break;
      case DC_GAS_COMP:
      case DC_GAS_H2O:
      case DC_GAS_CO2:
      case DC_GAS_H2:
      case DC_GAS_N2:    DCcon =  pm.Wx[xCH]*(pm.P*bar_to_Pa);
                          break;
      case DC_SOL_IDEAL:
      case DC_SOL_MINOR:
      case DC_SOL_MAJOR: DCcon =  pm.Wx[xCH];
                          break;
      case DC_SUR_GROUP:
      case DC_SSC_A0:
      case DC_SSC_A1:
      case DC_SSC_A2:
      case DC_SSC_A3:
      case DC_SSC_A4:
      case DC_WSC_A0:
      case DC_WSC_A1:
      case DC_WSC_A2:
      case DC_WSC_A3:
      case DC_WSC_A4:
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
      case DC_IESC_A:
      case DC_IEWC_B:     DCcon =  pm.Y_m[xCH];
                          break;
      case DC_PEL_CARRIER:
      case DC_SUR_MINAL:
      case DC_SUR_CARRIER: DCcon =  pm.Wx[xCH];
                           break;
      default:
          break; // error in DC class code
      }
   return DCcon;
  }

  // Retrieves the current (dual-thermodynamic) chemical potential of DC (xCH is DC DCH index)
  // directly from GEM IPM work structure, also for any DC not included into DATABR or having zero amount.
  // Parameter norm defines in wnich units the chemical potential value is returned:
  // false - in J/mol; true (default) - in mol/mol
  double TNode::DC_mu(const long int xCH, bool norm)
  {
	double muDC = pm.Fx[xCH];
 //  for(long ii=0; ii<CSD->nIC; ii++ )
 //    muDC += pm.A[  xCH * CSD->nIC + ii ] * (pm.U[ii]);
    if( norm )
        return muDC/pm.RT; // (R_CONSTANT * (CNode->TK));
    else
        return muDC;
  }

  // Retrieves the standard chemical potential of DC (xCH is DC DCH index) directly
  // from GEM IPM work structure at current pressure and temperature,
  // also for any DC not included into DATABR or having zero amount.
  // Parameter norm defines in which units the chemical potential value is returned:
  // false - in J/mol; true (default) - in mol/mol
  double TNode::DC_mu0(const long int xCH, bool norm)
  {
	return  DC_G0( xCH, CNode->P, CNode->TK, norm );
	/*double  G0 = pm.G0[xCH];
    if( norm )
      return G0;
    else
      return G0*pm.RT;
    */
  }

//---------------------------------------------------------//

void TNode::allocMemory()
{
// memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;

    // mem_set( CSD, 0, sizeof(DATACH) );
    datach_reset();
    // mem_set( CNode, 0, sizeof(DATABR) );
    databr_reset( CNode, 2 );

// internal structures
    multi = new TMulti();
//    pm = multi->GetPM();
//    profil = new TProfil( multi );
//    TProfil::pm = profil;
    
  //  multi->setProfil(profil);
}

void TNode::freeMemory()
{
   datach_free();
   // CSD = 0;
   delete CSD;
   CNode = databr_free( CNode );

  delete multi;
}

// Constructor of the class instance in memory for standalone GEMIPM2K or coupled program
TNode::TNode()
{
  CSD = 0;
  CNode = 0;
  allocMemory();
//  na = this;
  dbr_file_name = "dbr_file_name";
}



TNode::~TNode()
{
   freeMemory();
   // na = 0;
}

// Extracting and packing GEM IPM results into work DATABR structure
void TNode::packDataBr()
{
 long int ii;

// set default data to DataBr
//   CNode->NodeStatusCH = NEED_GEM_AIA;
   if( pm.pNP == 0 )
    CNode->NodeStatusCH = NEED_GEM_AIA;
  else
     CNode->NodeStatusCH = NEED_GEM_SIA;

   CNode->TK = pm.TCc+C_to_K; //25
   CNode->P = pm.Pc*bar_to_Pa; //1
//   CNode->IterDone = pm.IT;
   CNode->IterDone = pm.ITF+pm.IT;   // Now complete number of FIA and IPM iterations
// values
  CNode->Vs = pm.VXc*1.e-6; // from cm3 to m3
  CNode->Gs = pm.FX;
  CNode->Hs = pm.HXc;
  CNode->IC = pm.IC;
  CNode->pH = pm.pH;
  CNode->pe = pm.pe;
//  CNode->Eh = pm.FitVar[3];  Bugfix 19.12.2006  KD
  CNode->Eh = pm.Eh;
  CNode->Ms = pm.MBX;

  // arrays
   for( ii=0; ii<CSD->nPHb; ii++ )
   {  CNode->xPH[ii] = pm.XF[ CSD->xph[ii] ];
      if( CSD->nAalp >0 )
       CNode->aPH[ii] = pm.Aalp[ CSD->xph[ii] ]*kg_to_g;
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
   {   CNode->vPS[ii] = pm.FVOL[ CSD->xph[ii] ]/m3_to_cm3;
       CNode->mPS[ii] = pm.FWGT[ CSD->xph[ii] ]/kg_to_g;
       CNode->xPA[ii] = pm.XFA[ CSD->xph[ii] ];
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
   for(long int jj=0; jj<CSD->nICb; jj++ )
   { long int   new_ndx= (ii*CSD->nICb)+jj,
           mul_ndx = ( CSD->xph[ii]*CSD->nIC )+ CSD->xic[jj];
     CNode->bPS[new_ndx] = pm.BF[ mul_ndx ];
   }
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
      CNode->xDC[ii] = pm.X[ CSD->xdc[ii] ];
      CNode->gam[ii] = pm.Gamma[ CSD->xdc[ii] ];
      CNode->dul[ii] = pm.DUL[ CSD->xdc[ii] ];// always for GEM2MT init
      CNode->dll[ii] = pm.DLL[ CSD->xdc[ii] ];// always for GEM2MT init
   }
   for( ii=0; ii<CSD->nICb; ii++ )
   {  CNode->bIC[ii] = pm.B[ CSD->xic[ii] ]; // always for GEM2MT  init
      CNode->rMB[ii] = pm.C[ CSD->xic[ii] ];
      CNode->uIC[ii] = pm.U[ CSD->xic[ii] ];
   }
}

// Unpacking work DATABR structure into MULTI
//(GEM IPM work structure): uses DATACH
//  if uPrimalSol is true then the primal solution (vectors x, gamma, IC etc.)
//  will be unpacked - as an option for PIA mode with previous GEM solution from
//  the same node.
//  If uPrimalSol = false then the primal solution data will not be unpacked
//  into the MULTI structure (AIA mode or SIA mode with primal solution retained
//    in the MULTI structure from previous IPM calculation)
void TNode::unpackDataBr( bool uPrimalSol )
{
 long int ii;

 char buf[300];
 sprintf( buf, "Node:%ld:time:%lg:dt:%lg", CNode->NodeHandle, CNode->Tm, CNode->dt );
 strncpy( pm.stkey, buf, EQ_RKLEN );
  pm.TCc = CNode->TK-C_to_K;
  pm.Tc = CNode->TK;
  pm.Pc  = CNode->P/bar_to_Pa;
  pm.VXc = CNode->Vs/1.e-6; // from cm3 to m3
  // Obligatory arrays - always unpacked!
  for( ii=0; ii<CSD->nDCb; ii++ )
  {
    pm.DUL[ CSD->xdc[ii] ] = CNode->dul[ii];
    pm.DLL[ CSD->xdc[ii] ] = CNode->dll[ii];
    if( pm.DUL[ CSD->xdc[ii] ] < pm.DLL[ CSD->xdc[ii] ] )
    {
       char buf[300];
       sprintf(buf, "Upper kinetic restrictions smolest than lower for DC&RC %-6.6s",
                         pm.SM[CSD->xdc[ii]] );
       Error("unpackDataBr", buf );
    }
  }
  for( ii=0; ii<CSD->nICb; ii++ )
  {
      pm.B[ CSD->xic[ii] ] = CNode->bIC[ii];
      if( ii < CSD->nICb-1 && pm.B[ CSD->xic[ii] ] < pa.p.DB )
      {
         char buf[300];
         sprintf(buf, "Bulk mole amounts of IC  %-6.6s is %lg",
                           pm.SB[CSD->xic[ii]], pm.B[ CSD->xic[ii] ] );
          Error("unpackDataBr", buf );
      }

  }
  for( ii=0; ii<CSD->nPHb; ii++ )
  {
    if( CSD->nAalp >0 )
        pm.Aalp[ CSD->xph[ii] ] = CNode->aPH[ii]/kg_to_g;
  }

 if( !uPrimalSol )
 {    //  Using primal solution retained in the MULTI structure instead -
    ; // the primal solution data from the DATABR structure are not unpacked
//   pm.IT = 0;
 }
 else {   // Unpacking primal solution provided in the node DATABR structure
  pm.IT = 0;
  pm.MBX = CNode->Ms;
  pm.IC = CNode->IC;
  pm.Eh = CNode->Eh;
  for( ii=0; ii<CSD->nDCb; ii++ )
  /*    pm.X[ CSD->xdc[ii] ] = */
        pm.Y[ CSD->xdc[ii] ] = CNode->xDC[ii];
  for( ii=0; ii<CSD->nDCb; ii++ )
  {
     pm.lnGam[ CSD->xdc[ii] ] = log( CNode->gam[ii] );
  }
  for( ii=0; ii<CSD->nPSb; ii++ )
   pm.FVOL[ CSD->xph[ii] ] = CNode->vPS[ii]*m3_to_cm3;
  for( ii=0; ii<CSD->nPSb; ii++ )
   pm.FWGT[ CSD->xph[ii] ] = CNode->mPS[ii]*kg_to_g;

  for( ii=0; ii<CSD->nPHb; ii++ )
  {
    pm.XF[ CSD->xph[ii] ] =
    pm.YF[ CSD->xph[ii] ] = CNode->xPH[ii];
  }

  for( long int k=0; k<CSD->nPSb; k++ )
  for(long int i=0; i<CSD->nICb; i++ )
  { long int dbr_ndx= (k*CSD->nICb)+i,
          mul_ndx = ( CSD->xph[k]*CSD->nIC )+ CSD->xic[i];
    pm.BF[ mul_ndx ] = CNode->bPS[dbr_ndx];
  }

  for( ii=0; ii<CSD->nPSb; ii++ )
   pm.XFA[ CSD->xph[ii] ] = pm.YFA[ CSD->xph[ii] ] = CNode->xPA[ii];

  for( ii=0; ii<CSD->nICb; ii++ )
   pm.C[ CSD->xic[ii] ] = CNode->rMB[ii];
  for( ii=0; ii<CSD->nICb; ii++ )
   pm.U[ CSD->xic[ii] ] = CNode->uIC[ii];
 }
//  End
}


// (3) Writes the contents of the work instance of the DATABR structure into a disk file with path name  fname.
//   Parameters:
//   fname         null-terminated (C) string containing a full path to the DBR disk file to be written.
//                 NULL  - the disk file name path stored in the  dbr_file_name  field of the TNode class instance
//                 will be used, extended with ".out".  Usually the dbr_file_name field contains the path to the last input DBR file.
//   binary_f      defines if the file is to be written in binary format (true or 1, good for interruption of coupled modeling task
//                 if called in the loop for each node), or in text format (false or 0, default).
//   with_comments (text format only): defines the mode of output of comments written before each data tag and  content
//                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
//                 If   false (0), comments will not be written.
//  brief_mode     if true, tells that do not write data items,  that contain only default values in text format
void  TNode::GEM_write_dbr( const char* fname, bool binary_f, bool with_comments, bool brief_mode )
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
         databr_to_text_file(out_br, with_comments, brief_mode, str_file.c_str() );
      }
   }

// (4) Produces a formatted text file with detailed contents (scalars and arrays) of the GEM IPM work structure.
// This call is useful when GEM_run() returns with a NodeStatusCH value indicating a GEM calculation error
// (see  above).  Another use is for a detailed comparison of a test system calculation after the version upgrade of GEMIPM2K.
// Parameters: fname   null-terminated (C) string containing a full path to the disk file to be written.
//                     NULL  - the disk file name path stored in the  dbr_file_name  field of the TNode class instance will be used,
//                     extended with ".dump.out".  Usually the dbr_file_name field contains the path to the last input DBR file.
   void  TNode::GEM_print_ipm( const char* fname )
   {
     gstring str_file;
     if( fname == 0)
    	   str_file = dbr_file_name + ".Dump.out";
     else
           str_file = fname;

           outMultiTxt( str_file.c_str()  );
   }


// (9) Optional, for passing the current mass transport iteration information into the work
// DATABR structure (e.g. for using it in tracing/debugging or in writing DBR files for nodes)
void TNode::GEM_set_MT(
   double p_Tm,      // actual total simulation time, s          +       -      -
   double p_dt       // actual time step, s                      +       -      -
)
{
  CNode->Tm = p_Tm;
  CNode->dt = p_dt;
}

// (6) Passes (copies) the GEMIPM2K input data from the work instance of DATABR structure.
//  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
//   chemical data for all nodes used in the mass-transport model.
void TNode::GEM_restore_MT(
    long int  &p_NodeHandle,   // Node identification handle
    long int  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                      //                                    GEM input output  FMT control
    double &p_TK,      // Temperature T, Kelvin                       +       -      -
    double &p_P,      // Pressure P,  Pa                              +       -      -
    double &p_Vs,     // Volume V of reactive subsystem,  m3         (+)      -      +
    double &p_Ms,     // Mass of reactive subsystem, kg               -       -      +
    double *p_bIC,    // Bulk mole amounts of IC  [nICb]              +       -      -
    double *p_dul,    // Upper restrictions to amounts of DC [nDCb]   +       -      -
    double *p_dll,    // Lower restrictions to amounts of DC [nDCb]   +       -      -
    double *p_aPH     // Specific surface areas of phases,m2/kg[nPHb] +       -      -
   )
{
  long int ii;
  p_NodeHandle = CNode->NodeHandle;
  p_NodeStatusCH = CNode->NodeStatusCH;
  p_TK = CNode->TK;
  p_P = CNode->P;
  p_Vs = CNode->Vs;
  p_Ms = CNode->Ms;
// Checking if no-LPP IA is Ok
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

// (7)  Retrieves the GEMIPM2 chemical speciation calculation results from the work DATABR structure instance
//   into memory provided by the mass transport part. Dimensions and order of elements in the arrays must correspond
//   to those in currently existing DATACH memory structure.
void TNode::GEM_to_MT(
   long int &p_NodeHandle,    // Node identification handle
   long int &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
   long int &p_IterDone,      // Number of iterations performed in the last GEM IPM calculation
                         //                                                  GEM input output  FMT control
    // Chemical scalar variables
    double &p_Vs,    // Total volume V of reactive subsystem at given P,T, m3    -      -      +     +
    double &p_Ms,    // Total mass of the reactive subsystem, kg                 -      -      +     +
    double &p_Gs,    // Total Gibbs energy of the reactive subsystem, J          -      -      +     +
    double &p_Hs,    // Total enthalpy of reactive subsystem, J (reserved)       -      -      +     +
    double &p_IC,    // Effective aqueous ionic strength, molal                  -      -      +     +
    double &p_pH,    // pH of aqueous solution                                   -      -      +     +
    double &p_pe,    // pe of aqueous solution                                   -      -      +     +
    double &p_Eh,    // Eh of aqueous solution, V                                -      -      +     +
    // Dynamic data - dimensions see in DATACH.H structure
    double  *p_rMB,  // Mole balance residuals for Independent Components [nICb] -      -       +     +
    double  *p_uIC,  // Dual solution: IC chemical potentials, mol/mol [nICb]    -      -       +     +
    double  *p_xDC,  // Primal solution: DC mole amounts  [nDCb]                 -      -       +     +
    double  *p_gam,  // External activity coefficients of DC [nDCb]              -      -       +     +
    double  *p_xPH,  // Total mole amounts of all phases [nPHb]                  -      -       +     +
    double  *p_vPS,  // Total volumes of multicomponent phases, m3   [nPSb]      -      -       +     +
    double  *p_mPS,  // Total mass of multicomponent phase (carrier),kg [nPSb]   -      -       +     +
    double  *p_bPS,  // Bulk compositions of phases  [nPSb][nICb]                -      -       +     +
    double  *p_xPA,   //Amount of carrier in a multicomponent asymmetric phase[nPSb]-    -       +     +
    double  *p_aPH   //surface area for                                 phase[nPSb]-    -       +     +
 )
{
   long int ii;
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
  {
    p_xPH[ii] = CNode->xPH[ii];
    p_aPH[ii] = CNode->aPH[ii];
  }
  for( ii=0; ii<CSD->nPSb; ii++ )
  {
    p_vPS[ii] = CNode->vPS[ii];
    p_mPS[ii] = CNode->mPS[ii];
    p_xPA[ii] = CNode->xPA[ii];
  }
  for( ii=0; ii<CSD->nPSb*CSD->nICb; ii++ )
    p_bPS[ii] = CNode->bPS[ii];
}

// (8) Loads the GEMIPM2K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(
  long int  p_NodeHandle,   // Node identification handle
  long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                    //                                              GEM input output  FMT control
  double p_TK,     // Temperature T, Kelvin                            +       -      -
  double p_P,      // Pressure P, Pa                                   +       -      -
  double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
  double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
  double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
  double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
  double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
  double *p_aPH    // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
 )
 {
     long int ii;
     bool useSimplex = false;

     CNode->NodeHandle = p_NodeHandle;
     CNode->NodeStatusCH = p_NodeStatusCH;
     CNode->TK = p_TK;
     CNode->P = p_P;
     CNode->Vs = p_Vs;
     CNode->Ms = p_Ms;
   // Checking if no-LPP IA is Ok
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
      if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
      // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
}

//(8a) Loads the GEMIPM2K input data for a given mass-transport node into the work instance of DATABR structure.
//This overloaded variant uses the xDC speciation vector for setting the
// new bulk chemical composition to be used in the next GEM_run() calculation.
void TNode::GEM_from_MT(
 long int  p_NodeHandle,   // Node identification handle
 long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                  //                                              GEM input output  FMT control
 double p_TK,     // Temperature T, Kelvin                            +       -      -
 double p_P,      // Pressure P, Pa                                   +       -      -
 double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
 double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
 double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
 double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
 double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
 double *p_aPH,    // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
 double *p_xDC    // Mole amounts of DCs [nDCb] - will be convoluted
                  // and added to the bIC GEM input vector (if full speciation
                  // and not just increments then p_bIC vector must be zeroed off -
                  // it will be calculated from p_xDC and stoichiometry matrix A
)
{
  long int ii;
  bool useSimplex = false;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->TK = p_TK;
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
   if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
     CNode->NodeStatusCH = NEED_GEM_AIA;
   // Switch only if SIA is ordered, leave if simplex is ordered (KD)

   // Optional part - convolution of xDC vector into bIC vector
   if( p_xDC )
   {  long int jj;
      // Correction of bIC vector by convoluting the amounts of DCs
      for( jj=0; jj<CSD->nDCb; jj++ )
        if( p_xDC[jj] )
          for( ii=0; ii<CSD->nICb; ii++ )
            CNode->bIC[ii] += p_xDC[jj] * DCaJI( jj, ii );
   }
}

//(8b) Loads the GEMIPM2K input data for a given mass-transport node into the work instance of DATABR structure.
//In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
//  structure contains a full list of Dependent Components used in GEM IPM2 calculations.
void TNode::GEM_from_MT(
 long int  p_NodeHandle,   // Node identification handle
 long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                  //                                              GEM input output  FMT control
 double p_TK,     // Temperature T, Kelvin                            +       -      -
 double p_P,      // Pressure P, Pa                                   +       -      -
 double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
 double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
 double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
 double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
 double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
 double *p_aPH,   // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
 double *p_xDC,  // Mole amounts of DCs [nDCb] - old primal soln.     +      -      -
 double *p_gam   // DC activity coefficients [nDCb] - old primal s.   +      -      -
)
{
  long int ii;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->TK = p_TK;
  CNode->P = p_P;
  CNode->Vs = p_Vs;
  CNode->Ms = p_Ms;
// Checking if no-LPP IA is Ok
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
   else if( CNode->NodeStatusCH == NEED_GEM_SIA )
            CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal solution provided!

//  Discuss the policy!
//   if( p_xDC )
//   {  long int jj;
//      // Correction of bIC vector by convoluting the amounts of DCs
//      for( jj=0; jj<CSD->nDCb; jj++ )
//        if( p_xDC[jj] )
//          for( ii=0; ii<CSD->nICb; ii++ )
//            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
//   }

}

//-----------------------End of node.cpp--------------------------



