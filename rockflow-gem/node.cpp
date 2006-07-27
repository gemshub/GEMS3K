//-------------------------------------------------------------------
// $Id: node.cpp 684 2005-11-23 13:17:15Z gems $
//
// C/C++  interface between GEM IPM and FMT node array
// Working whith DATACH and DATABR structures
//
// Copyright (C) 2004-2005 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#include "node.h"
#include "gdatastream.h"
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

void  TNode::check_TP()
{
   bool ok = true;
   double T_ = CNode->T, P_ = CNode->P;

   if( CNode->T-C_to_K < CSD->Tval[0] )
   { ok = false;
     CNode->T = C_to_K + CSD->Tval[0];
   }
   if( CNode->T-C_to_K > CSD->Tval[CSD->nTp-1] )
   { ok = false;
     CNode->T = C_to_K + CSD->Tval[CSD->nTp-1];
   }

  if( !ok )
  {
    fstream f_log("ipmlog.txt", ios::out|ios::app );
    f_log << "In node "<< CNode->NodeHandle << "  Given T= "<<  T_ <<
             "  is beyond the range for thermodynamic data;" <<
             " set to T= " << CNode->T << endl;
  }

  ok = true;
  if( CNode->P < CSD->Pval[0] )
  { ok = false;
    CNode->P = CSD->Pval[0];
  }
  if( CNode->P > CSD->Pval[CSD->nPp-1] )
  { ok = false;
    CNode->P = CSD->Pval[CSD->nPp-1];
  }

  if( !ok )
  {
    fstream f_log("ipmlog.txt", ios::out|ios::app );
    f_log << "In node "<< CNode->NodeHandle << "  Given P= "<<  P_ <<
           "  is beyond the range for thermodynamic data;" <<
           " set to P= " << CNode->P << endl;
  }
}
//-------------------------------------------------------------------------
// GEM_run()
// GEM IPM calculation of equilibrium state for the current node.
// Mode - mode of GEMS calculation
//
//  Function returns: NodeStatus codes GEMS
//   ( OK; GEMIPM2K calculation error; system error )
//
//-------------------------------------------------------------------
int  TNode::GEM_run()
{
//  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
// f_log << " MAIF_CALC begin Mode= " << p_NodeStatusCH endl;
//---------------------------------------------
// Checking T and P
   check_TP();
// Unpacking work DATABR structure into MULTI (GEM IPM work structure): uses DATACH
   unpackDataBr();
// set up Mode
//   CNode->NodeStatusCH = (short)Mode;
// GEM IPM calculation of equilibrium state in MULTI
    TProfil::pm->calcMulti();
// Extracting and packing GEM IPM results into work DATABR structure
    packDataBr();

/**************************************************************
// only for testing output results for files
    GEM_printf( "calc_multi.ipm", "calculated_dbr.dat", "calculated.dbr" );
********************************************************* */

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
     f_log << "gems2: Unknown exception: program aborted" << endl;
       if( CNode->NodeStatusCH  == NEED_GEM_AIA )
         CNode->NodeStatusCH = ERR_GEM_AIA;
       else
         CNode->NodeStatusCH = ERR_GEM_PIA;
    }
   return CNode->NodeStatusCH;
}

 // reads work node (DATABR structure) from a  file
int  TNode::GEM_read_dbr( bool binary_f, char *fname )
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
  } catch(TError& err)
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
// reads in the data from MULTI, DATACH, DATABR files prepared
// using the GEMS-PSI RMT module
//  Parameters:
//  ipmfiles_lst_name - name of a text file that contains:
//    " -t/-b <MULTI file name> -t/-b <dataCH file name>,
//    <dataBR file name1>, ...  <dataBR file nameN> "
//    These files (one MULTI & dataCH file, at least one dataBR file) must exist in
//    the current directory; the dataBR files in the above list are indexed
//    as 1, 2, ... N (node types) and must contain valid initial chemical
//    systems (of the same structure described in the dataCH file) to set up
//    the initial state of the FMT node array. If -t flag is specified
//    then dataCH and dataBR files must be in text (ASCII) format;
//    if -b or nothing is specified then dataCH and dataBR files are
//    assumed to be binary (little-endian) files.
//  nodeTypes[nNodes] - array of node type (fortran) indexes of dataBR files in
//    the ipmfiles_lst_name list. This array, for each FMT node, specifies
//    from which dataBR file the initial chemical system should be taken.
//  Function returns:
//   0: OK; 1: GEMIPM read file error; -1: System error (e.g. memory allocation)
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
     bool binary_mult = true;
     bool binary_f = true;
     gstring chbr_in = ipmfiles_lst_name;

// Get path
      size_t pos = chbr_in.rfind("/");
      gstring Path = "";
      if( pos < npos )
       Path = chbr_in.substr(0, pos+1);

//  open file stream and reading name of MULTI file
//  -t/-b  "<MULTI file name>"
      fstream f_chbr(chbr_in.c_str(), ios::in );
      ErrorIf( !f_chbr.good() , chbr_in.c_str(), "Fileopen error");

      gstring datachbr_file;
      f_getline( f_chbr, datachbr_file, ' ');

//Testing flag "-t" or "-b" (by default "-b")   // use bynary or text files for Multi
      pos = datachbr_file.find( '-');
      if( pos != /*gstring::*/npos )
      {
         if( datachbr_file[pos+1] == 't' )
            binary_mult = false;
         f_getline( f_chbr, datachbr_file, ' ');
      }

    // Reading structure MULTI (GEM IPM work structure)
    gstring multu_in = Path + datachbr_file;

// Reading name of dataCH file and names of dataBR files
//  -t/-b  "<dataCH file name>" ,"<dataBR file1 name>", ..., "<dataBR fileN name>"
      f_getline( f_chbr, datachbr_file, ' ');

//Testing flag "-t" or "-b" (by default "-b")   // use bynary or text files as input
      pos = datachbr_file.find( '-');
      if( pos != /*gstring::*/npos )
      {
         if( datachbr_file[pos+1] == 't' )
            binary_f = false;
         f_getline( f_chbr, datachbr_file, ',');
      }

// Reading dataCH structure from file
     gstring dat_ch = Path + datachbr_file;
      if( binary_f )
      {  GemDataStream f_ch(dat_ch, ios::in|ios::binary);
         datach_from_file(f_ch);
       }
      else
      { fstream f_ch(dat_ch.c_str(), ios::in );
         ErrorIf( !f_ch.good() , dat_ch.c_str(), "DataCH Fileopen error");
         datach_from_text_file(f_ch);
      }

     i = 0;
     while( !f_chbr.eof() )  // For all dataBR files listed
     {

#ifndef IPMGEMPLUGIN
   pVisor->Message( 0, "GEM2MT node array",
      "Reading from disk a set of node array files to resume an interrupted RMT task.\n"
           "Please, wait...", i, nNodes() );
#endif

// Reading work dataBR structure from file
         f_getline( f_chbr, datachbr_file, ',');

         gstring dbr_file = Path + datachbr_file;
         if( binary_f )
         {
             GemDataStream in_br(dbr_file, ios::in|ios::binary);
             databr_from_file(in_br);
          }
         else
          {   fstream in_br(dbr_file.c_str(), ios::in );
                 ErrorIf( !in_br.good() , datachbr_file.c_str(),
                    "DataBR Fileopen error");
               databr_from_text_file(in_br);
          }

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
     }
#ifndef IPMGEMPLUGIN
   pVisor->CloseMessage();
#endif

    ErrorIf( i==0, datachbr_file.c_str(), "GEM_init() error: No dataBR files read!" );
    checkNodeArray( i, nodeTypes, datachbr_file.c_str()  );

// Reading structure MULTI (GEM IPM work structure)
if( binary_mult )
 {
   GemDataStream f_m(multu_in, ios::in|ios::binary);
#ifdef IPMGEMPLUGIN
    profil->readMulti(f_m);
#else
    TProfil::pm->readMulti(f_m);
#endif
  }
  else
  {
#ifdef IPMGEMPLUGIN
        profil->readMulti(multu_in.c_str());
#endif
  }
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

// Return DCH index of IC by Name or -1 if illegal name
int TNode::IC_name_to_x( const char *Name )
{
  uint len = strlen( Name );

  for(int ii = 0; ii<CSD->nIC; ii++ )
       if(!memcmp(Name, CSD->ICNL[ii], min(len,MaxICN)))
         return ii;
  return -1;
}

// Return DCH index of DC by Name or -1 if illegal name
int TNode::DC_name_to_x( const char *Name )
{
  uint len = strlen( Name );

  for(int ii = 0; ii<CSD->nDC; ii++ )
       if(!memcmp(Name, CSD->DCNL[ii], min(len,MaxDCN)))
         return ii;
  return -1;
}

// Return DCH index of Ph by Name or -1 if illegal name
int TNode::Ph_name_to_x( const char *Name )
{
  uint len = strlen( Name );

  for(int ii = 0; ii<CSD->nPH; ii++ )
       if(!memcmp(Name, CSD->PHNL[ii], min(len,MaxPHN)))
         return ii;
  return -1;
}

// Return for IComp DBR index from DCH index or -1 if not used in the data bridge
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

//---------------------------------------------------------//

void TNode::allocMemory()
{
// alloc memory for data bridge structures
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
  delete[] multi;
  delete[] profil;
#endif
}

#ifndef IPMGEMPLUGIN

// Make start DATACH and DATABR data from GEMS internal data (MULTI and other)
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
// realloc memory for     DATACH  *CSD;  and  DATABR  *CNode;

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

// These dimensionalities define sizes of dynamic data in DATABT structure!!!

  CSD->nICb = (short)selIC.GetCount();
  CSD->nDCb = (short)selDC.GetCount();
  CSD->nPHb = (short)selPH.GetCount();
  CSD->nPSb = 0;
  for( ii=0; ii< selPH.GetCount(); ii++, CSD->nPSb++ )
   if( selPH[ii] >= pmm->FIs )
       break;
  CSD->uRes2 = 0;
  CSD->dRes1 = 0.;
  CSD->dRes2 = 0.;

  CSD->Ttol = Ttol_;
  CSD->Ptol = Ptol_;

// realloc structures DataCh&DataBr

  datach_realloc();
  databr_realloc();

// set dynamic data to DataCH

  memcpy( CSD->nDCinPH, pmm->L1 , CSD->nPH*sizeof(short));
  for( ii=0; ii< selIC.GetCount(); ii++ )
    CSD->xIC[ii] = (short)selIC[ii];
  for( ii=0; ii< selDC.GetCount(); ii++ )
    CSD->xDC[ii] = (short)selDC[ii];
  for( ii=0; ii< selPH.GetCount(); ii++ )
    CSD->xPH[ii] = (short)selPH[ii];

  memcpy( CSD->A, pmm->A , CSD->nIC*CSD->nDC*sizeof(float));

  for( i1=0; i1< CSD->nIC; i1++ )
     CSD->ICmm[i1] = pmm->Awt[i1];

  memcpy( CSD->DCmm, pmm->MM , CSD->nDC*sizeof(double));
  memset( CSD->DD, 0, CSD->nDCs*sizeof(double));

  for( ii=0; ii<CSD->nIC; ii++ )
     memcpy( CSD->ICNL[ii], pmm->SB[ii] , MaxICN*sizeof(char));
  memcpy( CSD->DCNL, pmm->SM , MaxDCN*CSD->nDC*sizeof(char));
  for( ii=0; ii< CSD->nPH; ii++ )
    memcpy( CSD->PHNL[ii], pmm->SF[ii]+4 , MaxPHN*sizeof(char));


  memcpy( CSD->ccIC, pmm->ICC , CSD->nIC*sizeof(char));
  memcpy( CSD->ccDC, pmm->DCC , CSD->nDC*sizeof(char));
  memcpy( CSD->ccPH, pmm->PHC , CSD->nPH*sizeof(char));

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

   memset( &CNode->T, 0, 32*sizeof(double));
   CNode->T = pmm->Tc; //25
   CNode->P = pmm->Pc; //1
   CNode->Ms = pmm->MBX;

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

// set calculated&dynamic data to DataBR

   packDataBr();

// must be changed to matrix structure  ???????
// setted CSD->nPp*CSD->nTp = 1
   for( i1=0; i1<CSD->nTp; i1++ )
    CSD->Tval[i1] = Tai[i1];
   for( i1=0; i1<CSD->nPp; i1++ )
    CSD->Pval[i1] = Pai[i1];

   getG0_V0_H0_Cp0_matrix();

}

void TNode::getG0_V0_H0_Cp0_matrix()
{

  double cT, cP/*, cDC*/;
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
    cT = CSD->Tval[ii];
    for( int jj=0; jj<CSD->nPp; jj++)
    {
      cP = CSD->Pval[jj];
     // calc new G0, V0, H0, Cp0, S0
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
}

#else

TNode::TNode()
{
  CSD = 0;
  CNode = 0;
  allocMemory();
  na = this;
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
#else

 // numbers
  if( pmm->pNP == 0 )
    CNode->NodeStatusCH = OK_GEM_AIA;
  else
    CNode->NodeStatusCH = OK_GEM_PIA;

#endif

   CNode->T = pmm->Tc; //25
   CNode->P = pmm->Pc; //1
   CNode->IterDone = pmm->IT;

// values
  CNode->Vs = pmm->VXc;
  CNode->Gs = pmm->FX;
  CNode->Hs = pmm->HXc;
  CNode->IC = pmm->IC;
  CNode->pH = pmm->pH;
  CNode->pe = pmm->pe;
  CNode->Eh = pmm->FitVar[3];
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
void TNode::unpackDataBr()
{
 short ii;
 double Gamm;
// numbers

  if( CNode->NodeStatusCH >= NEED_GEM_PIA )
   pmm->pNP = 1;
  else
   pmm->pNP = 0; //  NEED_GEM_AIA;
  CNode->IterDone = 0;
  pmm->IT = 0;
// values
  pmm->Tc = CNode->T;
  pmm->Pc  = CNode->P;
  pmm->MBX = CNode->Ms;
  pmm->IC = CNode->IC;
  pmm->FitVar[3] = CNode->Eh;
// arrays
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
    pmm->DUL[ CSD->xDC[ii] ] = CNode->dul[ii];
    pmm->DLL[ CSD->xDC[ii] ] = CNode->dll[ii];
   }
   for( ii=0; ii<CSD->nICb; ii++ )
    pmm->B[ CSD->xIC[ii] ] = CNode->bIC[ii];

// added  to compare SD 15/07/04
   for( ii=0; ii<CSD->nDCb; ii++ )
    pmm->X[ CSD->xDC[ii] ] = pmm->Y[ CSD->xDC[ii] ] = CNode->xDC[ii];
   for( ii=0; ii<CSD->nDCb; ii++ )
   {
      Gamm = CNode->gam[ii];
      pmm->Gamma[ CSD->xDC[ii] ] = Gamm;
      pmm->lnGmo[ CSD->xDC[ii] ] = pmm->lnGam[ CSD->xDC[ii] ] = log(Gamm);
   }
   for( ii=0; ii<CSD->nPHb; ii++ )
   {
     pmm->XF[ CSD->xPH[ii] ] = pmm->YF[ CSD->xPH[ii] ] = CNode->xPH[ii];
     if( CSD->nAalp >0 )
          pmm->Aalp[ CSD->xPH[ii] ] = CNode->aPH[ii];
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
    pmm->FVOL[ CSD->xPH[ii] ] = CNode->vPS[ii];
   for( ii=0; ii<CSD->nPSb; ii++ )
    pmm->FWGT[ CSD->xPH[ii] ] = CNode->mPS[ii];

   for( ii=0; ii<CSD->nPSb; ii++ )
   for(short jj=0; jj<CSD->nICb; jj++ )
   { int new_ndx= (ii*CSD->nICb)+jj,
           mul_ndx = ( CSD->xPH[ii]*CSD->nIC )+ CSD->xIC[jj];
     pmm->BF[ mul_ndx ] = CNode->bPS[new_ndx];
   }
   for( ii=0; ii<CSD->nPSb; ii++ )
    pmm->XFA[ CSD->xPH[ii] ] = pmm->YFA[ CSD->xPH[ii] ] = CNode->xPA[ii];

   for( ii=0; ii<CSD->nICb; ii++ )
    pmm->C[ CSD->xIC[ii] ] = CNode->rMB[ii];
   for( ii=0; ii<CSD->nICb; ii++ )
    pmm->U[ CSD->xIC[ii] ] = CNode->uIC[ii];

}

void  TNode::GEM_printf( const char* multi_file,
                             const char* databr_text,
                             const char* databr_bin )
{
//**************************************************************
// only for testing output results for files
// binary DATABR
    gstring strr;
   if( databr_bin )
   {  strr = databr_bin;
      GemDataStream out_br(strr, ios::out|ios::binary);
      databr_to_file(out_br);
   }
// text DATABR
   if( databr_text )
   {  fstream out_br_t(databr_text, ios::out );
      ErrorIf( !out_br_t.good() , databr_text,
                "DataBR text file open error");
      databr_to_text_file(out_br_t);
   }
// output multy
    if( multi_file )
   {  strr = multi_file;
      GemDataStream o_m( strr, ios::out|ios::binary);
       TProfil::pm->outMulti(o_m, strr );
    }
//********************************************************* */

}

#ifdef IPMGEMPLUGIN

// calculation mode: passing input GEM data changed on previous FMT iteration
//                   into work DATABR structure
void TNode::GEM_from_MT(
   short  p_NodeHandle,   // Node identification handle
   short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double p_T,      // Temperature T, K                        +       -      -
   double p_P,      // Pressure P, bar                         +       -      -
   double p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
   double p_Ms,     // Mass of reactive subsystem, kg          -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
   double *p_dul,   // upper kinetic restrictions [nDCb]        +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]        +       -      -
   double *p_aPH  // Specific surface areas of phases (m2/g)      +       -     -

)
{
  int ii;
  bool useSimplex = false;

  CNode->NodeHandle = p_NodeHandle;
  CNode->NodeStatusCH = p_NodeStatusCH;
  CNode->T = p_T;
  CNode->P = p_P;
  CNode->Vs = p_Vs;
  CNode->Ms = p_Ms;
// Checking if no-simplex IA is Ok
   for( ii=0; ii<CSD->nICb; ii++ )
   {  //  Sveta 11/02/05 for test
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
   double &p_T,      // Temperature T, K                        +       -      -
   double &p_P,      // Pressure P, bar                         +       -      -
   double &p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
   double &p_Ms,     // Mass of reactive subsystem, kg          -       -      +
   double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
   double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
   double *p_aPH  // Specific surface areas of phases (m2/g)      +       -     -
   )
{
 int ii;
  p_NodeHandle = CNode->NodeHandle;
  p_NodeStatusCH = CNode->NodeStatusCH;
  p_T = CNode->T;
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

  memcpy( p_xDC, CNode->xDC, CSD->nDCb*sizeof(double) );
  memcpy( p_gam, CNode->gam, CSD->nDCb*sizeof(double) );
  memcpy( p_xPH, CNode->xPH, CSD->nPHb*sizeof(double) );
  memcpy( p_vPS, CNode->vPS, CSD->nPSb*sizeof(double) );
  memcpy( p_mPS, CNode->mPS, CSD->nPSb*sizeof(double) );
  memcpy( p_bPS, CNode->bPS, CSD->nPSb*CSD->nICb*sizeof(double) );
  memcpy( p_xPA, CNode->xPA, CSD->nPSb*sizeof(double) );
  memcpy( p_rMB, CNode->rMB, CSD->nICb*sizeof(double) );
  memcpy( p_uIC, CNode->uIC, CSD->nICb*sizeof(double) );
}

#endif

//-----------------------End of node.cpp--------------------------



