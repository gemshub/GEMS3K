//-------------------------------------------------------------------
// $Id: main.cpp 792 2006-09-19 08:10:41Z gems $
//
// Debugging version of a finite-difference 1D advection-diffusion
// mass transport model supplied by Dr. Frieder Enzmann (Uni Mainz)
// coupled with GEMIPM2K module for calculation of chemical equilibria
//
// Direct access to the TNodeArray class for storing all data for nodes
//
// Copyright (C) 2005,2007 S.Dmytriyeva, F.Enzmann, D.Kulik
//
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <iomanip>
#include "m_gem2mt.h"
#include "io_keyvalue.h"
#include "io_simdjson.h"

#include <dirent.h>
#include <sys/stat.h>

bool DirExists( const char* aPath )
{
    if ( aPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (aPath);

    if (pDir != NULL)
    {
        bExists = true;
        (void) closedir(pDir);
    }

    return bExists;
}


TGEM2MT* TGEM2MT::pm;

TGEM2MT::TGEM2MT( uint /*nrt*/ )
{
  mtp=&mt[0];
  set_def( 0 );
  ////mtp->PvMO =   S_ON;
  ////mtp->iStat =  AS_READY;
  na = 0;
  pa_mt = 0;
}

TGEM2MT::~TGEM2MT()
{
  mem_kill(0);
  if( na )
   delete na;
  if( pa_mt )
    delete pa_mt;
}

//=======================================================================================

//Calculate record
void TGEM2MT::RecCalc()
{
    try
    {
       bool iRet;

      if( mtp->PsVTK != S_OFF )
      {
         if( !DirExists( pathVTK.c_str() ) )
#ifndef  __unix
            mkdir( pathVTK.c_str() );
#else
             mkdir( pathVTK.c_str(), 0755 );
#endif
       // vfMakeDirectory(window(), pathVTK.c_str() );
      }

      if( mtp->iStat != AS_RUN  )
      {
        mtp->gStat = GS_GOING;
        mt_reset();
        mtp->gStat = GS_DONE;
      }

      // internal calc
      iRet = internalCalc();
      if(!iRet)
            mtp->iStat = AS_DONE;
      //else we have a stop point

     }
     catch( TError& xcpt )
     {
        mtp->gStat = GS_ERR;
        mtp->iStat = AS_INDEF;
        Error(  xcpt.title.c_str(), xcpt.mess.c_str() );
     }
}


// read TGEM2MT structure from file
int TGEM2MT::ReadTask( const char *gem2mt_in1, const char *vtk_dir )
{

 // read GEM2MT structure from file
  try
  {
   std::string gem2mt_in = gem2mt_in1;
   std::fstream ff(gem2mt_in, std::ios::in );
   ErrorIf( !ff.good() , gem2mt_in, "Fileopen error");

   if( gem2mt_in.rfind(".json") != std::string::npos )
   {
       io_formats::SimdJsonRead in_format( ff, "", "gem2mt" );
       from_text_file( in_format );
   }
   else
   {
       io_formats::KeyValueRead in_format( ff );
       from_text_file( in_format );
   }

   pathVTK = vtk_dir;
   if( !pathVTK.empty() )
   {
       if(  DirExists( pathVTK.c_str() ) )
          pathVTK += "/";
       else pathVTK =""; // this directory do not exist
   }
   return 0;
  }
  catch(TError& err)
  {
      std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
  }
  return 1;
}

// Write TGEM2MT structure to file
int TGEM2MT::WriteTask( const char *gem2mt_out1 )
{
 // write GEM2MT structure to file
  try
  {
   std::string gem2mt_out = gem2mt_out1;
   std::fstream ff(gem2mt_out, std::ios::out );
   ErrorIf( !ff.good() , gem2mt_out, "Fileopen error");

   if( gem2mt_out.rfind(".json") != std::string::npos )
   {
       io_formats::SimdJsonWrite out_format( ff, "", true );
       to_text_file( out_format, true, false );
   }
   else
   {
       io_formats::KeyValueWrite out_format( ff );
       to_text_file( out_format, true, false );
   }

   return 0;
  }
  catch(TError& err)
  {
      std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
      f_log << err.title.c_str() << "  : " << err.mess.c_str() <<std:: endl;
  }
  return 1;
}


// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up NodeArray and ParticleArray classes
int TGEM2MT::MassTransInit( const char *lst_f_name, const char *dbr_lst_f_name )
{
  int ii;

  // define name of vtk file
  std::string lst_in = lst_f_name;
  size_t pos = lst_in.rfind("\\");
  size_t pos2 = lst_in.rfind("/");
  if( pos == std::string::npos )
      pos = pos2;
  else
      if( pos2 < std::string::npos)
         pos = std::max(pos, pos2 );
  if( pos < std::string::npos )
  {
      if( pathVTK.empty() )
      {
          pathVTK = lst_in.substr(0, pos+1);
          pathVTK += "VTK/";
      }
      lst_in = lst_in.substr(pos+1);
  }
  pos = lst_in.find(".");
  lst_in = lst_in.substr(0, pos);
  pos = lst_in.find("-");
  lst_in = lst_in.substr(0, pos);
  nameVTK = lst_in;
  prefixVTK = lst_in;


  // The NodeArray must be allocated here
  TNodeArray::na = na = new TNodeArray( /* mtp->xC,mtp->yC,mtp->zC */ mtp->nC );

 // Prepare the array for initial conditions allocation
  long int* nodeType = new long int[mtp->nC];
  for( ii =0; ii<mtp->nC; ii++ )
         nodeType[ii] = mtp->DiCp[ii][0];

  // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
  // if mtp->iStat == AS_RUN we resume calculation
  if( na->GEM_init( lst_f_name, dbr_lst_f_name, nodeType, mtp->iStat == AS_RUN ) )
        return 1;  // error reading files

  CalcIPM( NEED_GEM_AIA, 0, mtp->nC, 0 ); //recalc all nodes ?

  for( ii=0; ii< mtp->nTai; ii++)
      mtp->Tval[ii] =  na->pCSD()->TKval[ii]-C_to_K;
  for( ii=0; ii< mtp->nPai; ii++)
      mtp->Pval[ii] =  na->pCSD()->Pval[ii]/bar_to_Pa;

  // use particles
  if( mtp->PsMode == RMT_MODE_W  )
  {
   na->SetGrid( mtp->sizeLc, mtp->grid );   // set up grid structure
   pa_mt = new TParticleArray( mtp->nPTypes, mtp->nProps,
         mtp->NPmean, mtp->ParTD, mtp->nPmin, mtp->nPmax, na );
   pa_mt->setUpCounters();
  }
     // put HydP
  if( mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B )
  {
      putHydP( na->pNodT0() );
      putHydP( na->pNodT1() );
}

  delete[] nodeType;
  return 0;
}


//==========================================================================================

// free dynamic memory in objects and values
void TGEM2MT::mem_kill(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
       "E05GTrem: Attempt to access corrupted dynamic memory.");

    //- if( mtp->lNam) delete[] mtp->lNam;
    //- if( mtp->lNamE) delete[] mtp->lNamE;
    //- if( mtp->tExpr) delete[] mtp->tExpr;
    //- if( mtp->gExpr) delete[] mtp->gExpr;
    if( mtp->sdref) delete[] mtp->sdref;
    if( mtp->sdval) delete[] mtp->sdval;
    if( mtp->DiCp) delete[] mtp->DiCp;
    if( mtp->FDLi) delete[] mtp->FDLi;
    //- if( mtp->PTVm) delete[] mtp->PTVm;
    //- if( mtp->StaP) delete[] mtp->StaP;
    if( mtp->xVTKfld) delete[] mtp->xVTKfld;
    //- if( mtp->xEt) delete[] mtp->xEt;
    //- if( mtp->yEt) delete[] mtp->yEt;
    //- if( mtp->Bn) delete[] mtp->Bn;
    if( mtp->HydP) delete[] mtp->HydP;
    //- if( mtp->qpi) delete[] mtp->qpi;
    //- if( mtp->qpc) delete[] mtp->qpc;
    //- if( mtp->xt) delete[] mtp->xt;
    //- if( mtp->yt) delete[] mtp->yt;
    //- if( mtp->CIb) delete[] mtp->CIb;
    //- if( mtp->CAb) delete[] mtp->CAb;
    if( mtp->FDLf) delete[] mtp->FDLf;
    if( mtp->PGT) delete[] mtp->PGT;
    if( mtp->Tval) delete[] mtp->Tval;
    if( mtp->Pval) delete[] mtp->Pval;
    if( mtp->nam_i) delete[] mtp->nam_i;
    //- if( mtp->for_i) delete[] mtp->for_i;
    //- if( mtp->stld) delete[] mtp->stld;
    //- if( mtp->CIclb) delete[] mtp->CIclb;
    //- if( mtp->AUcln) delete[] mtp->AUcln;
    if( mtp->FDLid) delete[] mtp->FDLid;
    if( mtp->FDLop) delete[] mtp->FDLop;
    if( mtp->FDLmp) delete[] mtp->FDLmp;
    if( mtp->MGPid) delete[] mtp->MGPid;
    if( mtp->UMGP) delete[] mtp->UMGP;
    //- if( mtp->SBM) delete[] mtp->SBM;
    if( mtp->BSF) delete[] mtp->BSF;
   if( mtp->MB) delete[] mtp->MB;
   if( mtp->dMB) delete[] mtp->dMB;
   if( mtp->DDc) delete[] mtp->DDc;
   if( mtp->DIc) delete[] mtp->DIc;
   if( mtp->DEl) delete[] mtp->DEl;
   if( mtp->for_e) delete[] mtp->for_e;
   //- if( mtp->xIC) delete[] mtp->xIC;
   //- if( mtp->xDC) delete[] mtp->xDC;
   //- if( mtp->xPH) delete[] mtp->xPH;
   if( mtp->grid) delete[] mtp->grid;
   if( mtp->NPmean) delete[] mtp->NPmean;
   if( mtp->nPmin) delete[] mtp->nPmin;
   if( mtp->nPmax) delete[] mtp->nPmax;
   if( mtp->ParTD) delete[] mtp->ParTD;
if( mtp->arr1) delete[] mtp->arr1;
if( mtp->arr2) delete[] mtp->arr2;
// work
   //- if( mtp->An) delete[] mtp->An;
   //- if( mtp->Ae) delete[] mtp->Ae;
   if( mtp->gfc) delete[] mtp->gfc;
   if( mtp->yfb) delete[] mtp->yfb;
   if( mtp->tt) delete[] mtp->tt;
   //- if( mtp->etext) delete[] mtp->etext;
   //- if( mtp->tprn) delete[] mtp->tprn;
   //- FreeNa();
   //- freeNodeWork();
}

// realloc dynamic memory
void TGEM2MT::mem_new(int q)
{
  ErrorIf( mtp!=&mt[q], GetName(),
      "E04GTrem: Attempt to access corrupted dynamic memory.");

 //- mtp->xIC = new long int[mtp->nICb];
 //- mtp->xDC = new long int[mtp->nDCb];
 //- mtp->xPH = new long int[mtp->nPHb];

 if( mtp->PvGrid == S_OFF )
   { if(mtp->grid) delete[] mtp->grid;
     mtp->grid = 0;
   }
 else
     mtp->grid = new double[ mtp->nC][3];

 if( mtp->PsMode == RMT_MODE_W  )
 {
   mtp->NPmean = new long int[ mtp->nPTypes];
   mtp->nPmin = new long int[ mtp->nPTypes];
   mtp->nPmax = new long int[ mtp->nPTypes];
   mtp->ParTD = new long int[mtp->nPTypes][6];
 }
 else
 {
  if(mtp->NPmean) delete[] mtp->NPmean;
  if(mtp->nPmin) delete[] mtp->nPmin;
  if(mtp->nPmax) delete[] mtp->nPmax;
  if(mtp->ParTD) delete[] mtp->ParTD;
  mtp->NPmean = 0;
  mtp->nPmin = 0;
  mtp->nPmax = 0;
  mtp->ParTD = 0;
 }
 mtp->nam_i= new char[ mtp->nIV][ MAXIDNAME ];
 //- mtp->PTVm = new double[ mtp->nIV][5];
 mtp->DiCp = new long int[ mtp->nC][2];
 //- mtp->StaP = new double[ mtp->nC ][4];

 if( mtp->PvnVTK == S_OFF )
    { if(mtp->xVTKfld) delete[] mtp->xVTKfld;
      mtp->xVTKfld = 0;
    }
 else
     mtp->xVTKfld = new long int[ mtp->nVTKfld][2];

 //- mtp->stld = new char[ mtp->nIV ][EQ_RKLEN];
 mtp->Tval  = new double[ mtp->nTai ];
 mtp->Pval  = new double[ mtp->nPai ];
 //- mtp->Bn = new double[ mtp->nIV][ mtp->Nb ];
 //- mtp->SBM = new char [ mtp->Nb][MAXICNAME+MAXSYMB];

 if( mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B )
      mtp->HydP = new double[ mtp->nC][SIZE_HYDP];
 else
    { if(mtp->HydP) delete[] mtp->HydP;
      mtp->HydP = 0;
    }

 //-if( mtp->PvICi == S_OFF )
 //-   {
 //-    if(mtp->CIb) delete[] mtp->CIb;
 //-    if(mtp->CIclb) delete[] mtp->CIclb;
 //-    mtp->CIb = 0;
 //-    mtp->CIclb = 0;
 //-   }
 //-   else
 //-   {
 //-    mtp->CIb = new double[ mtp->nIV][mtp->Nb];
 //-    mtp->CIclb =  new char[ mtp->Nb ];
 //-   }

 //- if( mtp->PvAUi == S_OFF )
 //-    {
 //-     if(mtp->CAb) delete[] mtp->CAb;
 //-      if(mtp->for_i) delete[] mtp->for_i;
 //-      if(mtp->AUcln) delete[] mtp->AUcln;
 //-      if(mtp->An) delete[] mtp->An;
 //-      mtp->CAb = 0;
 //-      mtp->for_i = 0;
 //-      mtp->AUcln = 0;
 //-      mtp->An = 0;
 //-      mtp->Lbi = 0;
 //-    }
 //-    else
 //-    {
 //-      mtp->CAb = new double[ mtp->nIV][ mtp->Lbi ];
 //-      mtp->for_i = new char[ mtp->Lbi][ MAXFORMUNITDT ];
 //-      mtp->AUcln = new char[ mtp->Lbi ];
 //-      mtp->An = new double[ mtp->Lbi][ mtp->Nb ];
 //-   }

 if( mtp->PvFDL == S_OFF )
   {
     if(mtp->FDLi) delete[] mtp->FDLi;
     if(mtp->FDLf) delete[] mtp->FDLf;
     if(mtp->FDLid) delete[] mtp->FDLid;
     if(mtp->FDLop) delete[] mtp->FDLop;
     if(mtp->FDLmp) delete[] mtp->FDLmp;
     mtp->FDLi = 0;
     mtp->FDLf = 0;
     mtp->FDLid = 0;
     mtp->FDLop = 0;
     mtp->FDLmp = 0;
     mtp->nFD = 0;
  }
   else
   {
      mtp->FDLi = new long int[ mtp->nFD][2];
      mtp->FDLf = new double[ mtp->nFD][4];
      mtp->FDLid= new char[ mtp->nFD][MAXSYMB];
      mtp->FDLop= new char[ mtp->nFD][MAXSYMB];
      mtp->FDLmp = new char[ mtp->nFD][MAXSYMB];
   }
 if( mtp->PvPGD == S_OFF )
   {
     if(mtp->PGT) delete[] mtp->PGT;
     if(mtp->MGPid) delete[] mtp->MGPid;
     if(mtp->UMGP) delete[] mtp->UMGP;
     mtp->PGT = 0;
     mtp->MGPid = 0;
     mtp->UMGP = 0;
     mtp->nPG = 0;
   }
   else
   {
     mtp->PGT  =  new double[ mtp->FIf*mtp->nPG ];
     mtp->MGPid = new char[ mtp->nPG ][ MAXSYMB ];
     mtp->UMGP = new char[ mtp->FIf ];
   }

 if( mtp->PvSFL == S_OFF )
  { if(mtp->BSF) delete[] mtp->BSF;
    mtp->BSF = 0;
  }
   else
      mtp->BSF = new double[ mtp->nSFD*mtp->Nf ];
 if( mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF )
   {
      mtp->MB =  new double[mtp->nC*mtp->Nf];
      mtp->dMB = new double[mtp->nC*mtp->Nf];
   }
   else
   {
     if(mtp->MB) delete[] mtp->MB;
     if(mtp->dMB) delete[] mtp->dMB;
      mtp->MB = 0;
      mtp->dMB = 0;
   }
   if( mtp->PvDDc == S_OFF )
   {
     if(mtp->DDc) delete[] mtp->DDc;
     mtp->DDc = 0;
   }
   else
     mtp->DDc = new double[mtp->Lsf];

   if( mtp->PvDIc == S_OFF )
   {
     if(mtp->DIc) delete[] mtp->DIc;
     mtp->DIc = 0;
   }
   else
       mtp->DIc = new double[ mtp->Nf ];

   if( mtp->nEl <= 0  )
   {
     if(mtp->DEl) delete[] mtp->DEl;
     if(mtp->for_e) delete[] mtp->for_e;
   //-   if(mtp->Ae) delete[] mtp->Ae;
     mtp->DEl = 0;
     mtp->for_e = 0;
   //-   mtp->Ae = 0;
     mtp->nEl = 0;
   }
   else
   {
       mtp->DEl = new double[ mtp->nEl ];
       mtp->for_e = new char[mtp->nEl][MAXFORMUNITDT];
   //-     mtp->Ae = new double[ mtp->nEl*mtp->Nb ];
   }

//----------------------------------------------------------------
   //- if( mtp->Nqpt > 0  )
   //-  mtp->qpi   = new double[mtp->Nqpt];
   //- else
   //- { if(mtp->qpi) delete[] mtp->qpi; mtp->qpi = 0;}

   //- if( mtp->Nqpg > 0  )
   //-  mtp->qpc   = new double[mtp->Nqpg];
   //- else
   //- { if(mtp->qpc) delete[] mtp->qpc; mtp->qpc = 0;}

   //- if( mtp->PvMSt == S_OFF )
   //- { if(mtp->tExpr) delete[] mtp->tExpr; mtp->tExpr = 0;}
   //- else
   //-    mtp->tExpr = new char[4096];

   //-if( mtp->PvMSg == S_OFF )
   //-   {
   //-    if(mtp->lNam) delete[] mtp->lNam;
   //-    if(mtp->gExpr) delete[] mtp->gExpr;
   //-    if(mtp->xt) delete[] mtp->xt;
   //-    if(mtp->yt) delete[] mtp->yt;
   //-    mtp->lNam = 0;
   //-    mtp->gExpr = 0;
   //-    mtp->xt = 0;
   //-    mtp->yt = 0;
   //-   }
   //-   else
   //-   {
   //-        mtp->lNam = new char[ mtp->nYS][ MAXGRNAME];
   //-        mtp->gExpr = new char[2048];
   //-        mtp->xt   = new double[ mtp->nC];
   //-        mtp->yt   = new double[ mtp->nC*mtp->nYS];
   //-   }

   //-if( mtp->PvEF == S_OFF )
   //-   {
   //-    if(mtp->lNamE) delete[] mtp->lNamE;
   //-    if(mtp->xEt) delete[] mtp->xEt;
   //-    if(mtp->yEt) delete[] mtp->yEt;
   //-    mtp->lNamE = 0;
   //-    mtp->xEt = 0;
   //-    mtp->yEt = 0;
   //-   }
   //-   else
   //-   {
   //-     mtp->lNamE = new char[ mtp->nYE][ MAXGRNAME];
   //-     mtp->xEt   = new double[ mtp->nE ];
   //-     mtp->yEt   = new double[ mtp->nE*mtp->nYE ];
   //-   }

    if( mtp->Nsd > 0 )
    {
        mtp->sdref = new char[ mtp->Nsd ][ V_SD_RKLEN ];
        mtp->sdval = new char[ mtp->Nsd ][ V_SD_VALEN ];
    }
    else
    {
      if(mtp->sdref) delete[] mtp->sdref;
      if(mtp->sdval) delete[] mtp->sdval;
      mtp->sdref = 0;
      mtp->sdval = 0;
    }
    //- mtp->etext = new char[4096];
    //- mtp->tprn = new char[2048];
    //mtp->gfc = (double *)aObj[ o_mtgfc].Free();
    //mtp->yfb = (double *)aObj[ o_mtyfb].Free();
    //mtp->tt = (double *)aObj[ o_mttt].Free();
}

//=============================================================

// Conversion of concentration units to moles
//
double TGEM2MT::Reduce_Conc( char UNITP, double Xe, double DCmw, double Vm,
    double R1, double Msys, double Mwat, double Vaq, double Maq, double Vsys )
{
    double Xincr = 0.;
    switch( UNITP )
    {  // Quantities
    case QUAN_MKMOL: /*'Y'*/
        Xincr = Xe / 1e6;
        goto FINISH;
    case QUAN_MMOL:  /*'h'*/
        Xincr = Xe / 1e3;
        goto FINISH;
    case QUAN_MOL:   /*'M'*/
        Xincr = Xe;
        goto FINISH;
    }
    if( DCmw > 1e-12 )
        switch( UNITP )
        {
        case QUAN_MGRAM: /*'y'*/
            Xincr = Xe / DCmw / 1e3;
            goto FINISH;
        case QUAN_GRAM:  /*'g'*/
            Xincr = Xe / DCmw;
            goto FINISH;
        case QUAN_KILO:  /*'G'*/
            Xincr = Xe * 1e3 / DCmw;
            goto FINISH;
        }
    /* Concentrations */
    if( fabs( R1 ) > 1e-12 )
        switch( UNITP )
        { // mole fractions relative to total moles in the system
        case CON_MOLFR:  /*'n'*/
            Xincr = Xe * R1;
            goto FINISH;
        case CON_MOLPROC:/*'N'*/
            Xincr = Xe / 100. * R1;
            goto FINISH;
        case CON_pMOLFR: /*'f'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow(10., -Xe )* R1;
            goto FINISH;
        }
    if( fabs( Vsys ) > 1e-12 && Vm > 1e-12 )
        switch( UNITP )   /* Volumes */
        {
        case CON_VOLFR:  /*'v'*/
            Xincr = Xe * Vsys * 1e3 / Vm;
            goto FINISH;
        case CON_VOLPROC:/*'V'*/
            Xincr = Xe * Vsys * 10. / Vm;
            goto FINISH;
        case CON_pVOLFR: /*'u'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vsys / Vm;
            goto FINISH;
        }
    if( fabs( Msys ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP ) // Mass fractions relative to mass of the system
        {
        case CON_WTFR:   /*'w'*/
            Xincr = Xe * Msys * 1e3 / DCmw;
            goto FINISH;
        case CON_WTPROC: /*'%'*/
            Xincr = Xe * Msys *10. / DCmw;
            goto FINISH;
        case CON_PPM:    /*'P'*/
            Xincr = Xe * Msys / 1e3 / DCmw;
            goto FINISH;
        }
    if( fabs( Mwat ) > 1e-12 )
        switch( UNITP ) /* Molalities */
        {
        case CON_MOLAL:  /*'m'*/
            Xincr = Xe * Mwat;
            goto FINISH;
        case CON_MMOLAL: /*'i'*/
            Xincr = Xe / 1e3 * Mwat;
            goto FINISH;
        case CON_pMOLAL: /*'p'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Mwat;
            goto FINISH;
        }
    if( fabs( Vaq ) > 1e-12 )
        switch( UNITP )  /* Molarities */
        {
        case CON_MOLAR:  /*'L'*/
            Xincr = Xe * Vaq;
            goto FINISH;
        case CON_MMOLAR: /*'j'*/
            Xincr = Xe * Vaq / 1e3;
            goto FINISH;
        case CON_pMOLAR: /*'q'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vaq;
            goto FINISH;
            /* g/l, mg/l, mkg/l */
        case CON_AQGPL: /*'d'*/
            Xincr = Xe * Vaq / DCmw;
            goto FINISH;
        case CON_AQMGPL: /*'e'*/
            Xincr = Xe * Vaq / DCmw / 1e3;
            goto FINISH;
        case CON_AQMKGPL: /*'b'*/
            Xincr = Xe * Vaq / DCmw / 1e6;
            goto FINISH;
        }
    if( fabs( Maq ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP )     /* Weight concentrations */
        {
        case CON_AQWFR:  /*'C'*/
            Xincr = Xe * Maq * 1e3 / DCmw;
            goto FINISH;
        case CON_AQWPROC:/*'c'*/
            Xincr = Xe * 10. * Maq / DCmw;
            goto FINISH;
        case CON_AQPPM:  /*'a'*/
            Xincr = Xe * Maq / 1e3 / DCmw;
            goto FINISH;
        }
    /* Error */
FINISH:
    return Xincr;
}

//---------------------------------------------------------------------------

