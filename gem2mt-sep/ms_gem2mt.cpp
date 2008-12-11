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
#include "io_arrays.h"
#include "ms_gem2mt.h"
#include "gstring.h"


TGEM2MT* TGEM2MT::pm;

TGEM2MT::TGEM2MT()
{
    mtp=&mt[0];

    // default data
    memset(mtp, 0, sizeof(GEM2MT) );
    mtp->PvMO =   S_ON;
    mtp->iStat =  AS_READY;

    na = 0;
    pa = 0;
}

TGEM2MT::~TGEM2MT()
{
  Free();
  if( na )
   delete na;
  if( pa )
    delete pa;
}

//Calculate record
void TGEM2MT::RecCalc( )
{
   // if( mtp->PsMode == 'A' || mtp->PsMode == 'D' || mtp->PsMode == 'T' ) 'W' 'V'
   {   // calculate start data
   
     bool iRet; 
     
     if( mtp->PsMode == GMT_MODE_F ) // Flux-box RMT scoping model
       iRet = CalcBoxModel( NEED_GEM_SIA );
     else
    	iRet =  Trans1D( NEED_GEM_SIA );

  }
}

// read TGEM2MT structure from file
int TGEM2MT::ReadTask( const char *unsp_in1 )
{
 // read GEM2MT structure from file
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
   fstream ff(unsp_in1, ios::in );
   ErrorIf( !ff.good() , unsp_in1, "Fileopen error");
   from_text_file( ff );
   return 0;
  }
  catch(TError& err)
  {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
  }
  return 1;
}

// Write TGEM2MT structure from file
int TGEM2MT::WriteTask( const char *unsp_in1 )
{
 // read GEM2MT structure from file
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
   fstream ff(unsp_in1, ios::out );
   ErrorIf( !ff.good() , unsp_in1, "Fileopen error");
   to_text_file( ff, true );
   gstring filename = unsp_in1;
           filename += ".res";
   //fstream ff1( filename.c_str(), ios::out );
   //ErrorIf( !ff1.good() , filename.c_str(), "Fileopen error");
   // result_to_text_file( ff1, true );
   return 0;
  }
  catch(TError& err)
  {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
  }
  return 1;
}

outField TGEM2MT_static_fields[21] =  {
 { "Mode", 1,0 },
 { "PvFDL", 0,0 },  // PvPGD, PvFDL, PvSFL
 { "PvGrid", 0,0 },
 { "Size" , 1,0 },
 { "nFD" , 0,0 },
 { "nPG" , 0,0 },
 { "nSFD" , 0,0 },
 { "MaxSteps", 1,0 },
 { "nPTypes", 0,0 },
 { "nProps", 0,0 },
 { "Tau", 1,0 },
 { "LSize", 0,0 },
 { "fVel", 1,0 },
 { "cLen", 1,0 },
 { "tf", 1,0 },
 { "cdv", 1,0 },
 { "cez", 1,0 },
 { "al_in", 1,0 },
 { "Dif_in", 1,0 },
 { "FIf", 0,0 },
 { "Nf", 0,0 }
  };

outField TGEM2MT_dynamic_fields[16] =  {
 { "DiCp", 1, 0 },
 { "HydP", 1, 0 },
 { "NPmean", 0, 0 },
 { "nPmin", 0, 0 },
 { "nPmax", 0, 0 },
 { "ParTD", 0, 0 },
 { "mGrid", 0, 0 },
 { "FDLi", 0, 0 },
 { "FDLf", 0, 0 },
 { "BSF", 0, 0 },
 { "PGT", 0, 0 },
 { "UMPG", 0, 0 },
 { "FDLid", 0, 0 },
 { "FDLop", 0, 0 },
 { "FDLmp", 0, 0 },
 { "MPGid", 0, 0 }
 };


// Reading TGEM2MT structure from text file
void TGEM2MT::from_text_file(fstream& ff)
{
// static arrays
 TReadArrays  rdar( 21, TGEM2MT_static_fields, ff);
 short nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
   case 0: rdar.readArray( "Mode", &mtp->PsMode, 1, 1);
           break;
   case 1: rdar.readArray( "PvFDL", &mtp->PvPGD, 3, 1);
           break;
   case 2: rdar.readArray( "PvGrid", &mtp->PvGrid, 1, 1);
           break;
   case 3: rdar.readArray( "Size", &mtp->nC, 1);
           //  mtp->nC = mtp->xC*mtp->yC*mtp->zC;
           break;
   case 4: rdar.readArray( "nFD", &mtp->nFD, 1);
           break;
   case 5: rdar.readArray( "nPG", &mtp->nPG, 1);
           break;
   case 6: rdar.readArray( "nSFD", &mtp->nSFD,  1);
           break;
   case 7: rdar.readArray( "MaxSteps", &mtp->ntM, 1);
           break;
   case 8: rdar.readArray( "nPTypes", &mtp->nPTypes,  1);
           break;
   case 9: rdar.readArray( "nProps", &mtp->nProps,  1);
           break;
   case 10: rdar.readArray( "Tau", mtp->Tau, 3);
           break;
   case 11: rdar.readArray( "LSize", mtp->sizeLc, 3);
           break;
   case 12: rdar.readArray( "fVel", &mtp->fVel, 1);
        break;
   case 13: rdar.readArray( "cLen", &mtp->cLen, 1);
         break;
   case 14: rdar.readArray( "tf", &mtp->tf, 1);
        break;
   case 15: rdar.readArray( "cdv", &mtp->cdv, 1);
        break;
   case 16: rdar.readArray( "cez", &mtp->cez, 1);
        break;
   case 17: rdar.readArray( "al_in", &mtp->al_in, 1);
        break;
   case 18: rdar.readArray( "Dif_in", &mtp->Dif_in, 1);
        break;
   case 19: rdar.readArray( "Nf", &mtp->Nf, 1);
        break;
   case 20: rdar.readArray( "FIf", &mtp->FIf, 1);
        break;
 }
   nfild = rdar.findNext();
 }

 if( mtp->PsMode == GMT_MODE_F )
 { 
    rdar.setAlws( 1 /*"PvFDL"*/);
    rdar.setAlws( 4 /*"nFD"*/);
    rdar.setAlws( 5 /*"nPG"*/);
    rdar.setAlws( 6 /*"nSFD"*/);
    rdar.setAlws( 19 /*"Nf"*/);
    rdar.setAlws( 20 /*"Fif"*/);
 }     

 if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
 {  
    rdar.setAlws( 2 /*"PvGrid"*/);
    rdar.setAlws( 8 /*"nPTypes"*/);
    rdar.setAlws( 9 /*"nProps"*/);
    rdar.setAlws( 11 /*"LSize"*/);
 }     
 // setup default DATACH readed after
  //mtp->Nf =  mtp->nICb = TNodeArray::na->pCSD()->nICb;
  //mtp->FIf = mtp->nPHb = TNodeArray::na->pCSD()->nPHb;
  // mtp->nDCb =TNodeArray::na->pCSD()->nDCb;
 
 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from TGEM2MT structure";
    Error( "Error", ret);
  }

  Alloc();

//dynamic data
 TReadArrays  rddar( 16, TGEM2MT_dynamic_fields, ff);

// Set up flags
 if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
 { 
   rddar.setAlws( 2 /*"NPmean"*/);
   rddar.setAlws( 3 /*"nPmin"*/);
   rddar.setAlws( 4 /*"nPmax"*/);
   rddar.setAlws( 5 /*"ParTD"*/);
   if( mtp->PvGrid != S_OFF )
	    rddar.setAlws( 6 /*"grid"*/);
 }     
 if( mtp->PsMode == GMT_MODE_F )
 {  
   if( mtp->PvFDL != S_OFF )
   {
	   rddar.setAlws( 7 /*"FDLi"*/);
	   rddar.setAlws( 8 /*"FDLf"*/);
	   rddar.setAlws( 12 /*"FDLid"*/);
	   rddar.setAlws( 13 /*"FDLop"*/);
	   rddar.setAlws( 14 /*"FDLmp"*/);
   }
  if( mtp->PvPGD != S_OFF )
  {
		   rddar.setAlws( 11 /*"UMPG"*/);
		   rddar.setAlws( 10 /*"PGT"*/);
		   rddar.setAlws( 15 /*"MPGid"*/);
  }
  if( mtp->PvSFL != S_OFF )
	   rddar.setAlws( 9 /*"BSF"*/);
 }     

 nfild = rddar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
   case 0: rddar.readArray( "DiCp", mtp->DiCp[0], mtp->nC*2 );
            break;
   case 1: rddar.readArray( "HydP", mtp->HydP[0], mtp->nC*SIZE_HYDP );
            break;
   case 2: rddar.readArray( "NPmean", mtp->NPmean, mtp->nPTypes );
            break;
   case 3: rddar.readArray( "nPmin", mtp->nPmin, mtp->nPTypes );
            break;
   case 4: rddar.readArray( "nPmax", mtp->nPmax, mtp->nPTypes );
            break;
   case 5: rddar.readArray( "ParTD", mtp->ParTD[0], mtp->nPTypes*6 );
             break;
   case 6: rddar.readArray(  "grid", mtp->grid[0], mtp->nC*3 );
             break;
   case 7: rddar.readArray( "FDLi", mtp->FDLi[0], mtp->nFD*2 );
             break;
   case 8: rddar.readArray( "FDLf", mtp->FDLf[0], mtp->nFD*4 );
             break;
   case 9: rddar.readArray(  "BSF", mtp->BSF, mtp->nSFD*mtp->Nf );
             break;
   case 10: rddar.readArray(  "PGT", mtp->PGT, mtp->FIf*mtp->nPG );
             break;
   case 11: rddar.readArray(  "UMPG", mtp->UMPG, mtp->FIf, 1 );
            break;
   case 12: rddar.readArray(  "FDLid", mtp->FDLid[0], mtp->nFD, MAXSYMB );
            break;
   case 13: rddar.readArray(  "FDLop", mtp->FDLop[0],  mtp->nFD, MAXSYMB  );
            break;
   case 14:rddar.readArray(  "FDLmp", mtp->FDLmp[0], mtp->nFD, MAXSYMB  );
            break;
   case 15:rddar.readArray(  "MPGid", mtp->MPGid[0], mtp->nPG, MAXSYMB );
            break;
  }
     nfild = rddar.findNext();
 }
 
 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from TGEM2MT structure";
    Error( "Error", ret);
  }

}


void TGEM2MT::Alloc()
{
  mtp->DiCp = new short[mtp->nC][2];
  mtp->HydP = new double[mtp->nC][SIZE_HYDP];
  if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
  {
    mtp->NPmean = new short[ mtp->nPTypes ];
    mtp->nPmin =  new short[ mtp->nPTypes ];
    mtp->nPmax =  new short[ mtp->nPTypes ];
    mtp->ParTD = new short[ mtp->nPTypes][6];
    if( mtp->PvGrid != S_OFF )
      mtp->grid = new float[mtp->nC][3];
  }
  if( mtp->PvFDL != S_OFF )
     {
        mtp->FDLi =  new short [mtp->nFD][2];
        mtp->FDLf =  new float [mtp->nFD][4];
        mtp->FDLid=  new char [mtp->nFD][MAXSYMB];
        mtp->FDLop=  new char [mtp->nFD][MAXSYMB];
        mtp->FDLmp = new char [mtp->nFD][MAXSYMB];
     }
   if( mtp->PvPGD != S_OFF )
     {
        mtp->PGT  =  new float[mtp->FIf*mtp->nPG];
        mtp->MPGid = new char [mtp->nPG][MAXSYMB];
        mtp->UMPG =  new char [ mtp->FIf ];
     }

   if( mtp->PvSFL != S_OFF )
       mtp->BSF =  new double[ mtp->nSFD*mtp->Nf];


   if( mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF )
     {
        mtp->MB =  new double[ mtp->nC * mtp->Nf];
        mtp->dMB = new double[ mtp->nC * mtp->Nf];
     }

}

void TGEM2MT::Free()
{
  if( mtp->DiCp  )
    delete[] mtp->DiCp;
  if( mtp->HydP  )
     delete[] mtp->HydP;
  if( mtp->NPmean  )
    delete[] mtp->NPmean;
  if( mtp->nPmin  )
    delete[] mtp->nPmin;
  if( mtp->nPmax  )
    delete[] mtp->nPmax;
  if( mtp->ParTD  )
    delete[] mtp->ParTD;
  if( mtp->grid  )
    delete[] mtp->grid;
  if(  mtp->FDLi )
	  delete[] mtp->FDLi;	  
  if( mtp->FDLf  )
    delete[] mtp->FDLf;
  if( mtp->FDLid  )
     delete[] mtp->FDLid;
  if( mtp->FDLop  )
    delete[] mtp->FDLop;
  if( mtp->FDLmp  )
    delete[] mtp->FDLmp;
  if( mtp->PGT  )
    delete[] mtp->PGT;
  if( mtp->MPGid  )
    delete[] mtp->MPGid;
  if( mtp->UMPG  )
    delete[] mtp->UMPG;
  if(  mtp->BSF )
	  delete[] mtp->BSF;
  if( mtp->MB  )
    delete[] mtp->MB;
  if( mtp->dMB  )
    delete[] mtp->dMB;
  if( mtp->gfc  )
    delete[] mtp->gfc;
  if( mtp->yfb  )
    delete[] mtp->yfb;
  if( tt  )
      delete[] tt;
}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up NodeArray and ParticleArray classes
int TGEM2MT::MassTransInit( const char *chbr_in1 )
{
  int ii;
  // The NodeArray must be allocated here
  TNodeArray::na = na = new TNodeArray( /* mtp->xC,mtp->yC,mtp->zC */ mtp->nC );

 // Prepare the array for initial conditions allocation
  long int* nodeType = new long int[mtp->nC];
  for( ii =0; ii<mtp->nC; ii++ )
         nodeType[ii] = mtp->DiCp[ii][0];

  // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
  if( na->GEM_init( chbr_in1, nodeType ) )
        return 1;  // error reading files

   // put HydP
  DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
  for( int jj=0; jj<mtp->nC; jj ++)
  {
     C0[jj]->NodeTypeHY = mtp->DiCp[jj][1];
     if( mtp->HydP )
     { C0[jj]->Vt = mtp->HydP[jj][0];
       C0[jj]->vp = mtp->HydP[jj][1];
       C0[jj]->eps = mtp->HydP[jj][2];
       C0[jj]->Km = mtp->HydP[jj][3];
       C0[jj]->al = mtp->HydP[jj][4];
       C0[jj]->Dif = mtp->HydP[jj][5];
       C0[jj]->hDl = C0[jj]->al*C0[jj]->vp+C0[jj]->Dif;
       C0[jj]->nto = mtp->HydP[jj][6];
     }
  }

 for ( ii=0; ii<mtp->nC; ii++)    // node iteration
  {
      na->CopyNodeFromTo( ii, mtp->nC, C0, na->pNodT1() );
  }  // ii    end of node iteration loop


  // use particles
  if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
  {
   na->SetGrid( mtp->sizeLc, mtp->grid );   // set up grid structure
   pa = new TParticleArray( mtp->nPTypes, mtp->nProps,
         mtp->NPmean, mtp->ParTD, mtp->nPmin, mtp->nPmax, na );
  }

  delete[] nodeType;
  return 0;
}

//---------------------------------------------------------------------------

