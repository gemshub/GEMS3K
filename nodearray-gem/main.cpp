//-------------------------------------------------------------------
// $Id: main.cpp 792 2006-09-19 08:10:41Z gems $
//
// Debugging version of a finite-difference 1D advection-diffusion
// mass transport model supplied by Dr. Frieder Enzmann (Uni Mainz)
// coupled with GEMIPM2K module for calculation of chemical equilibria
//
// Direct access to the TNodeArray class for storing all data for nodes
//
// Copyright (C) 2005 S.Dmytriyeva, F.Enzmann, D.Kulik
//
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include "io_arrays.h"
#include "ms_gem2mt.h"
#include "gstring.h"

istream&
f_getline(istream& is, gstring& str, char delim);

//---------------------------------------------------------------------------
// Test of 1D advection (finite difference method provided by Dr. F.Enzmann,
// Uni Mainz) coupled with GEMIPM2K kernel (PSI) using the TNodeArray class.
// Finite difference calculations split over independent components
// (through bulk composition of aqueous phase).
// Experiments with smoothing terms on assigning differences to bulk composition
// of nodes

int main( int argc, char* argv[] )
 {
     int       RetC = 0;
     gstring gem2mt_in1 = "gem2mt_init.txt";
     //gstring multu_in1 = "MgWBoundC.ipm";
     gstring chbr_in1 = "ipmfiles-dat.lst";

// from argv
      //if (argc >= 2 )
      //  multu_in1 = argv[1];
      if (argc >= 2 )
       chbr_in1 = argv[1];
      if (argc >= 3 )
        gem2mt_in1 = argv[2];

   try{

// The NodeArray must be allocated here
    TGEM2MT::pm = new TGEM2MT();

// Here we read the GEM2MT structure, prepared from GEMS or by hand
   if( TGEM2MT::pm->MassTransSetUp( gem2mt_in1.c_str() ))
     return 1;  // error reading files

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if( TGEM2MT::pm->MassTransInit( chbr_in1.c_str() ) )
      return 1;  // error reading files

// here we call the mass-transport finite-difference coupled routine
   RetC = TGEM2MT::pm->Trans1D( NEED_GEM_AIA );
   }
   catch(TError& err)
       {
        fstream f_log("ipmlog.txt", ios::out|ios::app );
        f_log << err.title.c_str() << ": " << err.mess.c_str() << endl;
        return 1;
       }

   return RetC;
}

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

outField TGEM2MT_static_fields[15] =  {
 { "Mode", 1,0 },
 { "Size" , 1,0 },
 { "MaxSteps", 1,0 },
 { "Tau", 1,0 },
 { "Grid", 1,0 },
 { "Types", 1,0 },
 { "Props", 1,0 },
 { "LSize", 1,0 },
 { "fVel", 1,0 },
 { "cLen", 1,0 },
 { "tf", 1,0 },
 { "cdv", 1,0 },
 { "cez", 1,0 },
 { "al_in", 1,0 },
 { "Dif_in", 1,0 }
 };

outField TGEM2MT_dynamic_fields[3] =  {
 { "DiCp", 1, 0 },
 { "NPmean", 1, 0 },
 { "mGrid", 1, 0 }
 };

// The mass transport start constant
int TGEM2MT::MassTransSetUp( const char *gem2mt_in1 )
{
/*
  mtp->PsMode = GMT_MODE_A;
  mtp->nC = 101;    // number of nodes (default 1500)
  mtp->ntM = 300;   // max number of time steps   10000

  mtp->cLen = 1.;      // length in m
  mtp->fVel = 1e-9;    // fluid velocity constant m/sec
  mtp->tf = 1.;     // time step reduce factor
  mtp->cdv = 1e-7;   // cutoff value for delta_T corrections for bulk compositions)
  mtp->cez = 1e-12;   // minimal allowed amount of element (except charge) in bulk composition

  mtp->Tau[START_] = 0.;
  mtp->Tau[STOP_] = 0.;
  mtp->Tau[STEP_] = 0.;
*/
 // read GEM2MT structure from file
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
   fstream ff(gem2mt_in1, ios::in );
   ErrorIf( !ff.good() , gem2mt_in1, "Fileopen error");

// static arrays
   TReadArrays  rdar( 15, TGEM2MT_static_fields, ff);
   short nfild = rdar.findNext();
   while( nfild >=0 )
   {
     switch( nfild )
     {
       case 0: rdar.readArray( "Mode", &mtp->PsMode, 1, 1);
            break;
       case 1: rdar.readArray( "Size", &mtp->xC, 3);
            mtp->nC = mtp->xC*mtp->yC*mtp->zC;
             break;
       case 2: rdar.readArray( "MaxSteps", &mtp->ntM, 1);
             break;
       case 3: rdar.readArray( "Tau", mtp->Tau, 3);
            break;
       case 4: rdar.readArray( "Grid", &mtp->PvGrid, 1, 1);
            break;
       case 5: rdar.readArray( "Types", &mtp->nPTypes, 1);
            break;
       case 6: rdar.readArray( "Props", &mtp->nProps, 1);
            break;
       case 7: rdar.readArray( "LSize", mtp->sizeLc, 3);
            break;
       case 8: rdar.readArray( "fVel", &mtp->fVel, 1);
            break;
       case 9: rdar.readArray( "cLen", &mtp->cLen, 1);
             break;
       case 10: rdar.readArray( "tf", &mtp->tf, 1);
            break;
       case 11: rdar.readArray( "cdv", &mtp->cdv, 1);
            break;
       case 12: rdar.readArray( "cez", &mtp->cez, 1);
            break;
       case 13: rdar.readArray( "al_in", &mtp->al_in, 1);
            break;
       case 14: rdar.readArray( "Dif_in", &mtp->Dif_in, 1);
            break;
     }
     nfild = rdar.findNext();
  }

   // testing read
   if( mtp->PsMode != GMT_MODE_W && mtp->PsMode != GMT_MODE_V )
   {
      rdar.setNoAlws( 4 /*"Grid"*/);
      rdar.setNoAlws( 5 /*"Types"*/);
      rdar.setNoAlws( 6 /*"Props"*/);
      rdar.setNoAlws( 7 /*"LSize"*/);
      mtp->PvGrid = '-';
   }
   gstring ret = rdar.testRead();
   if( !ret.empty() )
   { ret += " - fields must be readed from TGEM2MT structure";
     Error( "Error", ret);
   }

 // realloc memory
   Alloc();
 // read arrays
   int ii, jj, indx;

   TReadArrays  rddar( 3, TGEM2MT_dynamic_fields, ff);
  // Set up flags
   if( mtp->PsMode != GMT_MODE_W && mtp->PsMode != GMT_MODE_V )
      rddar.setNoAlws( 1 /*"NPmean"*/);
   if( mtp->PvGrid == S_OFF )
      rddar.setNoAlws( 2 /*"mGrid"*/);

   nfild = rddar.findNext();
   while( nfild >=0 )
   {
     switch( nfild )
     {
// line %5s "#", %6s "Init", %6s "Type",
//   %12s "Vt-m**3", %12s "vp-m/sec", %12s "porosity", %12s "Km-m**2",
//   %12s "al-m", %12s "hDl-m**2/s", %12s "nto"
// list #DiCp %5g index, %6g all #DiCp, %12.6g all #HydP
       case 0:  rddar.skipSpace();
               for( ii=0; ii<mtp->nC; ii++ )
               {
                 ff >> indx >> mtp->DiCp[ii][0] >> mtp->DiCp[ii][1];
                 for( jj=0; jj< SIZE_HYDP; jj++)
                    ff >> mtp->HydP[ii][jj];
               }
            break;
       case 1:  rddar.skipSpace();
           if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
// line ## prn=:mtPsfl[0] = "W" | mtPsfl[0] = "V"; ##
//    %5s "#", %6s "Mean", %6s "Min", %6s "Max",
//    %6s "ptype", %6s "mmode", %6s "tcode", %6s "ips", %6s "res", %6s "res"
// list #NPmean %5g index, %6g #NPmean, %6g #nPmin ,%6g #nPmax ,
//    %6g all #ParTD
              for( ii=0; ii<mtp->nPTypes; ii++ )
              {
                ff >> indx >> mtp->NPmean[ii];
                ff >> mtp->nPmin[ii] >> mtp->nPmax[ii];
                for( jj=0; jj< 6; jj++)
                   ff >> mtp->ParTD[ii][jj];
              }
            break;
       case 2:  rddar.skipSpace();
          if( mtp->PvGrid != S_OFF )
//     line ## prn=:mtPvfl[8] <> "-"; ##  %s "Grid"
//     list #mGrid  %12.6g all #mGrid
             for( ii=0; ii<mtp->nPTypes; ii++ )
               ff >> mtp->grid[ii][0]  >> mtp->grid[ii][1] >> mtp->grid[ii][2];
          break;
     }
     nfild = rddar.findNext();
   }
   ret = rddar.testRead();
   if( !ret.empty() )
   { ret += " - fields must be readed from TGEM2MT structure";
     Error( "Error", ret);
   }

   return 0;
  }
  catch(TError& err)
  {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
  }
  return 1;
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
}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up NodeArray and ParticleArray classes
int TGEM2MT::MassTransInit( const char *chbr_in1 )
{
  int ii;
  // The NodeArray must be allocated here
  TNodeArray::na = na = new TNodeArray( mtp->xC,mtp->yC,mtp->zC/*mtp->nC*/ );

 // Prepare the array for initial conditions allocation
  int* nodeType = new int[mtp->nC];
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

