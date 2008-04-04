//-------------------------------------------------------------------
// $Id: ms_multi_file.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of writing/reading IPM work data structure files
//
// Copyright (C) 2006-2007 S.Dmytriyeva
//
// This file is part of the GEM-Vizor library and GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

#include <math.h>

#include "io_arrays.h"
#include "m_param.h"
#include "gdatastream.h"

void TMulti::getLsModsum( int& LsModSum, int& LsIPxSum )
{  LsModSum = 0;
   LsIPxSum = 0;
   for(int i=0; i<pm.FIs; i++)
   {
     LsModSum += (pm.LsMod[i*3]*pm.LsMod[i*3+2]);
     LsIPxSum += (pm.LsMod[i*3]*pm.LsMod[i*3+1]);
   }
}


void TMulti::getLsMdcsum( int& LsMdcSum )
{  LsMdcSum = 0;
   for(int i=0; i<pm.FIs; i++)
     LsMdcSum += (pm.LsMdc[i]*pm.L1[i]);
}

//---------------------------------------------------------//

// writing MULTI to binary file
void TMulti::to_file( GemDataStream& ff, gstring& path  )
{
   if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( GetName(), "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;


#ifndef IPMGEMPLUGIN

   PAalp = syp->PAalp;
   PSigm = syp->PSigm;
   EpsW = TProfil::pm->tpp->EpsW;
   RoW = TProfil::pm->tpp->RoW;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
   EpsW = EpsW_;
   RoW = RoW_;
#endif


   ff.writeArray(pm.stkey, sizeof(char)*(EQ_RKLEN+5));
   ff.writeArray( &pm.N, 38);
   ff.writeArray(&pm.TC, 55);
   ff << PAalp;
   ff << PSigm;
   ff << EpsW;
   ff << RoW;

   //dynamic values

    // Part 1

    /* need  always to alloc vectors */
   ff.writeArray(pm.L1,  pm.FI);
   ff.writeArray(pm.muk, pm.FI);
   ff.writeArray(pm.mui, pm.N);
   ff.writeArray(pm.muj, pm.L);
   ff.writeArray(pm.DUL, pm.L);
   ff.writeArray(pm.DLL, pm.L);
   ff.writeArray(pm.Vol, pm.L);
   ff.writeArray(pm.Pparc, pm.L);
   ff.writeArray(pm.MM, pm.L);
   ff.writeArray(pm.Awt, pm.N);
   ff.writeArray(pm.A,  pm.N*pm.L);
   ff.writeArray(pm.XFs, pm.FI);
   ff.writeArray(pm.Falps, pm.FI);
   ff.writeArray(pm.G, pm.L);
   ff.writeArray(pm.G0, pm.L);
   ff.writeArray(pm.lnGam, pm.L);
   ff.writeArray(pm.lnGmo, pm.L);
   ff.writeArray(pm.B, pm.N);
   ff.writeArray(pm.U,  pm.N);
   ff.writeArray(pm.U_r, pm.N);
   ff.writeArray(pm.C, pm.N);
   ff.writeArray(pm.XF, pm.FI);
   ff.writeArray(pm.YF, pm.FI);
   ff.writeArray(pm.Falp, pm.FI);
   ff.writeArray(pm.X, pm.L);
   ff.writeArray(pm.Y, pm.L);
   ff.writeArray(pm.XY, pm.L);
   ff.writeArray(pm.MU, pm.L);
   ff.writeArray(pm.EMU,  pm.L);
   ff.writeArray(pm.NMU, pm.L);
   ff.writeArray(pm.W, pm.L);
   ff.writeArray(pm.F, pm.L);
   ff.writeArray(pm.F0, pm.L);
   ff.writeArray(pm.YOF, pm.FI);

   ff.writeArray((char*)pm.SB, (MAXICNAME+MAXSYMB)*pm.N);
   ff.writeArray((char*)pm.SB1, MAXICNAME * pm.N);
   ff.writeArray((char*)pm.SFs, (MAXPHNAME+MAXSYMB)*pm.FI);
   ff.writeArray((char*)pm.SM, MAXDCNAME * pm.L);
   ff.writeArray((char*)pm.SF, (MAXPHNAME+MAXSYMB)*pm.FI);
   ff.writeArray((char*)pm.SM2, MAXDCNAME * pm.Ls);
   ff.writeArray((char*)pm.SF2, (MAXPHNAME+MAXSYMB)*pm.FIs);

   ff.writeArray( pm.RLC, pm.L);
   ff.writeArray( pm.RSC, pm.L);
   ff.writeArray( pm.ICC, pm.N);
   ff.writeArray( pm.DCC, pm.L);
   ff.writeArray( pm.PHC, pm.FI);
   ff.writeArray( pm.DCCW, pm.L);

   ff.writeArray( pm.lnGmM, pm.L);
   ff.writeArray( pm.GEX,   pm.L);
   ff.writeArray( pm.FVOL, pm.FI);
   ff.writeArray( pm.FWGT, pm.FI);

    if( pm.L > 0 )
    {
      ff.writeArray(pm.Y_la, pm.L);
      ff.writeArray(pm.Y_w, pm.L);
      ff.writeArray(pm.Fx, pm.L);
      ff.writeArray(pm.Wx, pm.L);
      ff.writeArray(pm.VL, pm.L);
      ff.writeArray(pm.Gamma, pm.L);
      ff.writeArray(pm.lnGmf, pm.L);
//      ff.writeArray(pm.D, pm.L);
    }

   // Part 2  not requited arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
      ff.writeArray(pm.BF, pm.FIs*pm.N);
      ff.writeArray(pm.XFA, pm.FIs);
      ff.writeArray(pm.YFA, pm.FIs);
      ff.writeArray(pm.LsMod, pm.FIs*3);
      ff.writeArray(pm.LsMdc, pm.FIs);
      int LsModSum;
      int LsIPxSum;
      int LsMdcSum;
      getLsModsum( LsModSum, LsIPxSum );
      getLsMdcsum( LsMdcSum );
      ff.writeArray(pm.IPx, LsIPxSum);
      ff.writeArray(pm.PMc, LsModSum);
      ff.writeArray(pm.DMc, LsMdcSum);
      ff.writeArray(pm.PUL, pm.FIs);
      ff.writeArray(pm.PLL, pm.FIs);

      ff.writeArray((char*)pm.sMod, 6*pm.FIs);
      ff.writeArray( pm.RFLC, pm.FIs);
      ff.writeArray( pm.RFSC, pm.FIs);
    }

    if( pm.LO > 1 )
    {
      ff.writeArray(pm.Y_m, pm.L);
      ff.writeArray(pm.IC_m, pm.N);
      ff.writeArray(pm.IC_lm, pm.N);
      ff.writeArray(pm.IC_wm,  pm.N);
    }

    /* dispersion and sorbtion phases */
    if( PAalp != S_OFF )
    {
      ff.writeArray(pm.Aalp, pm.FI);
      ff.writeArray((float *)pm.Xr0h0, pm.FI*2);
    }

   if( PSigm != S_OFF )
      ff.writeArray(pm.Sigw, pm.FI);

    if( PSigm != S_OFF )
      ff.writeArray(pm.Sigg, pm.FI);

    if( pm.E )
    {
      ff.writeArray(pm.EZ,  pm.L);
      ff.writeArray(pm.Xcond, pm.FI);
      ff.writeArray(pm.Xeps,  pm.FI);
    }

    if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
      ff.writeArray((char*)pm.SCM, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.Nfsp, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.MASDT, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XcapA, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XcapB, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XcapD, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XcapF, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XdlA, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XdlB, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XdlD, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiB, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiD, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.XlamA, pm.FIs*pm.FIat);
      ff.writeArray((float*)pm.Xetaf, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XetaA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XetaB, pm.FIs*pm.FIat);
ff.writeArray((double*)pm.XetaD, pm.FIs*pm.FIat);     // added 12.09.05   KD
      ff.writeArray((double*)pm.XFTS, pm.FIs*pm.FIat);

ff.writeArray((short*)pm.SATX, pm.Lads*4);
ff.writeArray(pm.SATT, pm.Lads);
ff.writeArray((float*)pm.MASDJ, pm.Lads*DFCN);
//      ff.writeArray(pm.MASDJ, pm.Ls);
ff.writeArray( (double*)pm.lnSAC, pm.Lads*4 );
ff.writeArray((char*)pm.SM3, MAXDCNAME * pm.Lads);
ff.writeArray( pm.DCC3, pm.Lads);
ff.writeArray((double*)pm.D, MST*MST);
    }

    if( pm.PG > 0 )
    {
      ff.writeArray(pm.Fug,  pm.PG);
      ff.writeArray(pm.Fug_l,  pm.PG);
      ff.writeArray(pm.Ppg_l,  pm.PG);
    }

    // Part 3

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
      ff.writeArray(pm.Wb, pm.Ls);
      ff.writeArray(pm.Wabs, pm.Ls);
      ff.writeArray(pm.Rion, pm.Ls);

      ff.writeArray(pm.Qp, pm.FIs*QPSIZE);
      ff.writeArray(pm.Qd, pm.FIs*QDSIZE);
   }

//  Added 16.11.2004 by Sveta
//   if( pm.sitNcat*pm.sitNcat )
//     ff.writeArray( pm.sitE, pm.sitNcat*pm.sitNan );
//   if( pm.sitNcat )
//     ff.writeArray( pm.sitXcat, pm.sitNcat );
//   if( pm.sitNan )
//     ff.writeArray( pm.sitXan, pm.sitNan );

/*   gstring Path_ = path;
     gstring dir;
     gstring name;
     gstring ext;

     u_splitpath( Path_, dir, name, ext );
     Path_ = u_makepath( dir, name, "txt" );

     to_text_file( path.c_str() );
*/
}

// reading MULTI from binary file
void TMulti::from_file( GemDataStream& ff )
{
   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;

   ff.readArray(pm.stkey, sizeof(char)*(EQ_RKLEN+5));
   ff.readArray( &pm.N, 38);
   ff.readArray(&pm.TC, 55);
   ff >> PAalp;
   ff >> PSigm;
   ff >> EpsW;
   ff >> RoW;



#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
   EpsW_ = EpsW;
   RoW_ =  RoW;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#endif

   //dynamic values

    // Part 1

    /* need  always to alloc vectors */
   ff.readArray(pm.L1,  pm.FI);
   ff.readArray(pm.muk, pm.FI);
   ff.readArray(pm.mui, pm.N);
   ff.readArray(pm.muj, pm.L);
   ff.readArray(pm.DUL, pm.L);
   ff.readArray(pm.DLL, pm.L);
   ff.readArray(pm.Vol, pm.L);
   ff.readArray(pm.Pparc, pm.L);
   ff.readArray(pm.MM, pm.L);
   ff.readArray(pm.Awt, pm.N);
   ff.readArray(pm.A,  pm.N*pm.L);
   ff.readArray(pm.XFs, pm.FI);
   ff.readArray(pm.Falps, pm.FI);
   ff.readArray(pm.G, pm.L);
   ff.readArray(pm.G0, pm.L);
   ff.readArray(pm.lnGam, pm.L);
   ff.readArray(pm.lnGmo, pm.L);
   ff.readArray(pm.B, pm.N);
   ff.readArray(pm.U,  pm.N);
   ff.readArray(pm.U_r, pm.N);
   ff.readArray(pm.C, pm.N);
   ff.readArray(pm.XF, pm.FI);
   ff.readArray(pm.YF, pm.FI);
   ff.readArray(pm.Falp, pm.FI);
   ff.readArray(pm.X, pm.L);
   ff.readArray(pm.Y, pm.L);
   ff.readArray(pm.XY, pm.L);
   ff.readArray(pm.MU, pm.L);
   ff.readArray(pm.EMU,  pm.L);
   ff.readArray(pm.NMU, pm.L);
   ff.readArray(pm.W, pm.L);
   ff.readArray(pm.F, pm.L);
   ff.readArray(pm.F0, pm.L);
   ff.readArray(pm.YOF, pm.FI);

   ff.readArray((char*)pm.SB, (MAXICNAME+MAXSYMB)*pm.N);
   ff.readArray((char*)pm.SB1, MAXICNAME * pm.N);
   ff.readArray((char*)pm.SFs, (MAXPHNAME+MAXSYMB)*pm.FI);
   ff.readArray((char*)pm.SM, MAXDCNAME * pm.L);
   ff.readArray((char*)pm.SF, (MAXPHNAME+MAXSYMB)*pm.FI);
   ff.readArray((char*)pm.SM2, MAXDCNAME * pm.Ls);
   ff.readArray((char*)pm.SF2, (MAXPHNAME+MAXSYMB)*pm.FIs);

   ff.readArray( pm.RLC, pm.L);
   ff.readArray( pm.RSC, pm.L);
   ff.readArray( pm.ICC, pm.N);
   ff.readArray( pm.DCC, pm.L);
   ff.readArray( pm.PHC, pm.FI);
   ff.readArray( pm.DCCW, pm.L);

   ff.readArray( pm.lnGmM, pm.L);
   ff.readArray( pm.GEX,   pm.L);
   ff.readArray( pm.FVOL, pm.FI);
   ff.readArray( pm.FWGT, pm.FI);

    if( pm.L > 0 )
    {
      ff.readArray(pm.Y_la, pm.L);
      ff.readArray(pm.Y_w, pm.L);
      ff.readArray(pm.Fx, pm.L);
      ff.readArray(pm.Wx, pm.L);
      ff.readArray(pm.VL, pm.L);
      ff.readArray(pm.Gamma, pm.L);
      ff.readArray(pm.lnGmf, pm.L);
//      ff.readArray(pm.D, pm.L);
    }

   // Part 2  not requited arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
      ff.readArray(pm.BF, pm.FIs*pm.N);
      ff.readArray(pm.XFA, pm.FIs);
      ff.readArray(pm.YFA, pm.FIs);
      ff.readArray(pm.LsMod, pm.FIs*3);
      ff.readArray(pm.LsMdc, pm.FIs);
      int LsModSum;
      int LsIPxSum;
      int LsMdcSum;
      getLsModsum( LsModSum, LsIPxSum );
      getLsMdcsum( LsMdcSum );
      pm.IPx = new short[LsIPxSum];
      pm.PMc = new float[LsModSum];
      pm.DMc = new float[LsMdcSum];
      ff.readArray(pm.IPx, LsIPxSum);
      ff.readArray(pm.PMc, LsModSum);
      ff.readArray(pm.DMc, LsMdcSum);
      ff.readArray(pm.PUL, pm.FIs);
      ff.readArray(pm.PLL, pm.FIs);

      ff.readArray((char*)pm.sMod, 6*pm.FIs);
      ff.readArray( pm.RFLC, pm.FIs);
      ff.readArray( pm.RFSC, pm.FIs);
    }

    if( pm.LO > 1 )
    {
      ff.readArray(pm.Y_m, pm.L);
      ff.readArray(pm.IC_m, pm.N);
      ff.readArray(pm.IC_lm, pm.N);
      ff.readArray(pm.IC_wm,  pm.N);
    }

    /* dispersion and sorbtion phases */
    if( PAalp != S_OFF )
    {
      ff.readArray(pm.Aalp, pm.FI);
      ff.readArray((float *)pm.Xr0h0, pm.FI*2);
    }

   if( PSigm != S_OFF )
      ff.readArray(pm.Sigw, pm.FI);

    if( PSigm != S_OFF )
      ff.readArray(pm.Sigg, pm.FI);

    if( pm.E )
    {
      ff.readArray(pm.EZ,  pm.L);
      ff.readArray(pm.Xcond, pm.FI);
      ff.readArray(pm.Xeps,  pm.FI);
    }

    if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORBTION AND ION IXCHANDG */
      ff.readArray((char*)pm.SCM, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.Nfsp, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.MASDT, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XcapA, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XcapB, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XcapD, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XcapF, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XdlA, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XdlB, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XdlD, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiB, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiD, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.XlamA, pm.FIs*pm.FIat);
      ff.readArray((float*)pm.Xetaf, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XetaA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XetaB, pm.FIs*pm.FIat);
ff.readArray((double*)pm.XetaD, pm.FIs*pm.FIat);    // added 12.09.05  by KD
      ff.readArray((double*)pm.XFTS, pm.FIs*pm.FIat);

ff.readArray((short*)pm.SATX, pm.Lads*4);
ff.readArray(pm.SATT, pm.Lads);
ff.readArray((float*)pm.MASDJ, pm.Lads*DFCN);
//      ff.readArray(pm.MASDJ, pm.Ls);
ff.readArray((double*)pm.lnSAC, pm.Lads*4);
ff.readArray((char*)pm.SM3, MAXDCNAME * pm.Lads);
ff.readArray( pm.DCC3, pm.Lads);
ff.readArray((double*)pm.D, MST*MST);

    }

    if( pm.PG > 0 )
    {
      ff.readArray(pm.Fug,  pm.PG);
      ff.readArray(pm.Fug_l,  pm.PG);
      ff.readArray(pm.Ppg_l,  pm.PG);
    }

    // Part 3

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
      ff.readArray(pm.Wb, pm.Ls);
      ff.readArray(pm.Wabs, pm.Ls);
      ff.readArray(pm.Rion, pm.Ls);

      ff.readArray(pm.Qp, pm.FIs*QPSIZE);
      ff.readArray(pm.Qd, pm.FIs*QDSIZE);
   }
//  Added 16.11.2004 by Sveta
//   if( pm.sitNcat*pm.sitNcat )
//     ff.readArray( pm.sitE, pm.sitNcat*pm.sitNan );
//   if( pm.sitNcat )
//     ff.readArray( pm.sitXcat, pm.sitNcat );
//   if( pm.sitNan )
//     ff.readArray( pm.sitXan, pm.sitNan );
}


#ifdef IPMGEMPLUGIN

// realloc dynamic memory
void TMulti::multi_realloc( char PAalp, char PSigm )
{
  int ii;
   if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( GetName(), "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

    // Part 1
     // need  always to alloc vectors
 pm.L1 = new short[pm.FI];
 memset(pm.L1, 0, pm.FI*sizeof(short));
 pm.muk = new short[pm.FI];
 for( ii=0; ii<pm.FI; ii++)
   pm.muk[ii] = ii;
 pm.mui = new short[pm.N];
 for( ii=0; ii<pm.N; ii++)
   pm.mui[ii] = ii;
 pm.muj = new short[pm.L];
 for( ii=0; ii<pm.L; ii++)
   pm.muj[ii] = ii;

 pm.DUL = new double[pm.L];
 for( ii=0; ii<pm.L; ii++ )         // 28/11/2006
  pm.DUL[ii] = 1e6;
 pm.DLL = new double[pm.L];
 memset(pm.DLL, 0, pm.L*sizeof(double));  // 28/11/2006
 pm.Vol = new double[pm.L];
 memset(pm.Vol, 0, pm.L*sizeof(double));
 pm.Pparc = new double[pm.L];
 for( ii=0; ii<pm.L; ii++ )        // 28/11/2006
  pm.Pparc[ii] = 1.;
 pm.MM = new double[pm.L];
 memset(pm.MM, 0, pm.L*sizeof(double));
 pm.Awt = new float[pm.N];
 memset(pm.Awt, 0, pm.N*sizeof(float));
 pm.A = new float[pm.N*pm.L];
 memset(pm.A, 0, pm.N*pm.L*sizeof(float));
 pm.XFs = new float[pm.FI];
 memset(pm.XFs, 0, pm.FI*sizeof(float));
 pm.Falps = new float[pm.FI];
 memset(pm.Falps, 0, pm.FI*sizeof(float));
 pm.G = new double[pm.L];
 memset(pm.G, 0, pm.L*sizeof(double));
 pm.G0 = new double[pm.L];
 memset(pm.G0, 0, pm.L*sizeof(double));
 pm.lnGam = new double[pm.L];
 memset(pm.lnGam, 0, pm.L*sizeof(double));
 pm.lnGmo = new double[pm.L];
 memset(pm.lnGmo, 0, pm.L*sizeof(double));
 pm.B = new double[pm.N];
 memset(pm.B, 0, pm.N*sizeof(double));
 pm.U = new double[pm.N];
 memset(pm.U, 0, pm.N*sizeof(double));
 pm.U_r = new double[pm.N];
 memset(pm.U_r, 0, pm.N*sizeof(double));
 pm.C = new double[pm.N];
 memset(pm.C, 0, pm.N*sizeof(double));
 pm.XF = new double[pm.FI];
 memset(pm.XF, 0, pm.FI*sizeof(double));
 pm.YF = new double[pm.FI];
 memset(pm.YF, 0, pm.FI*sizeof(double));
 pm.Falp = new double[pm.FI];
 memset(pm.Falp, 0, pm.FI*sizeof(double));
 pm.X = new double[pm.L];
 memset(pm.X, 0, pm.L*sizeof(double));
 pm.Y = new double[pm.L];
 memset(pm.Y, 0, pm.L*sizeof(double));
 pm.XY = new double[pm.L];
 memset(pm.XY, 0, pm.L*sizeof(double));
 pm.MU = new double[pm.L];
 memset(pm.MU, 0, pm.L*sizeof(double));
 pm.EMU = new double[pm.L];
 memset(pm.EMU, 0, pm.L*sizeof(double));
 pm.NMU = new double[pm.L];
 memset(pm.NMU, 0, pm.L*sizeof(double));
 pm.W = new double[pm.L];
 memset(pm.W, 0, pm.L*sizeof(double));
 pm.F = new double[pm.L];
 memset(pm.F, 0, pm.L*sizeof(double));
 pm.F0 = new double[pm.L];
 memset(pm.F0, 0, pm.L*sizeof(double));
 pm.YOF = new double[pm.FI];
 memset(pm.YOF, 0, pm.FI*sizeof(double)); // 28/11/2006

    pm.SB = new char[pm.N][MAXICNAME+MAXSYMB];
    memset(pm.SB, 0, pm.N*(MAXICNAME+MAXSYMB)*sizeof(char));
    pm.SB1 = new char[pm.N][MAXICNAME];
    memset(pm.SB1, 0, pm.N*(MAXICNAME)*sizeof(char));
    pm.SFs = new char[pm.FI][MAXPHNAME+MAXSYMB];
    memset(pm.SFs, 0, pm.FI*(MAXPHNAME+MAXSYMB)*sizeof(char));
    pm.SM = new char[pm.L][MAXDCNAME];
    memset(pm.SM, 0, pm.L*(MAXDCNAME)*sizeof(char));
    pm.SF = new char[pm.FI][MAXPHNAME+MAXSYMB];
    memset(pm.SF, 0, pm.FI*(MAXPHNAME+MAXSYMB)*sizeof(char));
    pm.SM2 = new char[pm.Ls][MAXDCNAME];
    memset(pm.SM2, 0, pm.Ls*(MAXDCNAME)*sizeof(char));
    pm.SF2 = new char[pm.FIs][MAXPHNAME+MAXSYMB];
    memset(pm.SF2, 0, pm.FIs*(MAXPHNAME+MAXSYMB)*sizeof(char));
    pm.RLC = new char[pm.L];
    pm.RSC = new char[pm.L];
    for( ii=0; ii<pm.L; ii++ )        // 28/11/2006
    { pm.RLC[ii] = 'B';
      pm.RSC[ii] = 'M';
    }
    pm.ICC = new char[pm.N];
    memset(pm.ICC, 0, pm.N*sizeof(char));
    pm.DCC = new char[pm.L];
    memset(pm.DCC, 0, pm.L*sizeof(char));
    pm.PHC = new char[pm.FI];
    memset(pm.PHC, 0, pm.FI*sizeof(char));
    pm.DCCW = new char[pm.L];
    memset(pm.DCCW, 0, pm.L*sizeof(char));

 pm.lnGmM = new double[pm.L];
 memset(pm.lnGmM, 0, pm.L*sizeof(double));
 pm.GEX = new double[pm.L];
 memset(pm.GEX, 0, pm.L*sizeof(double)); // 28/11/2006
 pm.FVOL = new double[pm.FI];
 memset(pm.FVOL, 0, pm.FI*sizeof(double));
 pm.FWGT = new double[pm.FI];
 memset(pm.FWGT, 0, pm.FI*sizeof(double));

 if( pm.L > 0 )
 {
   pm.Y_la = new double[pm.L];
   memset(pm.Y_la, 0, pm.L*sizeof(double));
   pm.Y_w = new double[pm.L];
   memset(pm.Y_w, 0, pm.L*sizeof(double));
   pm.Fx = new double[pm.L];
   memset(pm.Fx, 0, pm.L*sizeof(double));
   pm.Wx = new double[pm.L];
   memset(pm.Wx, 0, pm.L*sizeof(double));
   pm.VL = new float[pm.L];
   memset(pm.VL, 0, pm.L*sizeof(float));
   pm.Gamma = new double[pm.L];
   memset(pm.Gamma, 0, pm.L*sizeof(double));
   pm.lnGmf = new double[pm.L];
   memset(pm.lnGmf, 0, pm.L*sizeof(double)); // 28/11/2006
//   pm.D = new double[pm.L];
 }
 else
 {
   pm.Y_la = 0;
   pm.Y_w = 0;
   pm.Fx = 0;
   pm.Wx = 0;
   pm.VL = 0;
   pm.Gamma = 0;
   pm.lnGmf = 0;
//   pm.D = 0;
 }

   // Part 2  not always required arrays

 if( pm.FIs > 0 && pm.Ls > 0 )
 {
   pm.BF = new double[pm.FIs*pm.N];
   memset(pm.BF, 0, pm.FIs*pm.N*sizeof(double));
   pm.BFC = new double[pm.N];
   memset(pm.BFC, 0, pm.N*sizeof(double));
   pm.XFA = new double[pm.FIs];
   memset(pm.XFA, 0, pm.FIs*sizeof(double));
   pm.YFA = new double[pm.FIs];
   memset(pm.YFA, 0, pm.FIs*sizeof(double));
   pm.LsMod = new short[pm.FIs*3];
   memset(pm.LsMod, 0, pm.FIs*3*sizeof(short));
   pm.LsMdc = new short[pm.FIs];
   memset(pm.LsMdc, 0, pm.FIs*sizeof(short));
   pm.IPx = 0;
   pm.PMc = 0;
   pm.DMc = 0;
   pm.PUL = new double[pm.FIs];
   for( ii=0; ii<pm.FIs; ii++ )        // 08/12/2007
    pm.PUL[ii] = 1e6;
   pm.PLL = new double[pm.FIs];
   memset(pm.PLL, 0, pm.FIs*sizeof(double));

   pm.sMod = new char[pm.FIs][6];
   memset(pm.sMod, 0, pm.FIs*6*sizeof(char));
   pm.RFLC = new char[pm.FIs];
   memset(pm.RFLC, 0, pm.FIs*sizeof(char));
   pm.RFSC = new char[pm.FIs];
   memset(pm.RFSC, 0, pm.FIs*sizeof(char));

 }
 else
 {
   pm.BF = 0;
   pm.BFC = 0;
   pm.XFA = 0;
   pm.YFA = 0;
   pm.LsMod = 0;
   pm.LsMdc = 0;
   pm.PMc = 0;
   pm.DMc = 0;
   pm.PUL = 0;
   pm.PLL = 0;

   pm.sMod = 0;
   pm.RFLC = 0;
   pm.RFSC = 0;
 }

 if( pm.LO > 1 )
 {
   pm.Y_m = new double[pm.L];
   memset(pm.Y_m, 0, pm.L*sizeof(double));
   pm.IC_m = new double[pm.N];
   memset(pm.IC_m, 0, pm.N*sizeof(double));
   pm.IC_lm = new double[pm.N];
   memset(pm.IC_lm, 0, pm.N*sizeof(double));
   pm.IC_wm = new double[pm.N];
   memset(pm.IC_wm, 0, pm.N*sizeof(double));
 }
 else
 {
   pm.Y_m = 0;
   pm.IC_m = 0;
   pm.IC_lm = 0;
   pm.IC_wm = 0;
 }

 // dispersion and sorbtion phases
 if( PAalp != S_OFF )
 {
   pm.Aalp = new float[pm.FI];
   memset(pm.Aalp, 0, pm.FI*sizeof(float)); // 28/11/2006
   pm.Xr0h0 = new float[pm.FI][2];
   memset(pm.Xr0h0, 0, pm.FI*2*sizeof(float));
 }
 else
 {
   pm.Aalp = 0;
   pm.Xr0h0 = 0;
 }

 if( PSigm != S_OFF )
 {   pm.Sigw = new float[pm.FI];
     pm.Sigg = new float[pm.FI];
     memset(pm.Sigw, 0, pm.FI*sizeof(float)); // 28/11/2006
     memset(pm.Sigg, 0, pm.FI*sizeof(float)); // 28/11/2006
 }
 else
 {   pm.Sigw = 0;
     pm.Sigg = 0;
 }

 if( pm.E )
 {
    pm.EZ = new double[pm.L];
    pm.Xcond = new float[pm.FI];
    pm.Xeps = new float[pm.FI];
    memset(pm.EZ, 0, pm.L*sizeof(double));
    memset(pm.Xcond, 0, pm.FI*sizeof(float));
    memset(pm.Xeps, 0, pm.FI*sizeof(float));
 }
 else
 {
    pm.EZ = 0;
    pm.Xcond = 0;
    pm.Xeps = 0;
 }

 if( pm.FIat > 0 /*&& pm.Lads > 0*/ && pm.FIs > 0 )
 { // ADSORBTION AND ION IXCHANDG
   pm.SATX = new short[pm.Lads][4];
   memset(pm.SATX, 0, pm.Lads*4*sizeof(short));
   pm.SCM  = new char[pm.FIs][MST];
   memset(pm.SCM, 0, pm.FIs*MST*sizeof(char));

    pm.Nfsp = new float[pm.FIs][MST];
    memset(pm.Nfsp, 0, pm.FIs*MST*sizeof(float));
    pm.MASDT = new float[pm.FIs][MST];
    memset(pm.MASDT, 0, pm.FIs*MST*sizeof(float));
    pm.XcapA = new float[pm.FIs][MST];
    memset(pm.XcapA, 0, pm.FIs*MST*sizeof(float));
    pm.XcapB = new float[pm.FIs][MST];
    memset(pm.XcapB, 0, pm.FIs*MST*sizeof(float));
    pm.XcapD = new float[pm.FIs][MST];
    memset(pm.XcapD, 0, pm.FIs*MST*sizeof(float));
    pm.XcapF = new float[pm.FIs][MST];
    memset(pm.XcapF, 0, pm.FIs*MST*sizeof(float));
    pm.XdlA = new float[pm.FIs][MST];
    memset(pm.XdlA, 0, pm.FIs*MST*sizeof(float));
    pm.XdlB = new float[pm.FIs][MST];
    memset(pm.XdlB, 0, pm.FIs*MST*sizeof(float));
    pm.XdlD = new float[pm.FIs][MST];
    memset(pm.XdlD, 0, pm.FIs*MST*sizeof(float));
    pm.XpsiA = new double[pm.FIs][MST];
    memset(pm.XpsiA, 0, pm.FIs*MST*sizeof(double));
    pm.XpsiB = new double[pm.FIs][MST];
    memset(pm.XpsiB, 0, pm.FIs*MST*sizeof(double));
    pm.XpsiD = new double[pm.FIs][MST];
    memset(pm.XpsiD, 0, pm.FIs*MST*sizeof(double));
    pm.XlamA = new float[pm.FIs][MST];
    memset(pm.XlamA, 0, pm.FIs*MST*sizeof(float));
    pm.Xetaf = new float[pm.FIs][MST];
    memset(pm.Xetaf, 0, pm.FIs*MST*sizeof(float));
    pm.XetaA = new double[pm.FIs][MST];
    memset(pm.XetaA, 0, pm.FIs*MST*sizeof(double));
    pm.XetaB = new double[pm.FIs][MST];
    memset(pm.XetaB, 0, pm.FIs*MST*sizeof(double));
    pm.XetaD = new double[pm.FIs][MST];
    memset(pm.XetaD, 0, pm.FIs*MST*sizeof(double));
    pm.MASDJ = new float[pm.Lads][DFCN];
    memset(pm.MASDJ, 0, pm.Lads*DFCN*sizeof(float));

//    pm.MASDJ = new float[pm.Ls];
pm.XFTS = new double[pm.FIs][MST];
memset(pm.XFTS, 0, pm.FIs*MST*sizeof(double));
pm.lnSAC = new double[pm.Lads][4];
memset(pm.lnSAC, 0, pm.Lads*4*sizeof(double));
pm.SATT = new char[pm.Lads];
memset(pm.SATT, 0, pm.Lads*sizeof(char));
pm.SM3 = new char[pm.Lads][MAXDCNAME];
memset(pm.SM3, 0, pm.Lads*MAXDCNAME*sizeof(char));
pm.DCC3 = new char[pm.Lads];
memset(pm.DCC3, 0, pm.Lads*sizeof(char));
pm.D = new double[MST][MST];
memset(pm.D, 0, MST*MST*sizeof(double));
 }
else
 { // ADSORPTION AND ION EXCHANGE
   pm.SCM  = 0;
    pm.Nfsp = 0;
    pm.MASDT = 0;
    pm.XcapA = 0;
    pm.XcapB = 0;
    pm.XcapD = 0;
    pm.XcapF = 0;
    pm.XdlA = 0;
    pm.XdlB = 0;
    pm.XdlD = 0;
    pm.XpsiA = 0;
    pm.XpsiB = 0;
    pm.XpsiD = 0;
    pm.XlamA = 0;
    pm.Xetaf = 0;
    pm.XetaA = 0;
    pm.XetaB = 0;
    pm.XetaD = 0;
    pm.MASDJ = 0;
    pm.XFTS = 0;
    pm.lnSAC = 0;
    pm.SATT = 0;
    pm.SM3 = 0;
    pm.DCC3 = 0;
    pm.D = 0;
 }

 if( pm.PG > 0 )
 {
  pm.Fug = new float[pm.PG];
  memset(pm.Fug, 0, pm.PG*sizeof(float));
  pm.Fug_l = new float[pm.PG];
  memset(pm.Fug_l, 0, pm.PG*sizeof(float));
  pm.Ppg_l = new float[pm.PG];
  memset(pm.Ppg_l, 0, pm.PG*sizeof(float));
 }
else
 {
  pm.Fug = 0;
  pm.Fug_l = 0;
  pm.Ppg_l = 0;
 }

   // Part 3
 if( pm.Ls > 1 && pm.FIs > 0 )
 {
    pm.Wb = new float[pm.Ls];
    memset(pm.Wb, 0, pm.Ls*sizeof(float));
    pm.Wabs = new float[pm.Ls];
    memset(pm.Wabs, 0, pm.Ls*sizeof(float));
    pm.Rion = new float[pm.Ls];
    memset(pm.Rion, 0, pm.Ls*sizeof(float));
    pm.Qp = new double[pm.FIs*QPSIZE];
    memset(pm.Qp, 0, pm.FIs*QPSIZE*sizeof(double));
    pm.Qd = new double[pm.FIs*QDSIZE];
    memset(pm.Qd, 0, pm.FIs*QPSIZE*sizeof(double));
 }
 else
 {
    pm.Wb = 0;
    pm.Wabs = 0;
    pm.Rion = 0;
    pm.Qp = 0;
    pm.Qd = 0;

 }
//  Added 16.11.2004 by Sveta
//    if( pm.sitNcat*pm.sitNcat )
//    { pm.sitE = new float[pm.sitNcat*pm.sitNan];
//      memset(pm.sitE, 0, pm.sitNcat*pm.sitNan*sizeof(float));
//    }
//    else
//       pm.sitE = 0;
//    if( pm.sitNcat )
//    {  pm.sitXcat = new short[pm.sitNcat];
//       memset(pm.sitXcat, 0, pm.sitNcat*sizeof(float));
//     }
//    else
//       pm.sitXcat = 0;
//    if( pm.sitNan )
//    {   pm.sitXan = new short[pm.sitNan];
//        memset(pm.sitXan, 0, pm.sitNan*sizeof(float));
//    }
//    else
//       pm.sitXan = 0;
}


// Reallocation of dynamic memory
void TMulti::multi_free()
{
    // Part 1
     // need  always to alloc vectors
 if( pm.L1) delete[] pm.L1;
 if( pm.muk) delete[] pm.muk;
 if( pm.mui) delete[] pm.mui;
 if( pm.muj) delete[] pm.muj;

 if( pm.DUL ) delete[] pm.DUL;
 if( pm.DLL ) delete[] pm.DLL;
 if( pm.Vol ) delete[] pm.Vol;
 if( pm.Pparc ) delete[] pm.Pparc;
 if( pm.MM ) delete[] pm.MM;
 if( pm.Awt ) delete[] pm.Awt;
 if( pm.A ) delete[] pm.A;
 if( pm.XFs ) delete[] pm.XFs;
 if( pm.Falps ) delete[] pm.Falps;
 if( pm.G ) delete[] pm.G;
 if( pm.G0 ) delete[] pm.G0 ;
 if( pm.lnGam ) delete[] pm.lnGam;
 if( pm.lnGmo ) delete[] pm.lnGmo;
 if( pm.B ) delete[] pm.B;
 if( pm.U ) delete[] pm.U;
 if( pm.U_r ) delete[] pm.U_r;
 if( pm.C ) delete[] pm.C;
 if( pm.XF ) delete[] pm.XF;
 if( pm.YF ) delete[] pm.YF;
 if( pm.Falp ) delete[] pm.Falp;
 if( pm.X ) delete[] pm.X;
 if( pm.Y ) delete[] pm.Y;
 if( pm.XY ) delete[] pm.XY;
 if( pm.MU ) delete[] pm.MU;
 if( pm.EMU ) delete[] pm.EMU;
 if( pm.NMU ) delete[] pm.NMU;
 if( pm.W ) delete[] pm.W;
 if( pm.F ) delete[] pm.F;
 if( pm.F0 ) delete[] pm.F0;
 if( pm.YOF ) delete[] pm.YOF;

 if(   pm.SB ) delete[] pm.SB;
 if(   pm.SB1 ) delete[] pm.SB1;
 if(   pm.SFs ) delete[] pm.SFs;
 if(   pm.SM ) delete[] pm.SM;
 if(   pm.SF ) delete[] pm.SF;
 if(   pm.SM2 ) delete[] pm.SM2;
 if(   pm.SF2 ) delete[] pm.SF2;
 if(   pm.RLC ) delete[] pm.RLC;
 if(   pm.RSC ) delete[] pm.RSC;
 if(   pm.ICC ) delete[] pm.ICC;
 if(   pm.DCC ) delete[] pm.DCC;
 if(   pm.PHC ) delete[] pm.PHC;
 if(   pm.DCCW ) delete[] pm.DCCW;
 if( pm.lnGmM ) delete[] pm.lnGmM;
 if( pm.GEX ) delete[] pm.GEX;
 if( pm.FVOL ) delete[] pm.FVOL;
 if( pm.FWGT ) delete[] pm.FWGT;

   if( pm.Y_la ) delete[] pm.Y_la;
   if( pm.Y_w ) delete[] pm.Y_w;
   if( pm.Fx ) delete[] pm.Fx;
   if( pm.Wx ) delete[] pm.Wx;
   if( pm.VL ) delete[] pm.VL;
   if( pm.Gamma ) delete[] pm.Gamma;
   if( pm.lnGmf ) delete[] pm.lnGmf;
//   if( pm.D ) delete[] pm.D;

   // Part 2  not requited arrays

   if( pm.BF ) delete[] pm.BF;
if( pm.BFC ) delete[] pm.BFC;
   if( pm.XFA ) delete[] pm.XFA;
   if( pm.YFA ) delete[] pm.YFA;
   if( pm.LsMod ) delete[] pm.LsMod;
   if( pm.LsMdc ) delete[] pm.LsMdc;
   if( pm.IPx ) delete[] pm.IPx;
   if( pm.PMc ) delete[] pm.PMc;
   if( pm.DMc ) delete[] pm.DMc;
   if( pm.PUL ) delete[] pm.PUL;
   if( pm.PLL ) delete[] pm.PLL;
   if( pm.sMod ) delete[] pm.sMod;
   if( pm.RFLC ) delete[] pm.RFLC;
   if( pm.RFSC ) delete[] pm.RFSC;

   if( pm.Y_m ) delete[] pm.Y_m;
   if( pm.IC_m ) delete[] pm.IC_m;
   if( pm.IC_lm ) delete[] pm.IC_lm;
   if( pm.IC_wm ) delete[] pm.IC_wm;

   if( pm.Aalp ) delete[] pm.Aalp;
   if( pm.Xr0h0 ) delete[] pm.Xr0h0;

   if( pm.Sigw ) delete[] pm.Sigw;
   if( pm.Sigg ) delete[] pm.Sigg;

   if( pm.EZ ) delete[] pm.EZ;
   if( pm.Xcond ) delete[] pm.Xcond;
   if( pm.Xeps ) delete[] pm.Xeps;


if( pm.SATX ) delete[] pm.SATX;
   if( pm.SCM ) delete[] pm.SCM;
   if( pm.Nfsp ) delete[] pm.Nfsp;
   if( pm.MASDT ) delete[] pm.MASDT;
   if( pm.XcapA ) delete[] pm.XcapA;
   if( pm.XcapB ) delete[] pm.XcapB;
   if( pm.XcapD ) delete[] pm.XcapD;
   if( pm.XcapF ) delete[] pm.XcapF;
   if( pm.XdlA ) delete[] pm.XdlA;
   if( pm.XdlB ) delete[] pm.XdlB;
   if( pm.XdlD ) delete[] pm.XdlD;
   if( pm.XpsiA ) delete[] pm.XpsiA;
   if( pm.XpsiB ) delete[] pm.XpsiB;
   if( pm.XpsiD ) delete[] pm.XpsiD;
   if( pm.XlamA ) delete[] pm.XlamA;
   if( pm.Xetaf ) delete[] pm.Xetaf;
   if( pm.XetaA ) delete[] pm.XetaA;
   if( pm.XetaB ) delete[] pm.XetaB;
if( pm.XetaD ) delete[] pm.XetaD;
if( pm.MASDJ ) delete[] pm.MASDJ;
   if( pm.XFTS ) delete[] pm.XFTS;
if( pm.lnSAC ) delete[] pm.lnSAC;
if( pm.SATT ) delete[] pm.SATT;
if( pm.SM3 ) delete[] pm.SM3;
if( pm.DCC3 ) delete[] pm.DCC3;
if( pm.D ) delete[] pm.D;


   if( pm.Fug ) delete[] pm.Fug;
   if( pm.Fug_l ) delete[] pm.Fug_l;
   if( pm.Ppg_l ) delete[] pm.Ppg_l;

   // Part 3
   if( pm.Wb ) delete[] pm.Wb;
   if( pm.Wabs ) delete[] pm.Wabs;
   if( pm.Rion ) delete[] pm.Rion;
   if( pm.Qp ) delete[] pm.Qp;
   if( pm.Qd ) delete[] pm.Qd;

//  Added 16.11.2004 by Sveta
//    if( pm.sitE )     delete[] pm.sitE;
//    if( pm.sitXcat )  delete[] pm.sitXcat;
//    if( pm.sitXan )    delete[] pm.sitXan;

    // optimization 08/02/2007
    Free_internal();
}

#endif

void TMulti::to_text_file( const char *path )
{
    //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;


#ifndef IPMGEMPLUGIN
   PAalp = syp->PAalp;
   PSigm = syp->PSigm;
   EpsW = TProfil::pm->tpp->EpsW;
   RoW = TProfil::pm->tpp->RoW;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
   EpsW = EpsW_;
   RoW = RoW_;
#endif

  fstream ff(path, ios::out );
  ErrorIf( !ff.good() , path, "Fileopen error");

  ff << pm.stkey << endl;

  TPrintArrays  prar(ff);

  prar.writeArray( "Short_Const",  &pm.N, 38 );
  prar.writeArray(  "Double_Const",  &pm.TC, 55 );
  prar.writeArray(  "EpsW", &EpsW, 1);
  prar.writeArray(  "RoW", &RoW, 1);

   //dynamic values

    // Part 1

    /* need  always to alloc vectors */
  prar.writeArray(  "L1", pm.L1,  pm.FI);
  prar.writeArray(  "muk", pm.muk, pm.FI);
  prar.writeArray(  "mui", pm.mui, pm.N);
  prar.writeArray(  "muj", pm.muj,  pm.L);
  prar.writeArray(  "DUL", pm.DUL,  pm.L);
  prar.writeArray(  "DLL", pm.DLL,  pm.L);
  prar.writeArray(  "Vol", pm.Vol,  pm.L);
  prar.writeArray(  "Pparc", pm.Pparc,  pm.L);
  prar.writeArray(  "MM", pm.MM,  pm.L);
  prar.writeArray(  "Awt", pm.Awt, pm.N);
  prar.writeArray(  "A", pm.A,  pm.N*pm.L);
  prar.writeArray(  "XFs", pm.XFs, pm.FI);
  prar.writeArray(  "Falps", pm.Falps,  pm.FI);
  prar.writeArray(  "G", pm.G,  pm.L);
  prar.writeArray(  "G0", pm.G0,  pm.L);
  prar.writeArray(  "lnGam", pm.lnGam,  pm.L);
  prar.writeArray(  "lnGmo", pm.lnGmo,  pm.L);
  prar.writeArray(  "B", pm.B,  pm.N);
  prar.writeArray(  "U", pm.U,  pm.N);
  prar.writeArray(  "U_r", pm.U_r,  pm.N);
  prar.writeArray(  "C", pm.C,  pm.N);
  prar.writeArray(  "XF", pm.XF,  pm.FI);
  prar.writeArray(  "YF", pm.YF,  pm.FI);
  prar.writeArray(  "Falp", pm.Falp,  pm.FI);
  prar.writeArray(  "X", pm.X,  pm.L);
  prar.writeArray(  "Y", pm.Y,  pm.L);
  prar.writeArray(  "XY", pm.XY,  pm.L);
  prar.writeArray(  "MU", pm.MU,  pm.L);
  prar.writeArray(  "EMU", pm.EMU,  pm.L);
  prar.writeArray(  "NMU", pm.NMU,  pm.L);
  prar.writeArray(  "W", pm.W,  pm.L);
  prar.writeArray(  "F", pm.F,  pm.L);
  prar.writeArray(  "F0", pm.F0,  pm.L);
  prar.writeArray(  "YOF", pm.YOF,  pm.FI);


  prar.writeArray(  "lnGmM", pm.lnGmM,  pm.L);
  prar.writeArray(  "GEX", pm.GEX,  pm.L);
  prar.writeArray(  "FVOL", pm.FVOL,  pm.FI);
  prar.writeArray(  "FWGT", pm.FWGT,  pm.FI);

    if( pm.L > 0 )
    {
     prar.writeArray(  "Y_la", pm.Y_la,  pm.L);
     prar.writeArray(  "Y_w", pm.Y_w,  pm.L);
     prar.writeArray(  "Fx", pm.Fx,  pm.L);
     prar.writeArray(  "Wx", pm.Wx,  pm.L);
     prar.writeArray(  "VL", pm.VL, pm.L);
     prar.writeArray(  "Gamma", pm.Gamma,  pm.L);
     prar.writeArray(  "lnGmf", pm.lnGmf,  pm.L);
//     prar.writeArray(  "D", pm.D,  pm.L);
    }

   // Part 2  not always required arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
     prar.writeArray(  "BF", pm.BF,  pm.FIs*pm.N);
     prar.writeArray(  "BFC", pm.BFC, pm.N);
     prar.writeArray(  "XFA", pm.XFA,  pm.FIs);
     prar.writeArray(  "YFA", pm.YFA,  pm.FIs);
     prar.writeArray(  "LsMod", pm.LsMod, pm.FIs*3);
     prar.writeArray(  "LsMdc", pm.LsMdc, pm.FIs);
     int LsModSum;
     int LsIPxSum;
     int LsMdcSum;
     getLsModsum( LsModSum, LsIPxSum );
     getLsMdcsum( LsMdcSum );
     prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
     prar.writeArray(  "PMc", pm.PMc,  LsModSum);
     prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);

     prar.writeArray(  "PUL", pm.PUL,  pm.FIs);
     prar.writeArray(  "PLL", pm.PLL,  pm.FIs);

    }

    if( pm.LO > 1 )
    {
     prar.writeArray(  "Y_m", pm.Y_m,  pm.L);
     prar.writeArray(  "IC_m", pm.IC_m,  pm.N);
     prar.writeArray(  "IC_lm", pm.IC_lm,  pm.N);
     prar.writeArray(  "IC_wm", pm.IC_wm,  pm.N);
    }

    // dispersed and sorption phases
    if( PAalp != S_OFF )
    {
     prar.writeArray(  "Aalp", pm.Aalp, pm.FI);
     prar.writeArray(  "Xr0h0", &pm.Xr0h0[0][0],  pm.FI*2);
    }

   if( PSigm != S_OFF )
     prar.writeArray(  "Sigw", pm.Sigw,  pm.FI);

    if( PSigm != S_OFF )
     prar.writeArray(  "Sigg", pm.Sigg,  pm.FI);

    if( pm.E )
    {
     prar.writeArray(  "EZ", pm.EZ,  pm.L);
     prar.writeArray(  "Xcond", pm.Xcond,  pm.FI);
     prar.writeArray(  "Xeps", pm.Xeps,  pm.FI);
    }

    if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
     prar.writeArray(  "Nfsp", &pm.Nfsp[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "MASDT", &pm.MASDT[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XcapA", &pm.XcapA[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XcapB", &pm.XcapB[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XcapD", &pm.XcapD[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XcapF", &pm.XcapF[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XdlA", &pm.XdlA[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XdlB", &pm.XdlB[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XdlD", &pm.XdlD[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XpsiA", &pm.XpsiA[0][0],  pm.FIs*pm.FIat);
     prar.writeArray(  "XpsiB", &pm.XpsiB[0][0],  pm.FIs*pm.FIat);
     prar.writeArray(  "XpsiD", &pm.XpsiD[0][0],  pm.FIs*pm.FIat);
     prar.writeArray(  "XlamA", &pm.XlamA[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "Xetaf", &pm.Xetaf[0][0], pm.FIs*pm.FIat);
     prar.writeArray(  "XetaA", &pm.XetaA[0][0],  pm.FIs*pm.FIat);
     prar.writeArray(  "XetaB", &pm.XetaB[0][0],  pm.FIs*pm.FIat);
prar.writeArray(  "XetaD", &pm.XetaD[0][0],  pm.FIs*pm.FIat);   // added 12.09.05 KD
     prar.writeArray(  "XFTS", &pm.XFTS[0][0],  pm.FIs*pm.FIat);

     prar.writeArray(  "SATX", &pm.SATX[0][0], pm.Lads*4);
//     prar.writeArray(  "MASDJ", pm.MASDJ, pm.Ls);
     prar.writeArray(  "MASDJ", &pm.MASDJ[0][0], pm.Lads*DFCN);
     prar.writeArray(  "lnSAC", &pm.lnSAC[0][0],  pm.Lads*4);
     prar.writeArray(  "D", &pm.D[0][0], MST*MST);
    }

    if( pm.PG > 0 )
    {
     prar.writeArray(  "Fug", pm.Fug, pm.PG);
     prar.writeArray(  "Fug_l", pm.Fug_l, pm.PG);
     prar.writeArray(  "Ppg_l", pm.Ppg_l, pm.PG);
    }

    // Part 3

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
     prar.writeArray(  "Wb", pm.Wb, pm.Ls);
     prar.writeArray(  "Wabs", pm.Wabs, pm.Ls);
     prar.writeArray(  "Rion", pm.Rion, pm.Ls);

     prar.writeArray(  "Qp", pm.Qp,  pm.FIs*QPSIZE);
     prar.writeArray(  "Qd", pm.Qd,  pm.FIs*QDSIZE);

    }

//  Added 16.11.2004 by Sveta
//   if( pm.sitNcat*pm.sitNcat )
//    prar.writeArray(  "sitE", pm.sitE, pm.sitNcat*pm.sitNan );
//   if( pm.sitNcat )
//    prar.writeArray(  "sitXcat", pm.sitXcat, pm.sitNcat );
//   if( pm.sitNan )
//     prar.writeArray(  "sitXan", pm.sitXan, pm.sitNan );
}



