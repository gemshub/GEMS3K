//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi_file.cpp
/// Implementation of writing/reading IPM I/O files of GEMS3K
//
// Copyright (c) 2006-2012 S.Dmytriyeva, D.Kulik
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

#include <cmath>

#include "io_arrays.h"
#include "m_param.h"
#include "node.h"
#include "gdatastream.h"

void TMulti::getLsModsum( long int& LsModSum, long int& LsIPxSum )
{  LsModSum = 0;
   LsIPxSum = 0;
   for(long int i=0; i<pm.FIs; i++)
   {
     LsModSum += (pm.LsMod[i*3]*pm.LsMod[i*3+2]);
     LsIPxSum += (pm.LsMod[i*3]*pm.LsMod[i*3+1]);
   }
}


void TMulti::getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum )
{  LsMdcSum = 0;
   LsMsnSum = 0;
   LsSitSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       LsMdcSum += (pm.LsMdc[i*3]*pm.L1[i]);
       LsMsnSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]*pm.L1[i]);
       LsSitSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]);
   }
 }

// dimensions from LsPhl array
void TMulti::getLsPhlsum( long int& PhLinSum,long int& lPhcSum )
{  PhLinSum = 0;
   lPhcSum = 0;

   for(long int i=0; i<pm.FI; i++)
   {
       PhLinSum += (pm.LsPhl[i*2]);
       lPhcSum += (/*pm.LsPhl[i*2]**/pm.LsPhl[i*2+1]);

   }
 }

// dimensions from LsMdc2 array
void TMulti::getLsMdc2sum( long int& DQFcSum,long int& rcpcSum )
{  DQFcSum = 0;
   rcpcSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       DQFcSum += (pm.LsMdc2[i*3]*pm.L1[i]);
//       rcpcSum += (pm.LsMdc2[i*3+1]*pm.L1[i]);
   }
 }

// dimensions from LsISmo array
void TMulti::getLsISmosum( long int& IsoCtSum,long int& IsoScSum, long int& IsoPcSum,long int& xSMdSum )
{  IsoCtSum = 0;
   IsoScSum = 0;
   IsoPcSum = 0;
   xSMdSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       IsoCtSum += (pm.LsISmo[i*4]*2);
       IsoScSum += (pm.LsISmo[i*4]*pm.LsISmo[i*4+1]);
       IsoPcSum += (pm.LsISmo[i*4+2]*pm.L1[i]);
       xSMdSum += (pm.LsISmo[i*4+3]*pm.L1[i]);
   }
 }

// dimensions from LsESmo array
void TMulti::getLsESmosum( long int& EImcSum,long int& mCDcSum )
{  EImcSum = 0;
   mCDcSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       mCDcSum += (pm.LsESmo[i*4+2]*pm.L1[i]);
       EImcSum += (pm.LsESmo[i*4]*pm.LsESmo[i*4+1]);
   }
 }

// dimensions from LsKin array
void TMulti::getLsKinsum( long int& xSKrCSum,long int& ocPRkC_feSArC_Sum,
              long int& rpConCSum,long int& apConCSum, long int& AscpCSum )
{  xSKrCSum = 0;
   ocPRkC_feSArC_Sum = 0;
   rpConCSum = 0;
   apConCSum = 0;
   AscpCSum = 0;

   for(long int i=0; i<pm.FI; i++)
   {
       xSKrCSum += (pm.LsKin[i*6+1]);
       ocPRkC_feSArC_Sum += (pm.LsKin[i*6]);
       rpConCSum += (pm.LsKin[i*6]*pm.LsKin[i*6+2]);
       apConCSum += (pm.LsKin[i*6]*pm.LsKin[i*6+1]*pm.LsKin[i*6+3]);
       AscpCSum += (pm.LsKin[i*6+4]);
   }
 }

// dimensions from LsUpt array
void TMulti::getLsUptsum(long int& UMpcSum, long int& xICuCSum )
{
   UMpcSum = 0;
   for(long int i=0; i<pm.FIs; i++)
   {
       UMpcSum += (pm.LsUpt[i*2]*pm.L1[i]);
   }
   xICuCSum = 0;
   for(long int i=0; i<pm.FIs; i++)
       xICuCSum += pm.LsUpt[i*2+1]; // pm.L1[i];
 }

void TMulti::setPa( TProfil *prof)
{
    paTProfil = &prof->pa;
}

#ifdef IPMGEMPLUGIN

/// Output to "ipmlog.txt" file Warnings
long int TMulti::testMulti( )
{
  if( pm.MK || pm.PZ )
  {
    if( paTProfil->p.PSM == 2 )
    {
      fstream f_log(node->ipmLogFile().c_str(), ios::out|ios::app );
      f_log << "Warning " << pm.stkey << ": " <<  pm.errorCode << ":" << endl;
      f_log << pm.errorBuf << endl;
    }
   return 1L;
  }
  return 0L	;
}

#else

long int TMulti::testMulti()
{
  //MULTI *pmp = multi->GetPM();
  if( pm.MK || pm.PZ )
  {
   if( paTProfil->p.PSM >= 2 )
   {
     fstream f_log(node->ipmLogFile().c_str(), ios::out|ios::app );
     f_log << "Warning " << pm.stkey << ": " <<  pm.errorCode << ":" << endl;
     f_log << pm.errorBuf << endl;
   }
   if( showMss )
   {
           addErrorMessage(" \nContinue?");
      switch( vfQuestion3(0, pm.errorCode, pm.errorBuf,
                           "&Yes", "&No", "&Yes to All" ))
       {
       case VF3_3:
           showMss=0l;
       case VF3_1:
           break;
       case VF3_2:
           Error(pmp->errorCode, pmp->errorBuf);
       }
   }

   return 1L;
  }

  return 0L	;
}

#endif

//---------------------------------------------------------//
/// Set default information
void TMulti::set_def( long int /*q*/)
{
    //mem_cpy( &pm.PunE, "jjbC", 4 );
    pm.PunE = 'j';         // Units of energy  { j;  J c C N reserved }
    pm.PunV = 'j';         // Units of volume  { j;  c L a reserved }
    pm.PunP = 'b';        // Units of pressure  { b;  B p P A reserved }
    pm.PunT = 'C';         // Units of temperature  { C; K F reserved }

    // mem_set( &pm.N, 0, 36*sizeof(long int));
    pm.N = 0;        	// N - number of IC in IPM problem
    pm.NR = 0;       	// NR - dimensions of R matrix
    pm.L = 0;        	// L -   number of DC in IPM problem
    pm.Ls = 0;       	// Ls -   total number of DC in multi-component phases
    pm.LO = 0;       	// LO -   index of water-solvent in IPM DC list
    pm.PG = 0;       	// PG -   number of DC in gas phase
    pm.PSOL = 0;     	// PSOL - number of DC in liquid hydrocarbon phase
    pm.Lads = 0;     	// Lads - number of DC in sorption phases
    pm.FI = 0;       	// FI -   number of phases in IPM problem
    pm.FIs = 0;      	// FIs -   number of multicomponent phases
    pm.FIa = 0;      	// FIa -   number of sorption phases
    pm.FI1 = 0;     // FI1 -   number of phases present in eqstate
    pm.FI1s = 0;    // FI1s -   number of multicomponent phases present in eqstate
    pm.FI1a = 0;    // FI1a -   number of sorption phases present in eqstate
    pm.IT = 0;      // It - number of completed IPM iterations
    pm.E = 0;       // PE - flag of electroneutrality constraint { 0 1 }
    pm.PD = 0;      // PD - mode of calling CalculateActivityCoefficients() { 0 1 2 3 4 }
    pm.PV = 0;      // PV - flag of system volume constraint { 0 1 }
    pm.PLIM = 0;    // PU - flag of activation of DC/phase restrictions { 0 1 }
    pm.Ec = 0;    // CalculateActivityCoefficients() return code: 0 (OK) or 1 (error)
    pm.K2 = 0;    // Number of Selekt2() loops
    pm.PZ = 0;    // Indicator of IPM-2 precision algorithm activation    funT = 0; sysT = 0;
    pm.pNP = 0; //Mode of FIA selection: 0- automatic-LPP = 0; 1- old eqstate = 0; -1-user's choice
    pm.pESU = 0;  // Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pm.pIPN = 0;  // State of IPN-arrays:  0-create; 1-available; -1 remake
    pm.pBAL = 0;  // State of reloading CSD:  1- BAL only; 0-whole CSD
    pm.tMin = G_TP;  // Type of thermodynamic potential to minimize
    pm.pTPD = 0;  // State of reloading thermod data: 0- all  1 - G0 only  2 - no
    pm.pULR = 0;  // Start recalc kinetic constraints (0-do not = 0; 1-do )internal
pm.pKMM = 0;
    pm.ITaia = 0;  // Number of IPM iterations completed in AIA mode (renamed from pRR1)
    pm.FIat = 0;   // max. number of surface site types
    pm.MK = 0;     // PM return code: 0 - continue;  1 - converged
    pm.W1 = 0;     // internal IPM-2 indicator
    pm.is = 0;     // is - index of IC for IPN equations ( CalculateActivityCoefficients() )
    pm.js = 0;     // js - index of DC for IPN equations ( CalculateActivityCoefficients() )
    pm.next = 0;
    pm.sitNcat = 0;    // SIT: number of cations
    pm.sitNan = 0;     // SIT: number of anions    
pm.ITau = -1;  // current time, s (kinetics)
pm.kTau = 0.;  // current time, s (kinetics)
pm.kdT = 0.;   // current time step, s (kinetics)

    // mem_set( &pm.TC, 0, 54*sizeof(double));
    pm.TC = pm.TCc = 0.; 	// Temperature T = 0.; min.-max. (0 = 0.;2000 C)
    pm.T = pm.Tc = 0.;   	// T = 0.; min.-max. K
    pm.P = pm.Pc = 0.;   	// Pressure P = 0.; min.-max.(0 = 0.;10000 bar)
    pm.VX_ = pm.VXc = 0.;    // V(X) - volume of the system = 0.; min.-max. = 0.; cm3
    pm.GX_ = pm.GXc = 0.;    // Gibbs potential of the system G(X) = 0.; min.-max. (J)
    pm.AX_ = pm.AXc = 0.;    // Helmholtz potential of the system F(X) = 0.; reserved
    pm.UX_ = pm.UXc = 0.;  	// Internal energy of the system U(X) = 0.; reserved
    pm.HX_ = pm.HXc = 0.; 	// Total enthalpy of the system H(X) = 0.; reserved
    pm.SX_ = pm.SXc = 0.; 	// Total entropy of the system S(X) = 0.; reserved
    pm.CpX_ = pm.CpXc = 0.;  // reserved
    pm.CvX_ = pm.CvXc = 0.;  // reserved
    pm.TMols = 0.;         // input total moles in b vector before rescaling
    pm.SMols = 0.;         // Standart total moles (upscaled) {10000}
    pm.MBX = 0.;        // Total mass of the system = 0.; kg
    pm.FX = 0.;    	// Current Gibbs potential of the system in IPM = 0.; moles
    pm.IC = 0.;         // Effective molal ionic strength of aqueous electrolyte
    pm.pH = 0.;         // pH of aqueous solution
    pm.pe = 0.;         // pe of aqueous solution
    pm.Eh = 0.;         // Eh of aqueous solution = 0.; V
    pm.DHBM = 0.;       // Adjusted balance precision criterion (IPM-2 )
    pm.DSM = 0.;        // min value phase DS (IPM-2)
    pm.GWAT = 0.;       // used in ipm_gamma()
    pm.YMET = 0.;       // reserved
    fillValue( pm.denW, 0., 5 );
    fillValue( pm.denWg, 0., 5 );
    fillValue( pm.epsW, 0., 5 );
    fillValue( pm.epsWg, 0., 5 );
    pm.PCI = 0.;        // Current value of Dikin criterion of IPM convergence DK>=DX
    pm.DXM = 0.;         // IPM convergence criterion threshold DX (1e-5)
    pm.lnP = 0.;        // log Ptotal
    pm.RT = 0.;         // RT: 8.31451*T (J/mole/K)
    pm.FRT = 0.;        // F/RT = 0.; F - Faraday constant = 96485.309 C/mol
    pm.Yw = 0.;         // Current number of moles of solvent in aqueous phase
    pm.ln5551 = 0.;     // ln(55.508373) = 4.0165339
    pm.aqsTail = 0.;    // v_j asymmetry correction factor for aqueous species
    pm.lowPosNum = 0.;  // Minimum physical DC amount (1.66e-24 mol)
    pm.logXw = 0.;      // work variable
    pm.logYFk = 0.;     // work variable
    pm.YFk = 0.;        // Current number of moles in a multicomponent phase
    pm.FitVar[0] =pm.FitVar[1] = pm.FitVar[2]= pm.FitVar[3]= pm.FitVar[4] = 0.;
    fillValue( pm.Tai, 0., 4 );
    fillValue( pm.Pai, 0., 4 );
    pm.SizeFactor = 1.; // using in TNode class

    // pointers
    pm.sitNcat = 0;
    pm.sitNan = 0;
    pm.L1    = 0;
    pm.LsMod = 0;
    pm.LsMdc = 0;
    pm.mui   = 0;
    pm.muk   = 0;
    pm.muj   = 0;
    pm.SATX =0;
    pm.DUL   = 0;
    pm.DLL   = 0;
    pm.fDQF   = 0;
    pm.PUL   = 0;
    pm.PLL   = 0;
    pm.YOF   = 0;
    pm.PMc   = 0;
    pm.DMc   = 0;
    pm.MoiSN  = 0;
    pm.SitFr  = 0;
    pm.Vol   = 0;
    pm.VL    = 0;
    pm.MM    = 0;
    pm.H0    = 0;
    pm.A0    = 0;
    pm.U0    = 0;
    pm.S0    = 0;
    pm.Cp0   = 0;
    pm.Pparc = 0;
    pm.Y_m   = 0;
    pm.Y_la  = 0;
    pm.Y_w   = 0;
    pm.Gamma = 0;
    pm.lnGmf = 0;
    pm.lnGmM = 0;
    pm.EZ    = 0;
    pm.Wb    = 0;
    pm.Wabs  = 0;
    pm.Rion  = 0;
    pm.Aalp  = 0;
    pm.Sigw  = 0;
    pm.Sigg  = 0;
    pm.Nfsp  = 0;
    pm.MASDT = 0;
    pm.FVOL  = 0;
    pm.FWGT  = 0;
    pm.XcapA = 0;
    pm.XcapB = 0;
    pm.XcapD = 0;
    pm.XdlA  = 0;
    pm.XdlB  = 0;
    pm.XdlD  = 0;
    pm.XpsiA = 0;
    pm.XpsiB = 0;
    pm.XpsiD = 0;
    pm.Xr0h0 = 0;
    pm.XlamA = 0;
    pm.Xetaf = 0;
    pm.Xcond = 0;
    pm.Xeps  = 0;
    pm.Awt   = 0;
    pm.A     = 0;
    pm.XFs   = 0;
        pm.Falps = 0;
pm.GamFs = 0;
        pm.Fug   = 0;
        pm.Fug_l = 0;
        pm.Ppg_l = 0;
        pm.XFTS  = 0;
        pm.MASDJ = 0;
        pm.G     = 0;
        pm.G0    = 0;
        pm.lnGam = 0;
        pm.lnGmo = 0;
//        pm.lnSAT = 0;
        pm.lnSAC = 0;
        pm.B     = 0;
        pm.U     = 0;
        pm.Uc     = 0;
        pm.Uefd     = 0;
        pm.U_r   = 0;
        pm.C     = 0;
        pm.IC_m  = 0;
        pm.IC_lm = 0;
        pm.IC_wm = 0;
        pm.BF    = 0;
        pm.BFC    = 0;
        pm.XF    = 0;
        pm.YF    = 0;
        pm.XFA   = 0;
        pm.YFA   = 0;
        pm.Falp  = 0;
        pm.XetaA = 0;
        pm.XetaB = 0;
        pm.XetaD = 0;
        pm.X     = 0;
        pm.Y     = 0;
        pm.XY    = 0;
        pm.XU    = 0;
        pm.Qp    = 0;
        pm.Qd    = 0;
        pm.MU    = 0;
        pm.EMU   = 0;
        pm.NMU   = 0;
        pm.W     = 0;
        pm.Fx    = 0;
        pm.Wx    = 0;
        pm.F     = 0;
        pm.F0    = 0;
        pm.D     = 0;
     //   pm.R     = 0;
     //   pm.R1    = 0;
        pm.sMod  = 0;
        pm.dcMod  = 0;
        pm.SB    = 0;
        pm.SB1    = 0;
        pm.SM    = 0;
        pm.SF    = 0;
        pm.SFs   = 0;
        pm.pbuf  = 0;
        pm.RLC   = 0;
        pm.RSC   = 0;
        pm.RFLC  = 0;
        pm.RFSC  = 0;
        pm.ICC   = 0;
        pm.DCC   = 0;
        pm.PHC   = 0;
        pm.SCM   = 0;
        pm.SATT  = 0;
        pm.DCCW  = 0;
        pm.XcapF = 0;
        pm.SM2    = 0;
        pm.SM3    = 0;
        pm.SF2    = 0;
        pm.DCC3   = 0;
        pm.IPx = 0;
        pm.ITF =  pm.ITG = 0;
        pm.VPh = 0;
        pm.GPh = 0;
        pm.HPh = 0;
        pm.SPh = 0;
        pm.CPh = 0;
        pm.APh = 0;
        pm.UPh = 0;


// New phase stuff 06/06/12
        pm.LsMdc2  = 0;
        pm.LsPhl   = 0;
        pm.PhLin   = 0;
// TSolMod stuff
        pm.lPhc   = 0;
        pm.DQFc   = 0;
//        pm.rcpc   = 0;
        pm.lnDQFt   = 0;
        pm.lnRcpt   = 0;
        pm.lnExet   = 0;
        pm.lnCnft   = 0;
//TSorpMod & TKinMet stuff
        pm.SorMc   = 0;
// TSorpMod stuff
        pm.LsESmo   = 0;
        pm.LsISmo   = 0;
        pm.xSMd   = 0;
        pm.EImc   = 0;
        pm.mCDc   = 0;
        pm.IsoPc   = 0;
        pm.IsoSc   = 0;
        pm.lnScalT   = 0;
        pm.lnSACT   = 0;
        pm.lnGammF   = 0;
        pm.CTerms   = 0;
        pm.IsoCt   = 0;
// TKinMet stuff
        pm.LsKin   = 0;
        pm.LsUpt   = 0;
        pm.xSKrC   = 0;
        pm.ocPRkC   = 0;
        pm.feSArC   = 0;
        pm.rpConC   = 0;
        pm.apConC   = 0;
        pm.AscpC   = 0;
        pm.UMpcC   = 0;
        pm.kMod   = 0;
        // new
        pm.PfFact  = 0;
        pm.PrT   = 0;
        pm.PkT   = 0;
        pm.PvT   = 0;
        pm.emRd   = 0;
        pm.emDf   = 0;
        pm.xICuC = 0;
}

//---------------------------------------------------------//
/// Writing MULTI to binary file
void TMulti::to_file( GemDataStream& ff  )
{
   if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( GetName(), "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

   //static values
   char PAalp;
   char PSigm;


#ifndef IPMGEMPLUGIN

   PAalp = TSyst::sm->GetSY()->PAalp;
   PSigm = TSyst::sm->GetSY()->PSigm;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
#endif


   ff.writeArray(pm.stkey, sizeof(char)*(EQ_RKLEN+5));
   ff.writeArray( &pm.N, 39);
   ff.writeArray(&pm.TC, 53);
   ff << PAalp;
   ff << PSigm;
   ff.writeArray( pm.denW, 5);
   ff.writeArray( pm.denWg, 5);
   ff.writeArray( pm.epsW, 5);
   ff.writeArray( pm.epsWg, 5);

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
   ff.writeArray((char*)pm.dcMod, 6*pm.L);

   ff.writeArray( pm.RLC, pm.L);
   ff.writeArray( pm.RSC, pm.L);
   ff.writeArray( pm.ICC, pm.N);
   ff.writeArray( pm.DCC, pm.L);
   ff.writeArray( pm.PHC, pm.FI);
   ff.writeArray( pm.DCCW, pm.L);

   ff.writeArray( pm.lnGmM, pm.L);
   ff.writeArray( pm.fDQF,  pm.L);
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
      ff.writeArray(pm.PUL, pm.FIs);
      ff.writeArray(pm.PLL, pm.FIs);
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
      ff.writeArray((double *)pm.Xr0h0, pm.FI*2);
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
      ff.writeArray((double*)pm.Nfsp, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.MASDT, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XcapA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XcapB, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XcapD, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XcapF, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XdlA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XdlB, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XdlD, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiB, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XpsiD, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XlamA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.Xetaf, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XetaA, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XetaB, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XetaD, pm.FIs*pm.FIat);
      ff.writeArray((double*)pm.XFTS, pm.FIs*pm.FIat);

ff.writeArray((long int*)pm.SATX, pm.Lads*4);
ff.writeArray(pm.SATT, pm.Lads);
ff.writeArray((double*)pm.MASDJ, pm.Lads*DFCN);
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

    // Part 3 phases
     if( pm.FIs > 0 && pm.Ls > 0 )
     {

       ff.writeArray((char*)pm.sMod, 8*pm.FIs);
       ff.writeArray(pm.LsMod, pm.FIs*3);
       ff.writeArray(pm.LsMdc, pm.FIs*3);
       ff.writeArray(pm.LsMdc2, pm.FIs*3);
       ff.writeArray(pm.LsPhl, pm.FI*2);
       long int LsModSum;
       long int LsIPxSum;
       long int LsMdcSum;
       long int LsMsnSum;
       long int LsSitSum;
       getLsModsum( LsModSum, LsIPxSum );
       getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
       ff.writeArray(pm.IPx, LsIPxSum);
       ff.writeArray(pm.PMc, LsModSum);
       ff.writeArray(pm.DMc, LsMdcSum);
       ff.writeArray(pm.MoiSN, LsMsnSum);
       ff.writeArray(pm.SitFr, LsSitSum);
       long int DQFcSum, rcpcSum;
       getLsMdc2sum( DQFcSum, rcpcSum );
       ff.writeArray(pm.DQFc, DQFcSum);
//       ff.writeArray(pm.rcpc, rcpcSum);
       long int PhLinSum, lPhcSum;
       getLsPhlsum( PhLinSum,lPhcSum );
       ff.writeArray((long int *)pm.PhLin, PhLinSum*2);
       ff.writeArray(pm.lPhc, lPhcSum);
       ff.writeArray(  pm.lnDQFt, pm.Ls);
       ff.writeArray(  pm.lnRcpt, pm.Ls);
       ff.writeArray(  pm.lnExet, pm.Ls);
       ff.writeArray(  pm.lnCnft, pm.Ls);
       ff.writeArray( pm.SorMc, pm.FIs*16 );

       //TSorpMod
       ff.writeArray(  pm.LsISmo, pm.FIs*4);
       ff.writeArray(  pm.LsESmo, pm.FIs*4);
       long int EImcSum, mCDcSum;
       getLsESmosum( EImcSum, mCDcSum );
       long int IsoCtSum, IsoScSum;
       long int IsoPcSum, xSMdSum;
       getLsISmosum( IsoCtSum,IsoScSum,IsoPcSum, xSMdSum );
        ff.writeArray(  pm.xSMd, xSMdSum);
        ff.writeArray(  pm.IsoPc,  IsoPcSum);
        ff.writeArray(  pm.IsoSc, IsoScSum);
        ff.writeArray(  pm.IsoCt,  IsoCtSum);
        ff.writeArray(  pm.EImc, EImcSum);
        ff.writeArray(  pm.mCDc,  mCDcSum);

        ff.writeArray(  pm.lnScalT, pm.Ls);
        ff.writeArray(  pm.lnSACT, pm.Ls);
        ff.writeArray(  pm.lnGammF, pm.Ls);
        ff.writeArray(  pm.CTerms, pm.Ls);

        // TKinMet stuff
        ff.writeArray(  pm.kMod[0], pm.FI*6 );
        ff.writeArray(  pm.LsKin, pm.FI*6);
        ff.writeArray(  pm.LsUpt, pm.FIs*2);

        long int UMpcSum, xICuCSum;
        getLsUptsum( UMpcSum, xICuCSum );
        long int xSKrCSum, ocPRkC_feSArC_Sum;
        long int rpConCSum, apConCSum, AscpCSum;
        getLsKinsum( xSKrCSum, ocPRkC_feSArC_Sum, rpConCSum, apConCSum, AscpCSum );

        ff.writeArray( pm.xSKrC, xSKrCSum);
        ff.writeArray( &pm.ocPRkC[0][0],  ocPRkC_feSArC_Sum*2);
        ff.writeArray( pm.feSArC, ocPRkC_feSArC_Sum);
        ff.writeArray( pm.rpConC,  rpConCSum);
        ff.writeArray( pm.apConC, apConCSum);
        ff.writeArray( pm.AscpC,  AscpCSum);
        ff.writeArray( pm.UMpcC, UMpcSum);

        ff.writeArray(  pm.PfFact, pm.FI );
        ff.writeArray(  pm.xICuC, xICuCSum );
        ff.writeArray( pm.PfFact, pm.FI);
        ff.writeArray( pm.PrT, pm.FI);
        ff.writeArray( pm.PkT, pm.FI);
        ff.writeArray( pm.PvT, pm.FI);
        ff.writeArray( pm.emRd, pm.Ls);
        ff.writeArray( pm.emDf, pm.Ls);
     }

    // Part 4

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
      ff.writeArray(pm.Wb, pm.Ls);
      ff.writeArray(pm.Wabs, pm.Ls);
      ff.writeArray(pm.Rion, pm.Ls);

      ff.writeArray(pm.Qp, pm.FIs*QPSIZE);
      ff.writeArray(pm.Qd, pm.FIs*QDSIZE);
   }
   	ff.writeArray( pm.H0, pm.L);
   	ff.writeArray( pm.A0, pm.L);
   	ff.writeArray( pm.U0, pm.L);
   	ff.writeArray( pm.S0, pm.L);
   	ff.writeArray( pm.Cp0, pm.L);
}

/// Reading MULTI from binary file
void TMulti::from_file( GemDataStream& ff )
{
   //static values
   char PAalp;
   char PSigm;

   ff.readArray(pm.stkey, sizeof(char)*(EQ_RKLEN+5));
   ff.readArray( &pm.N, 39);
   ff.readArray(&pm.TC, 53);
   ff >> PAalp;
   ff >> PSigm;
   ff.readArray( pm.denW, 5);
   ff.readArray( pm.denWg, 5);
   ff.readArray( pm.epsW, 5);
   ff.readArray( pm.epsWg, 5);

#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#else
   dyn_new();
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
   ff.readArray((char*)pm.dcMod, 6*pm.L);

   ff.readArray( pm.RLC, pm.L);
   ff.readArray( pm.RSC, pm.L);
   ff.readArray( pm.ICC, pm.N);
   ff.readArray( pm.DCC, pm.L);
   ff.readArray( pm.PHC, pm.FI);
   ff.readArray( pm.DCCW, pm.L);

   ff.readArray( pm.lnGmM, pm.L);
   ff.readArray( pm.fDQF, pm.L);
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
      ff.readArray(pm.PUL, pm.FIs);
      ff.readArray(pm.PLL, pm.FIs);
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
      ff.readArray((double *)pm.Xr0h0, pm.FI*2);
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
      ff.readArray((double*)pm.Nfsp, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.MASDT, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XcapA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XcapB, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XcapD, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XcapF, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XdlA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XdlB, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XdlD, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiB, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XpsiD, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XlamA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.Xetaf, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XetaA, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XetaB, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XetaD, pm.FIs*pm.FIat);
      ff.readArray((double*)pm.XFTS, pm.FIs*pm.FIat);

ff.readArray((long int*)pm.SATX, pm.Lads*4);
ff.readArray(pm.SATT, pm.Lads);
ff.readArray((double*)pm.MASDJ, pm.Lads*DFCN);
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

    // Part 2  not requited arrays
     if( pm.FIs > 0 && pm.Ls > 0 )
     {
       ff.readArray((char*)pm.sMod, 8*pm.FIs);
       ff.readArray(pm.LsMod, pm.FIs*3);
       ff.readArray(pm.LsMdc, pm.FIs*3);
       ff.readArray(pm.LsMdc2, pm.FIs*3);
       ff.readArray(pm.LsPhl, pm.FI*2);

       long int LsModSum;
       long int LsIPxSum;
       long int LsMdcSum;
       long int LsMsnSum;
       long int LsSitSum;
       getLsModsum( LsModSum, LsIPxSum );
       getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
       long int DQFcSum, rcpcSum;
       getLsMdc2sum( DQFcSum, rcpcSum );
       long int PhLinSum, lPhcSum;
       getLsPhlsum( PhLinSum,lPhcSum );

 #ifdef IPMGEMPLUGIN
       pm.IPx = new long int[LsIPxSum];
       pm.PMc = new double[LsModSum];
       pm.DMc = new double[LsMdcSum];
       pm.MoiSN = new double[LsMsnSum];
       pm.SitFr = new double[LsSitSum];
       pm.DQFc = new double[DQFcSum];
//       pm.rcpc = new double[rcpcSum];
       pm.PhLin = new long int[PhLinSum][2];
       pm.lPhc = new double[lPhcSum];

 #else
       pm.IPx = (long int *)aObj[ o_wi_ipxpm ].Alloc(LsIPxSum, 1, L_);
       pm.PMc = (double *)aObj[ o_wi_pmc].Alloc( LsModSum, 1, D_);
       pm.DMc = (double *)aObj[ o_wi_dmc].Alloc( LsMdcSum, 1, D_ );
       pm.MoiSN = (double *)aObj[ o_wi_moisn].Alloc( LsMsnSum, 1, D_ );
       pm.SitFr  = (double *)aObj[ o_wo_sitfr ].Alloc( LsSitSum, 1, D_ );
       pm.DQFc = (double *)aObj[ o_wi_dqfc].Alloc( DQFcSum, 1, D_ );
//       pm.rcpc  = (double *)aObj[ o_wi_rcpc ].Alloc( rcpcSum, 1, D_ );
       pm.PhLin = (long int (*)[2])aObj[ o_wi_phlin].Alloc( PhLinSum, 2, L_ );
       pm.lPhc  = (double *)aObj[ o_wi_lphc ].Alloc( lPhcSum, 1, D_ );
 #endif

       ff.readArray(pm.IPx, LsIPxSum);
       ff.readArray(pm.PMc, LsModSum);
       ff.readArray(pm.DMc, LsMdcSum);
       ff.readArray(pm.MoiSN, LsMsnSum);
       ff.readArray(pm.SitFr, LsSitSum);
       ff.readArray(pm.DQFc, DQFcSum);
//       ff.readArray(pm.rcpc, rcpcSum);
       ff.readArray((long int *)pm.PhLin, PhLinSum*2);
       ff.readArray(pm.lPhc, lPhcSum);

       ff.readArray(  pm.lnDQFt, pm.Ls);
       ff.readArray(  pm.lnRcpt, pm.Ls);
       ff.readArray(  pm.lnExet, pm.Ls);
       ff.readArray(  pm.lnCnft, pm.Ls);
       ff.readArray( pm.SorMc, pm.FIs*16 );

       //TSorpMod
       ff.readArray(  pm.LsISmo, pm.FIs*4);
       ff.readArray(  pm.LsESmo, pm.FIs*4);
       long int EImcSum, mCDcSum;
       getLsESmosum( EImcSum, mCDcSum );

       long int IsoCtSum, IsoScSum;
       long int IsoPcSum, xSMdSum;
       getLsISmosum( IsoCtSum,IsoScSum,IsoPcSum, xSMdSum );

#ifdef IPMGEMPLUGIN
       pm.xSMd = new long int[xSMdSum];
       pm.IsoPc = new double[IsoPcSum];
       pm.IsoSc = new double[IsoScSum];
       pm.IsoCt = new char[IsoCtSum];
       pm.EImc = new double[EImcSum];
       pm.mCDc = new double[mCDcSum];
#else
        pm.xSMd = (long int*)aObj[ o_wi_xsmd].Alloc( xSMdSum, 1, L_ );
        pm.IsoPc = (double*)aObj[ o_wi_isopc].Alloc( IsoPcSum, 1, D_ );
        pm.IsoSc = (double*)aObj[ o_wi_isosc].Alloc( IsoScSum, 1, D_ );
        pm.IsoCt = (char*)aObj[ o_wi_isoct].Alloc( IsoCtSum, 1, A_ );
        pm.EImc = (double*)aObj[ o_wi_eimc].Alloc( EImcSum, 1, D_ );
        pm.mCDc = (double*)aObj[ o_wi_mcdc].Alloc( mCDcSum, 1, D_ );
#endif
        ff.readArray(  pm.xSMd, xSMdSum);
        ff.readArray(  pm.IsoPc,  IsoPcSum);
        ff.readArray(  pm.IsoSc, IsoScSum);
        ff.readArray(  pm.IsoCt,  IsoCtSum);
        ff.readArray(  pm.EImc, EImcSum);
        ff.readArray(  pm.mCDc,  mCDcSum);

        ff.readArray(  pm.lnScalT, pm.Ls);
        ff.readArray(  pm.lnSACT, pm.Ls);
        ff.readArray(  pm.lnGammF, pm.Ls);
        ff.readArray(  pm.CTerms, pm.Ls);

        // TKinMet stuff
        ff.readArray(  pm.kMod[0], pm.FI*6 );
        ff.readArray(  pm.LsKin, pm.FI*6);
        ff.readArray(  pm.LsUpt, pm.FIs*2);

        long int UMpcSum, xICuCSum;
        getLsUptsum( UMpcSum, xICuCSum );
        long int xSKrCSum, ocPRkC_feSArC_Sum;
        long int rpConCSum, apConCSum, AscpCSum;
        getLsKinsum( xSKrCSum, ocPRkC_feSArC_Sum, rpConCSum, apConCSum, AscpCSum );

#ifdef IPMGEMPLUGIN
        pm.xSKrC = new long int[xSKrCSum];
        pm.ocPRkC = new long int[ocPRkC_feSArC_Sum][2];
        pm.feSArC = new double[ocPRkC_feSArC_Sum];
        pm.rpConC = new double[rpConCSum];
        pm.apConC = new double[apConCSum];
        pm.AscpC = new double[AscpCSum];
        pm.UMpcC = new double[UMpcSum];
        pm.xICuC = new long int[xICuCSum];
#else
        pm.xSKrC = (long int*)aObj[ o_wi_jcrdc].Alloc( xSKrCSum, 1, L_ );
        pm.ocPRkC = (long int(*)[2])aObj[ o_wi_ocprkc].Alloc( ocPRkC_feSArC_Sum, 2, L_ );
        pm.feSArC = (double*)aObj[ o_wi_fsac].Alloc( ocPRkC_feSArC_Sum, 1, D_ );
        pm.rpConC = (double*)aObj[ o_wi_krpc].Alloc( rpConCSum, 1, D_ );
        pm.apConC = (double*)aObj[ o_wi_apconc].Alloc( apConCSum, 1, D_ );
        pm.AscpC = (double*)aObj[ o_wi_ascpc].Alloc( AscpCSum, 1, D_ );
        pm.UMpcC = (double*)aObj[ o_wi_umpc].Alloc( UMpcSum, 1, D_ );
        pm.xICuC = (long int *)aObj[o_wi_xicuc ].Alloc( xICuCSum, 1, L_ );
#endif
        ff.readArray( pm.xSKrC, xSKrCSum);
        ff.readArray( &pm.ocPRkC[0][0],  ocPRkC_feSArC_Sum*2);
        ff.readArray( pm.feSArC, ocPRkC_feSArC_Sum);
        ff.readArray( pm.rpConC,  rpConCSum);
        ff.readArray( pm.apConC, apConCSum);
        ff.readArray( pm.AscpC,  AscpCSum);
        ff.readArray( pm.UMpcC, UMpcSum);

        ff.readArray(  pm.PfFact, pm.FI );
        ff.readArray(  pm.xICuC, xICuCSum );
        ff.readArray( pm.PfFact, pm.FI);
        ff.readArray( pm.PrT, pm.FI);
        ff.readArray( pm.PkT, pm.FI);
        ff.readArray( pm.PvT, pm.FI);
        ff.readArray( pm.emRd, pm.Ls);
        ff.readArray( pm.emDf, pm.Ls);
     }

     // Part 4

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
      ff.readArray(pm.Wb, pm.Ls);
      ff.readArray(pm.Wabs, pm.Ls);
      ff.readArray(pm.Rion, pm.Ls);

      ff.readArray(pm.Qp, pm.FIs*QPSIZE);
      ff.readArray(pm.Qd, pm.FIs*QDSIZE);
   }
   	ff.readArray( pm.H0, pm.L);
   	ff.readArray( pm.A0, pm.L);
   	ff.readArray( pm.U0, pm.L);
   	ff.readArray( pm.S0, pm.L);
   	ff.readArray( pm.Cp0, pm.L);
}


/// Writing structure MULTI ( free format file  )
void TMulti::to_text_file( const char *path, bool append )
{
    //static values
   char PAalp;
   char PSigm;

#ifndef IPMGEMPLUGIN
   PAalp = TSyst::sm->GetSY()->PAalp;
   PSigm = TSyst::sm->GetSY()->PSigm;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
#endif

   ios::openmode mod = ios::out;
    if( append )
     mod = ios::out|ios::app;
  fstream ff(path, mod );
  ErrorIf( !ff.good() , path, "Fileopen error");

  if( append )
   ff << "\nNext record" << endl;
  ff << pm.stkey << endl;
//  TProfil::pm->pa.p.write(ff);

  TPrintArrays  prar(0,0,ff);

  prar.writeArray( "Short_PARAM",  &paTProfil->p.PC, 10L );
  prar.writeArray( "Double_PARAM",  &paTProfil->p.DG, 28L );
  prar.writeArray( "Short_Const",  &pm.N, 39L );
  prar.writeArray(  "Double_Const",  &pm.TC, 53, 20 );
  prar.writeArray(  "EpsW", pm.epsW, 5);
  prar.writeArray(  "EpsWg", pm.epsWg, 5);
  prar.writeArray(  "DenW", pm.denW, 5);
  prar.writeArray(  "DenWg", pm.denWg, 5);
  ff << endl << "Error Code " << pm.errorCode << endl;
  ff << "Error Message" << pm.errorBuf << endl;

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
  prar.writeArray(  "Uc", &pm.Uc[0][0],  pm.N*2);
  prar.writeArray(  "Uefd", pm.Uefd,  pm.N);
  prar.writeArray(  "U_r", pm.U_r,  pm.N);
  prar.writeArray(  "C", pm.C,  pm.N);
  prar.writeArray(  "XF", pm.XF,  pm.FI);
  prar.writeArray(  "YF", pm.YF,  pm.FI);
  prar.writeArray(  "Falp", pm.Falp,  pm.FI);
  prar.writeArray(  "X", pm.X,  pm.L);
  prar.writeArray(  "Y", pm.Y,  pm.L);
  prar.writeArray(  "XY", pm.XY,  pm.L);
  prar.writeArray(  "XU", pm.XU,  pm.L);
  prar.writeArray(  "MU", pm.MU,  pm.L);
  prar.writeArray(  "EMU", pm.EMU,  pm.L);
  prar.writeArray(  "NMU", pm.NMU,  pm.L);
  prar.writeArray(  "W", pm.W,  pm.L);
  prar.writeArray(  "F", pm.F,  pm.L);
  prar.writeArray(  "F0", pm.F0,  pm.L);
  prar.writeArray(  "YOF", pm.YOF,  pm.FI);


  prar.writeArray(  "lnGmM", pm.lnGmM,  pm.L);
  prar.writeArray(  "fDQF", pm.fDQF,  pm.L);
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
     prar.writeArray(  "XetaD", &pm.XetaD[0][0],  pm.FIs*pm.FIat);
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

    // Part 3  new Phase definition
     if( pm.FIs > 0 && pm.Ls > 0 )
     {
      prar.writeArray(  "sMod", &pm.sMod[0][0], pm.FIs,8L);
      prar.writeArray(  "LsMod", pm.LsMod, pm.FIs*3);
      long int LsModSum;
      long int LsIPxSum;
      getLsModsum( LsModSum, LsIPxSum );
      prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
      prar.writeArray(  "PMc", pm.PMc,  LsModSum);
      long int LsMdcSum;
      long int LsMsnSum;
      long int LsSitSum;
      prar.writeArray(  "LsMdc", pm.LsMdc, pm.FIs*3);
      getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
      prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);
      prar.writeArray(  "MoiSN", pm.MoiSN,  LsMsnSum);
      prar.writeArray(  "SitFr", pm.SitFr,  LsSitSum);
      long int DQFcSum, rcpcSum;
      getLsMdc2sum( DQFcSum, rcpcSum );
      prar.writeArray(  "LsMdc2", pm.LsMdc2, pm.FIs*3);
      prar.writeArray(  "DQFc", pm.DQFc,  DQFcSum);
//      prar.writeArray(  "rcpc", pm.rcpc,  rcpcSum);
      long int PhLinSum, lPhcSum;
      getLsPhlsum( PhLinSum,lPhcSum );
      prar.writeArray(  "LsPhl", pm.LsPhl, pm.FI*2);
      prar.writeArray(  "PhLin", &pm.PhLin[0][0], PhLinSum*2);
      prar.writeArray(  "lPhc", pm.lPhc,  lPhcSum);

      prar.writeArray(  "lnDQFt", pm.lnDQFt, pm.Ls);
      prar.writeArray(  "lnRcpt", pm.lnRcpt, pm.Ls);
      prar.writeArray(  "lnExet", pm.lnExet, pm.Ls);
      prar.writeArray(  "lnCnft", pm.lnCnft, pm.Ls);

      prar.writeArray(  "SorMc", pm.SorMc, pm.FIs*16, 16L);

       // TSorpMod stuff
      long int IsoCtSum, IsoScSum;
      long int IsoPcSum, xSMdSum;
      getLsISmosum( IsoCtSum,IsoScSum,IsoPcSum, xSMdSum );
      prar.writeArray(  "LsISmo", pm.LsISmo, pm.FIs*4);
      prar.writeArray(  "xSMd", pm.xSMd, xSMdSum);
      prar.writeArray(  "IsoPc", pm.IsoPc,  IsoPcSum);
      prar.writeArray(  "IsoSc", pm.IsoSc, IsoScSum);
      prar.writeArray(  "IsoCt", pm.IsoCt,  IsoCtSum, 1L);
      long int EImcSum, mCDcSum;
      getLsESmosum( EImcSum, mCDcSum );
      prar.writeArray(  "LsESmo", pm.LsESmo, pm.FIs*4);
      prar.writeArray(  "EImc", pm.EImc, EImcSum);
      prar.writeArray(  "mCDc", pm.mCDc,  mCDcSum);

      prar.writeArray(  "lnScalT", pm.lnScalT, pm.Ls);
      prar.writeArray(  "lnSACT", pm.lnSACT, pm.Ls);
      prar.writeArray(  "lnGammF", pm.lnGammF, pm.Ls);
      prar.writeArray(  "CTerms", pm.CTerms, pm.Ls);

      // TKinMet stuff
      prar.writeArray(  "kMod", &pm.kMod[0][0], pm.FI, 6L);
      long int xSKrCSum, ocPRkC_feSArC_Sum;
      long int rpConCSum, apConCSum, AscpCSum;
      getLsKinsum( xSKrCSum, ocPRkC_feSArC_Sum, rpConCSum, apConCSum, AscpCSum );
      prar.writeArray(  "LsKin", pm.LsKin, pm.FI*6);
      prar.writeArray(  "xSKrC", pm.xSKrC, xSKrCSum);
      prar.writeArray(  "ocPRkC", &pm.ocPRkC[0][0],  ocPRkC_feSArC_Sum*2);
      prar.writeArray(  "feSArC", pm.feSArC, ocPRkC_feSArC_Sum);
      prar.writeArray(  "rpConC", pm.rpConC,  rpConCSum);
      prar.writeArray(  "apConC", pm.apConC, apConCSum);
      prar.writeArray(  "AscpC", pm.AscpC,  AscpCSum);
      long int UMpcSum, xICuCSum;
      getLsUptsum( UMpcSum, xICuCSum );
      prar.writeArray(  "LsUpt", pm.LsUpt, pm.FIs*2);
      prar.writeArray(  "UMpcC", pm.UMpcC, UMpcSum);

      prar.writeArray(  "PfFact", pm.PfFact, pm.FI);
      prar.writeArray(  "PrT", pm.PrT, pm.FI);
      prar.writeArray(  "PkT", pm.PkT, pm.FI);
      prar.writeArray(  "PvT", pm.PvT, pm.FI);
      prar.writeArray(  "emRd", pm.emRd, pm.Ls);
      prar.writeArray(  "emDf", pm.emDf, pm.Ls);
      if( pm.xICuC )
      {
        prar.writeArray(  "xICuC", pm.xICuC, xICuCSum);
      }
   }

    // Part 4

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
     prar.writeArray(  "Wb", pm.Wb, pm.Ls);
     prar.writeArray(  "Wabs", pm.Wabs, pm.Ls);
     prar.writeArray(  "Rion", pm.Rion, pm.Ls);

     prar.writeArray(  "Qp", pm.Qp,  pm.FIs*QPSIZE);
     prar.writeArray(  "Qd", pm.Qd,  pm.FIs*QDSIZE);

    }

    if(pm.H0)
    	prar.writeArray("H0",pm.H0, pm.L);
    if(pm.A0)
    	prar.writeArray("A0",pm.A0, pm.L);
    if(pm.U0)
    	prar.writeArray("U0",pm.U0, pm.L);
    if(pm.S0)
    	prar.writeArray("S0",pm.S0, pm.L);
    if(pm.Cp0)
    	prar.writeArray("Cp0",pm.Cp0, pm.L);

    prar.writeArray(  "VPh", &pm.VPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "GPh", &pm.GPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "HPh", &pm.HPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "SPh", &pm.SPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "CPh", &pm.CPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "APh", &pm.APh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "UPh", &pm.UPh[0][0], pm.FIs*MIXPHPROPS);

}

//--------------------- End of ms_multi_file.cpp ---------------------------


