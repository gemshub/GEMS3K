//-------------------------------------------------------------------
// $Id$
//
/// \file tsolmod_multi_file.cpp
/// Implementation of writing/reading IPM I/O files of GEMS3K
//
// Copyright (c) 2023 S.Dmytriyeva, D.Kulik
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

#include "v_service.h"
#include "io_template.h"
#include "io_keyvalue.h"
#include "gdatastream.h"
#include "jsonconfig.h"
#include "tsolmod_multi.h"

void TSolModMulti::getLsModsum( long int& LsModSum, long int& LsIPxSum )
{  LsModSum = 0;
   LsIPxSum = 0;
   for(long int i=0; i<pm.FIs; i++)
   {
     LsModSum += (pm.LsMod[i*3]*pm.LsMod[i*3+2]);
     LsIPxSum += (pm.LsMod[i*3]*pm.LsMod[i*3+1]);
   }
}


void TSolModMulti::getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum )
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
void TSolModMulti::getLsPhlsum( long int& PhLinSum,long int& lPhcSum )
{  PhLinSum = 0;
   lPhcSum = 0;

   for(long int i=0; i<pm.FI; i++)
   {
       PhLinSum += (pm.LsPhl[i*2]);
       lPhcSum += (/*pm.LsPhl[i*2]**/pm.LsPhl[i*2+1]);

   }
 }

// dimensions from LsMdc2 array
void TSolModMulti::getLsMdc2sum( long int& DQFcSum,long int& rcpcSum )
{  DQFcSum = 0;
   rcpcSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       DQFcSum += (pm.LsMdc2[i*3]*pm.L1[i]);
//       rcpcSum += (pm.LsMdc2[i*3+1]*pm.L1[i]);
   }
 }

// dimensions from LsISmo array
void TSolModMulti::getLsISmosum( long int& IsoCtSum,long int& IsoScSum, long int& IsoPcSum,long int& xSMdSum )
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
void TSolModMulti::getLsESmosum( long int& EImcSum,long int& mCDcSum )
{  EImcSum = 0;
   mCDcSum = 0;

   for(long int i=0; i<pm.FIs; i++)
   {
       mCDcSum += (pm.LsESmo[i*4+2]*pm.L1[i]);
       EImcSum += (pm.LsESmo[i*4]*pm.LsESmo[i*4+1]);
   }
 }

// dimensions from LsKin array
void TSolModMulti::getLsKinsum( long int& xSKrCSum,long int& ocPRkC_feSArC_Sum,
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
void TSolModMulti::getLsUptsum(long int& UMpcSum, long int& xICuCSum )
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

//---------------------------------------------------------//
/// Set default information
void TSolModMulti::set_def( int )
{
    //mem_cpy( &pm.PunE, "jjbC", 4 );
    fillValue( pm.stkey, '\0', EQ_RKLEN);
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
    pm.GWAT = 55.50837344;       // used in ipm_gamma()
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
    pm.L1    = nullptr;
    pm.LsMod = nullptr;
    pm.LsMdc = nullptr;
    pm.mui   = nullptr;
    pm.muk   = nullptr;
    pm.muj   = nullptr;
    pm.SATX = nullptr;
    pm.DUL   = nullptr;
    pm.DLL   = nullptr;
    pm.fDQF   = nullptr;
    pm.PUL   = nullptr;
    pm.PLL   = nullptr;
    pm.YOF   = nullptr;
    pm.PMc   = nullptr;
    pm.DMc   = nullptr;
    pm.MoiSN  = nullptr;
    pm.SitFr  = nullptr;
    pm.Vol   = nullptr;
    pm.VL    = nullptr;
    pm.MM    = nullptr;
    pm.H0    = nullptr;
    pm.A0    = nullptr;
    pm.U0    = nullptr;
    pm.S0    = nullptr;
    pm.Cp0   = nullptr;
    pm.Pparc = nullptr;
    pm.Y_m   = nullptr;
    pm.Y_la  = nullptr;
    pm.Y_w   = nullptr;
    pm.Gamma = nullptr;
    pm.lnGmf = nullptr;
    pm.lnGmM = nullptr;
    pm.EZ    = nullptr;
    pm.Wb    = nullptr;
    pm.Wabs  = nullptr;
    pm.Rion  = nullptr;
    pm.Aalp  = nullptr;
    pm.Sigw  = nullptr;
    pm.Sigg  = nullptr;
    pm.Nfsp  = nullptr;
    pm.MASDT = nullptr;
    pm.FVOL  = nullptr;
    pm.FWGT  = nullptr;
    pm.XcapA = nullptr;
    pm.XcapB = nullptr;
    pm.XcapD = nullptr;
    pm.XdlA  = nullptr;
    pm.XdlB  = nullptr;
    pm.XdlD  = nullptr;
    pm.XpsiA = nullptr;
    pm.XpsiB = nullptr;
    pm.XpsiD = nullptr;
    pm.Xr0h0 = nullptr;
    pm.XlamA = nullptr;
    pm.Xetaf = nullptr;
    pm.Xcond = nullptr;
    pm.Xeps  = nullptr;
    pm.Awt   = nullptr;
    pm.A     = nullptr;
    pm.XFs   = nullptr;
        pm.Falps = nullptr;
pm.GamFs = nullptr;
        pm.Fug   = nullptr;
        pm.Fug_l = nullptr;
        pm.Ppg_l = nullptr;
        pm.XFTS  = nullptr;
        pm.MASDJ = nullptr;
        pm.G     = nullptr;
        pm.G0    = nullptr;
        pm.lnGam = nullptr;
        pm.lnGmo = nullptr;
//        pm.lnSAT = nullptr;
        pm.lnSAC = nullptr;
        pm.B     = nullptr;
        pm.U     = nullptr;
        pm.Uc     = nullptr;
        pm.Uefd     = nullptr;
        pm.U_r   = nullptr;
        pm.C     = nullptr;
        pm.IC_m  = nullptr;
        pm.IC_lm = nullptr;
        pm.IC_wm = nullptr;
        pm.BF    = nullptr;
        pm.BFC    = nullptr;
        pm.XF    = nullptr;
        pm.YF    = nullptr;
        pm.XFA   = nullptr;
        pm.YFA   = nullptr;
        pm.Falp  = nullptr;
        pm.XetaA = nullptr;
        pm.XetaB = nullptr;
        pm.XetaD = nullptr;
        pm.X     = nullptr;
        pm.Y     = nullptr;
        pm.XY    = nullptr;
        pm.XU    = nullptr;
        pm.Qp    = nullptr;
        pm.Qd    = nullptr;
        pm.MU    = nullptr;
        pm.EMU   = nullptr;
        pm.NMU   = nullptr;
        pm.W     = nullptr;
        pm.Fx    = nullptr;
        pm.Wx    = nullptr;
        pm.F     = nullptr;
        pm.F0    = nullptr;
        pm.D     = nullptr;
     //   pm.R     = nullptr;
     //   pm.R1    = nullptr;
        pm.sMod  = nullptr;
        pm.dcMod  = nullptr;
        pm.SB    = nullptr;
        pm.SB1    = nullptr;
        pm.SM    = nullptr;
        pm.SF    = nullptr;
        pm.SFs   = nullptr;
        pm.pbuf  = nullptr;
        pm.RLC   = nullptr;
        pm.RSC   = nullptr;
        pm.RFLC  = nullptr;
        pm.RFSC  = nullptr;
        pm.ICC   = nullptr;
        pm.DCC   = nullptr;
        pm.PHC   = nullptr;
        pm.SCM   = nullptr;
        pm.SATT  = nullptr;
        pm.DCCW  = nullptr;
        pm.XcapF = nullptr;
        pm.SM2    = nullptr;
        pm.SM3    = nullptr;
        pm.SF2    = nullptr;
        pm.DCC3   = nullptr;
        pm.IPx = nullptr;
        pm.ITF =  pm.ITG = 0;
        pm.VPh = nullptr;
        pm.GPh = nullptr;
        pm.HPh = nullptr;
        pm.SPh = nullptr;
        pm.CPh = nullptr;
        pm.APh = nullptr;
        pm.UPh = nullptr;


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

/// Writing structure MULTI ( free format file  )
void TSolModMulti::to_text_file( const char *path, bool append )
{
    //static values
   char PAalp;
   char PSigm;
   get_PAalp_PSigm( PAalp, PSigm);

   std::ios::openmode mod = std::ios::out;
    if( append )
     mod = std::ios::out|std::ios::app;
  std::fstream ff(GemsSettings::with_directory(path), mod );
  ErrorIf( !ff.good() , path, "Fileopen error");

  io_formats::KeyValueWrite out_format( ff );
  out_format.put_head( "", "ipm");
  io_formats::TPrintArrays<io_formats::KeyValueWrite>  prar( 0, {}, out_format );

  if( append )
   prar.writeComment( true,"\nNext record" );
  prar.writeComment( true, char_array_to_string(pm.stkey, EQ_RKLEN)+"\n" );
  //  TProfil::pm->pa.p.write(ff);

  prar.writeArray( "Short_PARAM",  &base_param()->PC, 10L );
  prar.writeArray( "Double_PARAM",  &base_param()->DG, 28L );
  prar.writeArray( "Short_Const",  &pm.N, 39L );
  prar.writeArray(  "Double_Const",  &pm.TC, 53, 20 );
  prar.writeArray(  "Add_Double_Const",  &pm.XwMinM, 12, 20 );
  prar.writeArray(  "EpsW", pm.epsW, 5);
  prar.writeArray(  "EpsWg", pm.epsWg, 5);
  prar.writeArray(  "DenW", pm.denW, 5);
  prar.writeArray(  "DenWg", pm.denWg, 5);
  prar.writeComment( true, std::string("Error Code ")+ pm.errorCode);
  prar.writeComment( true, std::string("Error Message") + pm.errorBuf);

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

void TSolModMulti::get_PAalp_PSigm( char& PAalp, char& PSigm)
{
    PAalp = PAalp_;
    PSigm = PSigm_;
}

void TSolModMulti::alloc_IPx( long int LsIPxSum )
{
    if( pm.IPx ) delete[] pm.IPx;
    pm.IPx = new long int[ LsIPxSum];
}

void TSolModMulti::alloc_PMc( long int LsModSum )
{
    if( pm.PMc ) delete[] pm.PMc;
    pm.PMc = new double[LsModSum];
}

void TSolModMulti::alloc_DMc( long int LsMdcSum )
{
    if( pm.DMc ) delete[] pm.DMc;
    pm.DMc = new double[LsMdcSum];
}

void TSolModMulti::alloc_MoiSN( long int LsMsnSum )
{
    if(pm.MoiSN) delete[] pm.MoiSN;
    pm.MoiSN = new double[LsMsnSum];
}

void TSolModMulti::alloc_SitFr( long int LsSitSum )
{
    if(pm.SitFr) delete[] pm.SitFr;
    pm.SitFr = new double[LsSitSum];
}

void TSolModMulti::alloc_DQFc( long int DQFcSum )
{
    if(pm.DQFc) delete[] pm.DQFc;
    pm.DQFc = new double[DQFcSum];
}

void TSolModMulti::alloc_PhLin( long int PhLinSum )
{
    if(pm.PhLin) delete[] pm.PhLin;
    pm.PhLin = new long int[PhLinSum][2];
}

void TSolModMulti::alloc_lPhc( long int lPhcSum )
{
    if(pm.lPhc) delete[] pm.lPhc;
    pm.lPhc = new double[lPhcSum];
}

void TSolModMulti::alloc_xSMd( long int xSMdSum )
{
    if(pm.xSMd) delete[] pm.xSMd;
    pm.xSMd = new long int[xSMdSum];
}

void TSolModMulti::alloc_IsoPc( long int IsoPcSum )
{
    if(pm.IsoPc) delete[] pm.IsoPc;
    pm.IsoPc = new double[IsoPcSum];
}

void TSolModMulti::alloc_IsoSc( long int IsoScSum )
{
    if(pm.IsoSc) delete[] pm.IsoSc;
    pm.IsoSc = new double[IsoScSum];
}

void TSolModMulti::alloc_IsoCt( long int IsoCtSum )
{
    if(pm.IsoCt) delete[] pm.IsoCt;
    pm.IsoCt = new char[IsoCtSum];
}

void TSolModMulti::alloc_EImc( long int EImcSum )
{
    if(pm.EImc) delete[] pm.EImc;
    pm.EImc = new double[EImcSum];
}

void TSolModMulti::alloc_mCDc( long int mCDcSum )
{
    if(pm.mCDc) delete[] pm.mCDc;
    pm.mCDc = new double[mCDcSum];
}

void TSolModMulti::alloc_xSKrC( long int xSKrCSum )
{
    if(pm.xSKrC) delete[] pm.xSKrC;
    pm.xSKrC = new long int[xSKrCSum];
}

void TSolModMulti::alloc_ocPRkC( long int ocPRkC_feSArC_Sum )
{
    if(pm.ocPRkC) delete[] pm.ocPRkC;
    pm.ocPRkC = new long int[ocPRkC_feSArC_Sum][2];
}

void TSolModMulti::alloc_feSArC( long int ocPRkC_feSArC_Sum )
{
    if(pm.feSArC) delete[] pm.feSArC;
    pm.feSArC = new double[ocPRkC_feSArC_Sum];
}

void TSolModMulti::alloc_rpConC( long int rpConCSum )
{
    if(pm.rpConC) delete[] pm.rpConC;
    pm.rpConC = new double[rpConCSum];
}

void TSolModMulti::alloc_apConC( long int apConCSum )
{
    if(pm.apConC) delete[] pm.apConC;
    pm.apConC = new double[apConCSum];
}

void TSolModMulti::alloc_AscpC( long int AscpCSum )
{
    if(pm.AscpC) delete[] pm.AscpC;
    pm.AscpC = new double[AscpCSum];
}

void TSolModMulti::alloc_UMpcC( long int UMpcSum )
{
    if(pm.UMpcC) delete[] pm.UMpcC;
    pm.UMpcC = new double[UMpcSum];
}

void TSolModMulti::alloc_xICuC( long int xICuCSum )
{
    if(pm.xICuC) delete[] pm.xICuC;
    pm.xICuC = new long int[xICuCSum];

}

//--------------------- End of tsolmod_multi_file.cpp ---------------------------


