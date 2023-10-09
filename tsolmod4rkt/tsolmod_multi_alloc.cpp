//-------------------------------------------------------------------
/// \file ms_multi_alloc.cpp
/// Allocation of coping IPM internal structure
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

#include "tsolmod_multi.h"
#include "v_service.h"
#include "verror.h"

/// Realloc dynamic memory
void TSolModMulti::multi_realloc( char PAalp, char PSigm )
{
    long int ii, jj ;
    if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( "Multi", "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

    // Part 1
    // need  always to alloc vectors
    pm.L1 = new long int[pm.FI];
    for( ii=0; ii<pm.FI; ii++)
    {   pm.L1[ii] = 0;
    }
    pm.DUL = new double[pm.L];
    pm.DLL = new double[pm.L];
    pm.Vol = new double[pm.L];
    pm.Pparc = new double[pm.L];
    pm.MM = new double[pm.L];
    pm.G = new double[pm.L];
    pm.G0 = new double[pm.L];
    pm.lnGam = new double[pm.L];
    pm.lnGmo = new double[pm.L];
    pm.X = new double[pm.L];
    pm.Y = new double[pm.L];
    pm.XY = new double[pm.L];
    pm.MU = new double[pm.L];
    pm.EMU = new double[pm.L];
    pm.NMU = new double[pm.L];
    pm.W = new double[pm.L];
    pm.F = new double[pm.L];
    pm.F0 = new double[pm.L];
    pm.RLC = new char[pm.L];
    pm.RSC = new char[pm.L];
    pm.DCC = new char[pm.L];
    pm.DCCW = new char[pm.L];
    pm.lnGmM = new double[pm.L];
    pm.fDQF = new double[pm.L]; //24
    for( ii=0; ii<pm.L; ii++ )
    {
        pm.DUL[ii] = 1e6;
        pm.DLL[ii] = 0.0;
        pm.Vol[ii] = 0.0;
        pm.Pparc[ii] = 1.;
        pm.MM[ii] = 0.0;
        pm.G[ii] = 0.0;
        pm.G0[ii] = 0.0;
        pm.lnGam[ii] = 0.0;
        pm.lnGmo[ii] = 0.0;
        pm.X[ii] = 0.0;
        pm.Y[ii] = 0.0;
        pm.XY[ii] = 0.0;
        pm.MU[ii] = 0.0;
        pm.EMU[ii] = 0.0;
        pm.NMU[ii] = 0.0;
        pm.W[ii] = 0.0;
        pm.F[ii] = 0.0;
        pm.F0[ii] = 0.0;
        pm.RLC[ii] = 'B';
        pm.RSC[ii] = 'M';
        pm.DCC[ii] = 0;
        pm.DCCW[ii] = 0;
        pm.lnGmM[ii] = 0.0;
        pm.fDQF[ii] = 0.0;
    }

    pm.A = new double[pm.N*pm.L];
    for( ii=0; ii<pm.N*pm.L; ii++ )
        pm.A[ii] = 0.0;

    pm.Awt = new double[pm.N];
    pm.B = new double[pm.N];
    pm.U = new double[pm.N];
    pm.U_r = new double[pm.N];
    pm.C = new double[pm.N];
    pm.ICC = new char[pm.N];  //6
    for( ii=0; ii<pm.N; ii++ )
    {
        pm.Awt[ii] = 0.0;
        pm.B[ii] = 0.0;
        pm.U[ii] = 0.0;
        pm.U_r[ii] = 0.0;
        pm.C[ii] = 0.0;
        pm.ICC[ii] = 0;
    }

    pm.XFs = new double[pm.FI];
    pm.Falps = new double[pm.FI];
    pm.XF = new double[pm.FI];
    pm.YF = new double[pm.FI];
    pm.Falp = new double[pm.FI];
    pm.YOF = new double[pm.FI];
    pm.PHC = new char[pm.FI];
    pm.FVOL = new double[pm.FI];
    pm.FWGT = new double[pm.FI]; //9
    for( ii=0; ii<pm.FI; ii++ )
    {
        pm.XFs[ii] = 0.0;
        pm.Falps[ii] = 0.0;
        pm.XF[ii] = 0.0;
        pm.YF[ii] = 0.0;
        pm.Falp[ii] = 0.0;
        pm.YOF[ii] = 0.0;
        pm.PHC[ii] = 0;
        pm.FVOL[ii] = 0.0;
        pm.FWGT[ii] = 0.0;
    }

    pm.SB = new char[pm.N][MAXICNAME+MAXSYMB];
    for( ii=0; ii<pm.N; ii++)
    {
        fillValue( pm.SB[ii], '\0', MAXICNAME+MAXSYMB);
    }
    pm.SF = new char[pm.FI][MAXPHNAME+MAXSYMB];
    for( ii=0; ii<pm.FI; ii++)
    {
        fillValue( pm.SF[ii], '\0', MAXPHNAME+MAXSYMB);
    }
    pm.SM = new char[pm.L][MAXDCNAME];
    for( ii=0; ii<pm.L; ii++)
        fillValue( pm.SM[ii], '\0', MAXDCNAME);
    pm.dcMod = new char[pm.L][6];

    if( pm.L > 0 )
    {
        pm.Y_la = new double[pm.L];
        pm.Y_w = new double[pm.L];
        pm.Fx = new double[pm.L];
        pm.Wx = new double[pm.L];
        pm.VL = new double[pm.L];
        pm.Gamma = new double[pm.L];
        pm.lnGmf = new double[pm.L]; //7
        pm.GamFs = new double[pm.L];
        for( ii=0; ii<pm.L; ii++ )
        {
            pm.Y_la[ii] = 0.0;
            pm.Y_w[ii] = 0.0;
            pm.Fx[ii] = 0.0;
            pm.Wx[ii] = 0.0;
            pm.VL[ii] = 0.0;
            pm.Gamma[ii] = 0.0;
            pm.lnGmf[ii] = 0.0;
            pm.GamFs[ii] = 0.0;
        }
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
        pm.GamFs = 0;
    }

    // Part 2  not always required arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        pm.BF = new double[pm.FIs*pm.N];
        for( ii=0; ii<pm.FIs*pm.N; ii++ )
            pm.BF[ii] = 0.0;
        pm.BFC = new double[pm.N];
        for( ii=0; ii<pm.N; ii++ )
            pm.BFC[ii] = 0.0;

        pm.XFA = new double[pm.FIs];
        pm.YFA = new double[pm.FIs];
        for( ii=0; ii<pm.FIs; ii++ )
        {
            pm.XFA[ii] = 0.0;
            pm.YFA[ii] = 0.0;
        }
        pm.RFLC = new char[pm.FIs];
        pm.RFSC = new char[pm.FIs];
        for( ii=0; ii<pm.FIs; ii++)
        {
            pm.RFLC[ii] = 0;
            pm.RFSC[ii] = 0;
        }
    }
    else
    {
        pm.BF = 0;
        pm.BFC = 0;
        pm.XFA = 0;
        pm.YFA = 0;
        pm.RFLC = 0;
        pm.RFSC = 0;
    }

    if( pm.LO > 1 )
    {
        pm.Y_m = new double[pm.L];
        for( ii=0; ii<pm.L; ii++ )
            pm.Y_m[ii] = 0.0;
        pm.IC_m = new double[pm.N];
        pm.IC_lm = new double[pm.N];
        pm.IC_wm = new double[pm.N];
        for( ii=0; ii<pm.N; ii++ )
        {
            pm.IC_m[ii] = 0.0;
            pm.IC_lm[ii] = 0.0;
            pm.IC_wm[ii] = 0.0;
        }
    }
    else
    {
        pm.Y_m = 0;
        pm.IC_m = 0;
        pm.IC_lm = 0;
        pm.IC_wm = 0;
    }

    if( pm.E )
    {
        pm.EZ = new double[pm.L];
        for( ii=0; ii<pm.L; ii++ )
            pm.EZ[ii] = 0.0;
    }
    else
    {
        pm.EZ = 0;
    }

    // Part 3
    if( pm.Ls > 1 && pm.FIs > 0 )
    {
        pm.Wb = new double[pm.Ls];
        pm.Wabs = new double[pm.Ls];
        pm.Rion = new double[pm.Ls];
        for( ii=0; ii<pm.Ls; ii++ )
        {
            pm.Wb[ii] = 0.;
            pm.Wabs[ii] = 0.;
            pm.Rion[ii] = 0.;
        }
        pm.Qp = new double[pm.FIs*QPSIZE];
        pm.Qd = new double[pm.FIs*QDSIZE];
        for( ii=0; ii<pm.FIs*QPSIZE; ii++ )
            pm.Qp[ii] = 0.;
        for( ii=0; ii<pm.FIs*QDSIZE; ii++ )
            pm.Qd[ii] = 0.;
    }
    else
    {
        pm.Wb = 0;
        pm.Wabs = 0;
        pm.Rion = 0;
        pm.Qp = 0;
        pm.Qd = 0;

    }

    // added SD 03/02/2009
    pm.XU = new double[pm.L];
    for( ii=0; ii<pm.L; ii++ )
        pm.XU[ii] = 0.;
    pm.Uc = new double[pm.N][2];
    pm.Uefd = new double[pm.N];
    for( ii=0; ii<pm.N; ii++ )
    {
        pm.Uc[ii][0] = 0.;
        pm.Uc[ii][1] = 0.;
        pm.Uefd[ii] = 0.;
    }

    pm.Cp0   = new double[pm.L];
    pm.H0    = new double[pm.L];
    pm.U0    = new double[pm.L];
    pm.S0    = new double[pm.L];
    pm.A0    = new double[pm.L];
    for( ii=0; ii<pm.L; ii++ )
    {
        pm.Cp0[ii]   = 0.;
        pm.H0[ii]    = 0.;
        pm.U0[ii]    = 0.;
        pm.S0[ii]    = 0.;
        pm.A0[ii]    = 0.;

    }
    pm.VPh   = new double[pm.FIs][MIXPHPROPS];
    pm.GPh   = new double[pm.FIs][MIXPHPROPS];
    pm.HPh   = new double[pm.FIs][MIXPHPROPS];
    pm.SPh   = new double[pm.FIs][MIXPHPROPS];
    pm.CPh   = new double[pm.FIs][MIXPHPROPS];
    pm.APh   = new double[pm.FIs][MIXPHPROPS];
    pm.UPh   = new double[pm.FIs][MIXPHPROPS];
    for( ii=0; ii<pm.FIs; ii++ )
        for( jj=0; jj<MIXPHPROPS; jj++ )
        {
            pm.VPh[ii][jj]  = 0.;
            pm.GPh[ii][jj]  = 0.;
            pm.HPh[ii][jj]  = 0.;
            pm.SPh[ii][jj]  = 0.;
            pm.CPh[ii][jj]  = 0.;
            pm.APh[ii][jj]  = 0.;
            pm.UPh[ii][jj]  = 0.;
        }

    // NEW phase definition

    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        pm.IPx = 0;
        pm.PMc = 0;
        pm.DMc = 0;
        pm.MoiSN = 0;
        pm.SitFr = 0;
        pm.sMod = new char[pm.FIs][8];
        for( ii=0; ii<pm.FIs; ii++)
        {
            fillValue( pm.sMod[ii], ' ', 8);
        }
        pm.LsMod = new long int[pm.FIs*3];
        pm.LsMdc = new long int[pm.FIs*3];
        pm.LsMdc2 = new long int[pm.FIs*3];
        for( ii=0; ii<pm.FIs*3; ii++ )
        {
            pm.LsMod[ii] =0;
            pm.LsMdc[ii] = 0;
            pm.LsMdc2[ii] = 0;
        }
        pm.PhLin = 0;
        pm.lPhc  = 0;
        pm.LsPhl = new long int[pm.FI*2];
        for( ii=0; ii<pm.FI*2; ii++ )
            pm.LsPhl[ii] =0;
        // TSolMod stuff
        pm.lPhc   = 0;
        pm.DQFc   = 0;
        //    pm.rcpc   = 0;
        pm.lnDQFt   = new double[pm.Ls];
        pm.lnRcpt   = new double[pm.Ls];
        pm.lnExet   = new double[pm.Ls];
        pm.lnCnft   = new double[pm.Ls];
        pm.CTerms   = new double[pm.Ls];
        for( ii=0; ii<pm.Ls; ii++ )
        {    pm.lnDQFt[ii] =0.;
            pm.lnRcpt[ii] =0.;
            pm.lnExet[ii] =0.;
            pm.lnCnft[ii] =0.;
            pm.CTerms[ii] =0.;
        }
    }
    else
    {
        pm.LsMod = 0;
        pm.LsMdc = 0;
        pm.PMc = 0;
        pm.DMc = 0;
        pm.MoiSN = 0;
        pm.SitFr = 0;
        pm.sMod = 0;

        pm.LsMdc2  = 0;
        pm.LsPhl   = 0;
        pm.PhLin   = 0;
        // TSolMod stuff
        pm.lPhc   = 0;
        pm.DQFc   = 0;
        //    pm.rcpc   = 0;
        pm.lnDQFt   = 0;
        pm.lnRcpt   = 0;
        pm.lnExet   = 0;
        pm.lnCnft   = 0;
        pm.CTerms   = 0;
    }
}

/// Free of dynamic memory
void TSolModMulti::multi_kill()
{
    // Part 1
    // need  always to alloc vectors
    if( pm.L1) delete[] pm.L1;
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
    if(   pm.SM ) delete[] pm.SM;
    if(   pm.SF ) delete[] pm.SF;
    if(   pm.dcMod ) delete[] pm.dcMod;
    if(   pm.RLC ) delete[] pm.RLC;
    if(   pm.RSC ) delete[] pm.RSC;
    if(   pm.ICC ) delete[] pm.ICC;
    if(   pm.DCC ) delete[] pm.DCC;
    if(   pm.PHC ) delete[] pm.PHC;
    if(   pm.DCCW ) delete[] pm.DCCW;
    if( pm.lnGmM ) delete[] pm.lnGmM;
    if( pm.fDQF ) delete[] pm.fDQF;
    if( pm.FVOL ) delete[] pm.FVOL;
    if( pm.FWGT ) delete[] pm.FWGT;

    if( pm.Y_la ) delete[] pm.Y_la;
    if( pm.Y_w ) delete[] pm.Y_w;
    if( pm.Fx ) delete[] pm.Fx;
    if( pm.Wx ) delete[] pm.Wx;
    if( pm.VL ) delete[] pm.VL;
    if( pm.Gamma ) delete[] pm.Gamma;
    if( pm.lnGmf ) delete[] pm.lnGmf;
    if( pm.GamFs ) delete[] pm.GamFs;

    // Part 2  not requited arrays
    if( pm.BF ) delete[] pm.BF;
    if( pm.BFC ) delete[] pm.BFC;
    if( pm.XFA ) delete[] pm.XFA;
    if( pm.YFA ) delete[] pm.YFA;
    if( pm.RFLC ) delete[] pm.RFLC;
    if( pm.RFSC ) delete[] pm.RFSC;

    if( pm.Y_m ) delete[] pm.Y_m;
    if( pm.IC_m ) delete[] pm.IC_m;
    if( pm.IC_lm ) delete[] pm.IC_lm;
    if( pm.IC_wm ) delete[] pm.IC_wm;
    if( pm.EZ ) delete[] pm.EZ;

    // Part 3
    if( pm.Wb ) delete[] pm.Wb;
    if( pm.Wabs ) delete[] pm.Wabs;
    if( pm.Rion ) delete[] pm.Rion;
    if( pm.Qp ) delete[] pm.Qp;
    if( pm.Qd ) delete[] pm.Qd;

    // added SD 03/02/2009
    if( pm.XU ) delete[] pm.XU;
    if( pm.Uc ) delete[] pm.Uc;
    if( pm.Uefd ) delete[] pm.Uefd;

    if(pm.H0)  	delete[] pm.H0;
    if(pm.A0)  	delete[] pm.A0;
    if(pm.U0)  	delete[] pm.U0;
    if(pm.S0)  	delete[] pm.S0;
    if(pm.Cp0) 	delete[] pm.Cp0;

    if(pm.VPh)  	delete[] pm.VPh;
    if(pm.GPh)  	delete[] pm.GPh;
    if(pm.HPh)  	delete[] pm.HPh;
    if(pm.SPh)  	delete[] pm.SPh;
    if(pm.CPh)  	delete[] pm.CPh;
    if(pm.APh)  	delete[] pm.APh;
    if(pm.UPh)  	delete[] pm.UPh;

    if( pm.LsMod ) delete[] pm.LsMod;
    if( pm.LsMdc ) delete[] pm.LsMdc;
    if( pm.IPx ) delete[] pm.IPx;
    if( pm.PMc ) delete[] pm.PMc;
    if( pm.DMc ) delete[] pm.DMc;
    if(pm.MoiSN) delete[] pm.MoiSN;
    if(pm.SitFr) delete[] pm.SitFr;
    if( pm.sMod ) delete[] pm.sMod;

    if(pm.DQFc) delete[] pm.DQFc;
    //    if(pm.rcpc) delete[] pm.rcpc;
    if(pm.LsMdc2) delete[] pm.LsMdc2;
    if(pm.PhLin) delete[] pm.PhLin;
    if(pm.lPhc) delete[] pm.lPhc;
    if(pm.LsPhl) delete[] pm.LsPhl;

    // TSolMod stuff
    if(pm.lnDQFt) delete[] pm.lnDQFt;
    if(pm.lnRcpt) delete[] pm.lnRcpt;
    if(pm.lnExet) delete[] pm.lnExet;
    if(pm.lnCnft) delete[] pm.lnCnft;
    if(pm.CTerms) delete[] pm.CTerms;
}

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
    pm.DUL   = nullptr;
    pm.DLL   = nullptr;
    pm.fDQF   = nullptr;
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
    pm.FVOL  = nullptr;
    pm.FWGT  = nullptr;
    pm.Awt   = nullptr;
    pm.A     = nullptr;
    pm.XFs   = nullptr;
    pm.Falps = nullptr;
    pm.GamFs = nullptr;
    pm.G     = nullptr;
    pm.G0    = nullptr;
    pm.lnGam = nullptr;
    pm.lnGmo = nullptr;
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
    pm.sMod  = nullptr;
    pm.dcMod  = nullptr;
    pm.SB    = nullptr;
    pm.SM    = nullptr;
    pm.SF    = nullptr;
    pm.RLC   = nullptr;
    pm.RSC   = nullptr;
    pm.RFLC  = nullptr;
    pm.RFSC  = nullptr;
    pm.ICC   = nullptr;
    pm.DCC   = nullptr;
    pm.PHC   = nullptr;
    pm.DCCW  = nullptr;
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
    pm.lnDQFt   = 0;
    pm.lnRcpt   = 0;
    pm.lnExet   = 0;
    pm.lnCnft   = 0;
    pm.CTerms   = 0;
}

//--------------------- end of tsolmod_multi_alloc.cpp ---------------------------
