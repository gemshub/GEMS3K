//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi_copy.cpp
/// Implementation of coping IPM internal structure
//
// Copyright (c) 2017-2020 S.Dmytriyeva, D.Kulik
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

/// Realloc dynamic memory
void TSolModMulti::multi_realloc( char PAalp, char PSigm )
{
    long int ii, jj ;
    if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( "Multi", "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

    // Part 1
    // need  always to alloc vectors
    pm.L1 = new long int[pm.FI];
    pm.muk = new long int[pm.FI];
    for( ii=0; ii<pm.FI; ii++)
    {   pm.L1[ii] = 0;
        pm.muk[ii] = ii;
    }
    pm.mui = new long int[pm.N];
    for( ii=0; ii<pm.N; ii++)
        pm.mui[ii] = ii;
    pm.muj = new long int[pm.L];
    for( ii=0; ii<pm.L; ii++)
        pm.muj[ii] = ii;

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
    pm.SB1 = new char[pm.N][MAXICNAME];
    for( ii=0; ii<pm.N; ii++)
    {
        fillValue( pm.SB[ii], '\0', MAXICNAME+MAXSYMB);
        fillValue( pm.SB1[ii], '\0', MAXICNAME);
    }
    pm.SF = new char[pm.FI][MAXPHNAME+MAXSYMB];
    pm.SFs = new char[pm.FI][MAXPHNAME+MAXSYMB];
    for( ii=0; ii<pm.FI; ii++)
    {
        fillValue( pm.SF[ii], '\0', MAXPHNAME+MAXSYMB);
        fillValue( pm.SFs[ii], '\0',MAXPHNAME+MAXSYMB);
    }
    pm.SM = new char[pm.L][MAXDCNAME];
    for( ii=0; ii<pm.L; ii++)
        fillValue( pm.SM[ii], '\0', MAXDCNAME);
    pm.SM2 = new char[pm.Ls][MAXDCNAME];
    for( ii=0; ii<pm.Ls; ii++)
        fillValue( pm.SM2[ii], '\0', MAXDCNAME);
    pm.SF2 = new char[pm.FIs][MAXPHNAME+MAXSYMB];
    for( ii=0; ii<pm.FIs; ii++)
        fillValue( pm.SF2[ii], '\0', MAXPHNAME+MAXSYMB);
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
        pm.GamFs = 0;
        //   pm.D = 0;
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
        pm.PUL = new double[pm.FIs];
        pm.PLL = new double[pm.FIs]; //5
        for( ii=0; ii<pm.FIs; ii++ )
        {
            pm.XFA[ii] = 0.0;
            pm.YFA[ii] = 0.0;
            pm.PUL[ii] = 1e6;
            pm.PLL[ii] = 0.0;
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
        pm.PUL = 0;
        pm.PLL = 0;
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

    // dispersion and sorption phases
    if( PAalp != S_OFF )
    {
        pm.Aalp = new double[pm.FI];
        for( ii=0; ii<pm.FI; ii++ )
            pm.Aalp[ii] = 0.0;
        pm.Xr0h0 = new double[pm.FI][2];
        for( ii=0; ii<pm.FI; ii++ )
            pm.Xr0h0[ii][0] =  pm.Xr0h0[ii][1] = 0.0;
    }
    else
    {
        pm.Aalp = 0;
        pm.Xr0h0 = 0;
    }

    if( PSigm != S_OFF )
    {   pm.Sigw = new double[pm.FI];
        pm.Sigg = new double[pm.FI];
        for( ii=0; ii<pm.FI; ii++ )
        {
            pm.Sigw[ii] = 0.0;
            pm.Sigg[ii] = 0.0;
        }
    }
    else
    {   pm.Sigw = 0;
        pm.Sigg = 0;
    }

    if( pm.E )
    {
        pm.EZ = new double[pm.L];
        for( ii=0; ii<pm.L; ii++ )
            pm.EZ[ii] = 0.0;
        pm.Xcond = new double[pm.FI];
        pm.Xeps = new double[pm.FI];
        for( ii=0; ii<pm.FI; ii++ )
        {
            pm.Xcond[ii] = 0.0;
            pm.Xeps[ii] = 0.0;
        }
    }
    else
    {
        pm.EZ = 0;
        pm.Xcond = 0;
        pm.Xeps = 0;
    }

    if( pm.FIat > 0 /*&& pm.Lads > 0*/ && pm.FIs > 0 )
    { // ADSORBTION AND ION IXCHANDG
        pm.SATX = new long int[pm.Lads][4];
        pm.MASDJ = new double[pm.Lads][DFCN];
        pm.lnSAC = new double[pm.Lads][4];
        for( ii=0; ii<pm.Lads; ii++ )
        {
            pm.SATX[ii][0] = pm.SATX[ii][1] = pm.SATX[ii][2] = pm.SATX[ii][3] = 0;
            pm.lnSAC[ii][0] = pm.lnSAC[ii][1] = pm.lnSAC[ii][2] = pm.lnSAC[ii][3] = 0.0;
            for( jj=0; jj<MST; jj++ )
                pm.MASDJ[ii][jj] = 0.0;
        }

        pm.SCM  = new char[pm.FIs][MST];
        pm.Nfsp = new double[pm.FIs][MST];
        pm.MASDT = new double[pm.FIs][MST];
        pm.XcapA = new double[pm.FIs][MST];
        pm.XcapB = new double[pm.FIs][MST];
        pm.XcapD = new double[pm.FIs][MST];
        pm.XcapF = new double[pm.FIs][MST];
        pm.XdlA = new double[pm.FIs][MST];
        pm.XdlB = new double[pm.FIs][MST];
        pm.XdlD = new double[pm.FIs][MST];
        pm.XpsiA = new double[pm.FIs][MST];
        pm.XpsiB = new double[pm.FIs][MST];
        pm.XpsiD = new double[pm.FIs][MST];
        pm.XlamA = new double[pm.FIs][MST];
        pm.Xetaf = new double[pm.FIs][MST];
        pm.XetaA = new double[pm.FIs][MST];
        pm.XetaB = new double[pm.FIs][MST];
        pm.XetaD = new double[pm.FIs][MST];
        pm.XFTS = new double[pm.FIs][MST];  //19
        for( ii=0; ii<pm.FIs; ii++ )
            for( jj=0; jj<MST; jj++ )
            {
                pm.SCM[ii][jj]  = 0;
                pm.Nfsp[ii][jj] = 0.0;
                pm.MASDT[ii][jj] = 0.0;
                pm.XcapA[ii][jj] = 0.0;
                pm.XcapB[ii][jj] = 0.0;
                pm.XcapD[ii][jj] = 0.0;
                pm.XcapF[ii][jj] = 0.0;
                pm.XdlA[ii][jj] = 0.0;
                pm.XdlB[ii][jj] = 0.0;
                pm.XdlD[ii][jj] = 0.0;
                pm.XpsiA[ii][jj] = 0.0;
                pm.XpsiB[ii][jj] = 0.0;
                pm.XpsiD[ii][jj] = 0.0;
                pm.XlamA[ii][jj] = 0.0;
                pm.Xetaf[ii][jj] = 0.0;
                pm.XetaA[ii][jj] = 0.0;
                pm.XetaB[ii][jj] = 0.0;
                pm.XetaD[ii][jj] = 0.0;
                pm.XFTS[ii][jj] = 0.0;
            }

        pm.SATT = new char[pm.Lads];
        pm.SM3 = new char[pm.Lads][MAXDCNAME];
        pm.DCC3 = new char[pm.Lads];
        for( ii=0; ii<pm.Lads; ii++)
        {
            fillValue( pm.SM3[ii], '\0', MAXDCNAME);
            pm.SATT[ii] = 0;
            pm.DCC3[ii] = 0;
        }

        pm.D = new double[MST][MST];
        for( ii=0; ii<MST; ii++ )
            for( jj=0; jj<MST; jj++ )
                pm.D[ii][jj] = 0.0;

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
        pm.Fug = new double[pm.PG];
        pm.Fug_l = new double[pm.PG];
        pm.Ppg_l = new double[pm.PG];
        for( ii=0; ii<pm.PG; ii++ )
        {
            pm.Fug[ii] = 0.;
            pm.Fug_l[ii] = 0.;
            pm.Ppg_l[ii] = 0.;
        }
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
        for( ii=0; ii<pm.Ls; ii++ )
        {    pm.lnDQFt[ii] =0.;
            pm.lnRcpt[ii] =0.;
            pm.lnExet[ii] =0.;
            pm.lnCnft[ii] =0.;
        }
        //TSorpMod & TKinMet stuff
        pm.SorMc   = new double[pm.FIs*16];
        for( ii=0; ii<pm.FIs*16; ii++ )
            pm.SorMc[ii] =0.;
        // TSorpMod stuff
        pm.LsESmo   = new long int[pm.FIs*4];
        pm.LsISmo   = new long int[pm.FIs*4];
        for( ii=0; ii<pm.FIs*4; ii++ )
        {    pm.LsESmo[ii] =0;
            pm.LsISmo[ii] =0;
        }
        pm.xSMd   = 0;
        pm.EImc   = 0;
        pm.mCDc   = 0;
        pm.IsoPc   = 0;
        pm.IsoSc   = 0;
        pm.lnScalT   = new double[pm.Ls];
        pm.lnSACT   = new double[pm.Ls];
        pm.lnGammF   = new double[pm.Ls];
        pm.CTerms   = new double[pm.Ls];
        pm.IsoCt   = 0;
        for( ii=0; ii<pm.Ls; ii++ )
        {    pm.lnScalT[ii] =0.;
            pm.lnSACT[ii] =0.;
            pm.lnGammF[ii] =0.;
            pm.CTerms[ii] =0.;
        }
        // TKinMet stuff
        pm.LsKin   = new long int[pm.FI*6];
        for( ii=0; ii<pm.FI*6; ii++ )
            pm.LsKin[ii] =0;
        pm.LsUpt   = new long int[pm.FIs*2];
        for( ii=0; ii<pm.FIs*2; ii++ )
            pm.LsUpt[ii] =0;
        pm.xSKrC   = 0;
        pm.ocPRkC   = 0;
        pm.feSArC   = 0;
        pm.rpConC   = 0;
        pm.apConC   = 0;
        pm.AscpC   = 0;
        pm.UMpcC   = 0;
        pm.kMod   = new char[pm.FI][6];
        pm.PfFact  = new double[pm.FI];
        pm.PrT   = new double[pm.FI];
        pm.PkT   = new double[pm.FI];
        pm.PvT   = new double[pm.FI];
        for( ii=0; ii<pm.FI; ii++)
        {
            fillValue( pm.kMod[ii], 'N', 6);
            pm.PfFact[ii] =0.;
            pm.PrT[ii] =0.;
            pm.PkT[ii] =0.;
            pm.PvT[ii] =0.;
        }
        pm.emRd   = new double[pm.Ls];
        pm.emDf   = new double[pm.Ls];
        for( ii=0; ii<pm.Ls; ii++)
        {
            pm.emRd[ii] =0.;
            pm.emDf[ii] =0.;
        }
        pm.xICuC = 0;

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

    Alloc_TSolMod( pm.FIs );
    //Alloc_TSorpMod( pm.FIs );
    //Alloc_TKinMet( pm.FI );

}

/// Free of dynamic memory
void TSolModMulti::multi_kill()
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
    //   if( pm.D ) delete[] pm.D;

    // Part 2  not requited arrays

    if( pm.BF ) delete[] pm.BF;
    if( pm.BFC ) delete[] pm.BFC;
    if( pm.XFA ) delete[] pm.XFA;
    if( pm.YFA ) delete[] pm.YFA;
    if( pm.PUL ) delete[] pm.PUL;
    if( pm.PLL ) delete[] pm.PLL;
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

    //  Added 16.11.2004 by Sveta
    //    if( pm.sitE )     delete[] pm.sitE;
    //    if( pm.sitXcat )  delete[] pm.sitXcat;
    //    if( pm.sitXan )    delete[] pm.sitXan;

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
    //TSorpMod & TKinMet stuff
    if(pm.SorMc) delete[] pm.SorMc;
    // TSorpMod stuff
    if(pm.LsESmo) delete[] pm.LsESmo;
    if(pm.LsISmo) delete[] pm.LsISmo;
    if(pm.xSMd) delete[] pm.xSMd;
    if(pm.EImc) delete[] pm.EImc;
    if(pm.mCDc) delete[] pm.mCDc;
    if(pm.IsoPc) delete[] pm.IsoPc;
    if(pm.IsoSc) delete[] pm.IsoSc;
    if(pm.lnScalT) delete[] pm.lnScalT;
    if(pm.lnSACT) delete[] pm.lnSACT;
    if(pm.lnGammF) delete[] pm.lnGammF;
    if(pm.CTerms) delete[] pm.CTerms;
    if(pm.IsoCt) delete[] pm.IsoCt;
    // TKinMet stuff
    if(pm.LsKin) delete[] pm.LsKin;
    if(pm.LsUpt) delete[] pm.LsUpt;
    if(pm.xSKrC) delete[] pm.xSKrC;
    if(pm.ocPRkC) delete[] pm.ocPRkC;
    if(pm.feSArC) delete[] pm.feSArC;
    if(pm.rpConC) delete[] pm.rpConC;
    if(pm.apConC) delete[] pm.apConC;
    if(pm.AscpC) delete[] pm.AscpC;
    if(pm.UMpcC) delete[] pm.UMpcC;
    if(pm.kMod) delete[] pm.kMod;
    if(pm.PfFact) delete[] pm.PfFact;
    if(pm.PrT) delete[] pm.PrT;
    if(pm.PkT) delete[] pm.PkT;
    if(pm.PvT) delete[] pm.PvT;
    if(pm.emRd) delete[] pm.emRd;
    if(pm.emDf) delete[] pm.emDf;
    if(pm.xICuC) delete[] pm.xICuC;

    // optimization 08/02/2007
    Free_TSolMod();
    //Free_TSorpMod();
    //Free_TKinMet();
    Free_internal();
    Free_uDD();
}

/// Internal memory allocation for TSolMod performance optimization
/// (since version 2.3.0)
void TSolModMulti::Alloc_TSolMod( long int newFIs )
{
  if(  phSolMod && ( newFIs == sizeFIs) )
    return;

  Free_TSolMod();
  // alloc memory for all multicomponents phases
  phSolMod = new  TSolMod *[newFIs];
  sizeFIs = newFIs;
 for( long int ii=0; ii<newFIs; ii++ )
          phSolMod[ii] = nullptr;
}

void TSolModMulti::Free_TSolMod()
{
  long int kk;

  if( phSolMod )
  {
      for(  kk=0; kk<sizeFIs; kk++ )
        if( phSolMod[kk] )
           delete phSolMod[kk];
      delete[]  phSolMod;
  }
  phSolMod = nullptr;
  sizeFIs = 0;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Internal memory allocation for IPM performance optimization
/// (since version 2.2.0)
//
void TSolModMulti::Alloc_A_B( long int newN )
{
  if( AA && BB && (newN == sizeN) )
    return;
  Free_A_B();
  AA = new  double[newN*newN];
  BB = new  double[newN];
  sizeN = newN;
}

void TSolModMulti::Free_A_B()
{
  if(AA)
    { delete[] AA; AA = 0; }
  if(BB)
    { delete[] BB; BB = 0; }
  sizeN = 0;
}

#define  a(j,i) ((*(pm.A+(i)+(j)*pm.N)))
/// Building an index list of non-zero elements of the matrix pm.A
void TSolModMulti::Build_compressed_xAN()
{
 long int ii, jj, k;

 // Calculate number of non-zero elements in A matrix
 k = 0;
 for( jj=0; jj<pm.L; jj++ )
   for( ii=0; ii<pm.N; ii++ )
     if( fabs( a(jj,ii) ) > 1e-12 )
       k++;

   // Free old memory allocation
    Free_compressed_xAN();

   // Allocate memory
   arrL = new long int[pm.L+1];
   arrAN = new long int[k];

   // Set indexes in the index arrays
   k = 0;
   for( jj=0; jj<pm.L; jj++ )
   { arrL[jj] = k;
     for( ii=0; ii<pm.N; ii++ )
       if( fabs( a(jj,ii) ) > 1e-12 )
       {
        arrAN[k] = ii;
        k++;
       }
   }
   arrL[jj] = k;
}
#undef a

void TSolModMulti::Free_compressed_xAN()
{
  if( arrL  )
    { delete[] arrL; arrL = 0;  }
  if( arrAN )
    { delete[] arrAN; arrAN = 0;  }
}

void TSolModMulti::Free_internal()
{
  Free_compressed_xAN();
  Free_A_B();
 }

/// Internal memory allocation for IPM performance optimization
void TSolModMulti::Alloc_internal()
{
// optimization 08/02/2007
 Alloc_A_B( pm.N );
 Build_compressed_xAN();
}


/// Added for implementation of divergence detection in dual solution 06.05.2011 DK
void TSolModMulti::Alloc_uDD( long int newN )
{
    if( U_mean && U_M2 && U_CVo && U_CV && ICNud && (newN == nNu) )
      return;
    Free_uDD();
    U_mean = new  double[newN]; // w3 u mean values for r
    U_M2 = new  double[newN];   // w3 u mean values for r-1
    U_CVo = new  double[newN];  // w3 u mean difference for r-1
    U_CV = new  double[newN];   // w3 u mean difference for r
    ICNud = new long int[newN];
    nNu = newN;
}

void TSolModMulti::Free_uDD()
{
    if( U_mean  )
      { delete[] U_mean; U_mean = 0; }
    if( U_M2  )
      { delete[] U_M2; U_M2 = 0; }
    if( U_CVo  )
      { delete[] U_CVo; U_CVo = 0; }
    if( U_CV  )
      { delete[] U_CV; U_CV = 0; }
    if( ICNud )
      { delete[] ICNud; ICNud = 0; }
    nNu = 0;
}


//--------------------- End of tsolmod_multi_copy.cpp ---------------------------
