//-------------------------------------------------------------------
// $Id$
//
/// \file tsolmod_multi_file.cpp
/// Implementation of writing IPM trace file and allocation model arrays
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
#include "io_template.h"
#include "io_keyvalue.h"
#include "gdatastream.h"
#include "jsonconfig.h"

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

void TSolModMulti::get_PAalp_PSigm(char& PAalp, char& PSigm)
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

//--------------------- end of tsolmod_multi_file.cpp ---------------------------

