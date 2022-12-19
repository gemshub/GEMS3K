//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi_diff.cpp
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

#include <fstream>
#include "node.h"
#include "num_methods.h"
#include "v_service.h"


const BASE_PARAM pa_p_ = 
        {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     30000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
}; // BASE_PARAM


TMultiBase::TMultiBase(TNode *na_)
{
    pa_standalone.reset( new BASE_PARAM() );
    *pa_standalone = pa_p_;

    pmp = &pm;
    node1 = na_; // parent
    sizeN = 0;
    AA = nullptr;
    BB = nullptr;
    arrL = nullptr;
    arrAN = nullptr;

    U_mean = nullptr;
    U_M2 = nullptr;
    U_CVo = nullptr;
    U_CV = nullptr;
    ICNud = nullptr;
    pm.errorCode[0] ='\0';
    pm.errorBuf[0] ='\0';

    sizeFIs = 0;
    phSolMod = nullptr;
    sizeFIa = 0;
    phSorpMod = nullptr;
    sizeFI = 0;
    phKinMet = nullptr;

    pmp->Guns = nullptr;
    pmp->Vuns = nullptr;
    pmp->tpp_G = nullptr;
    pmp->tpp_S = nullptr;
    pmp->tpp_Vm = nullptr;
    set_def();
}

bool TMultiBase::testTSyst() const
{
    return true;
}

void TMultiBase::get_PAalp_PSigm( char& PAalp, char& PSigm)
{
    PAalp = PAalp_;
    PSigm = PSigm_;
}

void TMultiBase::STEP_POINT( const char* str)
{
    ipm_logger->debug( std::string("STEP_POINT ")+str);
}

void TMultiBase::alloc_IPx( long int LsIPxSum )
{
    if( pm.IPx ) delete[] pm.IPx;
    pm.IPx = new long int[ LsIPxSum];
}

void TMultiBase::alloc_PMc( long int LsModSum )
{
    if( pm.PMc ) delete[] pm.PMc;
    pm.PMc = new double[LsModSum];
}

void TMultiBase::alloc_DMc( long int LsMdcSum )
{
    if( pm.DMc ) delete[] pm.DMc;
    pm.DMc = new double[LsMdcSum];
}

void TMultiBase::alloc_MoiSN( long int LsMsnSum )
{
    if(pm.MoiSN) delete[] pm.MoiSN;
    pm.MoiSN = new double[LsMsnSum];
}

void TMultiBase::alloc_SitFr( long int LsSitSum )
{
    if(pm.SitFr) delete[] pm.SitFr;
    pm.SitFr = new double[LsSitSum];
}

void TMultiBase::alloc_DQFc( long int DQFcSum )
{
    if(pm.DQFc) delete[] pm.DQFc;
    pm.DQFc = new double[DQFcSum];
}

void TMultiBase::alloc_PhLin( long int PhLinSum )
{
    if(pm.PhLin) delete[] pm.PhLin;
    pm.PhLin = new long int[PhLinSum][2];
}

void TMultiBase::alloc_lPhc( long int lPhcSum )
{
    if(pm.lPhc) delete[] pm.lPhc;
    pm.lPhc = new double[lPhcSum];
}

void TMultiBase::alloc_xSMd( long int xSMdSum )
{
    if(pm.xSMd) delete[] pm.xSMd;
    pm.xSMd = new long int[xSMdSum];
}

void TMultiBase::alloc_IsoPc( long int IsoPcSum )
{
    if(pm.IsoPc) delete[] pm.IsoPc;
    pm.IsoPc = new double[IsoPcSum];
}

void TMultiBase::alloc_IsoSc( long int IsoScSum )
{
    if(pm.IsoSc) delete[] pm.IsoSc;
    pm.IsoSc = new double[IsoScSum];
}

void TMultiBase::alloc_IsoCt( long int IsoCtSum )
{
    if(pm.IsoCt) delete[] pm.IsoCt;
    pm.IsoCt = new char[IsoCtSum];
}

void TMultiBase::alloc_EImc( long int EImcSum )
{
    if(pm.EImc) delete[] pm.EImc;
    pm.EImc = new double[EImcSum];
}

void TMultiBase::alloc_mCDc( long int mCDcSum )
{
    if(pm.mCDc) delete[] pm.mCDc;
    pm.mCDc = new double[mCDcSum];
}

void TMultiBase::alloc_xSKrC( long int xSKrCSum )
{
    if(pm.xSKrC) delete[] pm.xSKrC;
    pm.xSKrC = new long int[xSKrCSum];
}

void TMultiBase::alloc_ocPRkC( long int ocPRkC_feSArC_Sum )
{
    if(pm.ocPRkC) delete[] pm.ocPRkC;
    pm.ocPRkC = new long int[ocPRkC_feSArC_Sum][2];
}

void TMultiBase::alloc_feSArC( long int ocPRkC_feSArC_Sum )
{
    if(pm.feSArC) delete[] pm.feSArC;
    pm.feSArC = new double[ocPRkC_feSArC_Sum];
}

void TMultiBase::alloc_rpConC( long int rpConCSum )
{
    if(pm.rpConC) delete[] pm.rpConC;
    pm.rpConC = new double[rpConCSum];
}

void TMultiBase::alloc_apConC( long int apConCSum )
{
    if(pm.apConC) delete[] pm.apConC;
    pm.apConC = new double[apConCSum];
}

void TMultiBase::alloc_AscpC( long int AscpCSum )
{
    if(pm.AscpC) delete[] pm.AscpC;
    pm.AscpC = new double[AscpCSum];
}

void TMultiBase::alloc_UMpcC( long int UMpcSum )
{
    if(pm.UMpcC) delete[] pm.UMpcC;
    pm.UMpcC = new double[UMpcSum];
}

void TMultiBase::alloc_xICuC( long int xICuCSum )
{
    if(pm.xICuC) delete[] pm.xICuC;
    pm.xICuC = new long int[xICuCSum];

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/// Output to "ipmlog.txt" file Warnings
long int TMultiBase::testMulti()
{
    if( pm.MK || pm.PZ )
    {
        if( base_param()->PSM >= 2 )
        {
            TNode::ipmlog_file->warn(" {} : {}:{}", char_array_to_string(pm.stkey, EQ_RKLEN), pm.errorCode, pm.errorBuf);
        }
        return 1L;
    }
    return 0L	;
}

bool TMultiBase::calculateActivityCoefficients_scripts( long int , long int , long int ,
                                                    long int , long int , long int , double )
{
    return true;
}

// Before Calculations
// Calculation by IPM (preparing for calculation, unpacking data) in GUI part
void TMultiBase::initalizeGEM_IPM_Data_GUI()
{
}

void TMultiBase::multiConstInit_PN()
{
    pm.PZ = base_param()->DW;  // in IPM
    //  pm.FitVar[0] = 0.0640000030398369;
}

void TMultiBase::GEM_IPM_Init_gui1()
{
}

void TMultiBase::GEM_IPM_Init_gui2()
{
}

/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
//
void TMultiBase::DC_LoadThermodynamicData(TNode* aNa ) // formerly CompG0Load()
{
    double TK, PPa;

    TNode* na = node1;
    if( aNa != nullptr )
        na = aNa;   // for reading GEMIPM files task
    ErrorIf( na == nullptr, "DCLoadThermodynamicData", "Could not be undefined node" );
    TK =  na->cTK();
    PPa = na->cP();

    // try generate thermodynamic data from ThermoEngine
    if( !na->load_all_thermodynamic_from_thermo( TK, PPa ))
    {
        load_all_thermodynamic_from_grid(na, TK, PPa );
    }
    pm.pTPD = 2;
}


/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void TMultiBase::load_all_thermodynamic_from_grid(TNode* aNa, double TK, double PPa )
{
    long int j, jj, k, xTP, jb, je=0;
    double Go, Gg=0., Ge=0., Vv, h0=0., S0 = 0., Cp0= 0., a0 = 0., u0 = 0.;
    double P = PPa/bar_to_Pa;
    DATACH  *dCH = aNa->pCSD();

    ipm_logger->info("Calc Lookup T: {}  P: {}", TK, PPa);
    if( dCH->nTp <1 || dCH->nPp <1 || aNa->check_TP( TK, PPa ) == false )
    {
        Error("load_all_thermodynamic_from_grid: ",
               std::string(" Temperature ")+std::to_string(TK)+" or pressure "+
               std::to_string(PPa)+" out of range, or no T/D data are provided" );
        return;
    }

    pm.T = pm.Tc = TK;
    pm.TC = pm.TCc = TK-C_to_K;
    pm.Pc = P;
    if( P < 1e-5 )
    { // Pressure at saturated H2O vapour at given temperature
        long int xT = aNa->check_grid_T(TK);
        if(xT>= 0)
            P = dCH->Psat[xT]/bar_to_Pa;
        else
            P =  LagranInterp( &PPa, dCH->TKval, dCH->Psat, PPa, TK, dCH->nTp, 1,6 )/bar_to_Pa;
    }
    pm.P = P;
    pm.RT = R_CONSTANT * pm.Tc;
    pm.FRT = F_CONSTANT/pm.RT;
    pm.lnP = log( P );

    xTP = aNa->check_grid_TP( TK, PPa );

    for( k=0; k<5; k++ )
    {
        jj =  k * aNa->gridTP();
        if( xTP >= 0 )
        {
            pm.denW[k] = dCH->denW[jj+xTP]/1e3;
            pm.epsW[k] = dCH->epsW[jj+xTP];
            pm.denWg[k] = dCH->denWg[jj+xTP]/1e3;
            pm.epsWg[k] = dCH->epsWg[jj+xTP];
        }
        else
        {
            pm.denW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,6 )/1e3;// from test denW enough
            pm.epsW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,5 );// from test epsW enough
            pm.denWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 )/1e3;
            pm.epsWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 );
        }
    }

#ifdef  USE_THERMO_LOG
    std::fstream f_log("thermodynamic-log-lookup.csv", std::ios::out/*|std::ios::app*/ );
    f_log << "\nCalc ThermoEngine;T;" << TK << ";P;" << PPa << "\n";
    f_log << "denW";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.denW[jj]);
    f_log << "\nepsW";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.epsW[jj]);
    f_log << "\ndenWg";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.denWg[jj]);
    f_log << "\nepsWg";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.epsWg[jj]);
#endif
    long int xVol =  getXvolume();

    for( k=0; k<pm.FI; k++ )
    {
        jb = je;
        je += pm.L1[k];
        // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
        // depending on the presence of these arrays in DATACH and Multi structures
        for( j=jb; j<je; j++ )
        {
            jj =  j * aNa->gridTP();
            if( xTP >= 0 )
            {
                Go = dCH->G0[ jj+xTP];
                Vv = dCH->V0[ jj+xTP]*1e5;
                if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
                if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
                if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
                if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
                if( dCH->U0 ) h0 = dCH->U0[ jj+xTP];
            }
            else
            {
                Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
                Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 5 )*1e5;
                if( dCH->S0 ) S0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->S0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp, 4 ); // from test S0[Ca+2] enough
                if( dCH->H0 ) h0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->H0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->Cp0 ) Cp0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->Cp0+jj,
                                                    PPa, TK, dCH->nTp, dCH->nPp, 3 ); // from test Cp0[Ca+2] not more
                if( dCH->A0 ) a0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->A0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->U0 ) u0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->U0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
            }
            if( pm.tpp_G )
                pm.tpp_G[j] = Go;
            if( pm.Guns )
                Gg = pm.Guns[j];
            else
                Gg = 0.;

            Ge = 0.;
            pm.G0[j] = ConvertGj_toUniformStandardState( Go+Gg+Ge, j, k ); // formerly Cj_init_calc()
            // Inside this function, pm.YOF[k] can be added!

                switch( pm.PV )
                { // put molar volumes of components into A matrix or into the vector of molar volumes
                // to be checked!
                case VOL_CONSTR:
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
                    if( xVol >= 0. )
                        pm.A[j*pm.N+xVol] = Vv;
                    [[fallthrough]];
                case VOL_CALC:
                case VOL_UNDEF:
                    if( pm.tpp_Vm )
                        pm.tpp_Vm[j] = Vv;
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
                    pm.Vol[j] = Vv  * 10.;
                    break;
                }
            if( pm.S0 ) pm.S0[j] = S0;
            if( pm.H0 ) pm.H0[j] = h0;
            if( pm.Cp0 ) pm.Cp0[j] = Cp0;
            if( pm.A0 ) pm.A0[j] = a0;
            if( pm.U0 ) pm.U0[j] = u0;

#ifdef  USE_THERMO_LOG
            f_log << "\n" << std::string(dCH->DCNL[j], 0, MaxDCN) << ";" << floating_point_to_string(Go)
                   << ";" << floating_point_to_string(pm.G0[j])
                   << ";" << floating_point_to_string(pm.Vol[j]);
            if( dCH->S0 ) f_log << ";" << floating_point_to_string(pm.S0[j]);
            if( dCH->H0 ) f_log << ";" << floating_point_to_string(pm.H0[j]);
            if( dCH->Cp0 ) f_log << ";" << floating_point_to_string(pm.Cp0[j]);
            if( dCH->A0 ) f_log << ";" << floating_point_to_string(pm.A0[j]);
            if( dCH->U0 ) f_log << ";" << floating_point_to_string(pm.U0[j]);
#endif
        }  // j
    } // k
}


//--------------------- End of ms_multi_diff.cpp ---------------------------
