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

#include "ms_multi.h"

TMultiBase::TMultiBase(TNode *na_)
{
    pa_standalone.reset( new BASE_PARAM() );
    *pa_standalone = pa_p_;

    pmp = &pm;
    node = na_; // parent
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

void TMultiBase::get_PAalp_PSigm( char& PAalp, char& PSigm)
{
   PAalp = PAalp_;
   PSigm = PSigm_;
}

void TMultiBase::STEP_POINT( const char* str)
{
   ipm_logger->info( std::string("STEP_POINT ")+str);
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

