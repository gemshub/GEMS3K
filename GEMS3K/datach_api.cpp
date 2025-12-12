//--------------------------------------------------------------------
// $Id$
//
/// \file datach_api.cpp
/// DATACH and DATABR structures allocations
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

#include <sstream>
#include "datach_api.h"
#include "gems3k_impex.h"
#include "gdatastream.h"
#include "io_keyvalue.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "v_service.h"

namespace  dbr_dch_api {

long int gridTP(const DATACH* pCSD)
{
    if( pCSD->mLook == 1L )
        return pCSD->nTp;
    else
        return (pCSD->nPp * pCSD->nTp);
}

std::string check_TP(const DATACH* CSD, double TK, double P)
{
    std::string error_msg;
    bool okT = true, okP = true;
    double T_=TK, P_=P;

    if( CSD->mLook == 1 )
    {
        for(long int  jj=0; jj<CSD->nPp; jj++)
            if( (fabs( P - CSD->Pval[jj] ) < CSD->Ptol ) && ( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol ) )
            {
                return error_msg;
            }
        Error( "check_TP: ", std::string("Temperature ")+std::to_string(TK)+
               " and pressure "+std::to_string(P)+" out of range");
        //return false;
    }
    else
    {
        if( TK <= CSD->TKval[0] - CSD->Ttol )
        { 				// Lower boundary of T interpolation interval
            okT = false;
            T_ = CSD->TKval[0] - CSD->Ttol;
        }
        if( TK >= CSD->TKval[CSD->nTp-1] + CSD->Ttol )
        {
            okT = false;
            T_ = CSD->TKval[CSD->nTp-1] + CSD->Ttol;
        }
        if( okT == false ) {
            error_msg += "Given TK=" + std::to_string(TK);
            error_msg += " is beyond the interpolation range for thermodynamic data near boundary T_= ";
            error_msg += std::to_string(T_);
        }

        if( P <= CSD->Pval[0] - CSD->Ptol )
        {
            okP = false;
            P_ = CSD->Pval[0] - CSD->Ptol;
        }
        if( P >= CSD->Pval[CSD->nPp-1] + CSD->Ptol )
        {
            okP = false;
            P_ = CSD->Pval[CSD->nPp-1] + CSD->Ptol;
        }
        if( !okP ) {
            error_msg += "Given P=" + std::to_string(P);
            error_msg += " is beyond the interpolation range for thermodynamic data near boundary P_= ";
            error_msg += std::to_string(P_);
        }
        return error_msg;
    }
    return error_msg;
}

long int check_grid_T(const DATACH* CSD, double TK)
{
    long int jj;
    for( jj=0; jj<CSD->nTp; jj++)
        if( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol )
            return jj;
    return -1;
}

long int check_grid_P(const DATACH* CSD, double P)
{
    long int jj;
    for( jj=0; jj<CSD->nPp; jj++)
        if( fabs( P - CSD->Pval[jj] ) < CSD->Ptol )
            return jj;
    return -1;
}

long int check_grid_TP(const DATACH* CSD, double TK, double P)
{
    long int xT, xP, ndx=-1;

    if( CSD->mLook == 1 )
    {
        for(long int  jj=0; jj<CSD->nPp; jj++)
            if( (fabs( P - CSD->Pval[jj] ) < CSD->Ptol ) && ( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol ) )
                return jj;
        Error( "check_grid_TP: " , std::string("Temperature ")+std::to_string(TK)+
               " and pressure "+std::to_string(P)+" out of grid" );
        //return -1;
    }
    else
    {
        xT = check_grid_T(CSD, TK);
        xP = check_grid_P(CSD, P);
        if( xT >=0 && xP>= 0 )
            ndx =  xP * CSD->nTp + xT;
        return ndx;
    }
    return ndx;
}

//-------------------------------------------------------------------------

// Writes CSD (DATACH structure) to a json string or key-value string
// \param brief_mode - Do not write data items that contain only default values
// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string datach_to_string(const std::string& current_set_name, const DATACH* pCSD, bool with_comments, bool brief_mode)
{
    std::stringstream ss;
    write_dch_format_stream(current_set_name, pCSD, ss, GEMS3KGenerator::default_type_f, with_comments, brief_mode );
    return ss.str();
}

// Reads CSD (DATACH structure) from a json string or key-value string
bool datach_from_string(const std::string& current_set_name, DATACH* pCSD, const std::string data)
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_dch_format_stream(current_set_name, pCSD, ss, GEMS3KGenerator::default_type_f );
    return true;
}

// Writes work node (DATABR structure) to a json string or key-value string
// \param brief_mode - Do not write data items that contain only default values
// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string databr_to_string(const std::string& current_set_name, const DATACH* pCSD, const DATABR* pCNode,
                             bool with_comments, bool brief_mode)
{
    std::stringstream ss;
    write_dbr_format_stream(current_set_name, pCSD, pCNode, ss, GEMS3KGenerator::default_type_f, with_comments, brief_mode );
    return ss.str();
}

// Reads work node (DATABR structure) from a json string or key-value string
bool databr_from_string(const std::string& current_set_name, const DATACH* pCSD, DATABR* pCNode, const std::string data)
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_dbr_format_stream(current_set_name, pCSD, pCNode, ss, GEMS3KGenerator::default_type_f );
    return true;
}


void write_dch_format_stream(const std::string& current_output_set_name, const DATACH* pCSD,
                             std::iostream& stream, GEMS3KGenerator::IOModes type_f, bool with_comments, bool brief_mode)
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonWrite out_format(stream, current_output_set_name);
        datach_to_text_file(pCSD, out_format, type_f==GEMS3KGenerator::f_thermofun, with_comments, brief_mode);
    }
#else
    {
        io_formats::SimdJsonWrite out_format(stream, current_output_set_name, with_comments);
        datach_to_text_file(pCSD, out_format, type_f==GEMS3KGenerator::f_thermofun, with_comments, brief_mode);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueWrite out_format( stream );
        datach_to_text_file(pCSD, out_format, type_f==GEMS3KGenerator::f_kv_thermofun, with_comments, brief_mode);
    }
        break;
    }
}

void read_dch_format_stream(const std::string& current_input_set_name, DATACH* pCSD,
                            std::iostream& stream, GEMS3KGenerator::IOModes  type_f)
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonRead in_format( stream, current_input_set_name, "dch" );
        datach_from_text_file(pCSD, in_format, type_f==GEMS3KGenerator::f_thermofun);
    }
#else
    {
        io_formats::SimdJsonRead in_format( stream, current_input_set_name, "dch" );
        datach_from_text_file(pCSD, in_format, type_f==GEMS3KGenerator::f_thermofun);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueRead in_format( stream );
        datach_from_text_file(pCSD, in_format, type_f==GEMS3KGenerator::f_kv_thermofun);
    }
        break;
    }
}

void write_dbr_format_stream(const std::string& current_output_set_name, const DATACH* pCSD, const DATABR* pCNode,
                             std::iostream& stream, GEMS3KGenerator::IOModes type_f,  bool with_comments, bool brief_mode)
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonWrite out_format( stream, current_output_set_name );
        databr_to_text_file(pCSD, pCNode, out_format, with_comments, brief_mode);
    }
#else
    {
        io_formats::SimdJsonWrite out_format( stream, current_output_set_name, with_comments );
        databr_to_text_file(pCSD, pCNode, out_format, with_comments, brief_mode);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueWrite out_format( stream );
        databr_to_text_file(pCSD, pCNode, out_format, with_comments, brief_mode);
    }
        break;
    }
}

void read_dbr_format_stream(const std::string& current_input_set_name, const DATACH* pCSD, DATABR* pCNode,
                            std::iostream& stream, GEMS3KGenerator::IOModes  type_f)
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonRead in_format( stream, current_input_set_name, "dbr" );
        databr_from_text_file(pCSD, pCNode, in_format);
    }
#else
    {
        io_formats::SimdJsonRead in_format( stream, current_input_set_name, "dbr" );
        databr_from_text_file(pCSD, pCNode, in_format);
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueRead in_format( stream );
        databr_from_text_file(pCSD, pCNode, in_format);
    }
        break;
    }
}

//==============================================================================

// allocating DataCH structure
void datach_realloc(DATACH* CSD)
{
    if( CSD->mLook == 1 &&  (CSD->nPp != CSD->nTp) )
        Error( "No-interpolation mode",
               "Different number of points for temperature and pressure ");

    datach_free(CSD);
    CSD->nDCinPH = new long int[CSD->nPH];

    if( CSD->nICb >0 )
        CSD->xic = new long int[CSD->nICb];
    else  CSD->xic = 0;
    if( CSD->nDCb >0 )
        CSD->xdc = new long int[CSD->nDCb];
    else  CSD->xdc = 0;
    if( CSD->nPHb >0 )
        CSD->xph = new long int[CSD->nPHb];
    else  CSD->xph = 0;

    CSD->A = new double[CSD->nIC*CSD->nDC];
    CSD->ICmm = new double[CSD->nIC];
    CSD->DCmm = new double[CSD->nDC];
    for(int  j=0; j<CSD->nDC; j++ )  {
        CSD->DCmm[j] = 0.0;
    }

    CSD->TKval = new double[CSD->nTp];
    CSD->Psat = new double[CSD->nTp];
    for(int  j=0; j<CSD->nTp; j++ )  {
        CSD->Psat[j] = 0.;
    }
    CSD->Pval = new double[CSD->nPp];

    CSD->denW = new double[ 5*gridTP(CSD)];
    CSD->denWg = new double[ 5*gridTP(CSD)];
    CSD->epsW = new double[ 5*gridTP(CSD)];
    CSD->epsWg = new double[ 5*gridTP(CSD)];
    for(int  j=0; j<5*gridTP(CSD); j++ )  {
        CSD->denW[j] = 0.;
        CSD->denWg[j] = 0.;
        CSD->epsW[j] = 0.;
        CSD->epsWg[j] = 0.;
    }

    CSD->G0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->V0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->H0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->S0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->Cp0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->A0 = new double[CSD->nDC*gridTP(CSD)];
    CSD->U0 = new double[CSD->nDC*gridTP(CSD)];

    if(  CSD->iGrd  )
        CSD->DD = new double[CSD->nDCs*gridTP(CSD)];
    else
        CSD->DD = 0;
    CSD->ICNL = new char[CSD->nIC][MaxICN];
    CSD->DCNL = new char[CSD->nDC][MaxDCN];
    CSD->PHNL = new char[CSD->nPH][MaxPHN];

    CSD->ccIC = new char[CSD->nIC];
    CSD->ccDC = new char[CSD->nDC];
    CSD->ccPH = new char[CSD->nPH];
}

// free dynamic memory
void datach_free(DATACH* CSD)
{
    if( CSD->nDCinPH )
    { delete[] CSD->nDCinPH;
        CSD->nDCinPH = 0;
    }
    if( CSD->xic )
    { delete[] CSD->xic;
        CSD->xic = 0;
    }
    if( CSD->xdc )
    { delete[] CSD->xdc;
        CSD->xdc = 0;
    }
    if( CSD->xph )
    { delete[] CSD->xph;
        CSD->xph = 0;
    }
    if( CSD->A )
    { delete[] CSD->A;
        CSD->A = 0;
    }
    if( CSD->ICmm )
    { delete[] CSD->ICmm;
        CSD->ICmm = 0;
    }
    if( CSD->DCmm )
    { delete[] CSD->DCmm;
        CSD->DCmm = 0;
    }

    if( CSD->TKval )
    { delete[] CSD->TKval;
        CSD->TKval = 0;
    }
    if( CSD->Psat )
    { delete[] CSD->Psat;
        CSD->Psat = 0;
    }
    if( CSD->Pval )
    { delete[] CSD->Pval;
        CSD->Pval = 0;
    }

    if( CSD->denW )
    { delete[] CSD->denW;
        CSD->denW = 0;
    }
    if( CSD->denWg )
    { delete[] CSD->denWg;
        CSD->denWg = 0;
    }
    if( CSD->epsW )
    { delete[] CSD->epsW;
        CSD->epsW = 0;
    }
    if( CSD->epsWg )
    { delete[] CSD->epsWg;
        CSD->epsWg = 0;
    }
    if( CSD->G0 )
    { delete[] CSD->G0;
        CSD->G0 = 0;
    }
    if( CSD->V0 )
    { delete[] CSD->V0;
        CSD->V0 = 0;
    }
    if( CSD->H0 )
    { delete[] CSD->H0;
        CSD->H0 = 0;
    }
    if( CSD->Cp0 )
    { delete[] CSD->Cp0;
        CSD->Cp0 = 0;
    }
    if( CSD->S0 )
    { delete[] CSD->S0;
        CSD->S0 = 0;
    }
    if( CSD->A0 )
    { delete[] CSD->A0;
        CSD->A0 = 0;
    }
    if( CSD->U0 )
    { delete[] CSD->U0;
        CSD->U0 = 0;
    }
    if( CSD->DD )
    { delete[] CSD->DD;
        CSD->DD = 0;
    }

    if( CSD->ICNL )
    { delete[] CSD->ICNL;
        CSD->ICNL = 0;
    }
    if( CSD->DCNL )
    { delete[] CSD->DCNL;
        CSD->DCNL = 0;
    }
    if( CSD->PHNL )
    { delete[] CSD->PHNL;
        CSD->PHNL = 0;
    }

    if( CSD->ccIC )
    { delete[] CSD->ccIC;
        CSD->ccIC = 0;
    }
    if( CSD->ccDC )
    { delete[] CSD->ccDC;
        CSD->ccDC = 0;
    }
    if( CSD->ccPH )
    { delete[] CSD->ccPH;
        CSD->ccPH = 0;
    }
    // delete[] CSD;
}

// Allocates DataBR structure
void databr_realloc(const DATACH* CSD, DATABR* CNode_)
{
    long int j,k;
    CNode_->bIC = new double[CSD->nICb];
    CNode_->rMB = new double[CSD->nICb];
    CNode_->uIC = new double[CSD->nICb];
    CNode_->bSP = new double[CSD->nICb];

    for(  j=0; j<CSD->nICb; j++ )
    {
        CNode_->rMB[j] = 0.;
        CNode_->uIC[j] = 0.;
        CNode_->bSP[j] = 0.;
    }

    CNode_->xDC = new double[CSD->nDCb];
    CNode_->gam = new double[CSD->nDCb];

    for(  j=0; j<CSD->nDCb; j++ )
    {
        CNode_->xDC[j] = 0.;
        CNode_->gam[j] = 1.;
    }

    //  default assignment
    CNode_->dul = new double[CSD->nDCb];
    for(  j=0; j<CSD->nDCb; j++ )
        CNode_->dul[j] = 1.0e6;            // default assignment
    CNode_->dll = new double[CSD->nDCb];
    for(  j=0; j<CSD->nDCb; j++ )
        CNode_->dll[j] = 0.0;              // default assignment

    if( CSD->nAalp >0 )
    {
        CNode_->aPH = new double[CSD->nPHb];
        for(  k=0; k<CSD->nPHb; k++ )
            CNode_->aPH[k] = 0.0;       // default assignment
    }
    else
        CNode_->aPH = 0;

    CNode_->xPH = new double[CSD->nPHb];
    CNode_->omPH = new double[CSD->nPHb];

    for(  k=0; k<CSD->nPHb; k++ )
    {
        CNode_->xPH[k] = 0.0;       // default assignment
        CNode_->omPH[k] = 0.0;
    }

    CNode_->vPS = new double[CSD->nPSb];
    CNode_->mPS = new double[CSD->nPSb];
    CNode_->bPS = new double[CSD->nPSb*CSD->nICb];
    CNode_->xPA = new double[CSD->nPSb];
    CNode_->amru = new double[CSD->nPSb];
    CNode_->amrl = new double[CSD->nPSb];

    for(  k=0; k<CSD->nPSb; k++ )
    {
        CNode_->vPS[k] = 0.0;
        CNode_->mPS[k] = 0.0;
        CNode_->xPA[k] = 0.0;
        for(  j=0; j<CSD->nICb; j++ )
            CNode_->bPS[k*CSD->nICb+j] = 0.0;
        CNode_->amru[k] = 1.0e6;
        CNode_->amrl[k] = 0.0;
    }
}

// free dynamic memory
void databr_free_internal(DATABR *CNode_)
{
    if( CNode_ == 0)
        return;

    if( CNode_->bIC )
    { delete[] CNode_->bIC;
        CNode_->bIC = 0;
    }
    if( CNode_->rMB )
    { delete[] CNode_->rMB;
        CNode_->rMB = 0;
    }
    if( CNode_->uIC )
    { delete[] CNode_->uIC;
        CNode_->uIC = 0;
    }

    if( CNode_->xDC )
    { delete[] CNode_->xDC;
        CNode_->xDC = 0;
    }
    if( CNode_->gam )
    { delete[] CNode_->gam;
        CNode_->gam = 0;
    }
    if( CNode_->dul )
    { delete[] CNode_->dul;
        CNode_->dul = 0;
    }
    if( CNode_->dll )
    { delete[] CNode_->dll;
        CNode_->dll = 0;
    }

    if( CNode_->aPH )
    { delete[] CNode_->aPH;
        CNode_->aPH = 0;
    }
    if( CNode_->xPH )
    { delete[] CNode_->xPH;
        CNode_->xPH = 0;
    }
    if( CNode_->omPH )
    { delete[] CNode_->omPH;
        CNode_->omPH = 0;
    }
    if( CNode_->vPS )
    { delete[] CNode_->vPS;
        CNode_->vPS = 0;
    }
    if( CNode_->mPS )
    { delete[] CNode_->mPS;
        CNode_->mPS = 0;
    }
    if( CNode_->bPS )
    { delete[] CNode_->bPS;
        CNode_->bPS = 0;
    }
    if( CNode_->xPA )
    { delete[] CNode_->xPA;
        CNode_->xPA = 0;
    }

    if( CNode_->bSP )
    { delete[] CNode_->bSP;
        CNode_->bSP = 0;
    }
    if( CNode_->amru )
    { delete[] CNode_->amru;
        CNode_->amru = 0;
    }
    if( CNode_->amrl )
    { delete[] CNode_->amrl;
        CNode_->amrl = 0;
    }
}

// set default values(zeros) for DATABR structure
void databr_reset(DATABR *CnNde1, long int level)
{
    //  FMT variables (units or dimensionsless) - to be used for storing them
    //  at the nodearray level = 0.; normally not used in the single-node FMT-GEM coupling
    CnNde1->Tm = 0.;
    CnNde1->dt = 0.;
#ifndef NO_NODEARRAYLEVEL
    CnNde1->Dif = 0.;
    CnNde1->Vt = 0.;
    CnNde1->vp = 0.;
    CnNde1->eps = 0.;
    CnNde1->Km = 0.;
    CnNde1->Kf = 0.;
    CnNde1->S = 0.;
    CnNde1->Tr = 0.;
    CnNde1->h = 0.;
    CnNde1->rho = 0.;
    CnNde1->al = 0.;
    CnNde1->at = 0.;
    CnNde1->av = 0.;
    CnNde1->hDl = 0.;
    CnNde1->hDt = 0.;
    CnNde1->hDv = 0.;
    CnNde1->nto = 0.; //19
#endif
    if(level <1 )
        return;

    CnNde1->NodeHandle = 0;
    CnNde1->NodeTypeHY = normal;
    CnNde1->NodeTypeMT = normal;
    CnNde1->NodeStatusFMT = Initial_RUN;
    CnNde1->NodeStatusCH = NEED_GEM_AIA;
    CnNde1->IterDone = 0;      //6

    // Chemical scalar variables
    CnNde1->TK = 0.;
    CnNde1->P = 0.;
    CnNde1->Vs = 0.;
    CnNde1->Vi = 0.;
    CnNde1->Ms = 0.;
    CnNde1->Mi = 0.;
    CnNde1->Gs = 0.;
    CnNde1->Hs = 0.;
    CnNde1->Hi = 0.;
    CnNde1->IC = 0.;
    CnNde1->pH = 0.;
    CnNde1->pe = 0.;
    CnNde1->Eh = 0.; //13

    if( level < 2 )
        return;

    // Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
    CnNde1->bIC = 0;
    CnNde1->rMB = 0;
    CnNde1->uIC = 0;
    CnNde1->xDC = 0;
    CnNde1->gam = 0;
    CnNde1->dul = 0;
    CnNde1->dll = 0;
    CnNde1->aPH = 0;
    CnNde1->xPH = 0;
    CnNde1->vPS = 0;
    CnNde1->mPS = 0;
    CnNde1->bPS = 0;
    CnNde1->xPA = 0;
    CnNde1->bSP = 0;
    CnNde1->amru = 0;
    CnNde1->amrl = 0;
    CnNde1->omPH = 0;
}

// set default values(zeros) for DATACH structure
void datach_reset(DATACH* CSD)
{
    CSD->nIC = 0;
    CSD->nDC = 0;
    CSD->nPH = 0;
    CSD->nPS = 0;
    CSD->nDCs = 0;
    CSD->nTp = 0;
    CSD->nPp = 0;
    CSD->iGrd = 0;
    CSD->nAalp = 0;
    CSD->nICb = 0;
    CSD->nDCb = 0;
    CSD->nPHb = 0;
    CSD->nPSb = 0;
    CSD->mLook = 0;
    // Lists = 0; vectors and matrices
    CSD->nDCinPH = 0;
    CSD->xic = 0;
    CSD->xdc = 0;
    CSD->xph = 0;  //18

    CSD->TKval = 0;
    CSD->Psat = 0;
    CSD->Pval = 0;
    CSD->A = 0;
    CSD->Ttol = 0.;
    CSD->Ptol = 0.;
    CSD->dRes1 = 0.;
    CSD->dRes2 = 0.;
    CSD->ICmm = 0;
    CSD->DCmm = 0;
    CSD->DD = 0;
    CSD->denW = 0;
    CSD->epsW = 0;
    CSD->denWg = 0;
    CSD->epsWg = 0;
    CSD->G0 = 0;
    CSD->V0 = 0;
    CSD->S0 = 0;
    CSD->H0 = 0;
    CSD->Cp0 = 0;
    CSD->A0 = 0;
    CSD->U0 = 0;
    CSD->ICNL = 0;
    CSD->DCNL = 0;
    CSD->PHNL = 0;
    CSD->ccIC = 0;
    CSD->ccDC = 0;
    CSD->ccPH = 0;
}

}  // dbr_dch_api
