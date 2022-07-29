//-------------------------------------------------------------------
// $Id$
//
/// \file node2.cpp
/// Implementation of TNode class functionality including initialization
/// and execution of the GEM IPM 3 kernel
/// Works with DATACH and DATABR structures
//
// Copyright (c) 2005-2012 S.Dmytriyeva, D.Kulik, G.Kosakowski, F.Hingerl
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

#include "node.h"
#include "num_methods.h"

double TNode::get_Ppa_sat( double Tk )
{
    long int i=0;
    for( i=0; i<CSD->nTp; i++ )
    {
        if( (CSD->TKval[i] + CSD->Ttol) > Tk && (CSD->TKval[i] - CSD->Ttol) < Tk )
        {
            if( CSD->Psat[i] > 1.1e-5 )
            {
                return CSD->Psat[i];
            }
            else
            {
                return 0;
            }
        }
    }
    return 0;
}

long int TNode::get_grid_index_Ppa_sat( double Tk )
{
    long int i=0;
    long int r=-1;
    for( i=0; i<CSD->nTp; i++ )
    {
        if( (CSD->TKval[i] + CSD->Ttol) > Tk && (CSD->TKval[i] - CSD->Ttol) < Tk )
        {
            if( CSD->Psat[i] > 1.1e-5 )
            {
                return i;
            }
            else
            {
                return r;
            }
        }
    }
    return r;
}

//return exact or interpolated saturated pressure value for the temperature Tk using the exported data from the *.dch file
double TNode::Get_Psat(double Tk)
{
    double psat = 0.;

    long int xT;

    xT = check_grid_T(Tk);

    if (xT >=0)
        psat = CSD->Psat[xT];
    else
        psat = LagranInterp1D(CSD->TKval, CSD->Psat,Tk, CSD->nTp, 3);

    return psat;
}

//-----------------------------------------------------------------
// work with lists

void *TNode::get_ptrTSolMod(int xPH) const
{
    return multi->pTSolMod(xPH);
}

//Returns DCH index of IC given the IC Name string (null-terminated)
// or -1 if no such name was found in the DATACH IC name list
long int TNode::IC_name_to_xCH( const char *Name ) const
{
    long int ii;
    size_t len = strlen( Name );
    len =  std::min<size_t>(len,MaxICN);

    for(ii = 0; ii<CSD->nIC; ii++ )
        if(!memcmp(Name, CSD->ICNL[ii], len ))
            if( len == MaxICN || CSD->ICNL[ii][len] == ' ' || CSD->ICNL[ii][len] == '\0' )
                return ii;
    return -1;
}

// Returns DCH index of DC given the DC Name string
// or -1 if no such name was found in the DATACH DC name list
long int TNode::DC_name_to_xCH( const char *Name ) const
{
    long int ii;
    size_t len = strlen( Name );
    len =  std::min<size_t>(len,MaxDCN);

    for( ii = 0; ii<CSD->nDC; ii++ )
        if(!memcmp(Name, CSD->DCNL[ii], std::min<size_t>(len,MaxDCN)))
            if( len == MaxDCN || CSD->DCNL[ii][len] == ' ' || CSD->DCNL[ii][len] == '\0' )
                return ii;
    return -1;
}

// Returns DCH index of Phase given the Phase Name string
// or -1 if no such name was found in the DATACH Phase name list
long int TNode::Ph_name_to_xCH( const char *Name ) const
{
    long int ii;
    size_t len = strlen( Name );
    len =  std::min<size_t>(len,MaxPHN);

    for( ii = 0; ii<CSD->nPH; ii++ )
        if(!memcmp(Name, CSD->PHNL[ii], std::min<size_t>(len,MaxPHN)))
            if( len == MaxPHN || CSD->PHNL[ii][len] == ' ' || CSD->PHNL[ii][len] == '\0' )
                return ii;
    return -1;
}

// Converts the IC DCH index into the IC DBR index
// or returns -1 if this IC is not used in the data bridge
long int TNode::IC_xCH_to_xDB( const long int xCH ) const
{
    for(long int ii = 0; ii<CSD->nICb; ii++ )
        if( CSD->xic[ii] == xCH )
            return ii;
    return -1;
}

// Converts the DC DCH index into the DC DBR index
// or returns -1 if this DC is not used in the data bridge
long int TNode::DC_xCH_to_xDB( const long int xCH ) const
{
    for(long int ii = 0; ii<CSD->nDCb; ii++ )
        if( CSD->xdc[ii] == xCH )
            return ii;
    return -1;
}

// Converts the Phase DCH index into the Phase DBR index
// or returns -1 if this Phase is not used in the data bridge
long int TNode::Ph_xCH_to_xDB( const long int xCH ) const
{
    for(long int ii = 0; ii<CSD->nPHb; ii++ )
        if( CSD->xph[ii] == xCH )
            return ii;
    return -1;
}

// Returns the DCH index of the first DC belonging to the phase with DCH index Phx
long int  TNode::Phx_to_DCx( const long int Phx ) const
{
    long int k, DCx = 0;
    for( k=0; k<CSD->nPHb; k++ )
    {
        if( k == Phx )
            break;
        DCx += CSD->nDCinPH[ k];
    }
    return DCx;
}

// Returns the DCH index of the Phase to which the Dependent Component with index xCH belongs
long int  TNode::DCtoPh_DCH( const long int xdc ) const
{
    long int k, DCx = 0;
    for( k=0; k<CSD->nPHb; k++ )
    {
        DCx += CSD->nDCinPH[ k];
        if( xdc < DCx )
            break;
    }
    return k;
}


// Returns the DCH index of the first DC belonging to the phase with DCH index Phx,
// plus returns through the nDCinPh (reference) parameter the number of DCs included into this phase
long int  TNode::PhtoDC_DCH( const long int Phx, long int& nDCinPh ) const
{
    long int k, DCx = 0;
    for( k=0; k<CSD->nPHb; k++ )
    {
        if( k == Phx )
            break;
        DCx += CSD->nDCinPH[ k];
    }
    nDCinPh = CSD->nDCinPH[k];
    return DCx;
}

// Returns the DBR index of the Phase to which the  Dependent Component with index xBR belongs
long int  TNode::DCtoPh_DBR( const long int xBR ) const
{
    long int DCxCH = DC_xDB_to_xCH( xBR );
    long int PhxCH = DCtoPh_DCH( DCxCH );
    return Ph_xCH_to_xDB(PhxCH);
}

// Returns the DBR index of the first DC belonging to the phase with DBR index Phx,
//plus returns through the nDCinPh (reference) parameter the number of DCs included into DBR for this phase
long int  TNode::PhtoDC_DBR( const long int Phx, long int& nDCinPh ) const
{
    long int ii, DCx, DCxCH, PhxCH, nDCinPhCH;

    PhxCH = Ph_xDB_to_xCH( Phx );
    DCxCH = PhtoDC_DCH( PhxCH, nDCinPhCH );

    DCx = -1;
    nDCinPh = 0;
    for( ii = 0; ii<CSD->nDCb; ii++ )
    {
        if( CSD->xdc[ii] >= DCxCH )
        {
            if( CSD->xdc[ii] >= DCxCH+nDCinPhCH  )
                break;
            nDCinPh++;
            if( DCx == -1)
                DCx = ii;
        }
    }
    return DCx;
}

// Test TK as lying in the vicinity of a grid point for the interpolation of thermodynamic data
// Return index of the node in lookup array or -1
long int  TNode::check_grid_T( double TK ) const
{
    long int jj;
    for( jj=0; jj<CSD->nTp; jj++)
        if( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol )
            return jj;
    return -1;
}

// Test P as lying in the vicinity of a grid point for the interpolation of thermodynamic data
// Return index of the node in lookup array or -1
long int  TNode::check_grid_P( double P ) const
{
    long int jj;
    for( jj=0; jj<CSD->nPp; jj++)
        if( fabs( P - CSD->Pval[jj] ) < CSD->Ptol )
            return jj;
    return -1;
}

// Tests TK and P as a grid point for the interpolation of thermodynamic data using DATACH
// lookup arrays. Returns -1L if interpolation is needed, or 1D index of the lookup array element
// if TK and P fit within the respective tolerances.
// For producing lookup arrays (in GEMS), we recommend using step for temperature less or equal to 10 degrees
// in order to assure good accuracy of interpolation especially for S0 and Cp0 of aqueous species.
long int  TNode::check_grid_TP(  double TK, double P ) const
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
        xT = check_grid_T( TK );
        xP = check_grid_P( P );
        if( xT >=0 && xP>= 0 )
            ndx =  xP * CSD->nTp + xT;
        return ndx;
    }
    return ndx;
}

// used in GEMSFIT only
//Sets new molar Gibbs energy G0(P,TK) value for Dependent Component
//in the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::Set_DC_G0(const long int xCH, const double P, const double TK, const double new_G0 )
{
    long int xTP, jj;

    if( check_TP( TK, P ) == false )
        return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
    {
        CSD->G0[ jj + xTP ]=new_G0;
        load_thermodynamic_data = false;
    }
    else {
        node_logger->error("ERROR: given P={} and TK={} pair is not provided in DATACH", P, TK);
    }
    return 0;
}

//Retrieves (interpolated) molar Gibbs energy G0(P,TK) value for Dependent Component
//from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
// Parameter norm defines in wnich units the value is returned: false - in J/mol; true (default) - in mol/mol
double TNode::DC_G0(const long int xCH, const double P, const double TK,  bool norm ) const
{
    long int xTP, jj;
    double G0;

    if( check_TP( TK, P ) == false )
        return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        G0 = CSD->G0[ jj + xTP ];
    else
        G0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->G0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 6 );

    if( norm )
        return G0/(R_CONSTANT * (TK));
    else
        return G0;
}

// Retrieves (interpolated, if necessary) molar volume V0(P,TK) value for Dependent Component (in J/Pa)
// from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_V0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double V0;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        V0 = CSD->V0[ jj + xTP ];
    else
        V0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->V0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 5 );
    return V0;
}


// Retrieves (interpolated) molar enthalpy H0(P,TK) value for Dependent Component (in J/mol)
// from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_H0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double H0;

    if( check_TP( TK, P ) == false )
        return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        H0 = CSD->H0[ jj + xTP ];
    else
        H0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->H0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 5 );
    return H0;
}

// Retrieves (interpolated) absolute molar enropy S0(P,TK) value for Dependent Component (in J/K/mol)
// from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_S0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double s0;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        s0 = CSD->S0[ jj + xTP ];
    else
        s0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->S0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 4 );
    return s0;
}

// Retrieves (interpolated) constant-pressure heat capacity Cp0(P,TK) value for Dependent Component (in J/K/mol)
// from the DATACH structure ( xCH is the DC DCH index) or 0.0, if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_Cp0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double cp0;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        cp0 = CSD->Cp0[ jj + xTP ];
    else
        cp0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->Cp0+jj,
                            P, TK, CSD->nTp, CSD->nPp, 3 );
    return cp0;
}

// Retrieves (interpolated) Helmholtz energy  of Dependent Component (in J/mol)
// from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_A0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double a0;

    if( check_TP( TK, P ) == false )
        return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        a0 = CSD->A0[ jj + xTP ];
    else
        a0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->A0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 5 );
    return a0;
}

// Retrieves (interpolated) Internal energy of  Dependent Component (in J/mol)
// from the DATACH structure ( xCH is the DC DCH index) or 7777777., if TK (temperature, Kelvin)
// or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DC_U0(const long int xCH, const double P, const double TK) const
{
    long int xTP, jj;
    double u0;

    if( check_TP( TK, P ) == false )
        return 7777777.;

    xTP = check_grid_TP( TK, P );
    jj =  xCH * gridTP();

    if( xTP >= 0 )
        u0 = CSD->U0[ jj + xTP ];
    else
        u0 = LagranInterp( CSD->Pval, CSD->TKval, CSD->U0+jj,
                           P, TK, CSD->nTp, CSD->nPp, 5 );
    return u0;
}

// Retrieves (interpolated) dielectric constant and its derivatives of liquid water at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
void TNode::EpsArrayH2Ow( const double P, const double TK, std::vector<double>& EpsAW )
{
    long int xTP;
    EpsAW.resize(5);
    long int nTP = CSD->nTp;

    if( check_TP( TK, P ) == false )
        return;

    xTP = check_grid_TP( TK, P );

    if( xTP >= 0 )
    {
        EpsAW[0] = CSD->epsW[ 0*nTP + xTP ];
        EpsAW[1] = CSD->epsW[ 1*nTP + xTP ];
        EpsAW[2] = CSD->epsW[ 2*nTP + xTP ];
        EpsAW[3] = CSD->epsW[ 3*nTP + xTP ];
        EpsAW[4] = CSD->epsW[ 4*nTP + xTP ];
    }
    else
    {
        EpsAW[0] = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+0*nTP,
                                 P, TK, CSD->nTp, CSD->nPp, 5 );
        EpsAW[1] = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+1*nTP,
                                 P, TK, CSD->nTp, CSD->nPp, 5 );
        EpsAW[2] = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+2*nTP,
                                 P, TK, CSD->nTp, CSD->nPp, 5 );
        EpsAW[3] = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+3*nTP,
                                 P, TK, CSD->nTp, CSD->nPp, 5 );
        EpsAW[4] = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+4*nTP,
                                 P, TK, CSD->nTp, CSD->nPp, 5 );
    }

}


// Retrieves (interpolated) density and its derivatives of liquid water at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
void TNode::DensArrayH2Ow( const double P, const double TK, std::vector<double>& DensAW )
{
    long int xTP;
    DensAW.resize(5);
    long int nTP = CSD->nTp;

    if( check_TP( TK, P ) == false )
        if( P > 1e-5 )
            return;


    xTP = check_grid_TP( TK, P );

    if( xTP < 0 )
        xTP = get_grid_index_Ppa_sat( TK );


    if( xTP >= 0 )
    {
        DensAW[0] = CSD->denW[ 0*nTP + xTP ];
        DensAW[1] = CSD->denW[ 1*nTP + xTP ];
        DensAW[2] = CSD->denW[ 2*nTP + xTP ];
        DensAW[3] = CSD->denW[ 3*nTP + xTP ];
        DensAW[4] = CSD->denW[ 4*nTP + xTP ];
    }
    else
    {
        DensAW[0] = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+0*nTP,
                                  P, TK, CSD->nTp, CSD->nPp, 5 );
        DensAW[1] = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+1*nTP,
                                  P, TK, CSD->nTp, CSD->nPp, 5 );
        DensAW[2] = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+2*nTP,
                                  P, TK, CSD->nTp, CSD->nPp, 5 );
        DensAW[3] = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+3*nTP,
                                  P, TK, CSD->nTp, CSD->nPp, 5 );
        DensAW[4] = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+4*nTP,
                                  P, TK, CSD->nTp, CSD->nPp, 5 );
    }
}


// Retrieves (interpolated) dielectric constant of liquid water at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::EpsH2Ow(const double P, const double TK)
{
    long int xTP, jj;
    double epsW;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 *gridTP();

    if( xTP >= 0 )
        epsW = CSD->epsW[ jj + xTP ];
    else
        epsW = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsW+jj,
                             P, TK, CSD->nTp, CSD->nPp, 5 );
    return epsW;
}

// Retrieves (interpolated) density of liquid water (in kg/m3) at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DenH2Ow(const double P, const double TK)
{
    long int xTP, jj;
    double denW;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * gridTP();

    if( xTP >= 0 )
        denW = CSD->denW[ jj + xTP ];
    else
        denW = LagranInterp( CSD->Pval, CSD->TKval, CSD->denW+jj,
                             P, TK, CSD->nTp, CSD->nPp, 5 );
    return denW;
}

// Retrieves (interpolated) dielectric constant of H2O vapor at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::EpsH2Og(const double P, const double TK)
{
    long int xTP, jj;
    double epsWg;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * gridTP();

    if( xTP >= 0 )
        epsWg = CSD->epsWg[ jj + xTP ];
    else
        epsWg = LagranInterp( CSD->Pval, CSD->TKval, CSD->epsWg+jj,
                              P, TK, CSD->nTp, CSD->nPp, 5 );
    return epsWg;
}

// Retrieves (interpolated) density of H2O vapor (in kg/m3) at (P,TK) from the DATACH structure or 0.0,
// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
double TNode::DenH2Og(const double P, const double TK)
{
    long int xTP, jj;
    double denWg;

    if( check_TP( TK, P ) == false )
        return 0.;

    xTP = check_grid_TP( TK, P );
    jj = 0; // 0 * gridTP();

    if( xTP >= 0 )
        denWg = CSD->denWg[ jj + xTP ];
    else
        denWg = LagranInterp( CSD->Pval, CSD->TKval, CSD->denWg+jj,
                              P, TK, CSD->nTp, CSD->nPp, 5 );
    return denWg;
}

//Retrieves the current phase volume in m3 ( xph is DBR phase index) in the reactive sub-system.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  TNode::Ph_Volume( const long int xBR ) const
{
    double vol;
    if( xBR < CSD->nPSb )
    {    // Phase-solution
        vol = CNode->vPS[xBR];
        // Perhaps not yet accounting for the volume of mixing!
    }
    else
    {
        long int xdc = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
        vol = DC_V0( xdc, CNode->P, CNode->TK );
        vol *= CNode->xDC[DC_xCH_to_xDB(xdc)];
    }
    return vol;
}

// Added on Oct 1, 2020
// Retrieves the current phase total Gibbs energy in J (xph is DBR phase index) in the reactive sub-system.
// Chemical potential of phase (J/mol) is Gibbs energy divided by moles of the phase.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  TNode::Ph_GibbsEnergy( const long int xph ) const
{
    long int xic;
    double mol = CNode->xPH[xph];
    if(mol < 1e-20)
        return 0.0;
    // Getting bulk elemental composition vector of the phase
    double *PhBI = new double[ CSD->nICb ];
    PhBI = Ph_BC( xph, PhBI );
    double phGn = 0.0, phG = 0.0;
    for(xic = 0L; xic < CSD->nICb; xic++ )
    {
        phGn += PhBI[xic] * (CNode->uIC[xic]);
    }
    delete[] PhBI;
    phG = phGn * R_CONSTANT * (CNode->TK); // phG in J/mol
    // Multiplying by moles amount of the phase
    return phG*mol;
}

// Retrieves the current phase enthalpy in J (xph is DBR phase index) in the reactive sub-system.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  TNode::Ph_Enthalpy( const long int xph ) const
{
    double ent, enth = 0.0;
    long int xdc, xdcb, xdce, nDCinPh, xch;

    double mol = CNode->xPH[xph];
    if(mol < 1e-20)
        return 0.0;
    // Getting the DBR index of the first DC belonging to the phase with DBR index xBR,
    // with nDCinPh being the number of DCs included into DBR for this phase
    xdcb = PhtoDC_DBR( xph, nDCinPh );
    xdce = xdcb + nDCinPh;
    node_logger->debug("Ph_Enthalpy xph: {}  xdcb: {}  xdce: {}", xph, xdcb, xdce);
    for(xdc = xdcb; xdc < xdce; xdc++ )
    {
        xch = DC_xDB_to_xCH( xdc ); // getting DCH index from DBR index of DC
        // Retrieves (interpolated) molar enthalpy H0(P,TK) value for Dependent Component (in J/mol)
        ent = DC_H0( xch, CNode->P, CNode->TK );
        if( ent < 7777777.0 )
            enth += ent * CNode->xDC[xdc];
        // else out of P or T range of interpolation
        node_logger->debug("Ph_Enthalpy xdc: {}  xch: {}  ent: {} enth: {}", xdc, xch, ent, enth);
    }
    // Not yet accounting for the enthalpy of mixing!  TBD

    return enth;
}

// Added on Oct 1, 2020
// Retrieves the current phase entropy in J/K (xph is DBR phase index) in the reactive sub-system.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  TNode::Ph_Entropy( const long int xph ) const
{
    double ent, entr = 0.0;
    long int xdc, xdcb, xdce, nDCinPh, xch;

    double mol = CNode->xPH[xph];
    if(mol < 1e-20)
        return 0.0;
    // Getting the DBR index of the first DC belonging to the phase with DBR index xBR,
    // with nDCinPh being the number of DCs included into DBR for this phase
    xdcb = PhtoDC_DBR( xph, nDCinPh );
    xdce = xdcb + nDCinPh;
    node_logger->debug("Ph_Entropy xph: {}  xdcb: {}  xdce: {}", xph, xdcb, xdce);
    for(xdc = xdcb; xdc < xdce; xdc++ )
    {
        xch = DC_xDB_to_xCH( xdc ); // getting DCH index from DBR index of DC
        // Retrieves (interpolated) molar entropy S0(P,TK) value for Dependent Component (in J/K/mol)
        ent = DC_S0( xch, CNode->P, CNode->TK );
        if( noZero(ent) )
            entr += ent * (CNode->xDC[xdc]);
        // else out of P or T range of interpolation
        node_logger->debug("Ph_Entropy xdc: {}  xch: {}  ent: {} entr: {}", xdc, xch, ent, entr);
    }
    // Not yet accounting for the enthalpy of mixing!  TBD
    return entr;  // in J/K
}

// Added on Oct 1, 2020
// Retrieves the current phase heat capacity Cp in J/K (xph is DBR phase index) in the reactive sub-system.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  TNode::Ph_HeatCapacityCp( const long int xph ) const
{
    double cap, capp = 0.0;
    long int xdc, xdcb, xdce, nDCinPh, xch;

    double mol = CNode->xPH[xph];
    if(mol < 1e-20)
        return 0.0;
    // Getting the DBR index of the first DC belonging to the phase with DBR index xBR,
    // with nDCinPh being the number of DCs included into DBR for this phase
    xdcb = PhtoDC_DBR( xph, nDCinPh );
    xdce = xdcb + nDCinPh;
    for(xdc = xdcb; xdc < xdce; xdc++ )
    {
        xch = DC_xDB_to_xCH( xdc ); // getting DCH index from DBR index of DC
        // Retrieves (interpolated) heat capacity Cp0(TK) value for Dependent Component (in J/K/mol)
        cap = DC_Cp0( xch, CNode->P, CNode->TK );
        if( noZero(cap) )
            capp += cap * (CNode->xDC[xdc]);
        // else out of P or T range of interpolation
        node_logger->debug("Ph_HeatCapacityCp xdc: {}  xch: {}  cap: {} capp: {}", xdc, xch, cap, capp);
    }
    // Not yet accounting for the enthalpy of mixing!  TBD
    return capp;  // in J/K
}

//Retrieves the current phase amount in moles (xph is DBR phase index) in the reactive sub-system.
double  TNode::Ph_Moles( const long int xBR ) const
{
    double mol;
    mol = CNode->xPH[xBR];

    return mol;
}

//Obsolete: Retrieves the current phase amount in moles ( xph is DBR phase index) in the reactive sub-system.
double  TNode::Ph_Mole( const long int xBR ) const
{
    double mol;
    mol = CNode->xPH[xBR];

    return mol;
}

// Retrieves the phase mass in kg (xph is DBR phase index).
// Works for multicomponent and for single-component phases. Returns 0.0 if phase amount is zero.
double  TNode::Ph_Mass( const long int xBR ) const
{
    double mass;
    if( xBR < CSD->nPSb )
        mass = CNode->mPS[xBR];
    else
    {
        long int xDC = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
        mass = CNode->xDC[ DC_xCH_to_xDB(xDC) ] * CSD->DCmm[xDC];
    }
    return mass;
}

// Retrieves the phase saturation index (xBR is DBR phase index).
// Works for multicomponent and for single-component phases.
double TNode::Ph_SatInd(const long int xBR ) const
{
    double SatX;
    SatX = CNode->omPH[xBR];
    return SatX;
    /*    double SatInd = 0.0;
    long int jj, dcx1, Ndc;
    dcx1 = PhtoDC_DBR( xph, Ndc );
    if( xph < CSD->nPSb )
    {
        for( jj=dcx1; jj<Ndc+dcx1; jj++)
            SatInd +=  Get_aDC( jj )/Get_gDC(jj);
    }
    else
      SatInd = Get_aDC( dcx1 );
    if( SatInd > 0.0 )
        SatInd = log10(SatInd);
    SatInd = pmm->Falp[xph]; // Falps[xph] contains zeros; // Fixed by DK on 8.10.2018 (temporarily)
    return SatInd;
*/
}

// Retrieval of the phase bulk composition ( xph is DBR phase index) into memory indicated by
// ARout (array of at least [dCH->nICb elements]). Returns pointer to ARout which may also be
// allocated inside of Ph_BC() in the case if parameter ARout = NULL is specified;
// to avoid a memory leak, you will have to free this memory wherever appropriate.
// This function works for multicomponent and for single-component phases
double* TNode::Ph_BC( const long int xBR, double* ARout ) const
{
    long int ii;
    if( !ARout )
        ARout = new double[ CSD->nICb ];   // Potential memory leak ! ! ! ! ! ! ! !

    if( xBR < CSD->nPSb )
        for( ii=0; ii<pCSD()->nICb; ii++ )
            ARout[ii] = CNode->bPS[ xBR * CSD->nICb + ii ];
    else
    {
        long int DCx = Phx_to_DCx( Ph_xDB_to_xCH(xBR) );
        for( ii=0; ii<pCSD()->nICb; ii++ )
        {
            ARout[ii] = CSD->A[ IC_xDB_to_xCH(ii) + DCx * CSD->nIC];
            ARout[ii] *= CNode->xDC[ DC_xCH_to_xDB(DCx) ];
        }
    }
    return ARout;
}

// Retrieval of (dual-thermodynamic) chemical potential of the DC (xdc is the DC DBR index).
// Parameter norm defines the scale: if true (1) then in mol/mol, otherwise in J/mol
double TNode::Get_muDC( const long int xdc, bool norm ) const
{	long int xCH, ii;
    double muDC = 0;

    xCH = DC_xDB_to_xCH(xdc);
    for( ii=0; ii<pCSD()->nICb; ii++ )
        muDC += CSD->A[  xCH * CSD->nIC + IC_xDB_to_xCH(ii) ] * CNode->uIC[ ii ];

    if( norm )
        return muDC;
    else
        return muDC*(R_CONSTANT * (CNode->TK));
}

//Retrieval of (dual-thermodynamic) activity of the DC (xdc is the DC DBR index)
//If parameter scale is true then activity is returned, if false then log10(activity)
double TNode::Get_aDC( const long int xdc, bool scale ) const
{
    double Mj  = Get_muDC( xdc, true );
    double Mj0 = DC_G0( DC_xDB_to_xCH(xdc), CNode->P, CNode->TK,  true );
    if( scale )
        return exp( Mj-Mj0 );
    else // decimal log
        return 0.4342944819 *( Mj-Mj0 );
    // return 	pow(10.0,pmm->Y_la[xCH]);
}

// Retrieves concentration of Dependent Component (xdc is the DC DBR index) in its phase
// in the respective concentration scale. For aqueous species, molality is returned;
// for gas species, mole fraction not partial pressure; for surface complexes - molality;
// for species in other phases - mole fraction. If DC has zero amount, the function returns 0.0.
double TNode::Get_cDC( const long int xdc ) const
{
    long int xph = DCtoPh_DBR( xdc);
    long int DCxCH = DC_xDB_to_xCH(xdc);
    double DCcon = 0.;

    switch( CSD->ccDC[DCxCH] )
    {
    case DC_SCP_CONDEN:

    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM:

    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR:
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
    case DC_SCM_SPECIES:

    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:
        if( noZero(CNode->xPH[xph]) )
            DCcon =  CNode->xDC[xdc]/CNode->xPH[xph];  //pmp->Wx[xCH];
        break;
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:   if( noZero(CNode->xPH[xph]) )
            DCcon =  CNode->xDC[xdc]/CNode->xPH[xph]; // *CNode->P;
        break;
    case DC_AQ_PROTON:
    case DC_AQ_SPECIES:
    case DC_AQ_SURCOMP:

    case DC_SUR_GROUP:
    case DC_SSC_A0:
    case DC_SSC_A1:
    case DC_SSC_A2:
    case DC_SSC_A3:
    case DC_SSC_A4:
    case DC_WSC_A0:
    case DC_WSC_A1:
    case DC_WSC_A2:
    case DC_WSC_A3:
    case DC_WSC_A4:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
    {
        double MMC = 0., Factor;
        long int PhxCH_aquel = Ph_xDB_to_xCH(0);
        long int H2Obr, jj;

        if(CSD->ccPH[PhxCH_aquel] != PH_AQUEL )
            break;                                  // no aquel phase in DATABR

        if( CNode->xPA[0] > 1e-19 )
        {
            for(jj=0; jj<CSD->nDCinPH[PhxCH_aquel]; jj++)
                if( CSD->ccDC[jj] == DC_AQ_SOLVENT || CSD->ccDC[jj] == DC_AQ_SOLVCOM )
                {  H2Obr = DC_xCH_to_xDB(jj);
                    if( H2Obr >= 0 )
                        MMC += CSD->DCmm[jj]*CNode->xDC[H2Obr]/CNode->xPA[0];
                }
        }
        else MMC=0.01801528; // Assuming water-solvent
        if( (CNode->xPA[0] > 1e-19) && (MMC > 1e-19) )
            Factor = 1./MMC/CNode->xPA[0]; // molality
        else Factor = 0.0;

        DCcon =  CNode->xDC[xdc]*Factor;  //pmp->Y_m[xCH];
    }
        break;
    default:
        break; // error in DC class code
    }
    return DCcon;
}

// Added 6.12.2011 DK
// Retrieves total dissolved aqueous molality of Independent Component with DBR index xIC
// or returns 0.0 if there is no water in the node or no aqueous phase in DATACH
double TNode::Get_mIC( const long xic ) const
{
    long int xaq = pmm->LO; // DC_name_to_xDB( "H2O@" ); //  24.03.21 "H2O(l)" index of H2O aq in DATABR
    double nAQ, nIC, scICinH2O, m_tot;
    if( (CNode->xPA && CNode->xPA[0] <= 0.) || xaq < 0L )  // check phase code in DATACH
        return 0.; // no water or zero amount of water
    nAQ = CNode->xPA[0];
    scICinH2O = DCaJI( xaq, xic ); // get stoich coeff of this IC in H2O
    nIC = CNode->bPS[ xic ]; // get amount of this IC in aq phase
    nIC -= scICinH2O * nAQ;
    m_tot = nIC * 55.5084350618 / nAQ;
    return m_tot;
}

// Added 31.01.2013 DM
// Retrives pH of the aqueous solution
double TNode::Get_pH( )
{
    double p_pH;
    p_pH = CNode->pH;
    return p_pH;
}

// Added 28.01.2014 DM
// Retrives pH of the aqueous solution
double TNode::Get_pe( )
{
    double p_pe;
    p_pe = CNode->pe;
    return p_pe;
}

// Added 12.06.2013 DM
// Retrives Eh of the aqueous solution
double TNode::Get_Eh( )
{
    double p_Eh;
    p_Eh = CNode->Eh;
    return p_Eh;
}

// Added 17.02.2014 DM
// Retrives effective molal ionic strength of aqueous solution
double TNode::Get_IC( )
{
    double p_IC;
    p_IC = CNode->IC;
    return p_IC;
}

// Access to equilibrium properties of phases and components using DATACH indexation

// Retrieves the current (dual-thermodynamic) activity of DC (xCH is DC DCH index)
// directly from GEM IPM work structure. Also activity of a DC not included into DATABR list
// can be retrieved. If DC has zero amount, its dual-thermodynamic activity is returned anyway.
// For single condensed phase component, this value has a meaning of the saturation index,
// also in the presence of metastability constraint(s).
double TNode::DC_a(const long int xCH)
{
    //double Mj  = DC_mu( xCH, true );
    //double Mj0 = DC_G0( xCH, CNode->P, CNode->TK,  true );
    //return (Mj-Mj0)/2.302585093;
    return 	pow(10.0,pmm->Y_la[xCH]);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Functions needed by GEMSFIT. Setting parameters for activity coefficient models.
//     aIPx     = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
//     aIPc     = pmp->PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
//     aDCc     = pmp->DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
//     NComp    = pmp->L1[k];         // Number of components in the phase
//     NPar     = pmp->LsMod[k*3];    // Number of interaction parameters
//     NPcoef   = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
//     MaxOrd   = pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)
//     NP_DC    = pmp->LsMdc[k*3];    // Number of non-ideality coeffs per one DC in multicomponent phase
//     NsSit    = pmp->LsMdc[k*3+1];  // Number of sublattices considered in a multisite mixing model (0 if no sublattices considered)
//     NsMoi    = pmp->LsMdc[k*3+2];  // Total number of moieties considered in sublattice phase model (0 if no sublattices considered)

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Retrieves indices of origin in TSolMod composite arrays for a phase of interest index_phase.
// Parameters IN: index_phase is the DCH index of phase of interest.
// Parameters OUT: ipaIPx, ipaIPc, ipaDCc are origin indices of this phase in aIPx, aIPc and aDCc arrays, respectively.
void TNode::Get_IPc_IPx_DCc_indices( long int &ipaIPx, long int &ipaIPc, long int &ipaDCc, const long int &index_phase )
{
    long int ip_IPx=0; long int ip_IPc=0; long int ip_DCc=0;

    for( int k=0; k < index_phase; k++ )
    {
        ip_IPx  += pmm->LsMod[k*3] * pmm->LsMod[k*3+1];
        ip_IPc  += pmm->LsMod[k*3] * pmm->LsMod[k*3+2];
        ip_DCc  += pmm->LsMdc[k] * pmm->L1[k];
    }
    ipaIPx = ip_IPx;
    ipaIPc = ip_IPc;
    ipaDCc = ip_DCc;
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Retrieves dimensions of TSolMod arrays for a phase of interest index_phase.
// Parameters IN: index_phase is the DCH index of phase of interest.
// Parameters OUT: NPar, NPcoef, MaxOrd, NComp, NP_DC, are number of interaction parameters, number of coefficients per parameter,
// maximum parameter order (i.e. row length in aIPx), number of components in the phase, and number of coefficients per component, respectively.
void TNode::Get_NPar_NPcoef_MaxOrd_NComp_NP_DC ( long int &NPar, long int &NPcoef, long int &MaxOrd,
                                                 long int &NComp, long int &NP_DC, const long int &index_phase )
{
    NPar   = pmm->LsMod[(index_phase)*3];
    NPcoef = pmm->LsMod[(index_phase)*3+2];
    MaxOrd = pmm->LsMod[(index_phase)*3+1];
    NComp  = pmm->L1[(index_phase)];
    NP_DC  = pmm->LsMdc[(index_phase)*3];
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Sets values of the aIPc array (of interaction parameter coefficients) for the solution phase of interest index_phase.
// Parameters IN: vaIPc - vector with the contents of the aIPc sub-array to be set; ipaIPc is the origin index (of the first element)
//    of the aIPc array; index_phase is the DCH index of phase of interest.
void TNode::Set_aIPc ( const std::vector<double> aIPc, const long int &ipaIPc, const long int &index_phase )
{
    long int rc, NPar, NPcoef;
    NPar = pmm->LsMod[ index_phase * 3 ];
    NPcoef =  pmm->LsMod[ index_phase * 3 + 2 ];
    if( aIPc.size() != (unsigned int)(NPar*NPcoef) )
    {
        node_logger->critical(" TNode::Set_aIPc() error: vector aIPc does not match the dimensions specified in the GEMS3K IPM file (NPar*NPcoef) !!!! \n"
                              " aIPc.size() = {}, NPar*NPcoef = {} bailing out now ... \n", aIPc.size(), NPar*NPcoef);
        exit(1);
    }
    for ( rc=0;rc<(NPar*NPcoef);rc++ )
    {
        (pmm->PMc[ ipaIPc + rc ]) = aIPc[ rc ];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
    }
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Gets values of the aIPc array (of interaction parameter coefficients) for the solution phase of interest index_phase.
// Parameters IN: ipaIPc is the origin index (of the first element) of the aIPc array; index_phase is the DCH index of phase of interest.
// Parameters OUT: returns vaIPc - vector with the contents of the aIPc sub-array.
void TNode::Get_aIPc ( std::vector<double> &aIPc, const long int &ipaIPc, const long int &index_phase )
{
    long int i, NPar, NPcoef;
    NPar   = pmm->LsMod[ index_phase * 3 ];
    NPcoef = pmm->LsMod[ index_phase * 3 + 2 ];
    aIPc.clear();
    aIPc.resize( (NPar*NPcoef) );
    i = 0;
    while (i<(NPar*NPcoef))
    {
        aIPc[ i ]   = pmm->PMc[ ipaIPc + i ];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
        i++;
    }
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Gets values of the aIPx list array (of indexes of interacting moieties or components) for the solution phase of interest index_phase.
// Parameters IN: ipaIPx is the origin index (of the first element) of the aIPx array; index_phase is the DCH index of phase of interest.
// Parameters OUT: returns vaIPx - vector with the contents of the aIPx sub-array.
void TNode::Get_aIPx ( std::vector<long int> &aIPx, const long int &ipaIPx, const long int &index_phase )
{
    long int i, NPar, MaxOrd;
    NPar   = pmm->LsMod[ index_phase * 3 ];
    MaxOrd = pmm->LsMod[ index_phase * 3 + 1 ];
    aIPx.clear();
    aIPx.resize( (NPar*MaxOrd) );
    i = 0;
    while (i<(NPar*MaxOrd))
    {
        aIPx[ i ]   = pmm->IPx[ ipaIPx + i];		// pointer to list of indices of interaction param coeffs, NPar * MaxOrd
        i++;
    }
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Sets values of the aDCc array (of components property coefficients) for the solution phase of interest index_phase.
// Parameters IN: vaDCc - vector with the contents of the aDCc sub-array to be set. ipaDCc is the origin index (of the first element)
//    of the aDCc array; index_phase is the DCH index of phase of interest.
void TNode::Set_aDCc( const std::vector<double> aDCc, const long int &ipaDCc, const long int &index_phase )
{
    long int rc, NComp, NP_DC;
    NComp = pmm->L1[ index_phase ];
    NP_DC = pmm->LsMdc[ index_phase ];
    if( aDCc.size() != (unsigned int)(NComp*NP_DC) )
    {
        node_logger->critical("TNode::Set_aDCc() error: vector aDCc does not match the dimensions specified in the GEMS3K IPM file (NComp*NP_DC) !!!! "
                              " aDCc.size() = {}, NComp*NP_DC = {} bailing out now ... \n", aDCc.size(), NComp*NP_DC);
        exit(1);
    }
    for ( rc=0;rc<(NComp*NP_DC);rc++ )
    {
        (pmm->DMc[ ipaDCc + rc ]) = aDCc[ rc ];		// end-member param coeffs, NComp * NP_DC
    }
}

// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class
// Gets values of the aDCc array (of components property coefficients) for the solution phase of interest index_phase.
// Parameters IN: ipaDCc is the origin index (of the first element) of the aDCc array; index_phase is the DCH index of phase of interest.
// Parameters OUT: returns vaDCc - vector with the contents of the aDCc sub-array.
void TNode::Get_aDCc( std::vector<double> &aDCc, const long int &index_phase_aDCc, const long int &index_phase )
{
    long int i, NComp, NP_DC;
    NComp = pmm->L1[ index_phase ];
    NP_DC = pmm->LsMdc[ index_phase ];
    aDCc.clear();
    aDCc.resize( (NComp*NP_DC) );
    i = 0;
    while (i<(NComp*NP_DC))
    {
        aDCc[ i ]   = pmm->DMc[ index_phase_aDCc + i ];  // pointer to list of indices of interaction param coeffs, NPar * MaxOrd
        i++;
    }
}

// direct access to set temperature in the current (work) node
void TNode::Set_Tk( const double &T_k )
{
    CNode->TK = T_k;
}

// direct access to set pressure (given in bar) in the current (work) node
void TNode::Set_Pb( const double &P_b )
{
    CNode->P = P_b * 1e5;  // in the node, pressure is given in Pa (see databr.h)!
}



// Retrieves the current concentration of Dependent Component (xCH is DC DCH index) in its
// phase directly from the GEM IPM work structure. Also activity of a DC not included into
// DATABR list can be retrieved. For aqueous species, molality is returned; for gas species,
// partial pressure; for surface complexes - density in mol/m2; for species in other phases -
// mole fraction. If DC has zero amount, the function returns 0.0.
double TNode::DC_c(const long int xCH)
{
    double DCcon = 0.;
    switch( pmm->DCC[xCH] )
    {
    case DC_SCP_CONDEN: DCcon =  pmm->Wx[xCH];
        break;
    case DC_AQ_PROTON:
    case DC_AQ_SPECIES:
    case DC_AQ_SURCOMP: DCcon =  pmm->Y_m[xCH];
        break;
    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM: DCcon =  pmm->Wx[xCH];
        break;
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:    DCcon =  pmm->Wx[xCH]*(pmm->P*bar_to_Pa);
        break;
    case DC_SOL_IDEAL:
    case DC_SCM_SPECIES:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR: DCcon =  pmm->Wx[xCH];
        break;
    case DC_SUR_GROUP:
    case DC_SSC_A0:
    case DC_SSC_A1:
    case DC_SSC_A2:
    case DC_SSC_A3:
    case DC_SSC_A4:
    case DC_WSC_A0:
    case DC_WSC_A1:
    case DC_WSC_A2:
    case DC_WSC_A3:
    case DC_WSC_A4:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:     DCcon =  pmm->Y_m[xCH];
        break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER: DCcon =  pmm->Wx[xCH];
        break;
    default:
        node_logger->warn(" error in DC class code {}", pmm->DCC[xCH]);
        break; // error in DC class code
    }
    return DCcon;
}

// Retrieves the current (dual-thermodynamic) chemical potential of DC (xCH is DC DCH index)
// directly from GEM IPM work structure, also for any DC not included into DATABR or having zero amount.
// Parameter norm defines in wnich units the chemical potential value is returned:
// false - in J/mol; true (default) - in mol/mol
double TNode::DC_mu(const long int xCH, bool norm)
{
    double muDC = pmm->Fx[xCH];
    //  for(long ii=0; ii<CSD->nIC; ii++ )
    //    muDC += pmm->A[  xCH * CSD->nIC + ii ] * (pmm->U[ii]);
    if( norm )
        return muDC/pmm->RT; // (R_CONSTANT * (CNode->TK));
    else
        return muDC;
}

// Retrieves the standard chemical potential of DC (xCH is DC DCH index) directly
// from GEM IPM work structure at current pressure and temperature,
// also for any DC not included into DATABR or having zero amount.
// Parameter norm defines in which units the chemical potential value is returned:
// false - in J/mol; true (default) - in mol/mol
double TNode::DC_mu0(const long int xCH, bool norm)
{
    return  DC_G0( xCH, CNode->P, CNode->TK, norm );
    /*double  G0 = pmm->G0[xCH];
    if( norm )
      return G0;
    else
      return G0*pmm->RT;
    */
}

//-----------------------End of node2.cpp--------------------------



