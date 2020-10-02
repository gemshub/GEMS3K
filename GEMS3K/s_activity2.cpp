//-------------------------------------------------------------------
// $Id: s_activity2.cpp 845 2013-06-20 15:58:57Z kulik $
//
/// \file ipm_chemical2.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, adsorption models etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, S.Dmitrieva, K.Chudnenko
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
//


#include <cmath>
#include "node.h"
#include "activities.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating bulk stoichiometry of a multicomponent phase
//
void TActivity::phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] )
{
    long int ii, i, j;
    double Xx;

    if( !A || !X || !BF )
        return;
    for( i=0; i<N; i++ )
          BF[i] = 0.;

    for( j=0; j<M; j++ )
    {
        Xx = X[j];
        if( fabs( Xx ) < 1e-16 )  // was 1e-12
            continue;
        for( ii=arrL[j+jb]; ii<arrL[j+jb+1]; ii++ )
        {  i = arrAN[ii];
            BF[i] += A[i+j*N] * Xx;
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Adds phase to total bulk stoichiometry of all solid phases in the system
// Done on request by TW in November 2006
//
void TActivity::phase_bfc( long int k, long int jj )
{
    long int ii, i, j;
    double Xx;

    if( act.PHC[k] == PH_AQUEL || act.PHC[k] == PH_GASMIX ||
        act.PHC[k] == PH_FLUID || act.PHC[k] == PH_PLASMA ||
        act.PHC[k] == PH_SIMELT || act.PHC[k] == PH_LIQUID )
        return;
    for( j=0; j<act.L1[k]; j++ )
    {
        Xx = act.X[j+jj];
        if( fabs( Xx ) < 1e-12 )
            continue;
        for( ii=arrL[j+jj]; ii<arrL[j+jj+1]; ii++ )
        {  i = arrAN[ii];
           act.BFC[i] += act.A[i+(jj+j)*act.N] * Xx;
        }
    }
}

/// Returns mass of all solid phases in grams (from the BFC vector)
double TActivity::bfc_mass( void )
{
   double TotalMass = 0.;
   for(long int i = 0; i<act.N; i++ )
     TotalMass += act.BFC[i]*act.Awt[i];
   return TotalMass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  a(j,i) ((*(act.A+(i)+(j)*act.N)))
/// Calculation of dual chemical potentials, activities, and primal
/// concentrations for DCs (indexed jb to je) in a k-th phase.
// Input arrays X, XF, XFA,  input factors: Factor, MMC
//
void TActivity::CalculateConcentrationsInPhase( double X[], double XF[], double /*XFA*/[],
              double Factor, double MMC, double /*Dsur*/, long int jb, long int je, long int k)
{
    long int j;
    double Muj, /* DsurT=0.0,*/ SPmol, lnFmol=4.016535;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    if( act.PHC[0] == PH_AQUEL )
    {  // mole fraction to molality conversion
        if( !k ) lnFmol = log(1000./MMC);  // aq species
        else lnFmol = log( H2O_mol_to_kg ); // 4.016535; 	   // other species
    }

    for( j=jb; j<je; j++ )
    { // loop over DC - with important bugfixes from 02.04.2003
        Muj = DC_DualChemicalPotential( act.U, act.A+j*act.N, act.NR, j );
        act.Fx[j] = Muj; // *pmp->RT; Fixed DK 12.03.2012

//        if( X[j] <= act.lowPosNum )
        if( X[j] <= act.DcMinM )
        { // zeroing off
            act.Wx[j] = 0.0;
//            act.VL[j] = log( lowPosNum );
//            act.Y_w[j] = 0.0;
            act.lnGam[j] = 0.0;
            if( act.PHC[0] == PH_AQUEL )
               act.Y_m[j] = 0.0;
            switch( act.DCC[j] ) // choice of expressions
            {                      // since 10.03.2008, changed the concept of DualTh activity
               case DC_SCP_CONDEN:
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
                    break;
               case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                    act.Y_la[j] = ln_to_lg*(Muj - act.G0[j] + lnFmol );
                    break;
               case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                    act.Y_la[j] = ln_to_lg* (Muj - act.G0[j] );
                    break;
               case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:   // gases
               case DC_GAS_H2: case DC_GAS_N2:
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
                    if( act.Pc > 1e-29 )
                        act.Y_la[j] += log10( act.Pc );
                    break;
               case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR: case DC_SOL_MINDEP: case DC_SOL_MAJDEP:
             case DC_SCM_SPECIES:
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
                    break;
               case DC_SUR_GROUP:
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] ); // + lnFmol ); corr. 06.10.10 DK
                    break;
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3:
               case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2:
               case DC_WSC_A3: case DC_WSC_A4: case DC_SUR_COMPLEX:
               case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] ); // + lnFmol ); corr. 06.10.10 DK
                    break;
               case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                    act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
                    break;
               default:
                    break; // error in DC class code
            }
            continue;
        }
        // calculation of the mole fraction
        act.Wx[j] = X[j]/XF[k];
//        if( X[j] > DcMinM )
//            act.VL[j] = log( act.Wx[j] );     // this is used nowhere except in some scripts. Remove?
//        else act.VL[j] = log( DcMinM );   // was lowPosNum
        act.Y_la[j] = 0.0;
        switch( act.DCC[j] ) // choice of expressions
        {
        case DC_SCP_CONDEN:
            act.Wx[j] = 1;
//            act.VL[j] = 0.0;
            if( act.LO )
            {   //  bugfix DK 08.02.10
                act.Y_m[j] = X[j]*Factor; // molality
            }
//            act.Y_w[j] = // mass % of the system
//                          1e2 * X[j] * act.MM[j] / act.MBX;
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
            act.FVOL[k] += act.Vol[j]*X[j];
            break;
        case DC_AQ_ELECTRON:
            act.Y_m[j] = 0.0;
            act.Y_la[j] = 0.0 - act.pe;
//            act.Y_w[j] = 0.0;
            break;
        case DC_AQ_PROTON:  // in molal scale!
            act.pH = -ln_to_lg*(Muj-act.G0[j] + lnFmol );
             [[fallthrough]];
        case DC_AQ_SPECIES:
            act.IC += 0.5* X[j]*Factor *(act.EZ[j]*act.EZ[j]); // increment to effective molal ionic strength
             [[fallthrough]];
        case DC_AQ_SURCOMP:
            SPmol = X[j]*Factor;  // molality
 //           act.IC += 0.5* SPmol *(act.EZ[j]*act.EZ[j]); // Bugfix DK 21.10.2011 'K' species don't count here!
            act.FVOL[k] += act.Vol[j]*X[j]; // fixed 04.02.03 KD
            act.Y_m[j] = SPmol;
            act.Y_la[j] = ln_to_lg*(Muj - act.G0[j] + lnFmol );
//            act.Y_w[j] = 1e6 * X[j] * act.MM[j] / act.FWGT[k];
//            if(  act.DCC[j] != DC_AQ_SURCOMP )
//            {  // Bugfix 21.10.2011  DK - excluding 'K' species from total IC molality
               //  Optimized for performance - calculation inline
//               for( i=arrL[j]; i<arrL[j+1]; i++ )
//               {  ii = arrAN[i];
//                  if( ii>= act.NR )
//                    continue;
//                  act.IC_m[ii] += SPmol* a(j,ii);  // total aqueous molality
//                  act.IC_wm[ii] += X[j]* a(j,ii);  // total aqueous mass concentration
//               }
//            }
            break;
        case DC_AQ_SOLVENT: // mole fractions of solvent
        case DC_AQ_SOLVCOM:
//            act.Y_m[j] = X[j]/XFA[k];  Replaced by DK 12.03.2012
            act.Y_m[j] = H2O_mol_to_kg;
//            act.Y_w[j] = 1e3*X[j]*act.MM[j]/act.FWGT[k];
            act.FVOL[k] += act.Vol[j]*X[j];
            act.Y_la[j] = ln_to_lg* (Muj - act.G0[j] );
            break;
        case DC_GAS_COMP:
        case DC_GAS_H2O:
        case DC_GAS_CO2:   // gases
        case DC_GAS_H2:
        case DC_GAS_N2:
            act.FVOL[k] += act.Vol[j]*X[j];
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
            if( act.Pc > 1e-9 )
                act.Y_la[j] += log10( act.Pc );
            break;
        case DC_SOL_IDEAL:
        case DC_SOL_MINOR:   //solution end member
        case DC_SOL_MAJOR:
        case DC_SOL_MINDEP:
        case DC_SOL_MAJDEP:
    case DC_SCM_SPECIES:
            act.FVOL[k] += act.Vol[j]*X[j];
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] );
            if( act.LO )
            {   //  bugfix DK 16.04.2012
                act.Y_m[j] = X[j]*Factor; // molality
            }
            break;
        case DC_SUR_GROUP: // adsorption:
            act.Y_m[j] = X[j]*Factor; // molality
//            act.Y_w[j] =  // mg/g sorbent
//                1e3 * X[j] * act.MM[j] / (MMC*XFA[k]);
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] ); // + lnFmol ); corr. 06.10.10
            act.FVOL[k] += act.Vol[j]*X[j]; // fixed 11.03.2008 KD
            break;
        case DC_SSC_A0:
        case DC_SSC_A1:
        case DC_SSC_A2:
        case DC_SSC_A3:
        case DC_SSC_A4:
        case DC_WSC_A0:
        case DC_WSC_A1:
        case DC_WSC_A2:
        case DC_WSC_A3:
        case DC_WSC_A4:  // case DC_SUR_GROUP:
        case DC_SUR_COMPLEX:
        case DC_SUR_IPAIR:
        case DC_IESC_A:
        case DC_IEWC_B:
            act.Y_m[j] = X[j]*Factor; // molality
//            act.Y_w[j] =  // mg/g sorbent
//                1e3 * X[j] * act.MM[j] / (MMC*XFA[k]);
//           DsurT = MMC * act.Aalp[k] * pa->p.DNS*1.66054e-6;
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] ); // + lnFmol - act.GEX[j] + Dsur + DsurT/( 1.0+DsurT )
            act.FVOL[k] += act.Vol[j]*X[j]; // fixed 11.03.2008   KD
            break;
        case DC_PEL_CARRIER:
        case DC_SUR_MINAL:
        case DC_SUR_CARRIER: // sorbent
            act.Y_m[j] = X[j]*Factor; // molality
//            act.Y_w[j] = 0.0;
            if( act.YF[0] >= act.DSM )
//              act.Y_w[j] = // mg of sorbent per kg aq solution
//                1e6 * X[j] * act.MM[j] / act.FWGT[0];
            act.Y_la[j] = ln_to_lg * ( Muj - act.G0[j] ); // - act.GEX[j] + Dsur - 1. + 1./(1.+Dsur) - DsurT + DsurT/(1+DsurT)
            act.FVOL[k] += act.Vol[j]*X[j];
            break;
        default:
            break; // error in DC class code
        }
        ;
    }   // j
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of derived values (concentrations etc.) on IPM iteration
///  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
// This function has to be rewritten using new set of built-in
// chemical functions.
//
void TActivity::CalculateConcentrations( double X[], double XF[], double XFA[])
{
    long int k, i, j, jj;
    double Factor=0.0, Dsur=0.0, MMC=0.0, VXc;//, YFk;
//    SPP_SETTING *pa = paTProfil;

//    if( act.Ls < 2 || !act.FIs )  Temporary disabled  09.03.2010 DK
//        return;

    for( i=0; i< act.N; i++ )
     act.BFC[i] = 0.;         // cleanup for phase_bfc() calculation

    for( j=0; j<act.Ls; j++ )
    {
        act.Wx[j] = 0.;
//        act.VL[j] = 0.;
    }
    j=0;
    VXc = 0.0;
    for( k=0; k<act.FI; k++ )
    { // cycle by phases
        i=j+act.L1[k];
        act.FWGT[k] = 0.0;
        act.FVOL[k] = 0.0;
        //   Dsur = 0.0;

        if( XF[k] > act.DSM && !( act.PHC[k] == PH_SORPTION && XFA[k] <= act.ScMinM ))
           phase_bfc( k, j );

        if( k >= act.FIs || act.L1[k] == 1 )
        { // this is a single- component phase
            act.Wx[j] = 1.0; // SD 04/05/2010
            if( XF[k] < act.DSM )
            {   // This phase to be zeroed off
                if( act.LO )
                    act.Y_m[j] = 0.0;
                act.Wx[j] = 0.0;
//                act.Y_w[j] = 0.0;
                act.Fx[j] = DC_DualChemicalPotential( act.U, act.A+j*act.N, act.NR, j );
                act.Y_la[j] = ln_to_lg * ( act.Fx[j] - act.G0[j] ); // -act.GEX[j]
//                act.Fx[j] *= act.RT;     // el-chem potential
                goto NEXT_PHASE;
            }
//            act.VL[j] = 0.0;
            if( act.LO && XFA[0] > 0 )
                act.Y_m[j] = X[j] * 1000./18.01528/XFA[0]; // molality
//            act.Y_w[j] = // mass % in the system
//                1e2 * X[j] * act.MM[j] / act.MBX;
            act.Fx[j] = DC_DualChemicalPotential( act.U, act.A+j*act.N, act.NR, j );
            act.Y_la[j] = ln_to_lg * ( act.Fx[j] - act.G0[j] ); // - act.GEX[j]
//            act.Fx[j] *= act.RT;     // el-chem potential
            act.FWGT[k] += X[j] * act.MM[j];
            act.FVOL[k] += X[j] * act.Vol[j];
            goto NEXT_PHASE;
        }
//        if( act.LO && !k )
//        {    // aqueous phase present
//            act.IC=0.;                   // fix 20.06.13 DK
//            for( ii=0; ii<act.N; ii++ )
//            {
//                act.IC_m[ii] = 0.0;
//                act.IC_lm[ii] = 0.0;
//                act.IC_wm[ii] = 0.0;
//            }
//        }
        if( XF[k] <= act.DSM ||
            (act.PHC[k] == PH_AQUEL && ( XFA[k] <= act.XwMinM || XF[k] <= act.DSM ) )
                || ( act.PHC[k] == PH_SORPTION && XFA[k] <= act.ScMinM ))
        {
            for( jj=0; jj<act.N; jj++)
             act.BF[k*act.N+jj] = 0.;

            for(jj=j; jj<i; jj++)   // Loop added 10.03.01  KD (GTDEMO)
            {
                act.Wx[j] = 0.0;
                if( act.LO )
                    act.Y_m[jj] = 0.0;
//                act.Y_w[jj] = 0.0;
                act.Fx[jj] = DC_DualChemicalPotential( act.U, act.A+jj*act.N, act.NR, jj );
                act.Y_la[jj] = ln_to_lg * ( act.Fx[jj] - act.G0[jj] );
                if(act.PHC[k] == PH_AQUEL ) // || act.PHC[k] == PH_SORPTION ) corr. 06.10.10 DK
                   act.Y_la[jj] += 1.74438;
                if(act.PHC[k] == PH_GASMIX || act.PHC[k] == PH_PLASMA )
                   act.Y_la[jj] += log10( act.Pc );
//                act.Fx[jj] *= act.RT;     // el-chem potential
                act.lnGam[jj] = 0.0;
            }
            goto NEXT_PHASE;
        }
        // calculate bulk stoichiometry of a multicomponent phase
        phase_bcs( act.N, act.L1[k], j, act.A+j*act.N, act.X+j, act.BF+k*act.N );

        switch( act.PHC[k] )
        {
        case PH_AQUEL:
            MMC = 0.0; // molar mass of carrier
//            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
//            if( XFA[k] > act.lowPosNum )
            if( XFA[k] > act.XwMinM )
            {
                for(jj=j; jj<i; jj++)
                    if( act.DCC[jj] == DC_AQ_SOLVENT ||
                            act.DCC[jj] == DC_AQ_SOLVCOM )
                        MMC += act.MM[jj]*X[jj]/XFA[k];
            }
            else MMC=18.01528; // Assuming water-solvent
//            if( (XFA[k] > act.lowPosNum) && (MMC > act.lowPosNum) )
            if( (XFA[k] > act.XwMinM) && (MMC > act.XwMinM ) )
                Factor = 1000./MMC/XFA[k]; // molality
            else Factor = 0.0;
//            act.IC=0.;
            act.pe = ln_to_lg* act.U[act.N-1];
            act.Eh = 0.000086 * act.U[act.N-1] * act.Tc;
             [[fallthrough]];
        case PH_GASMIX:
        case PH_FLUID:
        case PH_PLASMA:
        case PH_SIMELT:
        case PH_HCARBL:
        case PH_SINCOND:
        case PH_SINDIS:
        case PH_LIQUID:
//    case PH_IONEX:
//   case PH_ADSORPT:
            //YFk = XF[k];
            for(jj=j; jj<i; jj++)
            {
                if( X[jj] > act.DcMinM)      // fixed 30.08.2009 DK
                    act.FWGT[k] += X[jj]*act.MM[jj];
            }
            break;
/*
        case PH_POLYEL:
        case PH_SORPTION: // only sorbent end-members!
            YFk = XFA[k];
            MMC=0.0;

            for( ist=0; ist<act.FIat; ist++ )
                act.XFTS[k][ist] = 0.0;
//           if( XFA[k] < act.lowPosNum ) XFA[k] = act.lowPosNum;
            if( XFA[k] < act.ScMinM ) XFA[k] = act.ScMinM;
            for( jj=j; jj<i; jj++ )
            {
               jja = jj - ( act.Ls - act.Lads );
                if( act.DCC[jj] == DC_SUR_CARRIER ||
                        act.DCC[jj] == DC_SUR_MINAL ||
                        act.DCC[jj] == DC_PEL_CARRIER )
                {
                    MMC += act.MM[jj]*X[jj]/XFA[k];
                    // Only sorbent mass
                    act.FWGT[k] += X[jj]*act.MM[jj];
                }
                else
                {
/ *!!!!! * /           ist = act.SATX[jja][XL_ST];
                    act.XFTS[k][ist] += X[jj];
                }
            }
            act.logYFk = log(act.YFk);
//            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
//            if( Dsur <= -1.0 ) Dsur = -0.999999;
            break;
*/
        default:
             return; // Phase class code error!
        }
        // calculation of species concentrations in k-th phase
        CalculateConcentrationsInPhase( X, XF, XFA, Factor, MMC, Dsur, j, i, k );

NEXT_PHASE:
        VXc += act.FVOL[k];
/*
        if( act.PHC[k] == PH_AQUEL && XF[k] > DSM && XFA[k] > XwMinM )
            for( ii=0; ii<act.NR; ii++ )
            {
               if( act.LO  )
               { if( act.IC_m[ii] >= pa->p.DB )
                    act.IC_lm[ii] = ln_to_lg*log( act.IC_m[ii] );
                else
                    act.IC_lm[ii] = 0;
                if( act.FWGT[k] >= pa->p.DB )
                    act.IC_wm[ii] *= act.Awt[ii]*1000./act.FWGT[k];
                else
                    act.IC_wm[ii] = 0;
               }
            }
*/
        j = i;
    }  // k
}


//--------------------- End of s_activity2.cpp ---------------------------
