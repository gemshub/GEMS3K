//-------------------------------------------------------------------
// $Id: ipm_chemical2.cpp 825 2006-03-29 07:10:23Z gems $
//
// Copyright (C) 1992-2007  D.Kulik, S.Dmitrieva, K.Chudnenko, I.Karpov
//
// Implementation of chemistry-specific functions (concentrations,
// activity coefficients, adsorption models etc.)
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806); Kulik (2000), Geoch.Cosmoch.Acta 64,3161-3179
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry and of the
// standalone GEMIPM2K code (define IPMGEMPLUGIN).
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#include <math.h>
#include "m_param.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculating bulk stoichiometry of a multicomponent phase
//
void TMulti::phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] )
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
// Adds phase to total bulk stoichiometry of all solid phases in the system
// Done on request by TW in November 2006
//
void TMulti::phase_bfc( long int k, long int jj )
{
    long int ii, i, j;
    double Xx;

    if( pmp->PHC[k] == PH_AQUEL || pmp->PHC[k] == PH_GASMIX ||
        pmp->PHC[k] == PH_FLUID || pmp->PHC[k] == PH_PLASMA ||
        pmp->PHC[k] == PH_SIMELT || pmp->PHC[k] == PH_LIQUID )
        return;
    for( j=0; j<pmp->L1[k]; j++ )
    {
        Xx = pmp->X[j+jj];
        if( fabs( Xx ) < 1e-12 )
            continue;
        for( ii=arrL[j+jj]; ii<arrL[j+jj+1]; ii++ )
        {  i = arrAN[ii];
           pmp->BFC[i] += pmp->A[i+(jj+j)*pmp->N] * Xx;
        }
    }
}

// returns mass of all solid phases in grams (from the BFC vector)
double TMulti::bfc_mass( void )
{
   double TotalMass = 0.;
   for(long int i = 0; i<pmp->N; i++ )
     TotalMass += pmp->BFC[i]*pmp->Awt[i];
   return TotalMass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  a(j,i) ((*(pmp->A+(i)+(j)*pmp->N)))
// Calculation of dual chemical potentials, activities, and primal
// concentrations for DCs (indexed jb to je) in a k-th phase.
// Input arrays X, XF, XFA,  input factors: Factor, MMC
//
//  Do we need this all in GEMIPM ?
//
void TMulti::PH_CalculateConcentrations( double X[], double XF[], double XFA[],
              double Factor, double MMC, double /*Dsur*/, long int jb, long int je, long int k)
{
    long int j, ii, i;
    double Muj, /* DsurT=0.0,*/ SPmol, lnFmol=4.016535;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pmp->PHC[0] == PH_AQUEL )
    {  // mole fraction to molality conversion
        if( !k ) lnFmol = log(1000./MMC);  // aq species
        else lnFmol = log( H2O_mol_to_kg ); // 4.016535; 	   // other species
    }

    for( j=jb; j<je; j++ )
    { // loop over DC - with important bugfixes from 02.04.2003
        Muj = DC_DualChemicalPotential( pmp->U, pmp->A+j*pmp->N, pmp->NR, j );
        pmp->Fx[j] = Muj * pmp->RT;     // el-chem potential

//        if( X[j] <= pmp->lowPosNum )
        if( X[j] <= pmp->DcMinM )
        { // zeroing off
            pmp->Wx[j] = 0.0;
            pmp->VL[j] = log( pmp->lowPosNum );
            pmp->Y_w[j] = 0.0;
            pmp->lnGam[j] = 0.0;
            if( pmp->PHC[0] == PH_AQUEL )
               pmp->Y_m[j] = 0.0;
            switch( pmp->DCC[j] ) // choice of expressions
            {                      // since 10.03.2008, changed the concept of DualTh activity
               case DC_SCP_CONDEN:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
                    break;
               case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                    pmp->Y_la[j] = ln_to_lg*(Muj - pmp->G0[j] + lnFmol );
                    break;
               case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                    pmp->Y_la[j] = ln_to_lg* (Muj - pmp->G0[j] );
                    break;
               case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:   // gases
               case DC_GAS_H2: case DC_GAS_N2:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
                    if( pmp->Pc > 1e-29 )
                        pmp->Y_la[j] += log10( pmp->Pc );
                    break;
               case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
                    break;
               case DC_SUR_GROUP:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] + lnFmol );
                    break;
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3:
               case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2:
               case DC_WSC_A3: case DC_WSC_A4: case DC_SUR_COMPLEX:
               case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] + lnFmol );
                    break;
               case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
                    break;
               default:
                    break; // error in DC class code
            }
            continue;
        }
        // calculation of the mole fraction
        pmp->Wx[j] = X[j]/XF[k];
        if( X[j] > min( pmp->lowPosNum, pmp->DcMinM ) )
            pmp->VL[j] = log( pmp->Wx[j] );     // this is used nowhere except in some scripts. Remove?
        else pmp->VL[j] = log( pmp->lowPosNum );   // debugging 29.11.05 KD
        pmp->Y_la[j] = 0.0;
        switch( pmp->DCC[j] ) // choice of expressions
        {
        case DC_SCP_CONDEN:
            pmp->Wx[j] = 1;
            pmp->VL[j] = 0.0;
            if( pmp->LO )
            {   //  bugfix DK 08.02.10
                pmp->Y_m[j] = X[j]*Factor; // molality
            }
            pmp->Y_w[j] = // mass % of the system
                          1e2 * X[j] * pmp->MM[j] / pmp->MBX;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            break;
        case DC_AQ_ELECTRON:
            pmp->Y_m[j] = 0.0;
            pmp->Y_la[j] = 0.0 - pmp->pe;
            pmp->Y_w[j] = 0.0;
            break;
        case DC_AQ_PROTON:  // in molal scale!
            pmp->pH = -ln_to_lg*(Muj-pmp->G0[j] + lnFmol );
        case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
            SPmol = X[j]*Factor;  // molality
            pmp->IC += 0.5* SPmol *(pmp->EZ[j]*pmp->EZ[j]); // increment to effective molal ionic strength
            pmp->FVOL[k] += pmp->Vol[j]*X[j]; // fixed 04.02.03 KD
            pmp->Y_m[j] = SPmol;
            pmp->Y_la[j] = ln_to_lg*(Muj - pmp->G0[j] + lnFmol );
            pmp->Y_w[j] = 1e6 * X[j] * pmp->MM[j] / pmp->FWGT[k];
//  Optimized for performance - calculation inline
            for( i=arrL[j]; i<arrL[j+1]; i++ )
            {  ii = arrAN[i];
               if( ii>= pmp->NR )
                continue;
                pmp->IC_m[ii] += SPmol* a(j,ii);  // total aqueous molality
                pmp->IC_wm[ii] += X[j]* a(j,ii);  // total aqueous mass concentration
            }
            break;
        case DC_AQ_SOLVENT: // mole fractions of solvent
        case DC_AQ_SOLVCOM:
            pmp->Y_m[j] = X[j]/XFA[k];
            pmp->Y_w[j] = 1e3*X[j]*pmp->MM[j]/pmp->FWGT[k];
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg* (Muj - pmp->G0[j] );
            break;
        case DC_GAS_COMP:
        case DC_GAS_H2O:
        case DC_GAS_CO2:   // gases
        case DC_GAS_H2:
        case DC_GAS_N2:
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
            if( pmp->Pc > 1e-9 )
                pmp->Y_la[j] += log10( pmp->Pc );
            break;
        case DC_SOL_IDEAL:
        case DC_SOL_MINOR:   //solution end member
        case DC_SOL_MAJOR:
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] );
            break;
        case DC_SUR_GROUP: // adsorption:
            pmp->Y_m[j] = X[j]*Factor; // molality
            pmp->Y_w[j] =  // mg/g sorbent
                1e3 * X[j] * pmp->MM[j] / (MMC*XFA[k]);
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] + lnFmol );
            pmp->FVOL[k] += pmp->Vol[j]*X[j]; // fixed 11.03.2008 KD
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
            pmp->Y_m[j] = X[j]*Factor; // molality
            pmp->Y_w[j] =  // mg/g sorbent
                1e3 * X[j] * pmp->MM[j] / (MMC*XFA[k]);
//           DsurT = MMC * pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] + lnFmol ); // - pmp->GEX[j] + Dsur + DsurT/( 1.0+DsurT )
            pmp->FVOL[k] += pmp->Vol[j]*X[j]; // fixed 11.03.2008   KD
            break;
        case DC_PEL_CARRIER:
        case DC_SUR_MINAL:
        case DC_SUR_CARRIER: // sorbent
            pmp->Y_m[j] = X[j]*Factor; // molality
            pmp->Y_w[j] = 0.0;
            if( pmp->YF[0] >= pmp->DSM )
              pmp->Y_w[j] = // mg of sorbent per kg aq solution
                1e6 * X[j] * pmp->MM[j] / pmp->FWGT[0];
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] ); // - pmp->GEX[j] + Dsur - 1. + 1./(1.+Dsur) - DsurT + DsurT/(1+DsurT)
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            break;
        default:
            break; // error in DC class code
        }
        ;
    }   // j
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of derived values (concentrations etc.) on IPM iteration
//  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
// This function has to be rewritten using new set of built-in
// chemical functions.
//
void TMulti::CalculateConcentrations( double X[], double XF[], double XFA[])
{
    long int k, ii, i, j, ist, jj, jja;
    double Factor=0.0, Dsur=0.0, MMC=0.0;
    SPP_SETTING *pa = &TProfil::pm->pa;

//    if( pmp->Ls < 2 || !pmp->FIs )  Temporary disabled  09.03.2010 DK
//        return;

    for( i=0; i< pmp->N; i++ )
     pmp->BFC[i] = 0.;

    for( j=0; j<pmp->Ls; j++ )
    {
        pmp->Wx[j] = 0.;
        pmp->VL[j] = 0.;
    }
    j=0;
    pmp->VXc = 0.0;
    for( k=0; k<pmp->FI; k++ )
    { // cycle by phases
        i=j+pmp->L1[k];
        pmp->FWGT[k] = 0.0;
        pmp->FVOL[k] = 0.0;
        //   Dsur = 0.0;

        if( XF[k] > pmp->DSM && !( pmp->PHC[k] == PH_SORPTION && XFA[k] <= pa->p.ScMin ))
           phase_bfc( k, j );

        if( k >= pmp->FIs || pmp->L1[k] == 1 )
        { // this is a single- component phase
            pmp->Wx[j] = 1.0; // SD 04/05/2010
            if( XF[k] < pmp->DSM )
            {
                if( pmp->LO )
                    pmp->Y_m[j] = 0.0;
                pmp->Y_w[j] = 0.0;
                pmp->Fx[j] = DC_DualChemicalPotential( pmp->U, pmp->A+j*pmp->N, pmp->NR, j );
                pmp->Y_la[j] = ln_to_lg * ( pmp->Fx[j] - pmp->G0[j] ); // -pmp->GEX[j]
                pmp->Fx[j] *= pmp->RT;     // el-chem potential
                goto NEXT_PHASE;
            }
            //pmp->Wx[j] = 1.0;
            pmp->VL[j] = 0.0;
            if( pmp->LO && XFA[0] > 0 )
                pmp->Y_m[j] = X[j] * 1000./18.01528/XFA[0]; // molality
            pmp->Y_w[j] = // mass % in the system
                1e2 * X[j] * pmp->MM[j] / pmp->MBX;
            pmp->Fx[j] = DC_DualChemicalPotential( pmp->U, pmp->A+j*pmp->N, pmp->NR, j );
            pmp->Y_la[j] = ln_to_lg * ( pmp->Fx[j] - pmp->G0[j] ); // - pmp->GEX[j]
            pmp->Fx[j] *= pmp->RT;     // el-chem potential
            pmp->FWGT[k] += X[j] * pmp->MM[j];
            pmp->FVOL[k] += X[j] * pmp->Vol[j];
            goto NEXT_PHASE;
        }
        if( pmp->LO && !k )
            for( ii=0; ii<pmp->N; ii++ )
            {
                pmp->IC_m[ii] = 0.0;
                pmp->IC_lm[ii] = 0.0;
                pmp->IC_wm[ii] = 0.0;
            }

//        if( XF[k] <= pmp->DSM ||
//     (pmp->PHC[k] == PH_AQUEL && ( XFA[k] <= pmp->lowPosNum*1e3 || XF[k] <= pa->p.XwMin ) )
//                || ( pmp->PHC[k] == PH_SORPTION && XFA[k] <= pa->p.ScMin ))
        if( XF[k] <= pmp->DSM ||
            (pmp->PHC[k] == PH_AQUEL && ( XFA[k] <= pmp->XwMinM || XF[k] <= pmp->DSM ) )
                || ( pmp->PHC[k] == PH_SORPTION && XFA[k] <= pmp->ScMinM ))
        {
            for( jj=0; jj<pmp->N; jj++)
             pmp->BF[k*pmp->N+jj] = 0.;

            for(jj=j; jj<i; jj++)   // Loop added 10.03.01  KD (GTDEMO)
            {
                pmp->Wx[j] = 0.0;
                if( pmp->LO )
                    pmp->Y_m[jj] = 0.0;
                pmp->Y_w[jj] = 0.0;
                pmp->Fx[jj] = DC_DualChemicalPotential( pmp->U, pmp->A+jj*pmp->N, pmp->NR, jj );
                pmp->Y_la[jj] = ln_to_lg * ( pmp->Fx[jj] - pmp->G0[jj] );
                if(pmp->PHC[k] == PH_AQUEL || pmp->PHC[k] == PH_SORPTION )
                   pmp->Y_la[jj] += 1.74438;
                if(pmp->PHC[k] == PH_GASMIX || pmp->PHC[k] == PH_PLASMA )
                   pmp->Y_la[jj] += log10( pmp->Pc );
                pmp->Fx[jj] *= pmp->RT;     // el-chem potential
                pmp->lnGam[jj] = 0.0;
            }
            goto NEXT_PHASE;
        }
        // calculate bulk stoichiometry of a multicomponent phase
        phase_bcs( pmp->N, pmp->L1[k], j, pmp->A+j*pmp->N, pmp->X+j, pmp->BF+k*pmp->N );

        switch( pmp->PHC[k] )
        {
        case PH_AQUEL:
            MMC = 0.0; // molar mass of carrier
//            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
//            if( XFA[k] > pmp->lowPosNum )
            if( XFA[k] > pmp->XwMinM )
            {
                for(jj=j; jj<i; jj++)
                    if( pmp->DCC[jj] == DC_AQ_SOLVENT ||
                            pmp->DCC[jj] == DC_AQ_SOLVCOM )
                        MMC += pmp->MM[jj]*X[jj]/XFA[k];
            }
            else MMC=18.01528; // Assuming water-solvent
//            if( (XFA[k] > pmp->lowPosNum) && (MMC > pmp->lowPosNum) )
            if( (XFA[k] > pmp->XwMinM) && (MMC > pmp->lowPosNum) )
                Factor = 1000./MMC/XFA[k]; // molality
            else Factor = 0.0;
            pmp->IC=0.;
            pmp->pe = ln_to_lg* pmp->U[pmp->N-1];
            pmp->Eh = 0.000086 * pmp->U[pmp->N-1] * pmp->T;
        case PH_GASMIX:
        case PH_FLUID:
        case PH_PLASMA:
        case PH_SIMELT:
        case PH_HCARBL:
        case PH_SINCOND:
        case PH_SINDIS:
        case PH_LIQUID:
            pmp->YFk = XF[k];
            for(jj=j; jj<i; jj++)
            {
                if( X[jj] > pmp->DcMinM)      // fixed 30.08.2009 DK
                    pmp->FWGT[k] += X[jj]*pmp->MM[jj];
            }
            break;
        case PH_POLYEL:
        case PH_SORPTION: // only sorbent end-members!
            pmp->YFk = XFA[k];
            MMC=0.0;

            for( ist=0; ist<pmp->FIat; ist++ )
                pmp->XFTS[k][ist] = 0.0;
//           if( XFA[k] < pmp->lowPosNum ) XFA[k] = pmp->lowPosNum;
            if( XFA[k] < pmp->ScMinM ) XFA[k] = pmp->ScMinM;
            for( jj=j; jj<i; jj++ )
            {
               jja = jj - ( pmp->Ls - pmp->Lads );
                if( pmp->DCC[jj] == DC_SUR_CARRIER ||
                        pmp->DCC[jj] == DC_SUR_MINAL ||
                        pmp->DCC[jj] == DC_PEL_CARRIER )
                {
                    MMC += pmp->MM[jj]*X[jj]/XFA[k];
                    // Only sorbent mass
                    pmp->FWGT[k] += X[jj]*pmp->MM[jj];
                }
                else
                {
/*!!!!!*/           ist = pmp->SATX[jja][XL_ST];
                    pmp->XFTS[k][ist] += X[jj];
                }
            }
            pmp->logYFk = log(pmp->YFk);
//            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
//            if( Dsur <= -1.0 ) Dsur = -0.999999;
            break;
        default:
             return; // Phase class code error!
        }
        // calculation of species concentrations in k-th phase
        PH_CalculateConcentrations( X, XF, XFA, Factor, MMC, Dsur, j, i, k );

NEXT_PHASE:
        pmp->VXc += pmp->FVOL[k];
        if( pmp->PHC[k] == PH_AQUEL && XF[k] > pmp->DSM && XFA[k] > pmp->XwMinM )
            for( ii=0; ii<pmp->NR; ii++ )
            {
               if( pmp->LO  )
               { if( pmp->IC_m[ii] >= pa->p.DB )
                    pmp->IC_lm[ii] = ln_to_lg*log( pmp->IC_m[ii] );
                else
                    pmp->IC_lm[ii] = 0;
                if( pmp->FWGT[k] >= pa->p.DB )
                    pmp->IC_wm[ii] *= pmp->Awt[ii]*1000./pmp->FWGT[k];
                else
                    pmp->IC_wm[ii] = 0;
               }
            }
        j = i;
    }  // k
}

//--------------------------------------------------------------------------------
// Calculation of surface charge densities on multi-surface sorption phase
void TMulti::IS_EtaCalc()
{
    long int k, i, ist, isp, j=0, ja;
    double XetaS=0., XetaW=0.,  Ez, CD0, CDb;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    for( k=0; k<pmp->FIs; k++ )
    { // loop over phases
        i=j+pmp->L1[k];
        if( pmp->FIat > 0 )
            for( ist=0; ist<pmp->FIat; ist++ )
            {
                pmp->XetaA[k][ist] = 0.0;
                pmp->XetaB[k][ist] = 0.0;
                pmp->XetaD[k][ist] = 0.0;     // added 12.09.05  KD
            }

        if( pmp->XF[k] <= pmp->DSM ||
                (pmp->PHC[k] == PH_AQUEL && ( pmp->X[pmp->LO] <= pmp->XwMinM //  pa->p.XwMin
                 || pmp->XF[k] <= pmp->DHBM ) )
             || (pmp->PHC[k] == PH_SORPTION && pmp->XF[k] <= pmp->ScMinM ) ) //  pa->p.ScMin) )
            goto NEXT_PHASE;

        switch( pmp->PHC[k] )
        {  // initialization according to the phase class
        case PH_AQUEL:  // aqueous solution
            pmp->Yw = pmp->XFA[k];
            XetaW = 0.0;
            break;
        case PH_PLASMA:
        case PH_SIMELT:
            XetaS = 0.0;
            break;
        case PH_POLYEL:
        case PH_SORPTION: // reserved
            break;
        default:
            break;
        }
        for( ; j<i; j++ )
        { // loop over DC for calculating total phase charge
            if( pmp->X[j] <= pmp->lowPosNum*100. )
                continue; // Skipping too low concentrations
            ja = j - ( pmp->Ls - pmp->Lads );

            switch( pmp->DCC[j] ) // select expressions for species classes
            {
            case DC_AQ_ELECTRON:    case DC_AQ_PROTON:    case DC_AQ_SPECIES:  case DC_AQ_SURCOMP:
                XetaW += pmp->X[j]*pmp->EZ[j];
            case DC_AQ_SOLVCOM:    case DC_AQ_SOLVENT:
                break;
            case DC_PEL_CARRIER:  case DC_SUR_MINAL:
            case DC_SUR_CARRIER: // charge of carrier: ???????
                                 // pmp->XetaA[k] += pmp->X[j]*pmp->EZ[j];
                break;
                // surface species
            case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2:  case DC_SSC_A3:
            case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1:  case DC_WSC_A2:
            case DC_WSC_A3: case DC_WSC_A4:
            case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:
            case DC_IESC_A:
            case DC_IEWC_B: // Get ist - index of surface type
                            // and  isp - index of surface plane
                ist = pmp->SATX[ja][XL_ST]; // / MSPN;
                isp = pmp->SATX[ja][XL_SP]; // % MSPN;
                        // isp  index of outer surface charge allocation  (new)
                // Getting charge distribution information
                CD0 = pmp->MASDJ[ja][PI_CD0];
                    // species charge that goes into 0 plane
                CDb = pmp->MASDJ[ja][PI_CDB];
          // species charge that goes into 1, 2 or 3 plane according to isp value
                Ez = pmp->EZ[j];  // take formula charge as default
                if( !isp )
                { // This is the 0 (A) plane only - no charge distribution!
                    if( fabs( CD0 ) > 1e-20 ) // Only if 0-plane charge is given in the table
                       Ez = CD0;
                    pmp->XetaA[k][ist] += pmp->X[j]*Ez;
                }
                else
                { // The charge distribution (CD) is specified
                    if( pmp->SCM[k][ist] == SC_MTL )
                    {   // Modified TL: Robertson, 1997; also XTLM Kulik 2002
//                        if( fabs( CDb ) > 1e-20 )  // Doubtful...
//                           Ez = CDb;
                        pmp->XetaB[k][ist] += pmp->X[j]*CDb;
                    }
                    else if( pmp->SCM[k][ist] == SC_TLM )
                    {
// New CD version of TLM Hayes & Leckie, 1987  added 25.10.2004
                        pmp->XetaB[k][ist] += pmp->X[j] * CDb;
                        pmp->XetaA[k][ist] += pmp->X[j] * CD0;
                    }
                    else if( pmp->SCM[k][ist] == SC_3LM )
                    {
// CD 3-layer model (Hiemstra e.a. 1996) added 12.09.2005 by KD
                        if( isp == 1 )
                            pmp->XetaB[k][ist] += pmp->X[j] * CDb;
                        if( isp == 2 )
                            pmp->XetaD[k][ist] += pmp->X[j] * CDb;
                        pmp->XetaA[k][ist] += pmp->X[j] * CD0;
                    }
                    else if( pmp->SCM[k][ist] == SC_BSM )
                    { // Basic Stern model Christl & Kretzschmar, 1999
// New CD version of BSM  added 25.10.2004
                        pmp->XetaB[k][ist] += pmp->X[j] * CDb;
                        pmp->XetaA[k][ist] += pmp->X[j] * CD0;
                    }
                    else if( pmp->SCM[k][ist] == SC_MXC )
                    { // BSM for ion exchange on perm.charge surface
                        if( fabs( CDb ) > 1e-20 )  // Doubtful...
                           Ez = CDb;
                        pmp->XetaB[k][ist] += pmp->X[j]*Ez;
                        pmp->XetaA[k][ist] += pmp->X[j]*CD0;  // added for testing
                    }
                    else if( pmp->SCM[k][ist] == SC_CCM )
                    { // Added 25.07.03 to implement the extended CCM Nilsson ea 1996
// New CD version of BSM  added 25.10.2004
                           pmp->XetaB[k][ist] += pmp->X[j] * CDb;
                           pmp->XetaA[k][ist] += pmp->X[j] * CD0;
                    }
                 //    case DC_SUR_DL_ION:  XetaS += pmp->X[j]*pmp->EZ[j];
                }
                break;
            default:
                XetaS += pmp->X[j]*pmp->EZ[j];
                break;
            }
        }   // j
        // compare pmp->Xetaf[k]+pmp->XetaA[k]+pmp->XetaB[k] and XetaS
        // Test XetaW
NEXT_PHASE:
        j = i;
        if( pmp->LO && !k && pmp->FIat > 0 )
        {
            pmp->XetaA[k][0] = XetaW;
            pmp->XetaB[k][0] = XetaW;
            pmp->XetaD[k][0] = XetaW;
        }
        if( (pmp->PHC[k] == PH_PLASMA || pmp->PHC[k] == PH_SIMELT)
                && pmp->FIat)
            pmp->XetaA[k][0] = XetaS;
    }  // k
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of the surface potential pm[q].XpsiD[k] on diffuse
// layer plane on k-th sorption phase from total charge pmp->Xeta[k]
// ( in moles ) using Gouy-Chapman equation.
// Strictly valid at PSI < 30 mV. Modified by DAK 5 Jan 2000
// to add a Basic Stern EDL model.
//    Added 13.03.2008 by DK: returns int value showing (if not 0)
//    that some extreme values were reached for charge densities or
//    electric potentials (for detecting bad PIA), 0 otherwise
//
long int
TMulti::GouyChapman(  long int, long int, long int k )
{
    long int ist, status=0;
    double SigA=0., SigD=0., SigB=0., SigDDL=0.,
      XetaA[MST], XetaB[MST], XetaD[MST], f1, f3, A, Sig, F2RT, I, Cap;
    if( pmp->XF[k] < pmp->ScMinM ) // TProfil::pm->pa.p.ScMin )
        return status; // no sorbent

    // sorbent mass in grams
    pmp->YFk = pmp->FWGT[k];
    if( pmp->XF[k] < pmp->DSM )
       pmp->YFk = pmp->lowPosNum;

    for( ist=0; ist<pmp->FIat; ist++ )  // loop over surface types
    {
        double PsiD=0.0, PSIo=0.0, PsiA=0.0, PsiB=0.0; // Cleanup 05.12.2009 DK
        double ConvFactor = 1.;

        XetaA[ist] = XetaB[ist] = XetaD[ist] = 0.0;
        if( pmp->SCM[k][ist] == SC_NOT_USED || pmp->Nfsp[k][ist] < 1e-9  )
            continue;
        ConvFactor = F_CONSTANT / pmp->YFk / pmp->Aalp[k] / pmp->Nfsp[k][ist];
        // Calculation of charge densities (now limited to total charges > max. balance residual)
        if( fabs( pmp->XetaA[k][ist]) > pmp->DHBM ) // pmp->lowPosNum*100. )
            XetaA[ist] = pmp->XetaA[k][ist] * ConvFactor; // in C/m2
        if( fabs( pmp->XetaB[k][ist]) > pmp->DHBM ) // pmp->lowPosNum*100. )   // moles
            XetaB[ist] = pmp->XetaB[k][ist] * ConvFactor; // C/m2
        if( fabs( pmp->XetaD[k][ist]) > pmp->DHBM ) // pmp->lowPosNum*100. ) // moles
            XetaD[ist] = pmp->XetaD[k][ist] * ConvFactor; // C/m2

        // Limit maximum charge densities to prevent divergence
        if( fabs(XetaA[ist]) > 1.4 )
        {
// cout << "EDL charge density A " << XetaA[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaA[ist] = XetaA[ist] < 0.0 ? -1.4: 1.4;
            status = 60;
        }
        if( fabs(XetaB[ist]) > 2.0 )
        {
// cout << "EDL charge density B " << XetaB[ist] << " truncated to +- 1.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaB[ist] = XetaB[ist] < 0.0 ? -2.0: 2.0;
            status = 61;
        }
        if( fabs(XetaD[ist]) > 1.4 )
        {
// cout << "EDL charge density D " << XetaD[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaD[ist] = XetaD[ist] < 0.0 ? -1.4: 1.4;
            status = 62;
        }

        SigA = 0.;  SigB = 0.;   SigD = 0.;  SigDDL = 0.; // Cleanup 07.12.2009 DK
        pmp->XcapD[k][ist] = 0.0;  pmp->XpsiD[k][ist] = 0.0;
        if( fabs( XetaA[ist] ) < pmp->DHBM  && // pmp->lowPosNum*1e6 &&
               fabs( XetaB[ist] ) < pmp->DHBM && // pmp->lowPosNum*1e6 &&
               fabs( XetaD[ist] ) < pmp->DHBM ) // pmp->lowPosNum*1e6 )
            goto GEMU_CALC;  // skipping at near-zero charge
        // calculating charge density at diffuse layer
        switch( pmp->SCM[k][ist] )
        {
        case SC_CCM:  // Constant-Capacitance Model Schindler, extension Nilsson
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigDDL = -SigA - XetaB[ist];
            SigB = XetaB[ist];
            break;
        case SC_DDLM: // Generalized Double Layer Model [Dzombak and Morel, 1990]
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigDDL = -SigA;
            SigB = 0.0;
            break;
        case SC_TLM:  // Triple-Layer Model [Hayes and Leckie, 1987]
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_MTL:  // Modified Triple-Layer Model [Robertson, 1997]
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_BSM: // Basic Stern model: [Christl and Kretzschmar, 1999]
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_3LM: // Three-layer model: Hiemstra ea 1996; Tadanier & Eick 2002
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigD = XetaD[ist];
            SigDDL = -SigA - SigB -SigD;
            break;
        case SC_MXC:  // BSM for Ion-Exchange on permanent-charge surface
            SigA = pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_NNE:  // Non-Electrostatic
            break;
        default:
            continue;
        }
//        if( fabs( SigD ) > 1 )
//        {
//cout << "EDL charge density D " << SigD << " truncated to +- 1 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
//            SigD = SigD < 0.0 ? -1.0: 1.0;
//        }
        // Gouy-Chapman equation
        // parameters of diffuse layer using [Damaskin, 1987,p.192-195]
        A = 1e-9;
        F2RT = pmp->FRT / 2.;
        Sig = SigDDL; //  - XetaW[ist] ;
        I=pmp->IC;
        if( I > 1e-7 )
            // Aq solution density Ro included acc. to [Machesky ea., 1999]
            // Only for the basic Stern model for now
        {
            double Ro = 1.0;
            if( pmp->SCM[k][ist] == SC_BSM && pmp->FVOL[0] > 1e-16 )
                Ro = pmp->FWGT[0] / pmp->FVOL[0];
            A = sqrt( 2000. * 8.854e-12 * pmp->epsW[0] * pmp->RT * I * Ro );
        }
        Cap = F2RT * sqrt( 4.*A*A + Sig*Sig );

        // SD: workaround because of problems with log argument
        f3 =  sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A);
// std::cout<< f1  << ' ' << f3 << endl;
        if( f3 < 1 )
        {
            f1 = exp( -3. * F2RT );
            if( f3<f1) f3 = f1;
        }
        else
        {
            f1 = exp( 3. * F2RT );
            if( f3>f1 ) f3 = f1;
        }
        PSIo = log(f3)/F2RT;
//          PSIo = log( sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A) ) / F2RT;
//          Cap0 = fabs(Sig/PSIo);
//          Del = A*1e9/(2.*I*F)/cosh(PSIo*F2RT);
//          pmp->XdlD[k] = Del;
        pmp->XcapD[k][ist] = Cap;
        pmp->XpsiD[k][ist] = PSIo;
        PsiD = PSIo;
        // Truncating diffuse plane potential to avoid divergence
        if( fabs( PsiD ) > 0.4 )
        {
// cout << "All EDL models: PsiD = " << PsiD << " truncated to +- 0.4 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            PsiD = PsiD<0? -0.4: 0.4;
            status = 63;
        }
GEMU_CALC:
        // calculating potentials at EDL planes
        switch( pmp->SCM[k][ist] )
        {
        case SC_DDLM: // Generalized Double Layer Model [Dzombak & Morel 1990]
            pmp->XpsiA[k][ist] = PsiD;
            pmp->XpsiB[k][ist] = PsiD;
            break;
        case SC_CCM:  // Constant-Capacitance Model Schindler, ext. Nilsson
            if( pmp->XcapB[k][ist] > 0.001 )
            {  // Classic CCM Schindler with inner-sphere species only
               PsiA = SigA / pmp->XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 ) // truncated 0-plane potential
               {
                   PsiA = PsiA<0? -0.7: 0.7;
                   status = 64;
               }
               pmp->XpsiA[k][ist] = PsiA;
            }
            else { // Extended CCM model [Nilsson ea 1996] as TLM with PsiD = 0
               PsiB = - SigB / pmp->XcapB[k][ist];
               if( fabs( PsiB ) > 0.3 )  // truncated B-plane potential
               {
                   PsiB = PsiB<0? -0.3: 0.3;
                   status = 65;
               }
               PsiA = PsiB + SigA / pmp->XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 )
               {
                  PsiA = PsiA<0? -0.7: 0.7;
                  status = 66;
               }
               pmp->XpsiA[k][ist] = PsiA;
               pmp->XpsiB[k][ist] = PsiB;
            }
            break;
        case SC_MTL:  // Modified Triple Layer Model for X- Robertson | Kulik
// PsiD = 0.0; // test
            PsiB = PsiD - SigDDL / pmp->XcapB[k][ist];
            if( fabs( PsiB ) > 0.6)  // truncated B-plane potential
            {
// cout << "EDL (MTL) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 67;
            }
            PsiA = PsiB + SigA / pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0 plane potential
            {
// cout << "EDL (MTL) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 68;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_TLM:  // Triple-Layer Model   [Hayes 1987]
            PsiB = PsiD - SigDDL / pmp->XcapB[k][ist];
            if( fabs( PsiB ) > 0.6 )  // // truncated B-plane potential
            {
// cout << "EDL (TLM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 69;
            }
            PsiA = PsiB + SigA / pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0-plane potential
            {
// cout << "EDL (TLM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << k << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 70;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_3LM: // Three-Layer Model [Hiemstra & van Riemsdijk 1996]
//            PsiB = PsiD + SigD / pmp->XcapB[k][ist];
// cout << "EDL (3LM) PsiB(D) = " << PsiB << "  IT= " << pmp->IT << " k= "
// << k << " ist= " << ist << endl;
            PsiB = PsiD + ( SigA + SigB ) / pmp->XcapB[k][ist];  // Compare!
// cout << "EDL (3LM) PsiB(AB) = " << PsiB << "  IT= " << pmp->IT << " k= "
// << k << " ist= " << ist << endl;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
// cout << "EDL (3LM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 71;
            }
            PsiA = PsiB + SigA / pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )   // truncated 0-plane potential
            {
// cout << "EDL (3LM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 72;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_BSM: /* Basic Stern model, Christl & Kretzschmar, 1999 */
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
//cout << "EDL (BSM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 73;
            }
            PsiA = PsiB + SigA / pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0-plane potential
            {
//cout << "EDL (BSM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 74;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_MXC:  // BSM for permanent charge surfaces
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
// cout << "EDL (MXC) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 75;
            }
            PsiA = PsiB + SigA / pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 ) // truncated 0-plane potential
            {
// cout << "EDL (MXC) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 76;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_NNE:  // Non-Electrostatic
            pmp->XpsiA[k][ist] = 0.0;
            pmp->XpsiB[k][ist] = 0.0;
            break;
        default:
            break;
        }
    }  // ist   end of loop over surface types
    return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Calculation of new surface activity coefficient terms SACT (Kulik, 2004)
//
//  Revised by KD in April 2004 (PSI) to introduce new activity
//  coefficient terms SACT rigorously derived from Langmuir and QCA
//  isotherms (Kulik 2006, Radiochimica Acta).
//
//  SACT are placed into pmp->lnGam[j], as other activity coefficients except
//  relative surface potentials (Coulombic terms) kept separately.
//
//  Old (obsolete) SAT calculations (Kulik 2000, 2002) are retained.
//
//  pmp->lnSAC[*][3] vector is now used to keep original DUL[j] to restore
//  them after IPM-2 refinements for surface complexes.
//    Added 13.03.2008 by DK: returns int value showing (if true)
//    that some extreme values were obtained for some SACTs,
//    0 otherwise (for detecting bad PIA)
//
long int
TMulti::SurfaceActivityCoeff( long int jb, long int je, long int, long int, long int k )
{
	long int status = 0;
        long int i, ii, j, ja, ist=0, iss, dent, Cj, iSite[MST];
    double XS0,  xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, ISAT, XSs,
           SATst, xjn, q1, q2, aF, cN, eF, lnGamjo, lnDiff, lnFactor;
    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pmp->XF[k] <= pmp->DSM ) // No sorbent retained by the IPM - phase killed
        return status;
    if( pmp->XFA[k] <=  pa->p.ScMin )  // elimination of sorption phase
        return status;  // No surface species left

    for(i=0; i<MST; i++)
    {
        iSite[i] = -1;
        for( ii=0; ii<MST; ii++ )
           pmp->D[i][ii] = 0.0;  // cleaning the matrix for sites totals
    }
    // Extraction of site indices for the neutral >OH group
    for( j=jb; j<je; j++ )
    {
        ja = j - ( pmp->Ls - pmp->Lads );
        if( pmp->SATT[ja] == SAT_SOLV )
        {
           ist = pmp->SATX[ja][XL_ST]; // / MSPN;
           iSite[ist] = j;
        }
        // Counting current sites totals
        if( pmp->DCC[j] == DC_SUR_CARRIER ||
            pmp->DCC[j] == DC_SUR_MINAL ||
            pmp->DCC[j] == DC_PEL_CARRIER ||
            pmp->SATT[ja] == SAT_SOLV )
            continue;
        // Calculating ist (index of surface type)
        ist = pmp->SATX[ja][XL_ST];  // / MSPN;
        if( ist < 0 || ist >= MST )
            ist = 0;  // default: zero surface type
        // Calculate iss - index of site on surface type
        iss = pmp->SATX[ja][XL_SI];
        if( iss < 0 || iss >= MST )
            iss = 0;  // default: zero site is the weekest and the most abundant one
        pmp->D[iss][ist] += pmp->X[j]; // adding to total amount on a site
    }

    for( j=jb; j<je; j++ )
    { // Main loop for DCs - surface complexes
        lnGamjo = pmp->lnGmo[j];             // bugfix 16.03.2008 DK
        if( pmp->X[j] < min( pmp->DcMinM, pmp->lowPosNum ) )
            continue;  // This surface DC has been killed by IPM
//        OSAT = pmp->lnGmo[j]; // added 6.07.01 by KDA
        ja = j - ( pmp->Ls - pmp->Lads );
        rIEPS =  pa->p.IEPS;   // default 1e-3 (for old SAT - 1e-9)
//        dent = 1;  // default - monodentate
        switch( pmp->DCC[j] )  // code of species class
        {
        default: // pmp->lnGam[j] = 0.0;  not a surface species
            continue;
        case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
        case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
        case DC_SUR_GROUP: case DC_IEWC_B: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:
        case DC_IESC_A:
            // Calculating ist (index of surface type)
            ist = pmp->SATX[ja][XL_ST]; // / MSPN;
            // Calculating iss - index of site on surf.type
            iss = pmp->SATX[ja][XL_SI];
            if( iss < 0 || iss >= MST )
              iss = 0;  // zero site is the weekest and the most abundant one
            // Cj - index of carrier DC
            Cj = pmp->SATX[ja][XL_EM];
            if( Cj < 0 )
            {  // Assigned to the whole sorbent
                XVk = pmp->XFA[k];
                Mm = pmp->FWGT[k] / XVk;
            }
            else
            { // Assigned to one of the sorbent end-members
                XVk = pmp->X[Cj];
                if( XVk < pmp->ScMinM ) // pmp->DSM*0.1 )
                    continue;    // This end-member is zeroed off by IPM
                Mm = pmp->MM[Cj] * XVk/pmp->XFA[k];  // mol.mass
            }
            XSk = pmp->XFTS[k][ist];  // Total moles of sorbates on surface type
            XSs = pmp->D[iss][ist];   // Total moles of SC on site type
            xj = pmp->X[j];           // Current moles of this surface species

            // Extracting isotherm parameters
            if( pmp->MASDJ )
            {
               cN = pmp->MASDJ[ja][PI_P2];  // Frumkin/Pivovarov water coord. number
               if( cN > 0 )
                   dent = cN;   // denticity for L and QCA isoterms
               else dent = 1;   // Cleanup DK 07.12.2009
               aF = pmp->MASDJ[ja][PI_P1];  // Frumkin lateral interaction energy term
            //   bet = pmp->MASDJ[ja][PI_P3];   // BET beta parameter (reserved)
            }
            else {  // defaults
               cN = 0.0; aF = 0.0; dent = 1; // bet = 1.0;
            }
            switch( pmp->SATT[ja] ) // selection of the SACT model
            {
            case SAT_L_COMP: // Competitive monodentate Langmuir on a surface and site type
                XSkC = XSs / XVk / Mm * 1e6 /pmp->Nfsp[k][ist]/
                      pmp->Aalp[k]/1.66054;  // per nm2
                XS0 = (fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                        // max. density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;   // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )               // Setting limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( (pa->p.PC == 3 && !pmp->W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                pmp->lnGam[j] = ISAT;
                pmp->lnSAC[ja][0] = ISAT;
                break;
              // (Non)competitive QCA-L for 1 to 4 dentate species
//            case SAT_QCA4_NCOMP: dent++;     // code '4'
//            case SAT_QCA3_NCOMP: dent++;     // code '3'
//            case SAT_QCA2_NCOMP:             // code '2'
//            case SAT_QCA1_NCOMP:             // code '1'
            case SAT_QCA_NCOMP:  // dent++;   bidentate is default for QCA
                xj0 =
                 (fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pmp->Nfsp[k][ist] * 1e6     // xj
                     /pmp->Aalp[k]/1.66054; // Density per nm2 on site type iss
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/(double)dent)
                     xj = xj0/(double)dent - rIEPS;  // upper limit
//                ISAT = 0.0;
                q2 = xj0 - xj*dent;  // Computing differences in QCA gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
                pmp->lnGam[j] = ISAT;
                pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_FRUM_COMP: // Frumkin (FFG) isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pmp->Nfsp[k][ist] * 1e6  //  xj
                     /pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                q2 = xj0 - xj*dent;  // Computing differences in QCA gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
                // Calculation of the Frumkin exponential term
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {  // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * xj*dent / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
                }
                pmp->lnGam[j] = ISAT + eF;
                pmp->lnSAC[ja][0] = ISAT;
                pmp->lnSAC[ja][1] = eF;
                break;
            case SAT_FRUM_NCOMP: // (Non)competitive Frumkin isotherm
                                 // for permanent charge surfaces
                XSkC = xj / XVk / Mm / pmp->Nfsp[k][ist] * 1e6
                       / pmp->Aalp[k]/1.66054;  // per nm2
                XS0 = (pmp->MASDJ[ja][PI_DEN]/pmp->Aalp[k]/1.66054);
                         // max.dens.per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if(( pa->p.PC == 3 && !pmp->W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                // Calculation of the Frumkin exponential term (competitive)
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else { // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * XSkC / XS0;  // Fi = Fi'/(kT) Bockris p.938
                }
                pmp->lnGam[j] = ISAT + eF;
                pmp->lnSAC[ja][0] = ISAT;
                pmp->lnSAC[ja][1] = eF;
                break;
            case SAT_PIVO_NCOMP: // (Non)competitive Pivovarov isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pmp->Nfsp[k][ist] * 1e6
                     /pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                q2 = xj0 - xj*dent;  // Computing differences in gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
               // Calculation of the Frumkin exponential term
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {
                   double pivovar;
                   eF = cN * aF * xj / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
            // Calculation of the Pivovarov 98 exponential correction term
                   pivovar = xj0 / ( xj0 + xj * ( cN -1 ));
                   eF *= pivovar;
                }
                pmp->lnGam[j] = ISAT + eF;
                pmp->lnSAC[ja][0] = ISAT;
                pmp->lnSAC[ja][1] = eF;
                break;
            case SAT_VIR_NCOMP: // Non-Competitive virial isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pmp->Nfsp[k][ist] * 1e6
                     /pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                ISAT = 0.0;
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {   // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * xj / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
                }
                pmp->lnGam[j] = ISAT + eF;
                pmp->lnSAC[ja][0] = ISAT;
                pmp->lnSAC[ja][1] = eF;
                break;
            case SAT_BET_NCOMP: // Non-competitive BET for surface precipitation
                ISAT = 0.0;
// To be completed
//
//
                pmp->lnGam[j] = ISAT;
                pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_INDEF: // No SAT calculation whatsoever
                pmp->lnGam[j] = 0.0;
                pmp->lnSAC[ja][0] = 0;
                break;
            default:        // pmp->lnGam[j] = 0.0;
                break;
//  Obsolete old SAT calculations (retained from previous versions)
            case SAT_COMP: // Competitive SAT (obsolete) on a surface type
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pmp->X[iSite[ist]]; // neutral site does not compete!
                XS0 = pmp->MASDT[k][ist] * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist]; // expected total in moles
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                XSkC = XSk - xjn - xj; // occupied by the competing species;
                                     // this sorbate cannot compete to itself
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                xj0 = XS0 - XSkC;    // expected moles of this sorbate
                if(xj >= xj0)
                       xj = xj0 - rIEPS; // Limits
                if( xj * 2 <= xj0 )
                    ISAT = 0.0;      // ideal case
                else
                {
                   q1 = xj0 - xj;
                   q2 = rIEPS * XS0;
                   if( (pa->p.PC == 3 && !pmp->W1) || pa->p.PC != 3 )
                   {
                      if( q1 > q2 )
                        q2 = q1;
                   }
                   else {
                      q2 = q1;
                      if( q2 <= 1e-33 )
                         q2 = 1e-33;
                   }
                   ISAT = log( xj ) - log( q2 );
                }
                pmp->lnGam[j] = ISAT;
                pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_NCOMP: // Non-competitive truncated Langmuir SAT (obsolete)
                // rIEPS = pa->p.IEPS * 2;
                xj0 = fabs( pmp->MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist]; // in moles
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0;  // relative IEPS
                if(xj >= xj0)
                     xj = xj0 - rIEPS;  // upper limit
                if( xj * 2.0 <= xj0 )   // Linear adsorption - to improve !
                    ISAT = 0.0;
                else
                {
                    q1 = xj0 - xj;      // limits: rIEPS to 0.5*xj0
                    q2 = xj0 * rIEPS;
                    if( pa->p.PC == 3 && pmp->W1 )
                       ISAT = log( xj ) - log( q1 );
                    else {
                       if( q1 > q2 )
                          ISAT = log( xj ) - log( q1 );
                       else
                          ISAT = log( xj ) - log( q2 );
                    }
                 }
                pmp->lnGam[j] = ISAT;
                pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_SOLV:  // Neutral surface site (e.g. >O0.5H@ group)
                            // applies to the whole surface type!
                XSs = 0.0;  // calc total moles on all sites on surface type
                for( i=0; i<MST; i++ )
                   XSs += pmp->D[i][ist];
                XSkC = XSs / XVk / Mm * 1e6  // total non-solvent surf.species
                   /pmp->Nfsp[k][ist]/ pmp->Aalp[k]/1.66054;  // per nm2
                XS0 = (max( pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN] ));
                SATst = pa->p.DNS*1.66054*pmp->Aalp[k]/XS0;
                XS0 = XS0 / pmp->Aalp[k]/1.66054;
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( (pa->p.PC == 3 && !pmp->W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                pmp->lnGam[j] = ISAT + log( SATst );
                pmp->lnSAC[ja][0] = ISAT;
                pmp->lnSAC[ja][1] = log(SATst);
                break;
            }
        }

        if( lnGamjo > pmp->DHBM )
        {                               // Workaround DK 07.12.2009
            lnFactor = 0.2;
            lnDiff = pmp->lnGam[j] - lnGamjo;
            if( fabs( lnDiff ) > fabs( lnFactor ) ) // e times
            {
                if( fabs( lnDiff ) > 6.907755 )  // 1000 times
                    status = 101;   // the SACT has changed too much;
                            // the threshold pa->p.IEPS needs adjustment!
                // Smoothing (provisional)
                pmp->lnGam[j] = lnDiff > 0? lnGamjo + lnDiff - lnFactor:
                                lnGamjo + lnDiff + lnFactor;
            }
        }
    }  // j
   return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculating demo partial pressures of gases (works only in GEMS-PSI)
//
void TMulti::GasParcP()
{
#ifndef IPMGEMPLUGIN

        long int k,  i, jj=0;
    long int jb, je, j;

    if( !pmp->PG )
        return;

  char (*SMbuf)[MAXDCNAME] =
      (char (*)[MAXDCNAME])aObj[ o_w_tprn].Alloc( pmp->PG, 1, MAXDCNAME );
  pm.Fug = (double *)aObj[ o_wd_fug].Alloc( pm.PG, 1, D_ );
  pm.Fug_l = (double *)aObj[ o_wd_fugl].Alloc( pm.PG, 1, D_ );
  pm.Ppg_l = (double *)aObj[ o_wd_ppgl].Alloc( pm.PG, 1, D_ );

    for( k=0, je=0; k<pmp->FIs; k++ ) // phase
    {
        jb = je;
        je = jb+pmp->L1[k];
        if( pmp->PHC[k] == PH_GASMIX || pmp->PHC[k] == PH_PLASMA
           || pmp->PHC[k] == PH_FLUID )
        {
            for( j=jb; j<je; j++,jj++ )
            {  // fixed 02.03.98 DK

                copyValues(SMbuf[jj], pmp->SM[j], MAXDCNAME );
                pmp->Fug_l[jj] = -(pmp->G0[j]+pmp->GEX[j]);
                if( pmp->Pc > 1e-9 )
                    pmp->Fug_l[jj] += log(pmp->Pc);
                for( i=0; i<pmp->N; i++ )
                    pmp->Fug_l[jj] += *(pmp->A+j*pmp->N+i) * pmp->U[i];
                if( pmp->Fug_l[jj] > -37. && pmp->Fug_l[jj] < 16. )
                    pmp->Fug[jj] = exp( pmp->Fug_l[jj] );
                else  pmp->Fug[jj] = 0.0;
                // Partial pressure
                pmp->Ppg_l[jj] = pmp->Fug_l[jj] - pmp->lnGam[j];
                pmp->Fug_l[jj] *= .43429448;
                pmp->Ppg_l[jj] *= .43429448;
            }
            // break;
        }
    }
#endif
}

//--------------------- End of ipm_chemical2.cpp ---------------------------
