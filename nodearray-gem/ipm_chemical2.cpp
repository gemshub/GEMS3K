//-------------------------------------------------------------------
// $Id: ipm_gamma.cpp 690 2006-03-29 07:10:23Z gems $
//
// Copyright (C) 1992-2000  D.Kulik, S.Dmitrieva, K.Chudnenko, I.Karpov
//
// Implementation of chemistry-specific functions (concentrations,
// activity coefficients, adsorption models etc.)
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806); Kulik (2000), Geoch.Cosmoch.Acta 64,3161-3179
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//
#include <math.h>
#include "m_param.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// Calc  partial pressures of gases - to remake using dual potentials
void TMulti::GasParcP()
{
    short k,  i,  jj;
    int jb, je, j;

    if( !pmp->PG )
        return;
    for( k=0, je=0; k<pmp->FIs; k++ ) // phase
    {
        jb = je;
        je = jb+pmp->L1[k];
        if( pmp->PHC[k] == PH_GASMIX || pmp->PHC[k] == PH_PLASMA
           || pmp->PHC[k] == PH_FLUID )
        {
#ifndef IPMGEMPLUGIN
            char (*SMbuf)[MAXDCNAME] =
                (char (*)[MAXDCNAME])aObj[ o_w_tprn].Alloc( pmp->PG, 1, MAXDCNAME );
            memcpy(SMbuf,pmp->SM[jb], pmp->PG*MAXDCNAME*sizeof(char) );
            //     aObj[ o_w_tprn].SetPtr( pmp->SM[jb] );
            //     aObj[ o_w_tprn].SetDim( pmp->PG, 1 );
#endif
            for( j=jb,jj=0; j<je; j++,jj++ )
            {  /* fixed 02.03.98 DAK */
                pmp->Fug_l[jj] = -(float)(pmp->G0[j]+pmp->GEX[j]); /* + pmp->lnGam[j] */
                if( pmp->Pc > 1e-9 )
                    pmp->Fug_l[jj] += (float)log(pmp->Pc); /* fixed 02.03.98 DAK */
                for( i=0; i<pmp->N; i++ )
                    pmp->Fug_l[jj] += *(pmp->A+j*pmp->N+i) * (float)pmp->U[i];
                if( pmp->Fug_l[jj] > -37. && pmp->Fug_l[jj] < 16. )
                    pmp->Fug[jj] = exp( (double)pmp->Fug_l[jj] );
                else  pmp->Fug[jj] = /* F_EMPTY*/ 0.0;
                /* Partial pressure ! */
                pmp->Ppg_l[jj] = pmp->Fug_l[jj] - (float)pmp->lnGam[j]/* + log(pmp->Pc)*/;
                pmp->Fug_l[jj] *= (float).43429448;
                pmp->Ppg_l[jj] *= (float).43429448;
            }
            break;
        }
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of pH via activities of H2O and OH- */
/* Suggested by V.A.Sinitsyn, Apr 7, 1997 */
// Not use !!! SD
double TMulti::pH_via_hydroxyl( double x[], double Factor, int j)
{
    double lnaH;
    int jwa, jhy;
    jwa = j+1;
    jhy = j-1;                      /* Dangerous ! */
    lnaH = - pmp->G0[jhy] + 4.016534 + pmp->G0[jwa] + pmp->lnGam[jwa]
           - pmp->lnGam[jhy] - log( x[jhy]*Factor );
    x[j] = exp( lnaH - pmp->lnGam[j] ) / Factor;
    return (-lnaH *ln_to_lg);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* calculate bulk stoichiometry of a multicomponent phase
*/
void TMulti::phase_bcs( int N, int M, float *A, double X[], double BF[] )
{
    int i, j;
    double Xx;

    if( !A || !X || !BF )
        return;
    memset( BF, 0, N*sizeof( double ));
    for( j=0; j<M; j++ )
    {
        Xx = X[j];
        if( fabs( Xx ) < 1e-12 )
            continue;
        for( i=0; i<N; i++ )
            BF[i] += (double)A[i+j*N] * Xx;
    }
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* calculate total bulk stoichiometry of  solid phases
*/
void TMulti::phase_bfc( int k, int jj )
{
    int i, j;
    double Xx;

    if( pmp->PHC[k] == PH_AQUEL  || pmp->PHC[k] == PH_GASMIX ||
        pmp->PHC[k] == PH_FLUID  || pmp->PHC[k] == PH_PLASMA )
        return;
    for( j=0; j<pmp->L1[k]; j++ )
    {
        Xx = pmp->X[j+jj];
        if( fabs( Xx ) < 1e-12 )
            continue;
        for( i=0; i<pmp->N; i++ )
           pmp->BFC[i] += (double)pmp->A[i+(jj+j)*pmp->N] * Xx;
    }
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#define  a(j,i) ((double)(*(pmp->A+(i)+(j)*pmp->N)))
/* Calculation of dual chemical potentials and concentrations
* for species (indexed jb to je) in a k-th phase.
* Input arrays X, XF, XFA,  input factors: Factor, MMC
* This function has to be rewritten using new set of built-in
* chemical functions.
*/
void TMulti::ConCalcDC( double X[], double XF[], double XFA[],
              double Factor, double MMC, double Dsur, int jb, int je, int k)
{
    int j, ii;
    double Muj, DsurT, SPmol, lnFmol=4.016535;
    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pmp->PHC[0] == PH_AQUEL )
    {  /* mole fraction to molality conversion */
        if( !k ) lnFmol = log(1000./MMC);  // aq species
        else lnFmol = 4.016535; 	   // other species
    }

    for( j=jb; j<je; j++ )
    { /* cycle by DC - important bugfixes 02.04.2003 */
        Muj = DualChemPot( pmp->U, pmp->A+j*pmp->N, pmp->NR );
        pmp->Fx[j] = Muj * pmp->RT;     /* el-chem potential in J/mole */

        if( X[j] <= pmp->lowPosNum )
        { /* zeroing off */
            pmp->Wx[j] = pmp->lowPosNum; // 0.0;  debugging 29.11.05
            pmp->VL[j] = (float)log( pmp->lowPosNum );
            pmp->Y_w[j] = 0.0;
            pmp->lnGam[j] = 0.0;
            if( pmp->PHC[0] == PH_AQUEL )
               pmp->Y_m[j] = 0.0;
            switch( pmp->DCC[j] ) /* choice of expressions */
            {
               case DC_SCP_CONDEN:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] -pmp->GEX[j]);
                    break;
               case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES:
                    pmp->Y_la[j] = ln_to_lg*(Muj - pmp->G0[j] -pmp->GEX[j]
                                        + Dsur + lnFmol);
                    break;
               case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                    pmp->Y_la[j] = ln_to_lg* (Muj - pmp->G0[j] - pmp->GEX[j]
                                        + Dsur - 1. + 1. / ( 1.+Dsur ) );
                    break;
               case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:   /* gases */
               case DC_GAS_H2: case DC_GAS_N2:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]);
                    if( pmp->Pc > 1e-9 )
                        pmp->Y_la[j] += log10( pmp->Pc );
                    break;
               case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR:
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j] );
                    break;
               case DC_SUR_GROUP:
                    DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                                         + Dsur + DsurT/( 1.0+DsurT ) + lnFmol );
                    break;
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3:
               case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2:
               case DC_WSC_A3: case DC_WSC_A4: case DC_SUR_COMPLEX:
               case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                    DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                                         + Dsur + DsurT/( 1.0+DsurT ) + lnFmol );
                    break; // Coulombic term to be considered !!!!!!!!!!
               case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: /* sorbent */
                    DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
                    pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                       + Dsur - 1. + 1./(1.+Dsur) - DsurT + DsurT/(1+DsurT) );
                    break;
               default:
                    break; /* error in DC class code */
            }
            continue;
        }
        /* calculation of mole fraction */
        pmp->Wx[j] = X[j]/XF[k];
        if( pmp->Wx[j] > pmp->lowPosNum )
            pmp->VL[j] = (float)log( pmp->Wx[j] );
        else pmp->VL[j] = (float)log( pmp->lowPosNum );   // debugging 29.11.05 KD
        pmp->Y_la[j] = 0.0;
        switch( pmp->DCC[j] ) /* choice of expressions */
        {
        case DC_SCP_CONDEN:
            pmp->Wx[j] = 1;
            pmp->VL[j] = 0.0;
            if( pmp->LO )
                pmp->Y_m[j] = X[j]*Factor; /* molality */
            else pmp->Y_m[j] = 0.0;
            pmp->Y_w[j] = /* mass % of the system */
                1e2 * X[j] * pmp->MM[j] / pmp->MBX;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j]- pmp->GEX[j] );
            //  obsolete                pmp->Y_la[j] = ln_to_lg* pmp->lnGam[j];
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            break;
        case DC_AQ_ELECTRON: /* pE  SPmol = X[j]*Factor; */
            //                          pe_ = -ln_to_lg*( log( SPmol ) + pmp->lnGam[j]);
            pmp->Y_m[j] = 0.0;
            pmp->Y_la[j] = 0.0 - pmp->pe;  /* ln_a = -pe_ / .43429448; */
            pmp->Y_w[j] = 0.0;
            //                          pmp->IC_m[pmp->N-1] += SPmol;
            break;
        case DC_AQ_PROTON:
            pmp->pH = -ln_to_lg*(Muj-pmp->G0[j]-pmp->GEX[j]+Dsur+lnFmol);
            /* obsolete DAK 11.01.00    if( X[j] < 3e-8 )
                                           pmp->pH = pH_via_hydroxyl( X, Factor, j );
                                        else pmp->pH = -ln_to_lg*(log( X[j]*Factor ) +
                                                     pmp->lnGam[j] );    */
        case DC_AQ_SPECIES:
            SPmol = X[j]*Factor;  /* molality */
            pmp->IC +=  /* calc effective molal ionic strength */
                0.5* SPmol *(pmp->EZ[j]*pmp->EZ[j]);

//    pmp->FVOL[k] += pmp->Vol[j]*SPmol;  Error - found by B.Lothenbach 03.02.03
          pmp->FVOL[k] += pmp->Vol[j]*X[j]; // fixed 04.02.03 KD
            pmp->Y_m[j] = SPmol;
            pmp->Y_la[j] = ln_to_lg*(Muj - pmp->G0[j] -pmp->GEX[j]
                                        + Dsur + lnFmol); //    Variant: Without Dsur?
            // obsolete        pmp->Y_la[j] = ln_to_lg* (log( SPmol ) + pmp->lnGam[j] );
            pmp->Y_w[j] = 1e6 * X[j] * pmp->MM[j] / pmp->FWGT[k];
            for( ii=0; ii<pmp->NR; ii++ )
            {
                pmp->IC_m[ii] += SPmol* a(j,ii);
                pmp->IC_wm[ii] += X[j]* a(j,ii);  // moles of element in aq spec
            }
            break;
        case DC_AQ_SOLVENT: /* mol fractions in solvent */
        case DC_AQ_SOLVCOM:
            pmp->Y_m[j] = X[j]/XFA[k];
            pmp->Y_w[j] = 1e3*X[j]*pmp->MM[j]/pmp->FWGT[k];
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg* (Muj - pmp->G0[j] - pmp->GEX[j]
                                        + Dsur - 1. + 1. / ( 1.+Dsur ) );
            // Obsolete DAK 11.01.00   pmp->Y_la[j] = ln_to_lg* (log(pmp->Wx[j]) + pmp->lnGam[j] );
            break;
        case DC_GAS_COMP:
        case DC_GAS_H2O:
        case DC_GAS_CO2:   /* gases */
        case DC_GAS_H2:    /* volume */
        case DC_GAS_N2:
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j] );
            if( pmp->Pc > 1e-9 )
                pmp->Y_la[j] += log10( pmp->Pc );
            // Obsolete               pmp->Y_la[j] = ln_to_lg * (log(pmp->Wx[j]) + pmp->lnGam[j] );
            break;
        case DC_SOL_IDEAL:
        case DC_SOL_MINOR:   /* volume */
        case DC_SOL_MAJOR:
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j] );
            // Obsolete               pmp->Y_la[j] = ln_to_lg * (log(pmp->Wx[j]) + pmp->lnGam[j] );
            break;
            /* adsorption: Simplified by DAK 11.01.00 */
        case DC_SUR_GROUP:
            pmp->Y_m[j] = X[j]*Factor; /* molality */
            pmp->Y_w[j] =  /* mg/g sorbent */
                1e3 * X[j] * pmp->MM[j] / (MMC*XFA[k]);
            DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                                         + Dsur + DsurT/( 1.0+DsurT ) + lnFmol );
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
        case DC_WSC_A4:  /* case DC_SUR_GROUP: */
        case DC_SUR_COMPLEX:
        case DC_SUR_IPAIR:
        case DC_IESC_A:
        case DC_IEWC_B:
            pmp->Y_m[j] = X[j]*Factor; /* molality */
            pmp->Y_w[j] =  /* mg/g sorbent */
                1e3 * X[j] * pmp->MM[j] / (MMC*XFA[k]);
            DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                                         + Dsur + DsurT/( 1.0+DsurT ) + lnFmol );
            break;
        case DC_PEL_CARRIER:
        case DC_SUR_MINAL:
        case DC_SUR_CARRIER: /* sorbent */
            pmp->Y_m[j] = X[j]*Factor; /* molality */
            pmp->Y_w[j] = 0.0;
            if(pmp->FWGT[0]>pmp->lowPosNum)
              pmp->Y_w[j] = /* mg of sorbent per kg aq solution */
                1e6 * X[j] * pmp->MM[j] / pmp->FWGT[0];
            DsurT = MMC * (double)pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
            pmp->Y_la[j] = ln_to_lg * ( Muj - pmp->G0[j] - pmp->GEX[j]
                           + Dsur - 1. + 1./(1.+Dsur) - DsurT + DsurT/(1+DsurT) );
            // Obsolete         pmp->Y_la[j] = ln_to_lg * (log(pmp->Wx[j]) + pmp->lnGam[j] );
            pmp->FVOL[k] += pmp->Vol[j]*X[j];
            break;
        default:
            break; /* error in DC class code */
        }
        ;
    }   /* j */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of derived values (concentrations etc.) on IPM iteration
*  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
* This function has to be rewritten using new set of built-in
* chemical functions.
*/
void TMulti::ConCalc( double X[], double XF[], double XFA[])
{
    int k, ii;
    int i, j, ist/*=0*/, jj/*=0*/, jja/*=0*/;
    double Factor=0.0, Dsur=0.0, MMC=0.0;
    SPP_SETTING *pa = &TProfil::pm->pa;

   // Kostya: debug calculating x from dual solution
      if( pmp->Ls < 2 || !pmp->FIs )
        return;

    memset( pmp->BFC, 0, sizeof(double)*pmp->N );

    for( j=0; j<pmp->Ls; j++ )
    {
        pmp->Wx[j] = 0.;
        pmp->VL[j] = 0.;
    }
    j=0;
    pmp->VXc = 0.0;
    for( k=0; k<pmp->FI; k++ )
    { /* cycle by phases */
        i=j+pmp->L1[k];
        pmp->FWGT[k] = 0.0;
        pmp->FVOL[k] = 0.0;
        //   Dsur = 0.0;

      if( XF[k] > pmp->DSM &&
        !( pmp->PHC[k] == PH_SORPTION && XFA[k] <= pa->p.ScMin ))
       phase_bfc( k, j );

        if( k >= pmp->FIs || pmp->L1[k] == 1 )
        { /* this is a single- component phase */
            if( XF[k] < pmp->DSM )
            {
                //     pmp->Y_la[j] = 0.0;
                if( pmp->LO )
                    pmp->Y_m[j] = 0.0;
                pmp->Y_w[j] = 0.0;
                pmp->Fx[j] = DualChemPot( pmp->U, pmp->A+j*pmp->N, pmp->NR );
                pmp->Y_la[j] = ln_to_lg * ( pmp->Fx[j] - pmp->G0[j] -pmp->GEX[j]/* + Dsur */ );
                pmp->Fx[j] *= pmp->RT;     /* el-chem potential */
                goto NEXT_PHASE;
            }
            pmp->Wx[j] = 1.0;
            pmp->VL[j] = 0.0;
            if( pmp->LO && XFA[0] > 0 )
                pmp->Y_m[j] = X[j] * 1000./18.01528/XFA[0]; /* molality */
            pmp->Y_w[j] = /* mass % in the system */
                1e2 * X[j] * pmp->MM[j] / pmp->MBX;
            pmp->Fx[j] = DualChemPot( pmp->U, pmp->A+j*pmp->N, pmp->NR );
            pmp->Y_la[j] = ln_to_lg * ( pmp->Fx[j] - pmp->G0[j] - pmp->GEX[j] /* + Dsur */ );
            pmp->Fx[j] *= pmp->RT;     /* el-chem potential */
            //     pmp->Y_la[j] = ln_to_lg * pmp->lnGam[j];
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

        if( XF[k] <= pmp->DSM ||
     (pmp->PHC[k] == PH_AQUEL && ( XFA[k] <= pmp->lowPosNum*1e3 || XF[k] <= pa->p.XwMin ) )
                || ( pmp->PHC[k] == PH_SORPTION && XFA[k] <= pa->p.ScMin ))
        {
            memset( pmp->BF+k*pmp->N, 0, sizeof(double)*pmp->N );
            for(jj=j; jj<i; jj++)   // Loop added 10.03.01  KD (GTDEMO)
            {
                pmp->Wx[j] = 0.0;
                if( pmp->LO )
                    pmp->Y_m[jj] = 0.0;
                pmp->Y_w[jj] = 0.0;
                pmp->Fx[jj] = DualChemPot( pmp->U, pmp->A+jj*pmp->N, pmp->NR );
                pmp->Y_la[jj] = ln_to_lg * ( pmp->Fx[jj] - pmp->G0[jj] /* + Dsur */ );
                if(pmp->PHC[k] == PH_AQUEL || pmp->PHC[k] == PH_SORPTION )
                   pmp->Y_la[jj] += 1.74438;
                if(pmp->PHC[k] == PH_GASMIX || pmp->PHC[k] == PH_PLASMA )
                   pmp->Y_la[jj] += log10( pmp->Pc );
                pmp->Fx[jj] *= pmp->RT;     /* el-chem potential */
                pmp->lnGam[jj] = 0.0;
            }
            goto NEXT_PHASE;
        }
        /* calculate bulk stoichiometry of a multicomponent phase */
        phase_bcs( pmp->N, pmp->L1[k], pmp->A+j*pmp->N, pmp->X+j,
                   pmp->BF+k*pmp->N );

        switch( pmp->PHC[k] )
        {
        case PH_AQUEL:
            MMC = 0.0; /* mol mass of carrier */
            //                     Dsur = XF[k] - XFA[k];
            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
            if( XFA[k] > pmp->lowPosNum )
            {
                for(jj=j; jj<i; jj++)
                    if( pmp->DCC[jj] == DC_AQ_SOLVENT ||
                            pmp->DCC[jj] == DC_AQ_SOLVCOM )
                        MMC += pmp->MM[jj]*X[jj]/XFA[k];
            }
            else MMC=18.01528; // Assuming water-solvent
            if( (XFA[k] > pmp->lowPosNum) && (MMC > pmp->lowPosNum) )
                Factor = 1000./MMC/XFA[k]; /* molality */
            else Factor = 0.0;
            pmp->IC=0.; /* Factor = 0.5*55.508373/pmp->Yw; */
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
                if( X[jj] > pmp->lowPosNum)    // 1e-12 )
                    pmp->FWGT[k] += X[jj]*pmp->MM[jj];
            }
            break;
        case PH_POLYEL:
        case PH_SORPTION: /* only sorbent end-members! */
            pmp->YFk = XFA[k];
            MMC=0.0;

            for( ist=0; ist<pmp->FIat; ist++ )
                pmp->XFTS[k][ist] = 0.0;
            if( XFA[k] < pmp->lowPosNum ) XFA[k] = pmp->lowPosNum;
            for( jj=j; jj<i; jj++ )
            {
               jja = jj - ( pmp->Ls - pmp->Lads );
                if( pmp->DCC[jj] == DC_SUR_CARRIER ||
                        pmp->DCC[jj] == DC_SUR_MINAL ||
                        pmp->DCC[jj] == DC_PEL_CARRIER )
                {
                    MMC += pmp->MM[jj]*X[jj]/XFA[k];
                    /* Only sorbent mass */
                    pmp->FWGT[k] += X[jj]*pmp->MM[jj];
                }
                else
                {
/*!!!!!*/           ist = pmp->SATX[jja][XL_ST]; //  / MSPN;
                    pmp->XFTS[k][ist] += X[jj];
                }
            }
            pmp->logYFk = log(pmp->YFk);
            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
            if( Dsur <= -1.0 ) Dsur = -0.999999; // Debugging!!!!!
// double DsurTT=0.0;  //  Temporarily!
// DsurTT = XFA[k] * MMC * pmp->Aalp[k] * pa->p.DNS*1.66054e-6;
// Dsur = DsurTT / XF[k] - 1.0;
            break;
        default:

            return; /* Phase class code error! */
        }
        /* cycle for species concentrations in a phase */
        ConCalcDC( X, XF, XFA, Factor, MMC, Dsur, j, i, k );

NEXT_PHASE:
        pmp->VXc += pmp->FVOL[k];
        if( pmp->PHC[k] == PH_AQUEL && XF[k] > pa->p.XwMin && XFA[k] > pmp->lowPosNum*1e3 )
            for( ii=0; ii<pmp->NR; ii++ )
            {
               if( pmp->LO  ) // Sveta 30/08/01
               { if( pmp->IC_m[ii] >= pa->p.DB )
                    pmp->IC_lm[ii] = ln_to_lg*log( pmp->IC_m[ii] );
                else
                    pmp->IC_lm[ii] = 0;
                if( pmp->FWGT[k] >= pa->p.DB )
                    pmp->IC_wm[ii] *= (double)pmp->Awt[ii]*1000./pmp->FWGT[k];
                else
                    pmp->IC_wm[ii] = 0;
               }
            }
        j = i;
    }  /* k */

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Calculation of the surface potential pm[q].XpsiD[k] on diffuse
*  layer plane on k-th sorption phase from total charge pmp->Xeta[k]
* ( in moles ) using Gouy-Chapman equation.
* Strictly valid at PSI < 30 mV. Modified by DAK 5 Jan 2000
* to add a Basic Stern EDL model.
*/
void TMulti::GouyChapman(  int /*jb*/, int /*je*/, int k )
{
    int ist;
    double SigA=0., SigD=0., SigB=0., SigDDL=0.,
      XetaA[MST], XetaB[MST], XetaD[MST], f1, f3, A/*=1e-9*/, Sig, F2RT, I, Cap;
    /* Del, F=F_CONSTANT, Cap0;  double XetaW[MST]; */
    if( pmp->XF[k] <  TProfil::pm->pa.p.ScMin )
        return; /* no sorbent */

    /* sorbent mass in grams */
    pmp->YFk = pmp->FWGT[k];
    if(pmp->YFk < pmp->lowPosNum*100.)
       pmp->YFk = pmp->lowPosNum*100.;

    for( ist=0; ist<pmp->FIat; ist++ ) /* cycle by surface types */
    {
        double PsiD, PSIo, PsiA, PsiB;

        if( pmp->SCM[k][ist] == SC_NOT_USED || pmp->Nfsp[k][ist] < 1e-9  )
            continue;
        /* Calculation of charge densities */
        if( fabs( pmp->XetaA[k][ist]) > pmp->lowPosNum*100. )
            XetaA[ist] = pmp->XetaA[k][ist]*F_CONSTANT/pmp->YFk/
                (double)pmp->Aalp[k]/(double)pmp->Nfsp[k][ist]; /* C/m2 */
        else XetaA[ist] = 0.0;
        if( fabs( pmp->XetaB[k][ist]) > pmp->lowPosNum*100. )/* moles */
            XetaB[ist] = pmp->XetaB[k][ist] *F_CONSTANT/pmp->YFk/
                (double)pmp->Aalp[k]/(double)pmp->Nfsp[k][ist]; /* C/m2 */
        else XetaB[ist] = 0.0;
 if( fabs( pmp->XetaD[k][ist]) > pmp->lowPosNum*100. )/* moles */
     XetaD[ist] = pmp->XetaD[k][ist] *F_CONSTANT/pmp->YFk/
               (double)pmp->Aalp[k]/(double)pmp->Nfsp[k][ist]; /* C/m2 */
 else XetaD[ist] = 0.0;
// Trial variant
// if( fabs( pmp->XetaD[0][ist]) > pmp->lowPosNum*100. )/* moles */
//     XetaW[ist] = pmp->XetaD[0][ist] *F_CONSTANT/pmp->YFk/pmp->Aalp[k]
//                  /pmp->Nfsp[k][ist]; /* C/m2 */
// else XetaW[ist] = 0.0;
//
        /* Limit charge densities to 0.7 C/m2 */
        if( fabs(XetaA[ist]) > 1.4 )
        {
// cout << "EDL charge density A " << XetaA[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaA[ist] = XetaA[ist] < 0.0 ? -1.4: 1.4;
        }
        if( fabs(XetaB[ist]) > 2.0 )
        {
// cout << "EDL charge density B " << XetaB[ist] << " truncated to +- 1.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaB[ist] = XetaB[ist] < 0.0 ? -2.0: 2.0;
        }
        if( fabs(XetaD[ist]) > 1.4 )
        {
// cout << "EDL charge density D " << XetaD[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
            XetaD[ist] = XetaD[ist] < 0.0 ? -1.4: 1.4;
        }
        if( fabs( XetaA[ist] ) < pmp->lowPosNum*1e6 &&
               fabs( XetaB[ist] ) < pmp->lowPosNum*1e6 &&
               fabs( XetaD[ist] ) < pmp->lowPosNum*1e6 )
            goto GEMU_CALC;  /* skipping at near-zero charge */

        SigD = 0.;
        /* calculating charge density at diffuse layer */
        switch( pmp->SCM[k][ist] )
        {
        case SC_CCM:  /* Constant-Capacitance Model Schindler ext. Nilsson */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
//            SigD = -SigA;
            SigDDL = -SigA - XetaB[ist];
            SigB = XetaB[ist];
            break;
        case SC_DDLM: /* Diff. Double Layer Model Dzombak and Morel, 1990 */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
            SigDDL = -SigA;
            SigB = 0.0;
            break;
        case SC_TLM:  /* Triple-Layer Model Hayes and Leckie, 1987 */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_MTL:  /* Modified Triple-Layer Model Robertson, 1997 */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_BSM: /* Basic Stern model: Christl and Kretzschmar, 1999 */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
 case SC_3LM: /* Three-layer model: Hiemstra ea 1996; Tadanier & Eick 2002 */
      SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
      SigB = XetaB[ist];
      SigD = XetaD[ist];
      SigDDL = -SigA - SigB -SigD;
      break;
        case SC_MXC:  /* BSM for Ion-Exchange on perm.charge surface */
            SigA = (double)pmp->Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_NNE:  /* Non-Electrostatic Sorption */
            SigA = 0;
            SigB = 0;
            SigD = 0;
            SigDDL = 0;
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
        /* Gouy-Chapman equation */
        /* params of diffuse layer using Damaskin, 1987,p.192-195 */
        A = 1e-9;
        F2RT = pmp->FRT / 2.;
        Sig = SigDDL; //  - XetaW[ist] ;
        I=pmp->IC;
        if( I > 1e-7 )
            /* Aq solution density Ro included acc. to Machesky ea., 1999 */
            /* Only for the basic Stern model for now (Jan 5, 2000 DAK)! */
        {
            double Ro = 1.0;
            if( pmp->SCM[k][ist] == SC_BSM && pmp->FVOL[0] > 1e-16 )
                Ro = pmp->FWGT[0] / pmp->FVOL[0];
            A = sqrt( 2000. * 8.854e-12 * pmp->epsW * pmp->RT * I * Ro );
        }
        Cap = F2RT * sqrt( 4.*A*A + Sig*Sig );

        /* Sveta: workaround because of problems with log argument */
        f3 =  sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A);
        //std::cout<< f1  << ' ' << f3 << endl;
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
        pmp->XcapD[k][ist] = (float)Cap;
        pmp->XpsiD[k][ist] = PSIo;
        PsiD = PSIo; /* PsiA = PsiB = 0; */
        if( fabs( PsiD ) > 0.4 )
        {
// cout << "All EDL models: PsiD = " << PsiD << " truncated to +- 0.4 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiD = PsiD<0? -0.4: 0.4;
        }
GEMU_CALC:
        /* calculating potentials at EDL planes */
        switch( pmp->SCM[k][ist] )
        {
        case SC_DDLM: /* Diffuse Double Layer Model  Dzombak & Morel 1990 */
            pmp->XpsiA[k][ist] = PsiD;
            pmp->XpsiB[k][ist] = PsiD;
            break;
        case SC_CCM:  /* Constant-Capacitance Model Schindler, ext. Nilsson */
            if( pmp->XcapB[k][ist] > 0.001 )
            {  // Classic CCM Schindler with inner-sphere species only
               PsiA = SigA / (double)pmp->XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 )
               {
                   PsiA = PsiA<0? -0.7: 0.7;
               }
               pmp->XpsiA[k][ist] = PsiA;
            }
            else { // Extended CCM model Nilsson ea 1996 as TLM with PsiD = 0
               PsiB = - SigB / (double)pmp->XcapB[k][ist];
               if( fabs( PsiB ) > 0.3 )  /* Cutoff potential */
               {
                   PsiB = PsiB<0? -0.3: 0.3;
               }
               PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 )
               {
                  PsiA = PsiA<0? -0.7: 0.7;
               }
               pmp->XpsiA[k][ist] = PsiA;
               pmp->XpsiB[k][ist] = PsiB;
            }
            break;
        case SC_MTL:  /* Modified Triple Layer Model for X- Robertson | Kulik */
// PsiD = 0.0; // test
            PsiB = PsiD - SigDDL / (double)pmp->XcapB[k][ist];
            if( fabs( PsiB ) > 0.6)  /* Cutoff potential */
            {
// cout << "EDL (MTL) PsiB = " << PsiB << " truncated to +- 0.3 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
            }
            PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )
            {
// cout << "EDL (MTL) PsiA = " << PsiA << " truncated to +- 1.0 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_TLM:  /* Triple-Layer Model   Hayes 1987 */
            PsiB = PsiD - SigDDL / (double)pmp->XcapB[k][ist];
            if( fabs( PsiB ) > 0.6 )  /* Cutoff potential */
            {
// cout << "EDL (TLM) PsiB = " << PsiB << " truncated to +- 0.3 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
            }
            PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )
            {
// cout << "EDL (TLM) PsiA = " << PsiA << " truncated to +- 0.7 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_3LM:  // Three-Layer Model Hiemstra & van Riemsdijk 1996
//            PsiB = PsiD + SigD / pmp->XcapB[k][ist];
// cout << "EDL (3LM) PsiB(D) = " << PsiB << "  IT= " << pmp->IT << " k= "
// << k << " ist= " << ist << endl;
            PsiB = PsiD + ( SigA + SigB ) / (double)pmp->XcapB[k][ist];  // Compare!
// cout << "EDL (3LM) PsiB(AB) = " << PsiB << "  IT= " << pmp->IT << " k= "
// << k << " ist= " << ist << endl;
            if( fabs( PsiB ) > 0.6 )  /* Cutoff potential */
            {
// cout << "EDL (3LM) PsiB = " << PsiB << " truncated to +- 0.3 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
            }
            PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )
            {
// cout << "EDL (3LM) PsiA = " << PsiA << " truncated to +- 0.7 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_BSM: /* Basic Stern model, Christl & Kretzschmar, 1999 */
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  /* Cutoff potential */
            {
//cout << "EDL (BSM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
            }
            PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )
            {
//cout << "EDL (BSM) PsiA = " << PsiA << " truncated to +- 1.0 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_MXC:  /* BSM for const.charge surfaces */
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  /* Cutoff potential */
            {
// cout << "EDL (MXC) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
            }
            PsiA = PsiB + SigA / (double)pmp->XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )
            {
// cout << "EDL (MXC) PsiA = " << PsiA << " truncated to +- 1.0 V" <<
//      "  IT= " << pmp->IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
            }
            pmp->XpsiA[k][ist] = PsiA;
            pmp->XpsiB[k][ist] = PsiB;
            break;
        case SC_NNE:  /* Non-Electrostatic Sorption */
            pmp->XpsiA[k][ist] = 0.0;
            pmp->XpsiB[k][ist] = 0.0;
            break;
        default:
            break;
        }
    }  // ist    surf. types
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* With new calculation of surface activity coeff terms (2004)
*  Retains old (obsolete) SAT calculations (Kulik 2000, 2002)
*  Revised by KD in April 2004 (PSI) to introduce new activity
*  coefficient terms SACT rigorously derived from Langmuir and QCA
*  isotherms (Kulik 2004, to be submitted).
*  SACT are placed into pmp->lnGam[j], as any activity coeffs.
*  pmp->lnSAC[*][3] vector is now used to keep original DUL[j] to restore
*  them after IPM-2 refinements for surface complexes.
*/
void TMulti::SurfaceActivityCoeff( int jb, int je, int /*jpb*/, int /*jdb*/, int k )
{
    int i, ii, j, ja/*=0*/, ist/*=0*/, iss/*=0*/, dent/*=1*/, Cj/*=0*/, iSite[6];
    double XS0,  xj0, XVk, XSk, XSkC/*=0.*/, xj/*=0.*/, Mm/*=0.*/,
             rIEPS, ISAT, /*SAT,*/ XSs/*=0.*/,
             /*bet,*/ SATst/*=1.0*/, xjn/*=0.*/, q1, q2, aF/*=0.*/, cN/*=1.*/, eF/*=0.*/;
    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pmp->XF[k] <= pmp->DSM ) /* No sorbent retained by the IPM */
        return;
    if( pmp->XFA[k] <=  pa->p.ScMin )  // fixed 26.11.2006 by KD - elimination of sorption phase
        return;  /* No surface species left */

    for(i=0; i<MST; i++)
    {
        iSite[i] = -1;
        for( ii=0; ii<MST; ii++ )
           pmp->D[i][ii] = 0.0;  // cleaning the matrix for sites totals
    }
    /* Extraction of site indices for the neutral >OH group */
    for( j=jb; j<je; j++ )
    {
        ja = j - ( pmp->Ls - pmp->Lads );
        if( pmp->SATT[ja] == SAT_SOLV )
        {
/*!!!!!*/  ist = pmp->SATX[ja][XL_ST]; // / MSPN;
           iSite[ist] = j;
        }
        // Counting current sites totals
        if( pmp->DCC[j] == DC_SUR_CARRIER ||
            pmp->DCC[j] == DC_SUR_MINAL ||
            pmp->DCC[j] == DC_PEL_CARRIER ||
            pmp->SATT[ja] == SAT_SOLV )
            continue;
        // Calculate ist - index of surface type
/*!!!!!*/  ist = pmp->SATX[ja][XL_ST];  // / MSPN;
        if( ist < 0 || ist >= MST )
            ist = 0;  // default: zero surface type
        // Calculate iss - index of site on surf.type
        iss = pmp->SATX[ja][XL_SI];
        if( iss < 0 || iss >= MST )
            iss = 0;  // default: zero site is the weekest and the most abundant one
        pmp->D[iss][ist] += pmp->X[j]; // adding to total amount on a site
    }

    for( j=jb; j<je; j++ )
    { // Main loop for DC surface complexes
        if( pmp->X[j] <= pmp->lowPosNum /* *10. */ )
            continue;  /* This surface DC has been killed by IPM */
//        OSAT = pmp->lnGmo[j]; // added 6.07.01 by KDA
        ja = j - ( pmp->Ls - pmp->Lads );
        rIEPS =  pa->p.IEPS;   // between 1e-8 and 1e-10 default 1e-9
//        dent = 1;  // default - monodentate
        switch( pmp->DCC[j] )  /* code of species class */
        {
        default: /* pmp->lnGam[j] = 0.0; */
            continue;
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
        case DC_SUR_GROUP:
        case DC_IEWC_B:
        case DC_SUR_COMPLEX:
        case DC_SUR_IPAIR:
        case DC_IESC_A:
            /* Calculate ist - index of surface type */
/*!!!!!*/   ist = pmp->SATX[ja][XL_ST]; // / MSPN;
            /* Calculate iss - index of site on surf.type */
            iss = pmp->SATX[ja][XL_SI];
            if( iss < 0 || iss >= MST )
              iss = 0;  // zero site is the weekest and the most abundant one
            /* Cj - index of carrier DC */
            Cj = pmp->SATX[ja][XL_EM];
            if( Cj < 0 )
            {  /* Assigned to the whole sorbent */
                XVk = pmp->XFA[k];
                Mm = pmp->FWGT[k] / XVk;
            }
            else
            { /* Assigned to one of the sorbent end-members */
                XVk = pmp->X[Cj];
                if( XVk < pmp->DSM/10.0 )
                    continue; /* This end-member is zeroed off by IPM */
                Mm = pmp->MM[Cj] * XVk/pmp->XFA[k];  // mol.mass
            }
            XSk = pmp->XFTS[k][ist]; // Tot.moles of sorbates on surf.type
  XSs = pmp->D[iss][ist]; // tot.moles of SC on site type
            xj = pmp->X[j];  // Current moles of this surf.species

// isotherm parameters
            if( pmp->MASDJ )
            {
               cN = (double)pmp->MASDJ[ja][PI_P2];  // Frumkin/Pivovarov water coord. number
               dent = (int)cN;              // dentateness for L and QCA isoterms
               aF = (double)pmp->MASDJ[ja][PI_P1];  // Frumkin lateral interaction energy term
            //   bet = pmp->MASDJ[ja][PI_P3]; // BET beta parameter
            }
            else {
               cN = 0.0; aF = 0.0; /*bet = 1.0;*/ dent = 1;
            }
            switch( pmp->SATT[ja] ) // selection of the activity coeff. model
            {
            case SAT_L_COMP: // Competitive monodent. Langmuir on a surface and site type
                XSkC = XSs / XVk / Mm * 1e6 /(double)pmp->Nfsp[k][ist]/
                      (double)pmp->Aalp[k]/1.66054;  // per nm2        check this!!!!!!!
                XS0 = (double)(fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054); // max.dens.per nm2
/*            XS0 = pmp->MASDT[k][ist] * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist]; expected total sites in moles */
                if( pa->p.PC == 1 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( pa->p.PC == 2 && !pmp->W1 || pa->p.PC != 2 )
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
//            case SAT_QCA4_NCOMP: dent++;     /* code '4' */
//            case SAT_QCA3_NCOMP: dent++;     /* code '3' */
//            case SAT_QCA2_NCOMP:             /* code '2' */
            case SAT_QCA_NCOMP:  // dent++;   bidentate is default for QCA
//            case SAT_QCA1_NCOMP:             /* code '1' */
                xj0 =
                 (double)(fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
/*              xj0 = fabs(pmp->MASDJ[j]) * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist];     in moles */
                xj = XSs / XVk / Mm / (double)pmp->Nfsp[k][ist] * 1e6     // xj
                     /(double)pmp->Aalp[k]/1.66054; // Density per nm2 on site type iss
                if( pa->p.PC == 1 )
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
//                if( pa->p.PC == 2 && pmp->W1 )
//                      ISAT = log( xj ) - log( q1 );
                pmp->lnGam[j] = ISAT;
pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_FRUM_COMP:
                dent = 1; // monodentate for now
                xj0 =
                 (double)(fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / (double)pmp->Nfsp[k][ist] * 1e6  //  xj
                     /(double)pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC == 1 )
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
            case SAT_FRUM_NCOMP: /* (Non)competitive Frumkin isotherm for perm.charge surfaces */
                XSkC = xj / XVk / Mm / (double)pmp->Nfsp[k][ist] * 1e6
                       / (double)pmp->Aalp[k]/1.66054;  // per nm2
                XS0 = (double)(pmp->MASDJ[ja][PI_DEN]/pmp->Aalp[k]/1.66054);
                         /* max.dens.per nm2 */
/*            XS0 = pmp->MASDT[k][ist] * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist]; expected total sites in moles */
                if( pa->p.PC == 1 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( pa->p.PC == 2 && !pmp->W1 || pa->p.PC != 2 )
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
            case SAT_PIVO_NCOMP: /* (Non)competitive Pivovarov isotherm */
                dent = 1; // monodentate for now
                xj0 =
                 (double)(fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / (double)pmp->Nfsp[k][ist] * 1e6
                     /(double)pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC == 1 )
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
                else {
                   double pivovar/* = 1.0*/;
                   eF = cN * aF * xj / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
            // Calculation of the Pivovarov 98 exponential correction term
                   pivovar = xj0 / ( xj0 + xj * ( cN -1 ));
                   eF *= pivovar;
                }
                pmp->lnGam[j] = ISAT + eF;
pmp->lnSAC[ja][0] = ISAT;
pmp->lnSAC[ja][1] = eF;
                break;
            case SAT_VIR_NCOMP: /* Non-Competitive virial isotherm */
                dent = 1; // monodentate for now
                xj0 =
                 (double)(fabs(pmp->MASDJ[ja][PI_DEN])/pmp->Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / (double)pmp->Nfsp[k][ist] * 1e6  // xj / XVk ...
                     /(double)pmp->Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC == 1 )
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
            case SAT_BET_NCOMP: /* Non-competitive BET for surf.precipitation */
                ISAT = 0.0;
// To be completed
//
//
                pmp->lnGam[j] = ISAT;
pmp->lnSAC[ja][0] = ISAT;
                break;
            case SAT_INDEF: /* No SAT calculation */
                pmp->lnGam[j] = 0.0;
pmp->lnSAC[ja][0] = 0;
                break;
            default:        /* pmp->lnGam[j] = 0.0; */
                break;
//  Old obsolete SAT calculations
            case SAT_COMP: // Competitive SAT (old, obsolete) on a surface type
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pmp->X[iSite[ist]]; // neutral site does not compete!
                XS0 = (double)pmp->MASDT[k][ist] * XVk * Mm / 1e6
                      * (double)pmp->Nfsp[k][ist]; // expected total in moles
                if( pa->p.PC == 1 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                XSkC = XSk - xjn - xj; // occupied by the competing species;
                                     // this sorbate cannot compete to itself
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                xj0 = XS0 - XSkC;    // expected moles of this sorbate
                if(xj >= xj0)
                       xj = xj0 - rIEPS; // limits: 2 * rIEPS to XS0 - rIEPS
                if( xj * 2 <= xj0 )
                    ISAT = 0.0;      // ideal case
                else
                {
                   q1 = xj0 - xj;
                   q2 = rIEPS * XS0;
                   if( pa->p.PC == 2 && !pmp->W1 || pa->p.PC != 2 )
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
            case SAT_NCOMP: /* Non-competitive truncated L SAT (old, obsolete) */
// rIEPS = pa->p.IEPS * 2;
                xj0 = fabs( (double)pmp->MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * (double)pmp->Nfsp[k][ist]; /* in moles */
                if( pa->p.PC == 1 )
                    rIEPS = pa->p.IEPS * xj0;  // relative IEPS
                if(xj >= xj0)
                     xj = xj0 - rIEPS;  // upper limit
                if( xj * 2.0 <= xj0 )  // Linear adsorption - to improve !
                    ISAT = 0.0;
                else
                {
                    q1 = xj0 - xj;  // limits: rIEPS to 0.5*xj0
                    q2 = xj0 * rIEPS;
                    if( pa->p.PC == 2 && pmp->W1 )
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
                   /(double)pmp->Nfsp[k][ist]/ (double)pmp->Aalp[k]/1.66054;  // per nm2
                XS0 = (double)(max( pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN] ));
                SATst = pa->p.DNS*1.66054*(double)pmp->Aalp[k]/XS0;
                XS0 = XS0 / (double)pmp->Aalp[k]/1.66054;
                if( pa->p.PC == 1 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( pa->p.PC == 2 && !pmp->W1 || pa->p.PC != 2 )
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
/*                if( pmp->IT > 0 &&
                  ( SAT < 1e-2 || SAT > 1.0/rIEPS || fabs(log(SAT)+log(SATst)-OSAT) > 3 ) )
                {
cout << "IT=" << pmp->IT << ": NSite SAT for j= " << j << " old lnSAT= " << OSAT << " new lnSAT= " << ISAT << endl;
cout << "     x[jn]= " << pmp->X[j] << " XSk= " << XSk << " XSkC=" << XSkC << " XS0=" << XS0 << " SATst=" << SATst << endl;
//  ISAT = OSAT + 1;           // experiment
*/                break;
            }
        }
    }  /* j */
}

//--------------------- End of ipm_chemical2.cpp ---------------------------
