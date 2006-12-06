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
/* Calculation of max.moles of surface species for SAT stabilization
*  to improve IPM-2 convergence at high SAT values  KD 08.03.02
*  xj0 values are placed as upper kinetic constraints
*/
void TMulti::XmaxSAT_IPM2()
{
    int i, j, ja/*=0*/, k, jb, je=0, ist=0, Cj/*=0*/, iSite[6];
    double XS0/*=0.*/, xj0/*=0.*/, XVk, XSk/*=0.*/, XSkC/*=0.*/, xj,
           Mm, rIEPS, /*oDUL,*/ xjn/*=0.*/;

  if(!pmp->DUL )   // not possible to install upper kinetic constraints!
      return;

  for( k=0; k<pmp->FIs; k++ )
  { /* loop on phases */
     jb = je;
     je += pmp->L1[k];
     if( pmp->PHC[k] != PH_SORPTION )
          continue;

    if( pmp->XFA[k] < pmp->DSM ) // No sorbent retained by the IPM
        continue;
    if( pmp->XF[k]-pmp->XFA[k] < pmp->lowPosNum /* *10. */ )
        continue;  /* No surface species */

    for(i=0; i<6; i++)
        iSite[i] = -1;

    /* Extraction of site indices */
    for( j=jb; j<je; j++ )
    {
        ja = j - ( pmp->Ls - pmp->Lads );
        if( pmp->SATT[ja] != SAT_SOLV )
        {
            if( pmp->DCC[j] == DC_PEL_CARRIER || pmp->DCC[j] == DC_SUR_MINAL ||
                    pmp->DCC[j] == DC_SUR_CARRIER ) continue;
/*!!!!!!*/  ist = pmp->SATX[ja][XL_ST]; // / MSPN; MSPN = 2 - number of EDL planes
            continue;
        }
/*!!!!!!*/  ist = pmp->SATX[ja][XL_ST]; //  / MSPN;
        iSite[ist] = j;     // To be checked !!!
    }

    for( j=jb; j<je; j++ )
    { /* Loop for DC */
        if( pmp->X[j] <= pmp->lowPosNum /* *10. */ )
            continue;  /* This surface DC has been killed by IPM */
        rIEPS = TProfil::pm->pa.p.IEPS;
//        oDUL = pmp->DUL[j];
        ja = j - ( pmp->Ls - pmp->Lads );

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
/*!!!!!!*/  ist = pmp->SATX[ja][XL_ST]; // / MSPN;
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
            XSk = pmp->XFTS[k][ist]; /* Tot.moles of sorbates on surf.type */
            xj = pmp->X[j];  /* Current moles of this surf.species */
//            a=1.0;  Frumkin factor - reserved for extension to FFG isotherm
            switch( pmp->SATT[ja] )
            {
            case SAT_COMP: /* Competitive surface species on a surface type */
                /* a = fabs(pmp->MASDJ[j]); */
// rIEPS = TProfil::pm->pa.p.IEPS * 0.1;
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pmp->X[iSite[ist]]; // neutral site does not compete!
                XS0 = max(pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN]);
                XS0 = XS0 * XVk * Mm / 1e6 * pmp->Nfsp[k][ist]; /* expected total in moles */
// Experimental  rIEPS = TProfil::pm->pa.p.IEPS * XS0;
                XSkC = XSk - xjn - xj; /* occupied by the competing species */
                if( XSkC < 0.0 )
                    XSkC = rIEPS; // 0.0;   ??????????
//                if( XSkC > XS0 )
//                    XSkC = XS0 - 2.0 * rIEPS;
                xj0 = XS0 - XSkC;    /* expected moles of this sorbate */
                if( xj0 > pmp->lnSAC[ja][3] )
                    xj0 = pmp->lnSAC[ja][3];
                if( xj0 < rIEPS )
                   xj0 = rIEPS;  /* ensuring that it will not zero off */
                pmp->DUL[j] = xj0; // XS0*(1.0-TProfil::pm->pa.p.IEPS);  //pa.p.IEPS;
// pmp->DUL[j] = XS0 - rIEPS;
                // Compare with old DUL from previous iteration!
/*                if( pmp->W1 != 1 && pmp->IT > 0 && fabs( (pmp->DUL[j] - oDUL)/pmp->DUL[j] ) > 0.1 )
                {
cout << "XmaxSAT_IPM2 Comp. IT= " << pmp->IT << " j= " << j << " oDUL=" << oDUL << " DUL=" << pmp->DUL[j] << endl;
                }
*/                break;
            case SAT_L_COMP:
            case SAT_QCA_NCOMP:
            case SAT_QCA1_NCOMP:
            case SAT_QCA2_NCOMP:
            case SAT_QCA3_NCOMP:
            case SAT_QCA4_NCOMP:
            case SAT_BET_NCOMP:
            case SAT_FRUM_COMP:
            case SAT_FRUM_NCOMP:
            case SAT_PIVO_NCOMP:
            case SAT_VIR_NCOMP:
            case SAT_NCOMP: /* Non-competitive surface species */
// rIEPS = TProfil::pm->pa.p.IEPS * 2;
//                xj0 = fabs(pmp->MASDJ[j]) * XVk * Mm / 1e6
//                      * pmp->Nfsp[k][ist];  in moles
                 xj0 = fabs( (double)pmp->MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * (double)pmp->Nfsp[k][ist]; /* in moles */
                 pmp->DUL[j] = xj0 - rIEPS;
// Experimental  rIEPS = TProfil::pm->pa.p.IEPS * xj0;
//                 if( xj0 > 2.0*rIEPS )
//                   pmp->DUL[j] = xj0 - rIEPS; // xj0*(1.0-TProfil::pm->pa.p.IEPS); //pa.p.IEPS;
//                 else pmp->DUL[j] = rIEPS;
// Compare with old DUL from previous iteration!
/*                if( pmp->W1 != 1 && pmp->IT > 0 && fabs( (pmp->DUL[j] - oDUL)/pmp->DUL[j] ) > 0.1 )
                {
cout << "XmaxSAT_IPM2 Ncomp IT= " << pmp->IT << " j= " << j << " oDUL=" << oDUL << " DUL=" << pmp->DUL[j] << endl;
                }
*/                break;

            case SAT_SOLV:  /* Neutral surface site (e.g. >O0.5H@ group) */
rIEPS = TProfil::pm->pa.p.IEPS;
                XS0 = (double)(max( pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN] ));
                XS0 = XS0 * XVk * Mm / 1e6 * (double)pmp->Nfsp[k][ist]; /* in moles */

// Experimental   rIEPS = TProfil::pm->pa.p.IEPS * XS0;
                pmp->DUL[j] =  XS0 - rIEPS; // xj0*(1.0-TProfil::pm->pa.p.IEPS);  //pa.p.IEPS;
                if( pmp->DUL[j] <= rIEPS )
                   pmp->DUL[j] = rIEPS;
                break;
// New methods added by KD 13.04.04
            case SAT_INDEF: /* No SAT calculation */
            default:        /* pmp->lnGam[j] = 0.0; */
                break;
            }
        }
     }  /* j */
  } /* k */
}

// clearing pmp->DUL constraints!
void TMulti::XmaxSAT_IPM2_reset()
{
    int j, ja/*=0*/, k, jb, je=0;

  if(!pmp->DUL )   // no upper kinetic constraints!
      return;

  for( k=0; k<pmp->FIs; k++ )
  { /* loop on phases */
     jb = je;
     je += pmp->L1[k];
     if( pmp->PHC[k] != PH_SORPTION )
          continue;

    for( j=jb; j<je; j++ )
    { /* Loop for DC */
      ja = j - ( pmp->Ls - pmp->Lads );
      if( pmp->lnSAC && ja >= 0 && ja < pmp->Lads )
          pmp->DUL[j] = pmp->lnSAC[ja][3];  // temp. storing initial DUL constr.
    }  /* j */
  } /* k */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// calc value of dual chemical potencial
double TMulti::DualChemPot( double U[], float AL[], int N )
{
    double Nu = 0.0;
    for(int i=0; i<N; i++ )
        Nu += AL[i]? U[i]*(double)(AL[i]): 0.0;
    return Nu;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  This procedure sets kinetic limits according to a given
//  concentration units
//  Needs much more work and elaboration!
//
void TMulti::Set_DC_limits( int Mode )
{
    double XFL, XFU, XFS=0., XFM, MWXW, MXV, XL, XU;
    int jb, je, j,k, MpL;
    vstr tbuf(80);

    if( !pmp->PLIM )
        return;  /* no limits */
// ???????????????????????????????????????
    ConCalc( pmp->X, pmp->XF, pmp->XFA );

    for(k=0; k<pmp->FI; k++)
        XFS+=pmp->XF[k];  /* calc sum of mol all phases */

    jb=0;
    for( k=0; k<pmp->FI; k++ )
    { /* cycle by phases */
        je=jb+pmp->L1[k];
        /*XFM=0.;*/
        MWXW =0.;
        MXV = 0.;
        XFL = 0.;
        XFU = 1e6;
        if( Mode && pmp->XF[k] < pmp->DSM )
            goto NEXT_PHASE;
        XFM = pmp->FWGT[k]; /* Mass of a phase */
        if( Mode )
        {
            MWXW = XFM/pmp->XF[k]; /* current mol weight of phase */
            MXV = pmp->FVOL[k]/pmp->XF[k]; /*current mol volume of phase  */
        }
        /* Check codes for phase DC */
        MpL=0;
        for( j=jb; j<je; j++ )
            if( pmp->RLC[j] != NO_LIM )
                MpL = 1;
        if( pmp->RFLC[k] == NO_LIM && !MpL )
        { /* check type restrictions on phase */
            goto NEXT_PHASE;
        }
        switch( pmp->RFSC[k] )
        { /* check scale restrictions on phase in all system */
        case QUAN_MOL:
            XFL = Mode? pmp->XF[k]: pmp->PLL[k];
            XFU = Mode? pmp->XF[k]: pmp->PUL[k];
            break;
        case CON_MOLAL:
            XFL = Mode? pmp->XF[k] /* *pmp->X[pmp->LO]/55.508373 */:
                  pmp->PLL[k]*pmp->GWAT/55.508373;
            XFU = Mode? pmp->XF[k] /* *pmp->X[pmp->LO]/55.508373 */:
                  pmp->PUL[k]*pmp->GWAT/55.508373;
            break;
        case CON_MOLFR:
            XFL = Mode? pmp->XF[k]: pmp->PLL[k]*XFS;
            XFU = Mode? pmp->XF[k]: pmp->PUL[k]*XFS;
            break;
        case CON_WTFR:   if(MWXW < 1e-15) break;  // Temp.fix
            XFL = Mode? pmp->XF[k]: pmp->PLL[k]*pmp->MBX/MWXW;
            XFU = Mode? pmp->XF[k]: pmp->PUL[k]*pmp->MBX/MWXW;
            break;
        case CON_VOLFR:   if(MXV < 1e-15) break; // Temp.fix
            XFL = Mode? pmp->XF[k]: pmp->PLL[k]*pmp->VXc/MXV;
            XFU = Mode? pmp->XF[k]: pmp->PUL[k]*pmp->VXc/MXV;
            break;
        default:
            ; /* do more? */
        }
//        if( pmp->RFLC[k] == NO_LIM )
//        {                            Temporary!
            XFL = 0.0;
            XFU = 1e6;
//        }
        for( j=jb; j<je; j++ )
        { /* cycle by DC. */
            if( pmp->RLC[j] == NO_LIM )
                continue;

            if( Mode )
            {
                XU = pmp->DUL[j];
                XL = pmp->DLL[j];
            }
            else
                switch( pmp->RSC[j] ) /* get initial limits */
                {
                case QUAN_MOL:
                    XU = pmp->DUL[j];
                    XL = pmp->DLL[j];
                    break;
                case CON_MOLAL:
                    XU = pmp->DUL[j]*pmp->GWAT/55.508373;
                    XL = pmp->DLL[j]*pmp->GWAT/55.508373;
                    break;
                case CON_MOLFR:
                    XU = pmp->DUL[j]*XFU;
                    XL = pmp->DLL[j]*XFL;
                    break;
                case CON_WTFR:
//Ask Dima!!! 20/04/2002
#ifndef IPMGEMPLUGIN
                    XU = pmp->DUL[j]*XFU*MWXW /
         TProfil::pm->MolWeight(pmp->N, pmp->Awt, pmp->A+j*pmp->N );
                    XL = pmp->DLL[j]*XFL*MWXW /
         TProfil::pm->MolWeight(pmp->N, pmp->Awt, pmp->A+j*pmp->N );

#endif
                    break;
                case CON_VOLFR:
                    XU = pmp->DUL[j]*XFU*MXV/ pmp->Vol[j];
                    XL = pmp->DLL[j]*XFL*MXV/ pmp->Vol[j];
                    break;
                default:
                    ; /* do more */
                }
            /* check combine */
            if( XU < 0.0 ) XU = 0.0;
            if( XU > 1e6 ) XU = 1e6;
            if( XL < 0.0 ) XL = 0.0;
            if( XL > 1e6 ) XL = 1e6;
            if( XU > XFU )
            {
 //               JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent upper metastability limits j=%d k=%d XU=%g XFU=%g",
                         j, k, XU, XFU );
                Error( "E11IPM Set_DC_limits(): ",tbuf.p );
//                XU = XFU; // - pmp->lowPosNum;
            }
            if( XL < XFL )
            {
//                JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent lower metastability limits j=%d k=%d XL=%g XFL=%g",
                         j, k, XL, XFL );
                Error( "E12IPM Set_DC_limits(): ",tbuf.p );
//                XL = XFL; // - pmp->lowPosNum;
            }
            pmp->DUL[j]=XU;
            pmp->DLL[j]=XL;
        }   /* j */
NEXT_PHASE:
        jb = je;

    }  /* k */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// Calc sum count of phases
void TMulti::TotalPhases( double X[], double XF[], double XFA[] )
{
    int jj, j, i, k;
    double XFw, XFs, x;

    j=0;
    for( k=0; k< pmp->FI; k++ )
    { /* cycle by phases */
        i=j+pmp->L1[k];
        XFw = 0.0;
        XFs=0.0; /* calc  count mol of carrier */
        for(jj=j; jj<i; jj++)
        {
            x = X[jj];
            if( pmp->DCCW[jj] == DC_ASYM_CARRIER && pmp->FIs )
                XFw += x;
            else XFs += x;
        }
        XF[k] = XFw + XFs;
        if( k < pmp->FIs )
            XFA[k] = XFw;
        j=i;
    }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Corrections to prime chemical potentials F0[j]
*  of j-th species in k-th phase among IPM main iterations
*  Returns double - value of corrected chem. potential.
*  If error, returns +7777777 J/mole.
*  Last modif. 05 Jan 2000 by DAK to include BSM EDL model.
*/
double TMulti::Ej_init_calc( double, int j, int k)
{
    int ja=0, ist/*=0*/, isp/*=0*/, jc=-1;
    double F0=0.0, Fold, dF0, Mk=0.0, Ez, psiA, psiB, CD0, CDb, ObS;
    SPP_SETTING *pa = &TProfil::pm->pa;

    Fold = pmp->F0[j];
    if( pmp->FIat > 0 && j < pmp->Ls && j >= pmp->Ls - pmp->Lads )
    {
        ja = j - ( pmp->Ls - pmp->Lads );
        jc = pmp->SATX[ja][XL_EM];
    }
    if( k < pmp->FIs && pmp->XFA[k] > 1e-12)
    {
           if( jc < 0 ) /* phase (carrier) molar mass g/mkmol */
              Mk = pmp->FWGT[k]/pmp->XFA[k]*1e-6;
           else Mk = pmp->MM[jc]*(pmp->X[jc]/pmp->XFA[k])*1e-6;
        /* DC carrier molar mass g/mkmol */
    }
    switch( pmp->DCC[j] )
    { /* analyse species class code */
    case DC_SCP_CONDEN:
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:
    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM:
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR:
        F0 = pmp->lnGmM[j];
        break;
        /*adsorption */
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
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
        F0 = pmp->lnGmM[j]; /* + pmp->lnGam[j]; */
        /* get ist - index of surface type and isp - index of surface plane  */
/*!!!!!*/  ist = pmp->SATX[ja][XL_ST];  // / MSPN;
/*!!!!!*/  isp = pmp->SATX[ja][XL_SP]; // % MSPN;
        CD0 = (double)pmp->MASDJ[ja][PI_CD0];  // species charge that goes into 0 plane
        CDb = (double)pmp->MASDJ[ja][PI_CDB];  // species charge that goes into B plane
        ObS = (double)pmp->MASDJ[ja][PI_DEN];  // obsolete - the sign for outer-sphere charge
        if( ObS >= 0.0 )
            ObS = 1.0;
        else ObS = -1.0;
        psiA = pmp->XpsiA[k][ist];
        psiB = pmp->XpsiB[k][ist];
        Ez = double(pmp->EZ[j]);
        if( !isp )  // This is the 0 (A) plane species
        {
            if( fabs( CD0 ) > 1e-20 )  // Doubtful...
                Ez = CD0;
            F0 += psiA * Ez * pmp->FRT;
        }
        else  /* This is B plane */
        {
            if( pmp->SCM[k][ist] == SC_MTL || pmp->SCM[k][ist] == SC_MXC )
            { /* Modified TL: Robertson, 1997; also XTLM Kulik 2002 */
                  if( fabs( CDb ) > 1e-20 )  // Doubtful...
                      Ez = CDb;
                  F0 += psiB * Ez * pmp->FRT;
            }
            if( pmp->SCM[k][ist] == SC_TLM )
            {
// New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* pmp->FRT;
// see also Table 4 in Zachara & Westall, 1999
// Old version:  TLM Hayes & Leckie, 1987 uses the sign indicator at density
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* pmp->FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* pmp->FRT;
                  }
               }
            }
            else if( pmp->SCM[k][ist] == SC_BSM || pmp->SCM[k][ist] == SC_CCM )
            { /* Basic Stern model Christl & Kretzschmar, 1999 */
// New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* pmp->FRT;
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* pmp->FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* pmp->FRT;
                  }
               }
            }
        }
        if( Mk > 1e-9 )
        {
            if( pmp->SCM[k][ist] == SC_MXC || pmp->SCM[k][ist] == SC_NNE ||
                    pmp->SCM[k][ist] == SC_IEV )
                F0 -= log( Mk * (double)(pmp->Nfsp[k][ist]) *
                   (double)(pmp->Aalp[k]) * pa->p.DNS*1.66054 );
            else F0 -= log( Mk * (double)(pmp->Nfsp[k][ist]) *
                  (double)(pmp->Aalp[k]) * pa->p.DNS*1.66054 );
            F0 -= (double)(pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 /
                  ( 1.0 + (double)(pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 );
        }
        break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:  // const charge of carrier - not yet done !!!!!!!!!!
    case DC_SUR_CARRIER:
        F0 -= (double)(pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 /
              ( 1.0 + (double)(pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 );
        F0 += (double)(pmp->Aalp[k])*Mk*pa->p.DNS*1.66054;
        break;
    }
    F0 += pmp->lnGam[j];

    if( k >= pmp->FIs )   //Sveta ?????????????
        return F0;
    // Smoothing procedure for highly non-ideal phases
    if( pmp->sMod[k][SGM_MODE] != SM_IDEAL
            || pmp->sMod[k][SCM_TYPE] != SC_NNE )
    {
        dF0 = F0 - Fold;
        if( pmp->X[j]>pmp->lowPosNum && fabs( dF0 ) >= 1e-5 ) /* to be checked !!! */
            F0 = Fold + dF0 * pmp->FitVar[3];
    }  // FitVar[3] = TinkleSuppressFactor(); see GammaCalc()
    return F0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of prime chemical potential F (returned)
* from moles of DC Y[], total moles of phase YF[] and standard
* Gibbs energy gT (obtained from pmp->G[])
* On error returns +7777777.
*/
double TMulti::PrimeChemPot(
    double G,      /* gT0+gEx */
    double logY,   /* ln x */
    double logYF,  /* ln Xa */
    double asTail, /* asymmetry non-log term or 0 for symmetric phases */
    double logYw,  /* ln Xv */
    char DCCW      /* generalized species class code */
)
{
    double F;
    switch( DCCW )
    {
    case DC_SINGLE:
        F = G;
        break;
    case DC_ASYM_SPECIES:
        F = G + logY - logYw + asTail;
        break;
    case DC_ASYM_CARRIER:
        F = G + logY - logYF + asTail + 1.0 -
            1.0/(1.0 - asTail);
        break;
    case DC_SYMMETRIC:
        F = G + logY - logYF;
        break;
    default:
        F = 7777777.;
    }
    return F;
}


// Kernel functions of IPM - rewritten by DAK for adsorption

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// VJ - Update of prime chemical potentials
void TMulti::PrimeChemicalPotentials( double F[], double Y[], double YF[], double YFA[] )
{
    int i,j,k;
    double v, Yf; // v is debug variable

    j=0;
    for( k=0; k<pmp->FI; k++ )
    { /* cycle by phase */
        i=j+pmp->L1[k];
        if( YF[k] <= pmp->lowPosNum*100. || ( pmp->PHC[k] == PH_AQUEL &&
        ( YF[k] <= TProfil::pm->pa.p.XwMin || Y[pmp->LO] <= pmp->lowPosNum*1e3 )))
            goto NEXT_PHASE;

        pmp->YFk = 0.0;
        Yf= YF[k]; /* calculate number of moles of carrier */
        if( pmp->FIs && k<pmp->FIs )
            pmp->YFk = YFA[k];
        if( Yf >= 1e6 )
        {                 // error - will result in zerodivide!
           gstring pbuf(pmp->SF[k],0,20);
           char buf[200];
           sprintf( buf, "Broken IPM solution: Phase %s  Yf= %lg", pbuf.c_str(), Yf );
           Error( "E13IPM PrimeChemicalPotentials():", buf);
//           Yf = pmp->YFk;
        }
        if( pmp->YFk > pmp->lowPosNum*10. )
        {
            pmp->logXw = log(pmp->YFk);
            pmp->aqsTail = 1.- pmp->YFk / Yf;
        }
        if( pmp->L1[k] > 1 )
        {
            pmp->logYFk = log( Yf );
        }
        if( pmp->PHC[k] == PH_AQUEL)
            /* ln moles of solvent in aqueous phase */
            pmp->Yw = pmp->YFk;
        for( ; j<i; j++ )
        { /*cycle by DC. */
            if( Y[j] < pmp->lowPosNum )
                continue;  /* exception by minimum DC quantity */
            /* calculate chemical potential of j-th DC */
            v = PrimeChemPot( pmp->G[j], log(Y[j]), pmp->logYFk,
                              pmp->aqsTail, pmp->logXw, pmp->DCCW[j] );
            F[j] = v;
        }   /* j */
NEXT_PHASE:
        j = i;
    }  /* k */
    if( pmp->Yw > pmp->lowPosNum*1e3 )
        pmp->logXw = log(pmp->Yw);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of Karpov stability criteria for a DC*/
double TMulti::KarpovCriterionDC(
    double *dNuG,  /* Nu[j]-c[j] difference - is modified here */
    double logYF,  /* ln Xa */
    double asTail, /* asymmetry correction (0 for symmetric phases) */
    double logYw,  /* ln Xv */
    double Wx,     /* mole fraction */
    char DCCW      /* generalized DC class code */
)
{
    double Fj;     /* output phase stability criterion */

    if( logYF > -35. && Wx > 1e-18 )
        switch( DCCW ) /* expressions for fj */
        {
        default: /* error code !!! */
            *dNuG = 36.;
        case DC_SINGLE:
            Wx = 1.0;
        case DC_SYMMETRIC:
            break;
        case DC_ASYM_SPECIES:
            *dNuG += logYw - logYF - asTail;
            break;
        case DC_ASYM_CARRIER:
            *dNuG += 1.0/(1.0 - asTail) - asTail - 1.0;
        }
    if( fabs( *dNuG ) > 35.)
        Fj = ( *dNuG > 0 )? 1.5860135e15: 6.305117e-16;
    else Fj = exp( *dNuG );
    Fj -= Wx;
    return Fj;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// calculation of Karpov stability criteria for all phases
void TMulti::f_alpha()
{
    short k;
    int j, ii;
    double *EMU,*NMU, YF, Nu, dNuG, Wx, Yj, Fj;

    EMU = pmp->EMU;
    NMU = pmp->NMU;
    memset( EMU, 0, pmp->L*sizeof(double));
    memset( NMU, 0, pmp->L*sizeof(double));

    j=0;
    pmp->YMET = 0.0;
    for( k=0; k<pmp->FI; k++ )
    { /*phases */
        ii=j+pmp->L1[k];
        pmp->Falp[k] = pmp->YMET; /* metastability parameter */
        pmp->logXw = -35.;
        pmp->logYFk = -35.;

        pmp->YFk = 0.0;
        YF= pmp->YF[k]; /* moles of carrier */
        if( pmp->FIs && k<pmp->FIs )
            pmp->YFk = pmp->YFA[k];
        if( pmp->YFk > 6.305117e-16 )
        {
            pmp->logXw = log(pmp->YFk);
            pmp->aqsTail = 1.- pmp->YFk / YF;
        }
        else
        {
            pmp->logXw = -35.;
            pmp->aqsTail = 0.0;
        }

        if( pmp->L1[k] > 1 && YF > 6.305117e-16 )
            pmp->logYFk = log( YF );
        else pmp->logYFk = -35.;
        if( pmp->PHC[k] == PH_AQUEL) /* number of moles of solvent */
            pmp->Yw = pmp->YFk;

        for( ; j<ii; j++ )
        { /* DC */
            Nu = DualChemPot( pmp->U, pmp->A+j*pmp->N, pmp->NR );
            dNuG = Nu - pmp->G[j]; /* this is -s_j (6pot paper 1) */
            Wx = 0.0;
            Yj = pmp->Y[j];
            if( YF > 1e-12 && Yj > pmp->lowPosNum*10. )
                Wx = Yj / YF; /*calc mol fraction of DC */
            /* calc Karpov criteria of DC */
            Fj = KarpovCriterionDC( &dNuG, pmp->logYFk, pmp->aqsTail,
                                    pmp->logXw, Wx, pmp->DCCW[j] );
            NMU[j] = dNuG;
            EMU[j] = Fj;
            if( YF > 1e-12 && Yj > pmp->lowPosNum*10. )
                pmp->Falp[k] += EMU[j]; /* calc Karpov criteria of phase */
        }   /* j */
        j = ii;
    }  /* k */
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of a species increment to total Gibbs energy G(X)
*  of the system (returned). On error returns +7777777.
*/
double TMulti::FreeEnergyIncr(
    double G,      /* gT0+gEx */
    double x,      /* x - moles of species */
    double logXF,  /* ln Xa - moles of phase */
    double logXw,  /* ln Xv - moles of the solvent/sorbent */
    char DCCW      /* generalized species class code */
)
{
    double Gi;

    switch( DCCW )
    {
    case DC_ASYM_SPECIES:
        Gi = x * ( G + log(x) - logXw );
        break;
    case DC_ASYM_CARRIER:
    case DC_SYMMETRIC:
        Gi = x * ( G + log(x) - logXF );
        break;
    case DC_SINGLE:
        Gi = G * x;
        break;
    default:
        Gi = 7777777.;
    }
    return Gi;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Calculation of the total Gibbs energy of the system G(X).
*  Parameter LM is IPM step size for calculation of new
*  quantities of all species (vector X[]) using direction
*  of descent (MU[] vector). If LM==0, this function
*  just copies vector Y[] into X[].
*  Returns value of G(X) in moles.
*/
double TMulti::GX( double LM  )
{
    int i, j, k;
    double x, XF, XFw, FX, Gi /* debug variable */;

    if( LM<pmp->lowPosNum)     /* copy vector Y into X */
        for(i=0;i<pmp->L;i++)
            pmp->X[i]=pmp->Y[i];
    else  /*calculate new values of X */
        for(i=0;i<pmp->L;i++)
        {  /* vector pmp->MU - the direction of descent! */
            pmp->X[i]=pmp->Y[i]+LM*pmp->MU[i];
            if( pmp->X[i] <  pmp->lowPosNum )
                pmp->X[i]=0.;
        }
    /*calculate new total quantities of phases */
    TotalPhases( pmp->X, pmp->XF, pmp->XFA );

    /* calculate G(X) */
    FX=0.;
    j=0;
    for( k=0; k<pmp->FI; k++ )
    { /* loop for phases */
        i=j+pmp->L1[k];
        XFw = 0.0;  /* calc moles of solvent/sorbent */
        if( pmp->FIs && k<pmp->FIs )
            XFw = pmp->XFA[k];
        if( XFw > pmp->lowPosNum*10. )
            pmp->logXw = log( XFw );
        /*   */
        XF = pmp->XF[k];
        if( XF <= pmp->lowPosNum*1000. ||
                (pmp->PHC[k] == PH_AQUEL && (XF <= pmp->DHBM
                || XFw <= TProfil::pm->pa.p.XwMin) )
                || ( pmp->PHC[k] == PH_SORPTION && XFw <= TProfil::pm->pa.p.ScMin ))
            goto NEXT_PHASE;
        pmp->logYFk = log( XF );

        for( ; j<i; j++ )
        { /* Species */
            x = pmp->X[j];
            if( x < pmp->lowPosNum*10. )
                continue;
            /* calc increment of G(x) */
            Gi = FreeEnergyIncr( pmp->G[j], x, pmp->logYFk, pmp->logXw,
                                 pmp->DCCW[j] );
            FX += Gi;
        }   /* j */
NEXT_PHASE:
        j = i;
    }  /* k */
    return(FX);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// GX for ProbeUncert in statistic
// not use !! SD
double TMulti::pb_GX( double *Gxx  )
{
    int i, j;
    short k;
    double Gi, x, XF, XFw, FX;

    /* calc G(X) */
    FX=0.;
    j=0;
    for( k=0; k<pmp->FI; k++ )
    { /*phases*/
        i=j+pmp->L1[k];
        XFw = 0.0;  /* calc number of moles of carrier*/
        if( pmp->FIs && k<pmp->FIs )
            XFw = pmp->XFA[k];
        if( XFw > pmp->lowPosNum*10. )
            pmp->logXw = log( XFw );
        /* calc new quantitys of phase  */
        XF = pmp->XF[k];
        if( XF <= pmp->lowPosNum*1000. ||
           (pmp->PHC[k] == PH_AQUEL && (XF <= pmp->DHBM
                || XFw <= TProfil::pm->pa.p.XwMin) )
                || ( pmp->PHC[k] == PH_SORPTION && XFw <= TProfil::pm->pa.p.ScMin ))
            goto NEXT_PHASE;
        pmp->logYFk = log( XF );

        for( ; j<i; j++ )
        { /*DC */
            x = pmp->X[j];
            if( x < pmp->lowPosNum*10. )
                continue;
            /* calc increment of G(x) */
            Gi = FreeEnergyIncr( Gxx[j], x, pmp->logYFk, pmp->logXw,
                                 pmp->DCCW[j] );
            FX += Gi;
        }   /* j */
NEXT_PHASE:
        j = i;
    }  /* k */
    return(FX);
}

// Calc value cj (G0) by type of DC and standart value g(T,P)
// k - index of phase, j - index DC in phase
// if error code return 777777777.
double TMulti::Cj_init_calc( double g0, int j, int k )
{
    double G, YOF;

    G = g0/pmp->RT;
    /*  if( k < pmp->FIs )  */
    YOF = pmp->YOF[k];     // J/g:   check this!   04.12.2006  DK
    /* Calculation of concentration scaling increments */
    switch( pmp->DCC[j] )
    { /* Aqueous electrolyte */
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
        G += pmp->ln5551;
        /* calc mol weight of solvent !!!!!!!!!!!!!!!!!!!!!!!!!!! */
    case DC_AQ_SOLVCOM:
    case DC_AQ_SOLVENT:
#ifndef IPMGEMPLUGIN
        if( syp->PYOF != S_OFF )
//            pmp->GEX[j] += YOF;
        G += YOF;
#endif
        break;
    case DC_GAS_COMP: /* gases except H2O and CO2 */
    case DC_GAS_H2O: /* index to switch off? */
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR: // changed by DK on 4.12.2006
        if( pmp->PHC[k] == PH_GASMIX || pmp->PHC[k] == PH_FLUID
            || pmp->PHC[k] == PH_PLASMA )
        {
//        if( pmp->Pparc[j] != 1.0 && pmp->Pparc[j] > 1e-30 )
//           G += log( pmp->Pparc[j] ); // log partial pressure/fugacity
//        else
               G += log( pmp->Pc ); // log general pressure (changed 04.12.2006)
        }
        /* Solution of non-electrolyte */
    case DC_SCP_CONDEN: /*a single-component phase */
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:
    case DC_PEL_CARRIER:
#ifndef IPMGEMPLUGIN
        if( syp->PYOF != S_OFF )
//            pmp->GEX[j] += YOF;
           G += YOF;
#endif
        break;
        /* Sorption phases */
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
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
        G += pmp->ln5551;
        break;
    default: /* error code */
        return 7777777.;
    }
//    return G += pmp->GEX[j];
    return G;
}

//----------------------------------------------------------------
// Kostya: debug of calculation X(j) from Au

#define  a(j,i) ((double)(*(pmp->A+(i)+(j)*pmp->N)))
void TMulti::Mol_u( double Y[], double X[], double XF[], double XFA[] )
{
    int i,j,ja/*=0*/,jj,ii,jb,je,k;
    int isp/*=0*/, ist/*=0*/;
    double Ez, Psi;   // added  KD 23.11.01
//    double SPmol,SPmol1;
    double  Dsur, DsurT, MMC, *XU;
    XU = new double[pmp->L];
    ErrorIf( !XU, "Simplex", "Memory alloc error ");
    memset(XU, 0, sizeof(double)*(pmp->L) );

 //   ofstream ofs("c:/gems_b/x_u.txt",ios::out | ios::app);

  jb=0;
  for( k=0; k<pmp->FI; k++ )
  { /* cycle by phases */
      je=jb+pmp->L1[k];
      Dsur=0.0; DsurT=0.0;
      if( pmp->PHC[k] == PH_AQUEL && XF[k] > pmp->lowPosNum )
        Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr.
      if( (pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL)
            && XFA[k] > pmp->lowPosNum )
      {
         MMC = 0.0; // calculation of molar mass of the sorbent
         for( jj=jb; jj<je; jj++ )
         {
            if( pmp->DCC[jj] == DC_SUR_CARRIER ||
                pmp->DCC[jj] == DC_SUR_MINAL ||
                pmp->DCC[jj] == DC_PEL_CARRIER )
                    MMC += pmp->MM[jj]*X[jj]/XFA[k];
         }
         Dsur = XFA[k]/XF[k] - 1.0;
         DsurT = MMC * (double)(pmp->Aalp[k]) * TProfil::pm->pa.p.DNS*1.66054e-6;
      }

    for(j=jb;j<je;j++)
    {
      if( XF[k] >= pmp->DSM ) // pmp->lowPosNum ) fixed by KD 23.11.01
      {
 //        XU[j] = -pmp->G0[j] -pmp->lnGam[j]  changed 5.12.2006
         XU[j] = -pmp->G0[j] - pmp->lnGam[j] - pmp->GEX[j]
                  + DualChemPot( pmp->U, pmp->A+j*pmp->N, pmp->NR );
         if( pmp->PHC[k] == PH_AQUEL ) // pmp->LO && k==0)
         {
            if(j == pmp->LO)
                ;     //     disabled by KD 23.11.01
//                XU[j] += Dsur - 1. + 1. / ( 1.+ Dsur ) + log(XF[k]);
            else
                XU[j] += Dsur + log(XFA[k]);
         }
         else if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
         {
            if( pmp->DCC[j] == DC_PEL_CARRIER ||
                 pmp->DCC[j] == DC_SUR_CARRIER ||
                 pmp->DCC[j] == DC_SUR_MINAL )
                ;    //     disabled by KD 23.11.01
//                XU[j] += Dsur - 1.0 + 1.0 / ( 1.0 + Dsur )
//                      - DsurT + DsurT / ( 1.0 + DsurT ) + log(XF[k]);
            else  {    // rewritten by KD  23.11.01
               // Psi = 0.0;
               ja = j - ( pmp->Ls - pmp->Lads );
               Ez = pmp->EZ[j];
               /* Get ist - index of surface type */
/* !!!!! */    ist = pmp->SATX[ja][XL_ST]; // / MSPN;
               /* and isp - index of surface plane  */
//               isp = pmp->SATX[ja][XL_ST] % MSPN;
               isp = pmp->SATX[ja][XL_SP];
               if( !isp )
                   /* This is A (0) plane */
                   Psi = pmp->XpsiA[k][ist];
               else /* This is B or another plane */
                   Psi = pmp->XpsiB[k][ist];
               XU[j] += Dsur + DsurT/( 1.0 + DsurT ) + log(XFA[k])+
               log( DsurT * (double)(pmp->Nfsp[k][ist]) ) - pmp->FRT * Ez * Psi;
             }
         }
         else
           XU[j] += log(XF[k]);

         if( XU[j] > -42. && XU[j] < 42. )
             XU[j] = exp( XU[j] );
         else
             XU[j] = /* F_EMPTY*/ 0.0;
         if( XU[j] <= pmp->lowPosNum )
             XU[j]=0.;
      }
      else
          XU[j]=0.;
    }
    jb = je;
  }  /* k */
 /*   if(pmp->LO&& k==0)
     for( j=jb; j<je; j++ )
    {
            if(j==pmp->LO)
                   continue;
            SPmol = X[j]*Factor;  // molality
            SPmol1 = XU[j]*Factor;

            for( ii=0; ii<pmp->N; ii++ )
            {
                pmp->IC_m[ii] += SPmol* a(j,ii);
                MU[ii] += SPmol1* a(j,ii);
            }
     }
  */

 /*   ofs<<"---------------------"<<endl;
     ofs<<pmp->stkey<<endl;
    char elemb[20];
   if(pmp->LO)
    for( ii=0; ii<pmp->N; ii++ )
    {
        strncpy(elemb, pmp->SB[ii], 4);
        elemb[4]=0;
        ofs << elemb << pmp->IC_m[ii];
        ofs << "\t  "<< MU[ii];
        if(pmp->IC_m[ii]>0) ofs<<" \t "<<log10(pmp->IC_m[ii]);
        if(MU[ii]>0)  ofs<<" \t "<< log10(MU[ii])<<endl;
    }
    ofs<<endl;
    for( j=0; j<pmp->L; j++ )
    {   strncpy(elemb, pmp->SM[j], 16);
        elemb[16]=0;
        ofs << elemb  <<" \t"<<X[j]<<" \t "<< XU[j]<<endl;
    }
   }    */
    for( j=0; j<pmp->L; j++ )
    { ii=0;
      if(TProfil::pm->pa.p.PLLG)
      { for( i=0; i<pmp->N-pmp->E; i++ )
        if(a(i,j) && pmp->B[i] < pmp->DHBM*pow(10.,TProfil::pm->pa.p.DT))
        { ii=1; break; }
      }
      else
        if(Y[j]<pmp->DHBM*pow(10.,TProfil::pm->pa.p.DT))
          ii=1;
      if (ii)
        X[j]=XU[j];
      else
        X[j]=Y[j];
    }
    TotalPhases( X, XF, XFA );

    delete[] XU;
//   ofs.close();
}


// Convert class codes of DC into general codes of IPM
void TMulti::ConvertDCC()
{
    int i, j, k, iRet=0;
    char DCCW;

    j=0;
    for( k=0; k< pmp->FI; k++ )
    { /* cycle by phases  */
        i=j+pmp->L1[k];
        if( pmp->L1[k] == 1 )
        {
            pmp->DCCW[j] = DC_SINGLE;
            goto NEXT_PHASE;
        }
        for( ; j<i; j++ )
        { /* cykle by DC. */
            switch( pmp->DCC[j] ) /* selection necessary expressions  v_j */
            {
            case DC_SCP_CONDEN:
                DCCW = DC_SINGLE;
                break;
            case DC_GAS_COMP:
            case DC_GAS_H2O:
            case DC_GAS_CO2:
            case DC_GAS_H2:
            case DC_GAS_N2:
            case DC_SOL_IDEAL:
            case DC_SOL_MINOR:
            case DC_SOL_MAJOR:
                DCCW = DC_SYMMETRIC;
                break;
            case DC_AQ_PROTON:
            case DC_AQ_ELECTRON:
            case DC_AQ_SPECIES:
                DCCW = DC_ASYM_SPECIES;
                break;
            case DC_AQ_SOLVCOM:
            case DC_AQ_SOLVENT:
                DCCW = DC_ASYM_CARRIER;
                break;
            case DC_IESC_A:
            case DC_IEWC_B:
                DCCW = DC_ASYM_SPECIES;
                break;
                /* Remapping */
            case DC_SUR_GROUP:
            case DC_SUR_COMPLEX:
                DCCW = DC_ASYM_SPECIES;
                pmp->DCC[j] = DC_SSC_A0;
                break;
            case DC_SUR_IPAIR:
                DCCW = DC_ASYM_SPECIES;
                pmp->DCC[j] = DC_WSC_A0;
                break;
            case DC_SUR_MINAL:
            case DC_SUR_CARRIER:
            case DC_PEL_CARRIER:
                DCCW = DC_ASYM_CARRIER;
                break;
            default:
                if( isdigit( pmp->DCC[j] ))
                {
                    if( pmp->PHC[k] == PH_SORPTION ||
                            pmp->PHC[k] == PH_POLYEL )
                    {
                        DCCW = DC_ASYM_SPECIES;
                        break;
                    }
                }
                DCCW = DC_SINGLE;
                iRet++;  /* error in code  */
            }
            pmp->DCCW[j] = DCCW;
        }   /* j */
NEXT_PHASE:
        j = i;
    }  /* k */
    ErrorIf( iRet>0, "Multi", "Error in DCC code.");
}

//--------------------- End of ipm_chemical.cpp ---------------------------
