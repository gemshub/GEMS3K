//-------------------------------------------------------------------
// $Id: ipm_chemical.cpp 825 2006-03-29 07:10:23Z gems $
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
#include<iomanip>

#include "m_param.h"
#ifndef IPMGEMPLUGIN
#include "service.h"
#include "stepwise.h"
#endif

// #define GEMITERTRACE

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of max.moles of surface species for SACT stabilization
//  to improve IPM-2 convergence at high SACT values  KD 08.03.02
//  xj0 values are placed as upper kinetic constraints
//
void TMulti::XmaxSAT_IPM2()
{
    long int i, j, ja, k, jb, je=0, ist=0, Cj, iSite[6];
    double XS0, xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, xjn;

  if(!pmp->DUL )   // not possible to install upper kinetic constraints!
      return;

  for( k=0; k<pmp->FIs; k++ )
  { // loop over phases
     jb = je;
     je += pmp->L1[k];
     if( pmp->PHC[k] != PH_SORPTION )
          continue;

    if( pmp->XFA[k] < pmp->DSM ) // No sorbent retained by the IPM
        continue;
    if( pmp->XF[k]-pmp->XFA[k] < min( pmp->lowPosNum, pmp->DcMinM ) ) // may need qd_real in subtraction!
        continue;  // No surface species

    for(i=0; i<6; i++)
        iSite[i] = -1;

    // Extraction of site indices
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
    { // Loop over DCs
        if( pmp->X[j] <= min( pmp->lowPosNum, pmp->DcMinM ) )
            continue;  // This surface DC has been killed by the IPM
        rIEPS = TProfil::pm->pa.p.IEPS;
        ja = j - ( pmp->Ls - pmp->Lads );

        switch( pmp->DCC[j] )  // code of species class
        {
        default: // pmp->lnGam[j] = 0.0;
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
            // Calculate ist - index of surface type
/*!!!!!!*/  ist = pmp->SATX[ja][XL_ST]; // / MSPN;
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
                if( XVk < pmp->DSM/10.0 )
                    continue; // This end-member is zeroed off by IPM
                Mm = pmp->MM[Cj] * XVk/pmp->XFA[k];  // mol.mass
            }
            XSk = pmp->XFTS[k][ist]; // Tot.moles of sorbates on surf.type
            xj = pmp->X[j];  // Current moles of this surf.species
//            a=1.0;  Frumkin factor - reserved for extension to FFG isotherm
            switch( pmp->SATT[ja] )
            {
            case SAT_COMP: // Competitive surface species on a surface type
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pmp->X[iSite[ist]]; // neutral site does not compete!
                XS0 = max(pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN]);
                XS0 = XS0 * XVk * Mm / 1e6 * pmp->Nfsp[k][ist];
                            // expected total in moles
                XSkC = XSk - xjn - xj; // occupied by the competing species
                if( XSkC < 0.0 )
                    XSkC = rIEPS;
                xj0 = XS0 - XSkC;    // expected moles of this sorbate
                if( xj0 > pmp->lnSAC[ja][3] )
                    xj0 = pmp->lnSAC[ja][3];
                if( xj0 < rIEPS )
                   xj0 = rIEPS;  // ensuring that it will not zero off
                pmp->DUL[j] = xj0;
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
            case SAT_NCOMP: // Non-competitive surface species
                 xj0 = fabs( pmp->MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * pmp->Nfsp[k][ist]; // in moles
                 pmp->DUL[j] = xj0 - rIEPS;
                 // Compare with old DUL from previous iteration!
/*                if( pmp->W1 != 1 && pmp->IT > 0 && fabs( (pmp->DUL[j] - oDUL)/pmp->DUL[j] ) > 0.1 )
                {
cout << "XmaxSAT_IPM2 Ncomp IT= " << pmp->IT << " j= " << j << " oDUL=" << oDUL << " DUL=" << pmp->DUL[j] << endl;
                }
*/                break;

            case SAT_SOLV:  // Neutral surface site (e.g. >O0.5H@ group)
                rIEPS = TProfil::pm->pa.p.IEPS;
                XS0 = (max( pmp->MASDT[k][ist], pmp->MASDJ[ja][PI_DEN] ));
                XS0 = XS0 * XVk * Mm / 1e6 * pmp->Nfsp[k][ist]; // in moles

                pmp->DUL[j] =  XS0 - rIEPS;
                if( pmp->DUL[j] <= rIEPS )
                   pmp->DUL[j] = rIEPS;
                break;
// New methods added by KD 13.04.04
            case SAT_INDEF: // No SAT calculation
            default:        // pmp->lnGam[j] = 0.0;
                break;
            }
        }
     }  // j
  } // k
}

/* // clearing pmp->DUL constraints!
void TMulti::XmaxSAT_IPM2_reset()
{
    long int j, ja, k, jb, je=0;

  if(!pmp->DUL )   // no upper kinetic constraints!
      return;

  for( k=0; k<pmp->FIs; k++ )
  { // loop on phases
     jb = je;
     je += pmp->L1[k];
     if( pmp->PHC[k] != PH_SORPTION )
          continue;

    for( j=jb; j<je; j++ )
    { // Loop for DC
      ja = j - ( pmp->Ls - pmp->Lads );
      if( pmp->lnSAC && ja >= 0 && ja < pmp->Lads )
          pmp->DUL[j] = pmp->lnSAC[ja][3];  // temp. storing initial DUL constr.
    }  // j
  } // k
}
*/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculating value of dual chemical potential
//     performance optimized version  (February 2007)
//
double TMulti::DC_DualChemicalPotential( double U[], double AL[], long int N, long int j )
{
    long int i, ii;
    double Nu = 0.0;
//    for(long int i=; i<N; i++ )
//    Nu += AL[i]? U[i]*(AL[i]): 0.0;
   for( i=arrL[j]; i<arrL[j+1]; i++ )
   {  ii = arrAN[i];
      if( ii>= N )
       continue;
       Nu += U[ii]*(AL[ii]);
   }
   return Nu;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  This procedure sets kinetic constraints according to a given
//  concentration units
//  Needs much more work, elaboration, and performance optimization
//
void TMulti::Set_DC_limits( long int Mode )
{
    double XFL, XFU, XFS=0., XFM, MWXW, MXV, XL=0., XU=0.;
    long int jb, je, j,k, MpL;
    char tbuf[150];

    if( !pmp->PLIM )
        return;  // no metastability limits to be set
// ???????????????????????????????????????
    CalculateConcentrations( pmp->X, pmp->XF, pmp->XFA );

    for(k=0; k<pmp->FI; k++)
        XFS+=pmp->XF[k];  // calculate sum of moles in all phases

    jb=0;
    for( k=0; k<pmp->FI; k++ )
    { // cycle over phases
        je=jb+pmp->L1[k];
//        XFM=0.;
        MWXW =0.;
        MXV = 0.;
        XFL = 0.;
        XFU = 1e6;
        if( Mode && pmp->XF[k] < pmp->DSM )
            goto NEXT_PHASE;
        XFM = pmp->FWGT[k]; // Mass of a phase
        if( Mode )
        {
            MWXW = XFM/pmp->XF[k];         // current molar mass of phase
            MXV = pmp->FVOL[k]/pmp->XF[k]; // current molar volume of phase
        }
        // Check codes for phase DC
        MpL=0;
        for( j=jb; j<je; j++ )
            if( pmp->RLC[j] != NO_LIM )
                MpL = 1;

if( k < pmp->FIs )
{					// Temporary workaround - DK  13.12.2007
        if( pmp->RFLC[k] == NO_LIM && !MpL )
        { // check type restrictions on phase
            goto NEXT_PHASE;
        }
        switch( pmp->RFSC[k] )
        { // check scale restrictions on phase in all system
        case QUAN_MOL:
            XFL = Mode? pmp->XF[k]: pmp->PLL[k];
            XFU = Mode? pmp->XF[k]: pmp->PUL[k];
            break;
        case CON_MOLAL:
            XFL = Mode? pmp->XF[k]: pmp->PLL[k]*pmp->GWAT/H2O_mol_to_kg;
            XFU = Mode? pmp->XF[k]: pmp->PUL[k]*pmp->GWAT/H2O_mol_to_kg;
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
            ; // do more?
        }
//        if( pmp->RFLC[k] == NO_LIM )
//        {                            Temporary!
            XFL = 0.0;
            XFU = 1e6;
//        }
}
        for( j=jb; j<je; j++ )
        { // loop over DCs
            if( pmp->RLC[j] == NO_LIM )
                continue;

            if( Mode )
            {
                XU = pmp->DUL[j];
                XL = pmp->DLL[j];
            }
            else
                switch( pmp->RSC[j] ) // get initial limits
                {
                case QUAN_MOL:
                    XU = pmp->DUL[j];
                    XL = pmp->DLL[j];
                    break;
                case CON_MOLAL:
                    XU = pmp->DUL[j]*pmp->GWAT/H2O_mol_to_kg;
                    XL = pmp->DLL[j]*pmp->GWAT/H2O_mol_to_kg;
                    break;
                case CON_MOLFR:
                    XU = pmp->DUL[j]*XFU;
                    XL = pmp->DLL[j]*XFL;
                    break;
                case CON_WTFR:
//Ask DK! 20/04/2002
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
                    ; // do more
                }
            // check combine
            if( XU < 0.0 ) XU = 0.0;
            if( XU > 1e6 ) XU = 1e6;
            if( XL < 0.0 ) XL = 0.0;
            if( XL > 1e6 ) XL = 1e6;
            if( XU > XFU )
            {
 //               JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent upper DC metastability limits j=%ld k=%ld XU=%g XFU=%g",
                         j, k, XU, XFU );
                Error( "E11IPM: Set_DC_limits(): ",tbuf );
//                XU = XFU; // - pmp->lowPosNum;
            }
            if( XL < XFL )
            {
//                JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent lower DC metastability limits j=%ld k=%ld XL=%g XFL=%g",
                         j, k, XL, XFL );
                Error( "E12IPM: Set_DC_limits(): ",tbuf );
//                XL = XFL; // - pmp->lowPosNum;
            }
            pmp->DUL[j]=XU;
            pmp->DLL[j]=XL;
        }   // j
NEXT_PHASE:
        jb = je;
    }  // k
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculating total amounts of phases
//
void TMulti::TotalPhasesAmounts( double X[], double XF[], double XFA[] )
{
    long int jj, j, i, k;
    double XFw, XFs, x;

    j=0;
    for( k=0; k< pmp->FI; k++ )
    { // cycle by phases
        i=j+pmp->L1[k];
        XFw = 0.0;
        XFs=0.0; // calculating mole amount of carrier (solvent/sorbent)
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Corrections to primal chemical potentials F0[j]
//  of j-th species in k-th phase among IPM main iterations
//  Returns double value of corrected chem. potential.
//  If error, returns +7777777 J/mole.
//  Last modif. 05 Jan 2000 by DK to include BSM EDL model.
//
double TMulti::DC_PrimalChemicalPotentialUpdate( long int j, long int k )
{
    long int ja=0, ist, isp, jc=-1;
    double F0=0.0, Fold, dF0, Mk=0.0, Ez, psiA, psiB, CD0, CDb, ObS;
    double FactSur, FactSurT;
    SPP_SETTING *pa = &TProfil::pm->pa;

    Fold = pmp->F0[j];
    if( pmp->FIat > 0 && j < pmp->Ls && j >= pmp->Ls - pmp->Lads )
    {
        ja = j - ( pmp->Ls - pmp->Lads );
        jc = pmp->SATX[ja][XL_EM];
    }
    if( k < pmp->FIs && pmp->XFA[k] > 1e-12)
    {
           if( jc < 0 ) // phase (carrier) molar mass g/mkmol
              Mk = pmp->FWGT[k]/pmp->XFA[k]*1e-6;
           else Mk = pmp->MM[jc]*(pmp->X[jc]/pmp->XFA[k])*1e-6;
        // DC carrier molar mass g/mkmol
    }
    switch( pmp->DCC[j] )
    { // analyse species class code
    case DC_SCP_CONDEN:
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
case DC_AQ_SURCOMP:
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
        // adsorption
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
        if( !pmp->Lads || !pmp->FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
        F0 = pmp->lnGmM[j]; /* + pmp->lnGam[j]; */
        // get ist - index of surface type and isp - index of surface plane
/*!!!!!*/  ist = pmp->SATX[ja][XL_ST];  // / MSPN;
/*!!!!!*/  isp = pmp->SATX[ja][XL_SP]; // % MSPN;
        CD0 = pmp->MASDJ[ja][PI_CD0];  // species charge that goes into 0 plane
        CDb = pmp->MASDJ[ja][PI_CDB];  // species charge that goes into B plane
        ObS = pmp->MASDJ[ja][PI_DEN];  // obsolete - the sign for outer-sphere charge
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
        else  // This is B plane
        {
            if( pmp->SCM[k][ist] == SC_MTL || pmp->SCM[k][ist] == SC_MXC )
            { // Modified TL: Robertson, 1997; also XTLM Kulik 2002
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
            { // Basic Stern model, Christl & Kretzschmar, 1999
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
        if( Mk > 1e-9 )  // Mk is carrier molar mass in g/mkmol
        {   // Correction for standard density, surface area and surface type fraction
        	FactSur = Mk * (pmp->Aalp[k]) * pa->p.DNS*1.66054;
        	    // FactSur is adsorbed mole amount at st. surf. density per mole of solid carrier
        	FactSurT = FactSur * (pmp->Nfsp[k][ist]);
        	if( pmp->SCM[k][ist] == SC_MXC || pmp->SCM[k][ist] == SC_NNE ||
                    pmp->SCM[k][ist] == SC_IEV )
						// F0 -= log( Mk * (pmp->Nfsp[k][ist]) *
						// (pmp->Aalp[k]) * pa->p.DNS*1.66054 );
                  F0 -= log( FactSurT );
            else  F0 -= log( FactSurT );
				// F0 -= log( Mk * (pmp->Nfsp[k][ist]) *
				// (pmp->Aalp[k]) * pa->p.DNS*1.66054 );
            F0 -= FactSur / ( 1. + FactSur );
				// F0 -= (pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 /
				// ( 1.0 + (pmp->Aalp[k])*Mk*pa->p.DNS*1.66054 );
        }
        break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:  // constant charge of carrier - not completed
    case DC_SUR_CARRIER: // Mk is carrier molar mass in g/mkmol
        if( !pmp->Lads || !pmp->FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
       	FactSur = Mk * (pmp->Aalp[k]) * pa->p.DNS*1.66054;
        F0 -= FactSur / ( 1. + FactSur );
        F0 += FactSur;
        break;
    }
    F0 += pmp->lnGam[j];

    if( k >= pmp->FIs )
        return F0;
    // Smoothing procedure for highly non-ideal systems
    if( pmp->sMod[k][SGM_MODE] != SM_IDEAL )  // check this condition for sublattice SS models!
            // || pmp->sMod[k][SCM_TYPE] != SC_NNE )  // changed, 14.07.2009 (TW)
    {
        double SmoSensT = 1e-5;   // to be adjusted
        dF0 = F0 - Fold;
        if( pmp->X[j] > min( pmp->lowPosNum, pmp->DcMinM ) && fabs( dF0 ) >= SmoSensT )
       	    F0 = Fold + dF0 * SmoothingFactor();    // Changed 18.06.2008 DK
    }
    return F0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of DC primal chemical potential F (return value)
// from moles of DC Y[], total moles of phase YF[] and DC partial
// molar Gibbs energy gT (obtained from pmp->G[]) which includes
// activity coefficient terms.
// On error returns F = +7777777.
//
double TMulti::DC_PrimaChemicalPotential(
    double G,      // gT0+gEx
    double logY,   // ln x
    double logYF,  // ln Xa
    double asTail, // asymmetry non-log term or 0 for symmetric phases
    double logYw,  // ln Xv
    char DCCW      // generalized species class code
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


// Kernel functions of IPM - rewritten by DK for adsorption
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// VJ - Update of primal chemical potentials
//
void
TMulti::PrimalChemicalPotentials( double F[], double Y[], double YF[], double YFA[] )
{
    long int i,j,k;
    double v, Yf; // v is debug variable

    for( j=0; j<pmp->L; j++)
       F[j] =0;

    j=0;
    for( k=0; k<pmp->FI; k++ )
    { // loop over phases
        i=j+pmp->L1[k];
//        if( YF[k] <= pmp->lowPosNum*100. || ( pmp->PHC[k] == PH_AQUEL &&
//        ( YF[k] <= TProfil::pm->pa.p.XwMin || Y[pmp->LO] <= pmp->lowPosNum*1e3 )))
        if( pmp->L1[k] == 1L && YF[k] < pmp->PhMinM )
        	goto NEXT_PHASE;
        if( YF[k] <= pmp->DSM || ( pmp->PHC[k] == PH_AQUEL &&
            ( YF[k] <= pmp->DSM || Y[pmp->LO] <= pmp->XwMinM )))
            goto NEXT_PHASE;

        pmp->YFk = 0.0;
        Yf= YF[k]; // calculate number of moles of carrier
        if( pmp->FIs && k<pmp->FIs )
            pmp->YFk = YFA[k];
        if( Yf >= 1e6 )
        {                 // error - will result in zerodivide!
           gstring pbuf(pmp->SF[k],0,20);
           char buf[200];
           sprintf( buf, "Broken phase amount from primal approximation: Phase %s  Yf= %lg", pbuf.c_str(), Yf );
           Error( "E13IPM: PrimalChemicalPotentials():", buf);
//           Yf = pmp->YFk;
        }
//        if( pmp->YFk > pmp->lowPosNum*10. )
        if( ( pmp->PHC[k] == PH_AQUEL && pmp->YFk >= pmp->XwMinM )
        		|| ( pmp->PHC[k] == PH_SORPTION && pmp->YFk >= pmp->ScMinM )
        		|| ( pmp->PHC[k] == PH_POLYEL && pmp->YFk >= pmp->ScMinM ) )
        {
            pmp->logXw = log(pmp->YFk);
            pmp->aqsTail = 1.- pmp->YFk / Yf;
        }
        if( pmp->L1[k] > 1 )
        {
            pmp->logYFk = log( Yf );
        }
        if( pmp->PHC[k] == PH_AQUEL )
            // ln moles of solvent in aqueous phase
            pmp->Yw = pmp->YFk;
        for( ; j<i; j++ )
        { //  cycle by DC
            if( Y[j] < min( pmp->DcMinM, pmp->lowPosNum ))
                continue;  // exception by minimum DC quantity
                           // calculate chemical potential of j-th DC
            v = DC_PrimaChemicalPotential( pmp->G[j], log(Y[j]), pmp->logYFk,
                              pmp->aqsTail, pmp->logXw, pmp->DCCW[j] );
            F[j] = v;
       }   // j
NEXT_PHASE:
        j = i;
    }  // k
    if( pmp->Yw >= pmp->DSM ) // pmp->lowPosNum*1e3 )
        pmp->logXw = log(pmp->Yw);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of a species contribution to the total Gibbs energy G(X)
//  of the system (return value).
//  On error returns +7777777.
//
double TMulti::DC_GibbsEnergyContribution(
    double G,      // gT0+gEx
    double x,      // x - mole amount of species
    double logXF,  // ln Xa - mole amount of phase
    double logXw,  // ln Xv - mole amount of the solvent/sorbent
    char DCCW      // generalized species class code
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of the total Gibbs energy of the system G(X)
// and copying of Y, YF vectors into X,XF, respectively.
//  Parameter LM is the IPM step size for calculation of new
//  quantities of all species (vector X[]) using the direction
//  of descent (MU[] vector). If LM == 0, this function
//  just copies vector Y[] into X[].
//  Returns value of G(X) in moles.
//
double TMulti::GX( double LM  )
{
    long int i, j, k;
    double x, XF, XFw, FX, Gi; // debug variable
//    double const1= pmp->lowPosNum*10.,
//           const2 = pmp->lowPosNum*1000.;

    if( LM < pmp->lowPosNum )     // copy vector Y into X
        for(i=0;i<pmp->L;i++)
            pmp->X[i]=pmp->Y[i];
    else  // calculate new values of X
        for(i=0;i<pmp->L;i++ )
        {  // gradient vector pmp->MU - the direction of descent!
            pmp->X[i]=pmp->Y[i]+LM*pmp->MU[i];
//            if( pmp->X[i] <  pmp->lowPosNum )   // this is the Ls set cutoff !!!!!!!!!!
            if( pmp->X[i] <  pmp->DcMinM )
            	pmp->X[i]=0.;
        }
    // calculate new total quantities of phases
    TotalPhasesAmounts( pmp->X, pmp->XF, pmp->XFA );

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<pmp->FI; k++ )
    { // loop for phases
        i=j+pmp->L1[k];
        pmp->logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( pmp->FIs && k<pmp->FIs )
            XFw = pmp->XFA[k];
 //       if( XFw > const1 )
        if( ( pmp->PHC[k] == PH_AQUEL && XFw >= pmp->XwMinM )
        		|| ( pmp->PHC[k] == PH_SORPTION && XFw >= pmp->ScMinM )
        		|| ( pmp->PHC[k] == PH_POLYEL && XFw >= pmp->ScMinM ) )
             pmp->logXw = log( XFw );
        /*   */
        XF = pmp->XF[k];
        if( !(pmp->FIs && k < pmp->FIs) )
        {
        	if( XF < pmp->PhMinM )
        		goto NEXT_PHASE;
        }
        else if( XF < pmp->DSM && pmp->logXw < -100. )
        	goto NEXT_PHASE;

        pmp->logYFk = log( XF );

        for( ; j<i; j++ )
        { // DCs (species)
            x = pmp->X[j];
            if( x < pmp->DcMinM )
                continue;
            // calculating increment of G(x)
            // Gi = DC_GibbsEnergyContribution( pmp->G[j], x, pmp->logYFk, pmp->logXw,
            //                     pmp->DCCW[j] );
            // call replaced here by inline variant for higher performance
            switch( pmp->DCCW[j] )
            {
             case DC_ASYM_SPECIES:
                    Gi = x * ( pmp->G[j] + log(x) - pmp->logXw );
                    break;
            case DC_ASYM_CARRIER:
            case DC_SYMMETRIC:
                   Gi = x * ( pmp->G[j] + log(x) - pmp->logYFk );
                   break;
            case DC_SINGLE:
                   Gi = pmp->G[j] * x;
                   break;
           default:
                    Gi = 7777777.;
           }
          FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
//    cout << "GX  " << setprecision(16) << scientific <<  FX << endl;
    return(FX);
}

#ifdef Use_qd_real
// Experimental variant with summation of G(X) with qd_real precision
// Added by SD, DK on 27.08.2009
qd_real TMulti::qdGX( double LM  )
{
    long int i, j, k;
    double x, XF, XFw, Gi; // debug variable
    qd_real FX;
//    double const1= pmp->lowPosNum*10.,
//           const2 = pmp->lowPosNum*1000.;
    if( LM < pmp->lowPosNum )     // copy vector Y into X
        for(i=0;i<pmp->L;i++)
            pmp->X[i]=pmp->Y[i];
    else  // calculate new values of X
        for(i=0;i<pmp->L;i++ )
        {  // vector pmp->MU - the direction of descent!
            pmp->X[i]=pmp->Y[i] + LM*pmp->MU[i];
//            if( pmp->X[i] <  pmp->lowPosNum )   // this is the Ls set cutoff !!!!!!!!!!
            if( pmp->X[i] <  pmp->DcMinM )
            	pmp->X[i]=0.;
        }
    // calculate new total quantities of phases
    TotalPhasesAmounts( pmp->X, pmp->XF, pmp->XFA );

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<pmp->FI; k++ )
    { // loop for phases
        i=j+pmp->L1[k];
        pmp->logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( pmp->FIs && k<pmp->FIs )
            XFw = pmp->XFA[k];
 //       if( XFw > const1 )
        if( ( pmp->PHC[k] == PH_AQUEL && XFw >= pmp->XwMinM )
        		|| ( pmp->PHC[k] == PH_SORPTION && XFw >= pmp->ScMinM )
        		|| ( pmp->PHC[k] == PH_POLYEL && XFw >= pmp->ScMinM ) )
             pmp->logXw = log( XFw );
        /*   */
        XF = pmp->XF[k];
        if( !(pmp->FIs && k < pmp->FIs) )
        {
        	if( XF < pmp->PhMinM )
        		goto NEXT_PHASE;
        }
        else if( XF < pmp->DSM && pmp->logXw < -100. )
        	goto NEXT_PHASE;
//        if( XF <= const2 ||
//                (pmp->PHC[k] == PH_AQUEL && (XF <= pmp->DHBM
//                || XFw <= TProfil::pm->pa.p.XwMin) )
//                || ( pmp->PHC[k] == PH_SORPTION && XFw <= TProfil::pm->pa.p.ScMin ))
//            goto NEXT_PHASE;
        pmp->logYFk = log( XF );

        for( ; j<i; j++ )
        { // DCs (species)
            x = pmp->X[j];
            if( x < pmp->DcMinM )
                continue;
            // calculating increment of G(x)
            // Gi = DC_GibbsEnergyContribution( pmp->G[j], x, pmp->logYFk, pmp->logXw,
            //                     pmp->DCCW[j] );
            // call replaced here by inline variant for higher performance
            switch( pmp->DCCW[j] )
            {
             case DC_ASYM_SPECIES:
                    Gi = x * ( pmp->G[j] + log(x) - pmp->logXw );
                    break;
            case DC_ASYM_CARRIER:
            case DC_SYMMETRIC:
                   Gi = x * ( pmp->G[j] + log(x) - pmp->logYFk );
                   break;
            case DC_SINGLE:
                   Gi = pmp->G[j] * x;
                   break;
           default:
                    Gi = 7777777.;
           }
          FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
//    pmp->qdFX = FX;
//    cout << "qdGX " << setprecision(20) << scientific << FX << endl;
    return(FX);
}
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Variant of GX() function for use in the UnSpace module (non-optimized)
// Should not be called from within GEMIPM!
//
double TMulti::pb_GX( double *Gxx  )
{
    long int i, j, k;
    double Gi, x, XF, XFw, FX;
    SPP_SETTING *pa = &TProfil::pm->pa;

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<pmp->FI; k++ )
    { // phase loop
        i=j+pmp->L1[k];
        pmp->logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( pmp->FIs && k<pmp->FIs )
            XFw = pmp->XFA[k];
 //       if( XFw > const1 )
        if( ( pmp->PHC[k] == PH_AQUEL && XFw >= pa->p.XwMin )
                        || ( pmp->PHC[k] == PH_SORPTION && XFw >= pa->p.ScMin )
                        || ( pmp->PHC[k] == PH_POLYEL && XFw >= pa->p.ScMin ) )
             pmp->logXw = log( XFw );
        /*   */
        XF = pmp->XF[k];
        if( !(pmp->FIs && k < pmp->FIs) )
        {
                if( XF < pa->p.PhMin )
        		goto NEXT_PHASE;
        }
        else if( XF < pa->p.DS && pmp->logXw < 100. )
        	goto NEXT_PHASE;
        pmp->logYFk = log( XF );

        for( ; j<i; j++ )
        { // DC loop
            x = pmp->X[j];
            if( x < pa->p.DcMin )
                continue;
            // calculating DC increment to G(x)
            Gi = DC_GibbsEnergyContribution( Gxx[j], x, pmp->logYFk, pmp->logXw,
                                 pmp->DCCW[j] );
            FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    return(FX);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Conversion of g(T,P) value for DCs into the uniform cj scale
// k - index of phase, j - index DC in phase
// if error code, returns 777777777.
//
double TMulti:: DC_ConvertGj_toUniform_cj( double g0, long int j, long int k )
{
    double G, YOF=0;

    G = g0/pmp->RT;
    if( pmp->YOF )
        YOF = pmp->YOF[k];     // J/g:   check this!   04.12.2006  DK
    // Calculation of standard concentration scaling terms
    switch( pmp->DCC[j] )
    { // Aqueous electrolyte
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
case DC_AQ_SURCOMP:
        G += pmp->ln5551;
        // calculate molar mass of solvent
    case DC_AQ_SOLVCOM:
    case DC_AQ_SOLVENT:
#ifndef IPMGEMPLUGIN
        if( syp->PYOF != S_OFF )
#endif
          if( YOF != 0.0 )
        	G += YOF;  // In GEMIPM2K, YOF[k] is the only way to influence G directly

        break;
    case DC_GAS_COMP: // gases except H2O and CO2
    case DC_GAS_H2O: // index to switch off?
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
        // non-electrolyte condensed mixtures
    case DC_SCP_CONDEN: // single-component phase
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:
    case DC_PEL_CARRIER:
#ifndef IPMGEMPLUGIN
        if( syp->PYOF != S_OFF )
#endif
          if( YOF != 0.0 )
        	 G += YOF;
        break;
        // Sorption phases
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
    default: // error - returning 7777777
        return 7777777.;
    }
    return G;
}

// Converting DC class codes into generic internal codes of IPM
//
void TMulti::ConvertDCC()
{
    long int i, j, k, iRet=0;
    char DCCW;

    j=0;
    for( k=0; k< pmp->FI; k++ )
    { // phase loop
        i=j+pmp->L1[k];
        if( pmp->L1[k] == 1 )
        {
            pmp->DCCW[j] = DC_SINGLE;
            goto NEXT_PHASE;
        }
        for( ; j<i; j++ )
        { // DC loop
            switch( pmp->DCC[j] ) // select v_j expression
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
            case DC_AQ_SURCOMP:
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
                // Remapping
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
                iRet++;  // error the class code
            }
            pmp->DCCW[j] = DCCW;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    ErrorIf( iRet>0, "E19IPM: ConvertDCC()", "Invalid DC class code. Memory corruption?");
}

// get the index of volume IC ("Vv") for the volume balance constraint
long int TMulti::getXvolume()
{
 long int ii, ret = -1;
 for( ii = pmp->N-1; ii>=0; ii--)
 {
  if( pmp->ICC[ii] == IC_VOLUME )
  { ret = ii; break; }
 }
 return ret;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of Karpov stability criteria for a DC
// Modified for kinetic constraints 05.11.2007 by DK
//
double TMulti::KarpovCriterionDC(
    double *dNuG,  // Nu[j]-c[j] difference - is modified here
    double logYF,  // ln Xa   (Xa is mole amount of the whole phase)
    double asTail, // asymmetry correction (0 for symmetric phases)
    double logYw,  // ln Xw   (Xw is mole amount of solvent)
    double Wx,     // mole fraction of this DC
    char DCCW      // Generic class code of DC
)
{
    double Fj=0.0;  // output phase stability criterion

    if( logYF > -35. && Wx > 1e-18 )    // Check thresholds!
        switch( DCCW ) // expressions for fj
        {
        default: // error code would be needed here !!!
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
    Fj -= Wx;                    // If Wx = 0 then this DC is not in L_S set

    return Fj;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// calculation of Karpov stability criteria for all phases
//
void TMulti::KarpovCriterionPH()
{
    bool KinConstr;
    long int k, j, ii;
    double *EMU,*NMU, YF, Nu, dNuG, Wx, Yj, Fj;
    SPP_SETTING *pa = &TProfil::pm->pa;

    EMU = pmp->EMU;
    NMU = pmp->NMU;
    for(ii=0; ii<pmp->L; ii++ )
        EMU[ii] = NMU[ii]=0.0;
    j=0;
    pmp->YMET = 0.0;
    for( k=0; k<pmp->FI; k++ )
    { // phases
        ii=j+pmp->L1[k];
        pmp->Falp[k] = pmp->YMET; // metastability parameter
        pmp->logXw = -35.;
        pmp->logYFk = -35.;

        pmp->YFk = 0.0;
        YF= pmp->YF[k]; // moles of carrier
        if( pmp->FIs && k<pmp->FIs )
            pmp->YFk = pmp->YFA[k];
        if( pmp->YFk > 6.305117e-16 )   // check threshold!
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
        if( pmp->PHC[k] == PH_AQUEL) // number of moles of solvent
            pmp->Yw = pmp->YFk;

// The code below was re-arranged by DK on 2.11.2007
        if(pmp->L1[k] == 1 )
        {   // This is a single-component phase - always included in L_S set
            KinConstr = false;
            Wx = 1.0;
            Yj = pmp->Y[j];
            Nu = DC_DualChemicalPotential( pmp->U, pmp->A+j*pmp->N, pmp->NR, j );
            dNuG = Nu - pmp->G[j]; // this is -s_j (6pot paper 1)
            if( // pmp->DUL[j] < pa->p.DKIN  ||    // DKIN (1e-6) is used here as tolerance
                 ( pmp->DUL[j] < 1e6 && Yj >= ( pmp->DUL[j] - pa->p.DKIN ) )
                || ( pmp->DLL[j] > 0 && Yj <= ( pmp->DLL[j] + pa->p.DKIN ) ) )
                KinConstr = true; // Avoiding phase with the amount lying on the non-trivial kinetic constraint
            Fj = KarpovCriterionDC( &dNuG, pmp->logYFk, pmp->aqsTail,
                            pmp->logXw, Wx, pmp->DCCW[j] );
            NMU[j] = dNuG;
            EMU[j] = Fj;
            if( KinConstr == false )
                pmp->Falp[k] = Fj;   // Karpov criterion of pure phase
        }
        else {
        // This is a multi-component phase
            for( ; j<ii; j++ )
            {
                KinConstr = false;
                Nu = DC_DualChemicalPotential( pmp->U, pmp->A+j*pmp->N, pmp->NR, j );
                dNuG = Nu - pmp->G[j]; // this is -s_j (6pot paper 1)
                Wx = 0.0;
                Yj = pmp->Y[j];
//                if( YF > pa->p.DS && Yj > pmp->lowPosNum )
                if( YF > pmp->DSM && Yj > pmp->DcMinM )
                    Wx = Yj / YF; // calculating mole fraction of DC
                if( // pmp->DUL[j] < pa->p.DKIN ||      // DKIN (1e-6) is used here as tolerance
                        ( pmp->DUL[j] < 1e6 && Yj >= ( pmp->DUL[j] - pa->p.DKIN ) )
                    || ( pmp->DLL[j] > 0 && Yj <= ( pmp->DLL[j] + pa->p.DKIN ) ) )
                    KinConstr = true; // Avoiding DC with the amount lying on the non-=trivial kinetic constraint
                // calculating Karpov stability criteria for DCs
                Fj = KarpovCriterionDC( &dNuG, pmp->logYFk, pmp->aqsTail,
                         pmp->logXw, Wx, pmp->DCCW[j] );
                NMU[j] = dNuG;  // dNuG is stored for all DCs, not only those in L_S set
// Experimental option - checking zeroed off DCs in multicomponent phases
// during the first PhaseSelect() run in the PIA mode of GEMIPM  (DK 11.01.2008)
if( pmp->pNP && !pmp->K2 )
{  // Checking L_S set and potentially stable zero DCs and phases
   if( YF >= pmp->DSM ) // pa->p.DS )
   {                            // phase is there
           if( Yj >= min(pmp->lowPosNum, pmp->DcMinM ) )
       {	                    // DC is there
                   if( KinConstr == false)
              pmp->Falp[k] += Fj; // incrementing Karpov stability criterion (only positive)
       }
           else {                  // DC is zeroed off
                   if( Fj > pa->p.DF*100. )
                  pmp->Falp[k] += Fj;    // we are interested only in a potentially stable DC with Fj >> DF ?
           }
       EMU[j] = Fj;
   }
   else {  // The whole phase is absent
           if( Fj > pa->p.DF*100. )
          pmp->Falp[k] += Fj;    // we are interested only in a potentially stable DC with Fj >> DF ?
       EMU[j] = Fj;              // To check values (remove this line later on)
   }
}
else {   // Standard checking in the L_S set only
//       if( YF >= pa->p.DS && Yj > pmp->lowPosNum )  // Checking L_S set
         if( YF >= pmp->DSM && Yj >= min( pmp->lowPosNum, pmp->DcMinM ) )  // Checking L_S set
         {
             if( KinConstr == false )
                 pmp->Falp[k] += Fj; // incrementing Karpov stability criterion for the phase
                    EMU[j] = Fj;
             }
             else
                    EMU[j] = 0;   // This DC is not in L_S set: e.g. metastability constrained
}
             }   // j
        }
        j = ii;
    }  // k
}

//===================================================================
// Parameters:
//     AmountCorrectionThreshold - the maximum DC amount correction that can be cleaned ( 1e-5 )
//     MjuDiffCutoff - normalized chem.pot. difference threshold (dMu = ln a - ln a,dual)
//
// Returns 0 if no subsequent refinement of mass balance is needed;
//         1 if species amounts were cleaned up to more than requested overall mass balance accuracy
//        -1 if cleanup has been done and the degeneration of the chemical system occurred
//         2 solution is seriously distorted and full PhaseSelect3() loop is necessary
//
long int TMulti::CleanupSpeciation( double AmountCorrectionThreshold, double MjuDiffCutoff )
{
    long int NeedToImproveMassBalance = 0, L1k, L1kZeroDCs, k, j, jb = 0;
    double MjuPrimal, MjuDual, MjuDiff, Yj, YjDiff=0., YjCleaned;
    double CutoffDistortionMBR = 0.1 * pmp->DHBM;
    bool KinConstr, Degenerated = false;
    SPP_SETTING *pa = &TProfil::pm->pa;

    PrimalChemicalPotentials( pmp->F, pmp->Y, pmp->YF, pmp->YFA );
    jb=0;
    for(k=0;k<pmp->FI;k++)
    {
       if( ( pmp->YF[k] >= pmp->DcMinM ) ) // Only in phase present in mass balance!
       {                            // (acc. to definition of the L_S set)
            L1k = pmp->L1[k]; // Number of components in the phase
            L1kZeroDCs = 0;
            for(j=jb; j<jb+L1k; j++)
            {
               Yj = YjCleaned = pmp->Y[j];
               KinConstr = false;
                // Detecting the DC having the non-trivial kinetic constraint
               // Fixing a very small component constrained from below
               if( pmp->DUL[j] < 1e6 && Yj >= ( pmp->DUL[j] - pa->p.DKIN ) )
               {
                   pmp->Y[j] = pmp->DUL[j];
                   KinConstr = true;
               }
               if( pmp->DLL[j] > 0 && Yj <= ( pmp->DLL[j] + pa->p.DKIN ) )
               { // Fixing a small component constrained from above
                   pmp->Y[j] = pmp->DLL[j];
                   KinConstr = true;
               }
               if( KinConstr == true )
               {
                   YjDiff = fabs( pmp->Y[j] -Yj );
                   if( YjDiff > CutoffDistortionMBR )
                       NeedToImproveMassBalance = 1;
                   if( YjDiff > AmountCorrectionThreshold )
                       NeedToImproveMassBalance = 2;
//                   continue;   // skipping if non-trivial metastability constraint
               }
               else if( Yj >= pmp->DcMinM )
               {   // we check in the Ls set only, except metastability constraints
                  MjuPrimal = pmp->F[j];   // normalized
                  MjuDual = pmp->Fx[j]/pmp->RT;
                  MjuDiff = MjuPrimal - MjuDual;
                  if( fabs( MjuDiff ) > MjuDiffCutoff )
                  {
                      if( L1k == 1 && MjuDiff > 0. )
                      {  // Pure phase
                         YjCleaned = 0.0; // Cleaning out a "phantom" pure phase
                      }
                      if(L1k > 1)
                      {  // Component of a solution phase
                         // The species is present in a larger or smaller amount than necessary
                           YjCleaned = Yj / exp( MjuDiff );
                      }
                      YjDiff = YjCleaned - Yj;
                      if( fabs( YjDiff ) > CutoffDistortionMBR )
                      {
                          NeedToImproveMassBalance = 1;
                          if( fabs( YjDiff ) > AmountCorrectionThreshold )
                          {   // Correction was too large
                              NeedToImproveMassBalance = 2;
                              // Temporary: only correction of the size of threshold
                              if( YjDiff > 0. )
                                pmp->Y[j] += AmountCorrectionThreshold;
                              else
                                pmp->Y[j] -= AmountCorrectionThreshold;
                          }
                          else {  // Reasonable correction
                              pmp->Y[j] = YjCleaned;
                          }
                      }
                      else {  // Correction that does not affect the mass balance
                          pmp->Y[j] = YjCleaned;
                      }
                      if( pmp->Y[j] < pmp->DcMinM )
                      {  // Corrected amount is too small - zeroed off
                         pmp->Y[j] = 0.;
                         L1kZeroDCs++;
                      }
                  }
               }
               else if( Yj < pmp->DcMinM )
               {  // already zero
                   L1kZeroDCs++;
               }
            }  // for j
            if(( pmp->L1[k] - L1kZeroDCs <= 1 && k < pmp->FIs )
                ||( pmp->L1[k] - L1kZeroDCs == 0 && k >= pmp->FIs ))
            {   Degenerated = true;
                NeedToImproveMassBalance = 1;
            }
       }
       jb+=pmp->L1[k];
   }
   if( NeedToImproveMassBalance )
   { // diagnostics to be implemented
      if( Degenerated && NeedToImproveMassBalance == 1 )
          NeedToImproveMassBalance = -1;
      if( Degenerated && NeedToImproveMassBalance == 2 )
      {  // Diagnostic output here
          NeedToImproveMassBalance = -2;
      }
   }
   return NeedToImproveMassBalance;
}

//====================================================================================
// New simplified PSSC() algorithm   DK 01.05.2010
// PhaseSelection() part only looks for phases to be inserted, also checks if some
// solution phases are unstable. Removal of unstable phases is done afterwards in
// CleanupSpeciation() function.
// As phase stability criterion, uses (log) phase stability (saturation) index
// computed from DualTh activities of components and activity coefficients
// returns 1L if Ok; 0 if one more IPM loop should be done;
//        -1L if 3 loops did not fix the problem
//
long int TMulti::PhaseSelectionSpeciationCleanup( long int &kfr, long int &kur, long int CleanupStatus )
{
    double logSI, PhaseAmount = 0., AmThExp, AmountThreshold = 0.;
    double Yj, YjDiff, YjCleaned=0., MjuPrimal, MjuDual, MjuDiff;
    double CutoffDistortionMBR = 0.1 * pmp->DHBM;
    bool KinConstrDC, KinConstrPh;
    bool MassBalanceViolation = false;
    bool NeedToImproveMassBalance = false;
    long int L1k, L1kZeroDCs, k, j, jb = 0, status,
        DCinserted = 0, DCremoved = 0, PHinserted = 0, PHremoved = 0;
    double MjuDiffCutoff = 1e-3; // InsValue;
    SPP_SETTING *pa = &TProfil::pm->pa;
    if( pa->p.GAS > 1e-6 )
         MjuDiffCutoff = pa->p.GAS;
    AmThExp = (double)abs( pa->p.PRD );
    if( AmThExp && AmThExp < 4.)
    {
        AmThExp = 4.;
    }
    AmountThreshold = pow(10.,-AmThExp);

    kfr = -1; kur = -1;
    (pmp->K2)++;
    for( j=0; j<pmp->L; j++ )
        pmp->XY[j]=pmp->Y[j];    // Storing a copy of the new speciation vector

    PrimalChemicalPotentials( pmp->F, pmp->Y, pmp->YF, pmp->YFA );
    StabilityIndexes( ); // Calculation of phase stability criteria

    for(k=0;k<pmp->FI;k++)
    {
       L1k = pmp->L1[k]; // Number of components in the phase
       KinConstrPh = false;
       for(j=jb; j<jb+L1k; j++)
       {  // Checking if a DC in phase is under kinetic control
//          Yj = pmp->Y[j];
          KinConstrDC = false;
          // Detecting the DC having the non-trivial kinetic constraint
          if( pmp->DUL[j] < 1e6 ) // && Yj >= ( pmp->DUL[j] - pa->p.DKIN ) )
              KinConstrDC = true;
          if( pmp->DLL[j] > 0.0 ) // && Yj <= ( pmp->DLL[j] + pa->p.DKIN ) )
              KinConstrDC = true;
          if( KinConstrDC == true ) // Bug fixed 21.07.2010 DK
              KinConstrPh = true;
       } //  j
       if( pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL
           || KinConstrPh == true )
           goto NextPhase;  // Temporary workaround
       PhaseAmount = pmp->XF[k];
       logSI = pmp->Falp[k];
       if( logSI > -pa->p.DFM && logSI < pa->p.DF && PhaseAmount < pmp->DSM )
       {  // Phase is stable and present in zero or less than DS amount - zeroing off
          bool RemFlagDC = false;
          for(j=jb; j<jb+L1k; j++)
          {
             if( pmp->Y[j] )
             {
                DCremoved++; RemFlagDC = true;
             }
             pmp->Y[j] = 0.;
          }
          if( RemFlagDC == true )
             PHremoved++;
          goto NextPhase;
       }
       if( logSI >= pa->p.DF )  // 2 - INSERTION CASE
       {  // this phase is stable or over-stable
           if( PhaseAmount < pmp->DSM ) // pmp->DFYsM )
           {  // phase appears to be lost - insertion of all components of the phase
               if( L1k > 1 )
                  DC_RaiseZeroedOff( jb, jb+L1k, k );
               else
                  pmp->Y[jb] = pmp->DFYsM; // Spec. value for pure phase insertion
               DCinserted += L1k;
               PHinserted++;
               kfr = k;
               MassBalanceViolation = true;
           } // otherwise (if present), the phase is cleaned up
           goto NextPhase;
       }
       if( logSI <= -pa->p.DFM )  // 3 - ELIMINATION CASE
       {
//         bool Incomplete = false;
           if( PhaseAmount >= pmp->DcMinM )
           {  // this phase is present - checking elimination if unstable
             kur = k;
             if( PhaseAmount <= AmountThreshold )
             {   // can be zeroed off
                for(j=jb; j<jb+L1k; j++)
                {
                    if( pmp->Y[j] >= pmp->DcMinM )
                        DCremoved++;
                    pmp->Y[j] = 0.;
                }
             }
             else { // Phase amount too high - elimination may break the mass balance
                MassBalanceViolation = true;
                for(j=jb; j<jb+L1k; j++)
                {
                    if( pmp->Y[j] >= pmp->DcMinM )
                        DCremoved++;
                    pmp->Y[j] = 0.;
                }
             }
             PHremoved++;
          }
          goto NextPhase;
       }
     NextPhase: jb+=pmp->L1[k];
   } // k
   // First loop over phases finished
   if( pa->p.PRD && MassBalanceViolation == false )  // Cleanup mode (PRD: amount threshold exponent)
   {  //  not done if insertion or elimination of some phases requires another IPM loop
     jb = 0;
     for(k=0;k<pmp->FI;k++)  // Cleanup loop on phases
     {
       L1k = pmp->L1[k]; // Number of components in the phase
       KinConstrPh = false;
       for(j=jb; j<jb+L1k; j++)
       {  // Checking if a DC in phase is under kinetic control
          Yj = pmp->Y[j];
          KinConstrDC = false;
          // Detecting the DC having the non-trivial kinetic constraint
          if( pmp->DUL[j] < 1e6 ) // && Yj >= ( pmp->DUL[j] - pa->p.DKIN ) )
          {
              if( Yj > pmp->DUL[j] )
                  // Fixing a very small component constrained from above
                  pmp->Y[j] = pmp->DUL[j];
              KinConstrDC = true;
          }
          if( pmp->DLL[j] > 0 ) // && Yj <= ( pmp->DLL[j] + pa->p.DKIN ) )
          {
              if( Yj < pmp->DLL[j] )
              // Fixing a small component constrained from below
                 pmp->Y[j] = pmp->DLL[j];
              KinConstrDC = true;
          }
          if( KinConstrDC == true )
          {
             YjDiff = fabs( pmp->Y[j] -Yj );
             if( YjDiff > CutoffDistortionMBR )
                  NeedToImproveMassBalance = true;
             if( YjDiff > AmountThreshold )
                  MassBalanceViolation = true;
             KinConstrPh = true;
          }
       } //  j
       if( KinConstrPh == true ) // || pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL )
          goto NextPhaseC;  // Temporary workaround
       PhaseAmount = pmp->XF[k];                           // bugfix 04.04.2011 DK
       logSI = pmp->Falp[k];  // phase stability criterion
       if( (logSI >= pa->p.DF || logSI <= -pa->p.DFM )
           && !(( pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL) && PhaseAmount < pmp->DSM ))
           goto NextPhaseC;  // Already done in previous part
//       PhaseAmount = pmp->XF[k];                         // bugfix 04.04.2011 DK
       if( logSI > -pa->p.DFM && PhaseAmount >= pmp->DSM )
       { // Cleaning up a phase which is (over)stable and present in mass balance
          bool Degenerated = false;
          L1kZeroDCs = 0;
          for(j=jb; j<jb+L1k; j++)
          {
//                if( pmp->DUL[j] < 1e6 || pmp->DLL[j] > 0 )
//                    continue; // metastability-controlled DCs are already fixed
            Yj = YjCleaned = pmp->Y[j];
            if( Yj >= pmp->DcMinM )
            {  // we check here in the Ls set only, except metastability constraints
               MjuPrimal = pmp->F[j];   // normalized
               MjuDual = pmp->Fx[j]/pmp->RT;
               MjuDiff = MjuPrimal - MjuDual;
               if( fabs( MjuDiff ) > MjuDiffCutoff )
               {
                  YjCleaned = Yj / exp( MjuDiff ); // also applies to a DC in a solution phase
                  if( L1k == 1 )
                  {  // Pure phase
                      if( logSI <= -0.4343*MjuDiffCutoff && YjCleaned < AmountThreshold )
                          YjCleaned = 0.;
                      if( logSI >= pa->p.DF && YjCleaned < pmp->DFYsM )
                      {   // over-stable phase in too small amount - insertion and next IPM loop (experimental)
                          YjCleaned = pmp->DFYsM;
                          kfr = k;
                          MassBalanceViolation = true;
                      }
                  }
               }
             }
             else if( !(pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL)  )
             {  // component is out of the Ls set - primal chemical potential not defined
                  if( L1k == 1 )
                     ; // L1kZeroDCs++; No way to improve for pure substance
                  else  // possibly lost DC in a solution phase (only estimated mole fraction available)
                     YjCleaned = pmp->EMU[j]*pmp->YF[k]; // does not work for adsorption yet
             }
             YjDiff = YjCleaned - Yj;
             if( fabs( YjDiff ) > CutoffDistortionMBR )
             {
                    NeedToImproveMassBalance = true;
                    if( fabs( YjDiff ) > AmountThreshold )
                    {   // Correction was too large - next IPM loop required
                       kur = k;
                       MassBalanceViolation = true;
                      // Provisional: only correction up to the amount threshold
                       if( YjDiff > 0. )
                          pmp->Y[j] += AmountThreshold;
                       else
                          pmp->Y[j] -= AmountThreshold;
                    }
                    else {  // Reasonable correction - no balance violation expected
                      pmp->Y[j] = YjCleaned;
                    }
              }
              else {  // Correction that does not affect the mass balance at all
                  pmp->Y[j] = YjCleaned;
              }
              if( pmp->Y[j] < pmp->DcMinM )
              {  // Corrected amount is too small - DC amount is zeroed off
                  pmp->Y[j] = 0.;
                  DCremoved++;
                  L1kZeroDCs++;
              }

        }  // for j
        if( L1k - L1kZeroDCs <= 1 && L1k > 1 )
        {
           if( L1k - L1kZeroDCs )
              Degenerated = true;
           else
              PHremoved++;
           NeedToImproveMassBalance = true;
        }
        if( L1k - L1kZeroDCs == 0 && L1k == 1 )
        {
           PHremoved++;
           NeedToImproveMassBalance = true;
        }
//        goto NextPhaseC;
     }
     NextPhaseC: jb+=pmp->L1[k];
   } // k
}

#ifndef IPMGEMPLUGIN
    STEP_POINT("PSSC()");
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelectionSpeciationCleanup()"
#endif
#endif
    // Analysis of phase selection and cleanup status
// PZ    // Indicator of PhaseSelection() status (since r1594):
//            0 untouched, 1 phase(s) inserted, 2 insertion done after 5 major IPM loops
// W1     // Indicator of CleanupSpeciation() status (since r1594) 0 untouched,
//           -1 phase(s) removed, 1 some DCs inserted
// K2     // Number of IPM loops performed ( >1 up to 3 because of PhaseSelection() )
    status = 1L;
    CleanupStatus = 0;
    if( !PHinserted )
    {  // No phases were inserted back to mass balance - only cleanup
       status = 1L;
       if( DCinserted )
       {  // some DC in multicomponent phases were restored
          CleanupStatus = 1L;
       }
       else if( NeedToImproveMassBalance )
       {
          CleanupStatus = -1L;
       }
       if( MassBalanceViolation )
       {
           if( pmp->K2 < 6 )
              status = 0L;  // attempt will be done to improve by doing one more IPM loop
           else
              status = -1L;  // five loops done, violation persistent - bail out
       }
       if( PHremoved || DCremoved )
           CleanupStatus = -1L;
    }
    else { // phases were inserted - go to another IPM loop
       if( pmp->K2 < 6 )
          status = 0L;
       else
          status = -1L;
    }
    if( status == -1L )
    {  // changes in Y vector are not accepted - restore (too many IPM loops)
       for(j=0;j<pmp->L;j++)
          pmp->Y[j]=pmp->XY[j];
    }
    return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// calculation of (logarithmic) stability indexes logSI for all phases
//
void TMulti::StabilityIndexes( void )
{
    long int L1k, k, j, jb = 0;
    double ln_ax_dual, gamma_primal, x_estimate, StabIndex, logSI;
    double lnFmol = log( H2O_mol_to_kg );  // may not work with mixed-solvent electrolyte
    double lnPc = 0., Xw = 1., lnXw = 0., lnFugPur=0.; // AqsTail = 0.;

    if( pmp->Pc > 1e-29 )
       lnPc = log( pmp->Pc );
    if( pmp->PHC[0] == PH_AQUEL && pmp->YFA[0] >= pmp->XwMinM  ) // number of moles of solvent
    {
        Xw = pmp->YFA[0] / pmp->YF[0];
//        AqsTail = 1. - Xw;
        lnXw = log( Xw );
    }
    jb=0;
    for(k=0;k<pmp->FI;k++)
    {
       L1k = pmp->L1[k]; // Number of components in the phase
       StabIndex = 0.;
       for(j=jb; j<jb+L1k; j++)
       {  // calculation for all components in all phases
          gamma_primal = pmp->Gamma[j];  // primal (external) activity coefficient
          if( gamma_primal < 1e-33 || gamma_primal > 1e33 )
              gamma_primal = 1.;
          ln_ax_dual = lg_to_ln * pmp->Y_la[j];  // DualTh activity
          if( ln_ax_dual < -777. )
              ln_ax_dual = -777.;
          lnFugPur = pmp->GEX[j];  // Pure fugacity or DQF parameter

          switch( pmp->DCC[j] ) // choice of corrections for estimated mole fractions
          {
             case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                  ln_ax_dual -= lnFmol - lnXw;
                  break;
             case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                  break;
             case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2: case DC_GAS_H2: case DC_GAS_N2:
                  ln_ax_dual -= lnPc + lnFugPur;
                  break;
             case DC_SCP_CONDEN: case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR:
                  ln_ax_dual -= lnFugPur;
                  break;
             case DC_SUR_GROUP:
                  ln_ax_dual -= lnFmol;   // maybe more correction is needed for surface species
                  break;
             case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
             case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
             case DC_SUR_COMPLEX: case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                  ln_ax_dual -= lnFmol;
  //                gamma_primal = exp( pmp->F0[j] );
                 break;
             case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                  ln_ax_dual -= lnFugPur;
  //                gamma_primal = exp( pmp->F0[j] );
                  break;
             default:
                  break; // error in DC class code
          }
          x_estimate = exp( ln_ax_dual )/gamma_primal;   // estimate of DC concentration
          StabIndex += x_estimate;  // Increment to stability index
          pmp->NMU[j] = log( x_estimate );  // may be used for something more constructive
          pmp->EMU[j] = x_estimate;         // stored the estimated mole fraction of phase component
       }  // for j
       logSI = log10( StabIndex );
       if( logSI < -333. )
           logSI = -333.;
       if( fabs( logSI ) < 1e-9 )
           logSI = 0.;
       pmp->Falp[k] = logSI; // NormDoubleRound( logSI, 3 );
       jb += pmp->L1[k];
    }  // for k
}

//=================================================================== old ===========================
// Checking Karpov phase stability criteria Fa for phases and DCs
//  using Selekt2() algorithm by Karpov & Chudnenko (1989)
//  modified by DK in 1995 and in 2007
//  Returns 0, if some phases were inserted and a new IPM loop is needed
//             (up to 3 loops possible);
//          1, if the IPM solution is final and consistent, no phases were inserted
//          -1, if the IPM solution is inconsistent after 3 Selekt2() loops
//  In this case, the index of most problematic phase is passed through kf or
//  ku parameter (parameter value -1 means that no problematic phases were found)
//
long int TMulti::PhaseSelect( long int &kfr, long int &kur, long int CleanupStatus ) //  rLoop )
{
    long int k, j, jb, kf, ku;
    double F1, F2, *F0; // , sfactor;
    SPP_SETTING *pa = &TProfil::pm->pa;
int rLoop = CleanupStatus;
rLoop = -1;
//    sfactor = calcSfactor();
    KarpovCriterionPH( );  // calculation of Karpov phase stability criteria (in pmp->Falp)
    F0 = pmp->Falp;

    (pmp->K2)++;
    kf = -1; ku = -1;  // Index for phase diagnostics
    F1 = pa->p.DF;  // Fixed 29.10.2007  DK
    F2 = -pa->p.DFM;  // Meaning of DFM changed 02.11.2007

    for(k=0;k<pmp->FI;k++)
    {
        if( F0[k] > F1 && pmp->YF[k] < pmp->DSM ) //  pa->p.DS )  // < pmp->lowPosNum?
        {            // stable phase not in mass balance - to be inserted
            F1=F0[k];
            kf=k;
        }
        if( F0[k] < F2 && pmp->YF[k] >= pmp->DSM ) // pa->p.DS )  // Fixed 2.11.2007
        {            // unstable phase in mass balance - to be excluded
            F2=F0[k];
            ku=k;
        }
    }
kfr = kf;
kur = ku;
    if( kfr < 0 && kur < 0 )
    {    // No phases to insert/exclude or no Fa distortions found
          // Successful end of iterations of SELEKT2()
        return 1L;
    }

    if( (F2 < -pa->p.DFM ) && ( ku >= 0 ) )
    {
        if( pmp->K2 > 4 ) // Three Selekt2() loops have already been done!
                return -1L;   // Persistent presence of unstable phase(s) - bad system!

        // Excluding problematic phases
        do
        {  // excluding all phases with  F2 < DF*sfactor
            for( jb=0, k=0; k < ku; k++ )
                 jb += pmp->L1[k];

            DC_ZeroOff( jb, jb+pmp->L1[ku], ku ); // Zeroing the phase off
            pmp->FI1--;
            // find a new phase to exclude, if any exists
            F2= -pa->p.DFM;
            ku = -1;
            for( k=0; k<pmp->FI; k++ )
                if( F0[k] < F2 && pmp->YF[k] >= pmp->DSM )
                {
                    F2=F0[k];
                    ku=k;
                }
        }
        while( ( F2 <= -pa->p.DFM ) && ( ku >= 0 ) );
    } //if ku

    // Inserting problematic phases
    if( F1 > pa->p.DF && kf >= 0 )
    {
        if( pmp->K2 > 4 )
           return -1L;   // Persistent absence of stable phase(s) - bad system!

        // There is a phase for which DF*sfactor threshold is exceeded
        do
        {   // insert this phase and set Y[j] for its components
            // with account for asymmetry and non-ideality
            for( jb=0, k=0; k < kf; k++ )
                 jb += pmp->L1[k];

            DC_RaiseZeroedOff( jb, jb+pmp->L1[kf], kf );

            pmp->FI1++;  // check phase rule

            if( pmp->FI1 >= pmp->NR+1 )
               break;   // No more phases can be inserted

            // find a new phase to insert, if any exists
            F1= pmp->lowPosNum; // was about 1e-16
            kf = -1;
            for( k=0; k<pmp->FI; k++ )
                if( F0[k] > F1 && pmp->YF[k] < pmp->DSM )
                {
                    F1=F0[k];
                    kf=k;
                }
        }
        while( F1 > pa->p.DF && kf >= 0 );
        // end of insertion cycle
    } // if kf changed SD 03/02/2009
// Raise zeros in DC amounts in phases-solutions in the case if some phases
// were inserted or excluded                          - experimental option!
//        double RaiseZeroVal = pmp->DHBM;  // Added 29.10.07  by DK
        jb=0;
        for(k=0;k<pmp->FIs;k++)
        {
            if( ( pmp->YF[k] >= pmp->DSM ) || ( pmp->pNP && rLoop < 0 ) ) // Only in phase present in mass balance!
            {                            // (acc. to definition of L_S set) PIA only if initial!
                 pmp->YF[k]=0.;
                 for(j=jb;j<jb+pmp->L1[k];j++)
                 {
                    if( pmp->Y[j] < min( pmp->lowPosNum, pmp->DcMinM ) )  // fixed 30.08.2009
                        pmp->Y[j] = RaiseDC_Value( j ); // bugfix 29.10.07
                    pmp->YF[k] += pmp->Y[j]; // calculate new amounts of phases
                 }
            }
            jb+=pmp->L1[k];
        }
    // } SD 03/02/2009
#ifndef IPMGEMPLUGIN
    STEP_POINT("Select2()");
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelection"
#endif
#endif
    if( pmp->K2 > 1 )
    { // more then the first step - but the IPM solution has not been improved
//      double RaiseZeroVal = pmp->DHBM*0.1;   // experimental
       for(j=0;j<pmp->L;j++)
          if( fabs(pmp->Y[j]- pmp->XY[j]) > RaiseDC_Value( j ) ) //
               goto S6;
       // pmp->PZ=2; // No significant change has been done by Selekt2()
       return 1L;
    }
S6: // copy of X vector has been changed by Selekt2() algorithm - store
    for(j=0;j<pmp->L;j++)
        pmp->XY[j]=pmp->Y[j];

    return 0L;  // Another loop is needed
}

// New function to improve on raising zero values in PhaseSelect() and after SolveSimplexLPP()
double TMulti::RaiseDC_Value( const long int j )
{
        double RaiseZeroVal = pmp->DFYsM;

        switch(pmp->DCC[j] )
        {
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
    case DC_AQ_SURCOMP:	RaiseZeroVal = pmp->DFYaqM;
                                                break;
    case DC_SOL_IDEAL:
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:	RaiseZeroVal = pmp->DFYidM;
                                                break;
    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM: RaiseZeroVal = pmp->DFYwM;
                                                break;
    case DC_SOL_MINOR:	RaiseZeroVal = pmp->DFYhM;
                                                break;
    case DC_SOL_MAJOR: RaiseZeroVal = pmp->DFYrM;
                                                break;
        // adsorption
    case DC_SSC_A0:    case DC_SSC_A1:    case DC_SSC_A2:    case DC_SSC_A3:    case DC_SSC_A4: // obsolete
    case DC_WSC_A0:	   case DC_WSC_A1:    case DC_WSC_A2:    case DC_WSC_A3:    case DC_WSC_A4: // obsolete
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:		RaiseZeroVal = pmp->DFYaqM;
                                                break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER: RaiseZeroVal = pmp->DFYrM;
                                                 break;
    case DC_SCP_CONDEN:	 RaiseZeroVal = pmp->DFYcM;
                                                 break;
    default:
                 break;
        }
        return RaiseZeroVal;
}

//--------------------- End of ipm_chemical.cpp ---------------------------
