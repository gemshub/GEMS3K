//-------------------------------------------------------------------
// $Id: ipm_chemical3.cpp 690 2006-03-29 07:10:23Z gems $
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
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//
#include <math.h>
#include "m_param.h"
#include "s_fgl.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of surface charge densities on multi-surface sorption phase
//
void TMulti::IS_EtaCalc()
{
    int k, i, ist, isp, j=0, ja;
    double XetaS=0., XetaW=0.,  Ez, CD0, CDb;
    SPP_SETTING *pa = &TProfil::pm->pa;

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
                (pmp->PHC[k] == PH_AQUEL && ( pmp->X[pmp->LO] <= pa->p.XwMin
                 || pmp->XF[k] <= pmp->DHBM ) )
             || (pmp->PHC[k] == PH_SORPTION && pmp->XF[k] <= pa->p.ScMin) )
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
            case DC_AQ_ELECTRON:    case DC_AQ_PROTON:    case DC_AQ_SPECIES:
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
                CD0 = (double)pmp->MASDJ[ja][PI_CD0];
                    // species charge that goes into 0 plane
                CDb = (double)pmp->MASDJ[ja][PI_CDB];
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
// Linking DOD for executing Phase mixing model scripts
//
void TMulti::pm_GC_ods_link( int k, int jb, int jpb, int jdb, int ipb )
{

#ifndef IPMGEMPLUGIN
    ErrorIf( k < 0 || k >= pmp->FIs , "GammaCalc", "Invalid link: k=0||>FIs" );
    aObj[ o_nsmod].SetPtr( pmp->sMod[k] );
    aObj[ o_nncp].SetPtr( pmp->LsMod+k*3 );
    aObj[ o_nncd].SetPtr( pmp->LsMdc+k );
    aObj[ o_ndc].SetPtr(  pmp->L1+k );
    aObj[ o_nez].SetPtr( pmp->EZ+jb );
    aObj[o_nez].SetN(  pmp->L1[k]);
    aObj[ o_npcv].SetPtr( pmp->PMc+jpb );
    aObj[o_npcv].SetDim( pmp->LsMod[k*3], pmp->LsMod[k*3+2]);
    //  Object for indexation of interaction parameters
    aObj[ o_nu].SetPtr( pmp->IPx+ipb ); // added 07.12.2006  KD
    aObj[o_nu].SetDim( pmp->LsMod[k*3], pmp->LsMod[k*3+1]);
    //
    aObj[ o_ndcm].SetPtr( pmp->DMc+jdb );
    aObj[o_ndcm].SetDim( pmp->L1[k], pmp->LsMdc[k] );
    aObj[ o_nmvol].SetPtr( pmp->Vol+jb );
    aObj[o_nmvol].SetN( pmp->L1[k]);
    aObj[ o_nppar].SetPtr(pmp->VL+jb );
    aObj[o_nppar].SetN(  pmp->L1[k]);
//    aObj[ o_ngtn].SetPtr( pmp->G0+jb );
 aObj[ o_ngtn].SetPtr( pmp->GEX+jb );     // changed 05.12.2006 by DK
//  aObj[ o_ngtn].SetPtr( pmp->lnGmM+jb );     // experimental: 21.02.2007 by DK   
    aObj[o_ngtn].SetN( pmp->L1[k] );
    aObj[ o_ngam].SetPtr( pmp->Gamma+jb ); // Gamma calculated
    aObj[o_ngam].SetN( pmp->L1[k] );
    aObj[ o_nlngam].SetPtr( pmp->lnGam+jb ); // ln Gamma calculated
    aObj[o_nlngam].SetN( pmp->L1[k]);
    aObj[ o_nas].SetPtr(  pmp->A+pmp->N*jb );
    aObj[o_nas].SetDim(  pmp->L1[k], pmp->N );
    aObj[ o_nxa].SetPtr(  pmp->XF+k );
    aObj[ o_nxaa].SetPtr(  pmp->XFA+k );
    if( pmp->FIat > 0 )
    {
        aObj[ o_nxast].SetPtr( pmp->XFTS[k] );
        aObj[ o_nxcec].SetPtr( pmp->MASDT[k] );
    }
    else
    {
        aObj[ o_nxast].SetPtr( 0 );
        aObj[ o_nxcec].SetPtr( 0 );
    }
    /* */
    aObj[ o_nbmol].SetPtr( pmp->FVOL+k );  // phase volume
    aObj[ o_nxx].SetPtr(  pmp->X+jb );
    aObj[o_nxx].SetN( pmp->L1[k]);
    aObj[ o_nwx].SetPtr(  pmp->Wx+jb );
    aObj[o_nwx].SetN( pmp->L1[k]);
    aObj[ o_nmju].SetPtr( pmp->Fx+jb );
    aObj[o_nmju].SetN( pmp->L1[k]);
    aObj[ o_nqp].SetPtr( pmp->Qp+k*QPSIZE );
    aObj[ o_nqd].SetPtr( pmp->Qd+k*QDSIZE );   // Fixed 7.12.04 by KD
#endif
}

// Calculation of smoothing factor for high non-ideality systems
// Returns current value of this factor
//
double TMulti::TinkleSupressFactor( double ag, int ir)
{
    double TF, dg, rg=0.0, irf;
    int Level, itqF, itq;

    Level = pmp->pRR1;
    irf = (double)ir;
    dg = TProfil::pm->pa.p.DGC;
    if( dg > -0.1 )
    {
        if(ag>1) ag=1;
        if(ag<0.1) ag=0.1;
        if(dg>0.15) dg=0.15;
        if( dg <= 0.0 )
          TF = ag;
        else
          TF = ag * ( 1 - pow(1-exp(-dg*irf),60.));
        Level =  (int)rg;
    }
    else
    {  // Old smoothing for solid solutions: -1.0 < pa.p.DGC < 0
        itq = ir/60;
        dg = fabs( dg );
        itqF = ir/(60/(itq+1))-itq;  // 0,1,2,4,5,6...
        if( itqF < Level )
            itqF = Level;
        TF = ag * pow( dg, itqF );
        Level = itqF;
    }
    pmp->pRR1 = (short)Level;
    return TF;
}

static double ICold=0.;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Main call point for calculation of DC activity coefficients (lnGam vector)
// Controls various built-in models, as well as generic Phase script calculation
//
// LinkMode is a parameter indicating the status of Gamma calculations:
//   LINK_TP_MODE  - calculation of equations depending on TP only;
//   LINK_UX_MODE  - calculation of equations depending on current
//                   IPM approximation of the equilibrium state;
//   LINK_FIA_MODE - calculation of Gammas on the initial approximation (FIA).
//    Added 13.03.2008 by DK: returns int value showing (if not 0) 
//    that some extreme values were obtained for some SACTs, PSIs,
//    or activity coefficients (for detecting bad PIA case). 0 if Ok.
//
int 
TMulti::GammaCalc( int LinkMode  )
{
    int k, j, jb, je=0, jpb, jpe=0, jdb, jde=0, ipb, ipe=0;
    char *sMod;  
    int statusGam=0, statusGC=0, statusSACT=0;
    double LnGam, pmpXFk;
    SPP_SETTING *pa = &TProfil::pm->pa;

//  high-precision IPM-2 debugging
//   if(pmp->PZ && pmp->W1 > 1 )
//     goto END_LOOP;
//   pmp->FitVar[3] = TinkleSupressFactor( pmp->FitVar[4], pmp->IT );

    // calculating concentrations of species in multi-component phases
    switch( LinkMode )
    {
    case LINK_FIA_MODE: // Initial approximation mode
        if( pmp->LO )
        {
            pmp->GWAT = 55.51;
            pmp->XF[0] = pmp->GWAT;
            for( j=0; j<pmp->L; j++ )
                if( pmp->X[j] < pmp->lowPosNum )
                    pmp->X[j] = pa->p.DFYaq;
            ConCalc( pmp->X, pmp->XF, pmp->XFA );
// if( pmp->E )
//    IS_EtaCalc( );
//          pmp->IC = max( pmp->MOL, pmp->IC );
            pmp->IC = 0.0;  // Important for the simplex FIA reproducibility
            if( pmp->E && pmp->FIat > 0 )
            {
               for( k=0; k<pmp->FIs; k++ )
               {
                 int ist;
                 if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
                   for( ist=0; ist<pmp->FIat; ist++ ) // loop over surface types
                   {
                     pmp->XetaA[k][ist] = 0.0;
                     pmp->XetaB[k][ist] = 0.0;
                     pmp->XpsiA[k][ist] = 0.0;
                     pmp->XpsiB[k][ist] = 0.0;
                     pmp->XpsiD[k][ist] = 0.0;
                     pmp->XcapD[k][ist] = 0.0;
                   }  // ist
               }  // k
            } // FIat
        }   // LO
        // Cleaning vectors of  activity coefficients
        for( j=0; j<pmp->L; j++ )
        {
            if( pmp->lnGmf[j] )
                pmp->lnGam[j] = pmp->lnGmf[j]; // setting up fixed act.coeff. for simplex IA
            else pmp->lnGam[j] = 0.;   // + pmp->lnGmm[j];
            pmp->Gamma[j] = 1.;
        }
        break;      // added as bugfix 04.03.2008  DK
    case LINK_TP_MODE:  // Built-in functions depending on T,P only
         pmp->FitVar[3] = 1.0;  // resetting the IPM smoothing factor
         
         for( k=0; k<pmp->FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += pmp->L1[k];
            if( pmp->L1[k] == 1 )
                continue;
            // Indexes for extracting data from IPx, PMc and DMc arrays
            ipb = ipe;
            ipe += pmp->LsMod[k*3]*pmp->LsMod[k*3+1];
            jpb = jpe;
            jpe += pmp->LsMod[k*3]*pmp->LsMod[k*3+2];
            jdb = jde;
            jde += pmp->LsMdc[k]*pmp->L1[k];
            sMod = pmp->sMod[k];

   		    double nxk = 1./pmp->L1[k];
            for( j= jb; j<je; j++ )
    		{  
if(pmp->XF[k] < pmp->lowPosNum )   // workaround 10.03.2008 DK 
   pmp->Wx[j] = nxk;  // need this eventually to avoid problems with zero mole fractions
    		   pmp->GEX[j] =0.0;  // cleaning GEX in TP mode!
    		   pmp->lnGmo[j] = pmp->lnGam[j]; // saving activity coefficients in TP mode
       	    }	  
            if( sMod[SGM_MODE] != SM_STNGAM )
                continue;  // The switch below is for built-in functions only!
            
            // the following section probably needs to be re-written to allow more flexible
            // combinations of fluid models for pure gases with gEX mixing models, 
            // scheme should probably be the same as in LINK_UX_MODE, 03.06.2008 (TW)
            switch( pmp->PHC[k] )
            {
               case PH_LIQUID:
               case PH_SINCOND:
               case PH_SINDIS:
               case PH_HCARBL:
               case PH_SIMELT:
                       SolModParPT( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] ); // new solution models (TW, DK 2007)
                    break;
               case PH_GASMIX:
               case PH_PLASMA:
               case PH_FLUID:
                     if( sMod[SPHAS_TYP] == SM_CGFLUID )
                     {
                       CGofPureGases( jb, je, jpb, jdb, k ); // CG2004 pure gas
                       break;
                     }
                     if( sMod[SPHAS_TYP] == SM_PRFLUID )
                       PRSVofPureGases( jb, je, jpb, jdb, k, ipb ); // PRSV pure gas
                     break;
               default: break;
            }
        } // k
        break;
    case LINK_UX_MODE: 
pmp->FitVar[3] = TinkleSupressFactor( pmp->FitVar[4], pmp->IT );  // Getting actual smoothing parameter    	
    	// calculating DC concentrations after this IPM iteration
        ConCalc( pmp->X, pmp->XF, pmp->XFA );
        // cleaning activity coefficients
        for( j=0; j<pmp->L; j++ )
        {
            pmp->lnGam[j] = 0.;
            pmp->Gamma[j] = 1.;
        }
        if( pmp->E && pmp->LO ) // checking electrostatics
        {
          IS_EtaCalc();  //  calculating charges and charge densities 
          if( pmp->FIat > 0 )
             for( k=0; k<pmp->FIs; k++ )
             {
               if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
               {  int ist;
                  for( ist=0; ist<pmp->FIat; ist++ ) // loop over surface types
                  {
                     pmp->XpsiA[k][ist] = 0.0;        // cleaning Psi before GouyChapman() 
                     pmp->XpsiB[k][ist] = 0.0;
                     pmp->XpsiD[k][ist] = 0.0;
                  }  // ist
                }
             }  // k
         } // pmp->E
        break;
    default:
        Error("GammaCalc","Invalid Link Mode.");
    }

    jpe=0; jde=0; ipe=0; je=0;
    for( k=0; k<pmp->FI; k++ )
    { // loop on phases
        jb = je;
        je += pmp->L1[k];
        if( pmp->L1[k] == 1 )
            goto END_LOOP;
        sMod = pmp->sMod[k];
        if( sMod[SGM_MODE] == SM_IDEAL )
            goto END_LOOP;
        pmpXFk = 0.;  // Added 07.01.05 by KD
        for( j = jb; j < je; j++ )
            pmpXFk += pmp->X[j];
        if( pmp->XF[k] < pmp->DSM ) // Bugfix by KD 09.08.2005 (bug report Th.Matschei)
            pmp->XF[k] = pmpXFk;
        // Indexes for extracting data from IPx, PMc and DMc arrays
        ipb = ipe;                  // added 07.12.2006 by KD
        ipe += pmp->LsMod[k*3]*pmp->LsMod[k*3+1];
        jpb = jpe;
        jpe += pmp->LsMod[k*3]*pmp->LsMod[k*3+2];  // Changed 07.12.2006  by KD
        jdb = jde;
        jde += pmp->LsMdc[k]*pmp->L1[k];
      
   if( LinkMode == LINK_UX_MODE && sMod[SGM_MODE] == SM_STNGAM )
   {                             // check that SGM_MODE for adsorption is not SM_IDEAL in Phase records!
        switch( pmp->PHC[k] )
        {  // calculating activity coefficients using built-in functions
          case PH_AQUEL:   // DH III variant consistent with HKF
             if( pmpXFk > pa->p.XwMin && pmp->X[pmp->LO] > pmp->lowPosNum*1e3 && pmp->IC > pa->p.ICmin )
             {
                switch( sMod[SPHAS_TYP] )
                {
                  case SM_AQDH3:
                       DebyeHueckel3Karp( jb, je, jpb, jdb, k );
                          break;
                  case SM_AQDH2:
                       DebyeHueckel2Kjel( jb, je, jpb, jdb, k );
                          break;
                  case SM_AQDH1:
                       DebyeHueckel1LL( jb, je, k );
                          break;
                  case SM_AQDHH:
                       DebyeHueckel3Hel( jb, je, jpb, jdb, k );
                          break;
                  case SM_AQDAV:
                       Davies03temp( jb, je, k );
                          break;
                  case SM_AQSIT:  // SIT - under construction
                       SIT_aqac_PSI( jb, je, jpb, jdb, k, ipb );
                          break;
                  default:
                          break;
                }
                ICold = pmp->IC;
             }
             goto END_LOOP;
             break;
          case PH_GASMIX:
          case PH_PLASMA:
          case PH_FLUID:
             if( pmpXFk > pmp->DSM && pmp->XF[k] > pa->p.PhMin )
             {
                 if( sMod[SPHAS_TYP] == SM_CGFLUID )
                    ChurakovFluid( jb, je, jpb, jdb, k );
                 if( sMod[SPHAS_TYP] == SM_PRFLUID )
                    PRSVFluid( jb, je, jpb, jdb, k, ipb );
                    // Added by Th.Wagner and DK on 20.07.06
             }
             goto END_LOOP;             
             break;
         case PH_LIQUID:
         case PH_SIMELT:
         case PH_SINCOND:
         case PH_SINDIS:
         case PH_HCARBL:  // solid and liquid mixtures
             if( pmpXFk > pmp->DSM )
             {
                switch( sMod[SPHAS_TYP] )
                {
                  case SM_REDKIS:  // left for compatibility with old projects
                       RedlichKister( jb, je, jpb, jdb, k );
                          break;
                  case SM_MARGB:   // left for compatibility with old projects
                       MargulesBinary( jb, je, jpb, jdb, k );
                          break;
                  case SM_MARGT:   // left for compatibility with old projects
                       MargulesTernary( jb, je, jpb, jdb, k );
                          break;
                  case SM_GUGGENM: // Redlich-Kister model (multicomponent), 2007 (TW)
                       SolModActCoeff( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] );
                          break;
                  case SM_VANLAAR: // VanLaar model (multicomponent), 2007 (TW)
                       SolModActCoeff( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] );
                          break;
                  case SM_REGULAR: // Regular model (multicomponent), 2007 (TW)
                       SolModActCoeff( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] );
                          break;
                  case SM_NRTLLIQ: // NRTL model (multicomponent), 03.06.2007 (TW)
                       SolModActCoeff( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] );
                          break;
                  case SM_WILSLIQ: // Wilson model (multicomponent), 09.06.2007 (TW)
                       SolModActCoeff( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP] );
                          break;              	  	
                  default:
                          break;
                }
             }
             goto END_LOOP;
             break;
        case PH_POLYEL:  // PoissonBoltzmann( q, jb, je, k ); break;
        case PH_SORPTION: // electrostatic potenials from Gouy-Chapman eqn
            if( pmp->PHC[0] == PH_AQUEL && pmpXFk > pmp->DSM
                && (pmp->XFA[0] > pmp->lowPosNum && pmp->XF[0] > pa->p.XwMin ))
            {
//              ConCalc( pmp->X, pmp->XF, pmp->XFA  );  Debugging
                    if( pmp->E )
                    {
//                       IS_EtaCalc();
                       statusGC = GouyChapman( jb, je, k );
                    // PoissonBoltzmann( q, jb, je, k )
                    }
        // Calculating surface activity coefficient terms
                    statusSACT = SurfaceActivityCoeff(  jb, je, jpb, jdb, k );
            }
//            if( sMod[SGM_MODE] == SM_IDEAL )
//                goto END_LOOP;
            break;
         default:
            goto END_LOOP;
       } // end switch 
   }  // end if LinkMode == LINK_UX_MODE 


#ifndef IPMGEMPLUGIN
// This part running Phase math scripts is not used in standalone GEMIPM2K
        // Link DOD and set sizes of work arrays
        pm_GC_ods_link( k, jb, jpb, jdb, ipb );
        pmp->is=0;
        pmp->js=0;
        pmp->next=1;

        switch( LinkMode )
        { // check the calculation mode
        case LINK_TP_MODE: // running TP-dependent scripts
            if(( sMod[SPHAS_DEP] == SM_TPDEP || sMod[SPHAS_DEP] == SM_UXDEP ) && qEp[k].nEquat() )
            {	// Changed on 26.02.2008 to try TW DQF scripts - DK
        	      qEp[k].CalcEquat();
            }
        	if((sMod[DCOMP_DEP] == SM_TPDEP || sMod[DCOMP_DEP] == SM_UXDEP) && qEd[k].nEquat() )
            {
                switch( sMod[DCE_LINK] )
                {
                case SM_PUBLIC:  // one script for all species
                    for( pmp->js=0, pmp->is=0; pmp->js<pmp->L1[k]; pmp->js++ )
                        qEd[k].CalcEquat();
                    break;
                case SM_PRIVATE_: // separate group of equations per species
                    qEd[k].CalcEquat();
                    break;
                }
            }           
        	break;
        case LINK_FIA_MODE: // cold start approximation
            goto END_LOOP;
        	break;
        case LINK_UX_MODE:  // the model is dependent on current concentrations on IPM iteration    	
            switch( pmp->PHC[k] )
            {  //
              case PH_AQUEL:  
                 if(!(pmpXFk > pa->p.XwMin && pmp->X[pmp->LO] > pmp->lowPosNum*1e3 && pmp->IC > pa->p.ICmin ))
                	 goto END_LOOP; 
                 break;
              case PH_GASMIX:
              case PH_PLASMA:
              case PH_FLUID:
                 if( !(pmpXFk > pmp->DSM && pmp->XF[k] > pa->p.PhMin))
                     goto END_LOOP;             
                 break;
             case PH_LIQUID:
             case PH_SIMELT:
             case PH_SINCOND:
             case PH_SINDIS:
             case PH_HCARBL:  // solid and liquid mixtures
                 if( !(pmpXFk > pmp->DSM) )
                     goto END_LOOP;
                 break;
            case PH_POLYEL:  // PoissonBoltzmann( q, jb, je, k ); break;
            case PH_SORPTION: // electrostatic potenials from Gouy-Chapman eqn
                if( !(pmp->PHC[0] == PH_AQUEL && pmpXFk > pmp->DSM
                    && (pmp->XFA[0] > pmp->lowPosNum && pmp->XF[0] > pa->p.XwMin )))
                    goto END_LOOP;
                break;
             default:
                goto END_LOOP;
           } // end switch         	
        	
            if( sMod[SPHAS_DEP] == SM_UXDEP && qEp[k].nEquat() )
                // Equations for the whole phase
                qEp[k].CalcEquat();
            if( sMod[DCOMP_DEP] == SM_UXDEP && qEd[k].nEquat() )
            {  // Equations for species
                switch( sMod[DCE_LINK] )
                {
                case SM_PUBLIC:  // one script for all species
                    for( pmp->js=0, pmp->is=0; pmp->js<pmp->L1[k]; pmp->js++ )
                        qEd[k].CalcEquat();
                    break;
                case SM_PRIVATE_:  // separate group of equations for each species
                    qEd[k].CalcEquat();
                    break;
                }
            }
            if( pmp->PHC[k] == PH_AQUEL )
                ICold = pmp->IC;
            break;
        default:
            Error("GammaCalc","Invalid LinkMode 2");
        } // end switch
#endif

END_LOOP:
        if( LinkMode == LINK_TP_MODE )  // TP mode - added 04.03.2008 by DK
        {
        	for( j=jb; j<je; j++ )
        	{
if( pmp->XF[k] < pmp->lowPosNum )   // workaround 10.03.2008 DK
	pmp->Wx[j] = 0.0;               //
        	   LnGam = pmp->lnGmo[j];
               pmp->lnGam[j] = LnGam;
               if( fabs( LnGam ) > 1e-9 && fabs( LnGam ) < 84. )  
                    pmp->Gamma[j] = exp( LnGam );
               else pmp->Gamma[j] = 1.0;
               pmp->F0[j] = Ej_init_calc( 0.0, j, k );
               pmp->G[j] = pmp->G0[j] + pmp->GEX[j] + pmp->F0[j];       		
        	}	
        }
        else 
        { // Real mode for activity coefficients
          for( j=jb; j<je; j++ )
          {
            if( fabs( 1.0-pmp->Gamma[j] ) > 1e-9 
            		&& pmp->Gamma[j] > 3.3e-37 && pmp->Gamma[j] < 3.03e+36 )   // > 1e-35 before 26.02.08
                pmp->lnGam[j] += log( pmp->Gamma[j] );
            LnGam = pmp->lnGam[j];
            pmp->lnGmo[j] = LnGam;
            if( fabs( LnGam ) > 1e-9 && fabs( LnGam ) < 84. )   // before 26.02.08: < 42.
                pmp->Gamma[j] = exp( LnGam );
            else pmp->Gamma[j] = 1.0;
            pmp->F0[j] = Ej_init_calc( 0.0, j, k );
            pmp->G[j] = pmp->G0[j] + pmp->GEX[j] + pmp->F0[j];
          }
        }   
    }  // k - end loop over phases
    //  if( wn[W_EQCALC].status )
    //  aSubMod[MD_EQCALC]->ModUpdate("PM_ipms   EqCalc() converged");
    if( statusGC )
        	return statusGC;     
    if( statusSACT )
    	return statusSACT; 
    return statusGam; 
}

// ----------------------------------------------------------------------------
// Built-in functions for activity coefficients
//
// aqueous electrolyte
// SIT NEA PSI (not yet official)
//
void
TMulti::SIT_aqac_PSI( int jb, int je, int jpb, int jdb, int k, int ipb )
{

    int j, icat, ian, /*ic, ia,*/  index1, index2, ip, NComp, NPar, NPcoef, MaxOrd;
    short *aIPx;
    float *aIPc, *aDCc;
    double T, A, B, I, sqI, bgi=0, Z2, lgGam, SumSIT;
//    float nPolicy;

    I= pmp->IC;
    if( I <  TProfil::pm->pa.p.ICmin )
        return;
    T = pmp->Tc;
    A = 1.82483e6 * sqrt( pmp->denW ) / pow( T*pmp->epsW, 1.5 );
    B = 50.2916 * sqrt( pmp->denW ) / sqrt( T*pmp->epsW );

//    molt = ( pmp->XF[0]-pmp->XFA[0] )*1000./18.01528/pmp->XFA[0]; // tot.molality
    sqI = sqrt( I );

    ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "SIT",
        "Error: A,B were not calculated - no values of RoW and EpsW !" );

    NComp = pmp->L1[k];          // Number of components in the phase
    NPar = pmp->LsMod[k*3];      // Number of interaction parameters
    NPcoef = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    MaxOrd =  pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)

    // These pointers provide direct access to parts of MULTI arrays related to this phase!
    aIPx = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
                          // -> NPar x MaxOrd   added 07.12.2006   KD
    aIPc = pmp->PMc+jpb;    // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    aDCc = pmp->DMc+jdb;    // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    // aWx = pmp->Wx+jb;       // End member mole fractions
    // alnGam = pmp->lnGam+jb; // End member ln activity coeffs


    // Calculation of EDH equation
//  bgi = bg;
    ian= -1;
    icat = -1;
    for( j=jb; j<je; j++ )
    {
// Determining the index of cation or anion
      if( pmp->EZ[j] < 0 )
          ian = j-jb; // ian++;
      else if( pmp->EZ[j] > 0 )
          icat = j-jb; // icat++;
//      else ;

      if( pmp->EZ[j] )
      {    // Charged species : calculation of the DH part
           Z2 = pmp->EZ[j]*pmp->EZ[j];
           lgGam = ( -A * sqI * Z2 ) / ( 1. + 1.5 * sqI );  // B * 4.562 = 1.5 at 25 C

// Calculation of SIT sum - new variant
           SumSIT = 0.;
           if( pmp->EZ[j] > 0 )  // cation
           {
              for( ip=0; ip<NPar; ip++ )
              {
                 index1 = aIPx[ip*MaxOrd];
                 if( index1 != icat )
                    continue;
                 index2 = aIPx[ip*MaxOrd+1];
                 SumSIT += (double)aIPc[ip*NPcoef]   // epsilon
                        * pmp->Y_m[jb+index2];
              }
           }
           else {   // anion
              for( ip=0; ip<NPar; ip++ )
              {
                 index2 = aIPx[ip*MaxOrd+1];
                 if( index2 != ian )
                    continue;
                 index1 = aIPx[ip*MaxOrd];  // index of cation
                 SumSIT += (double)aIPc[ip*NPcoef]  // epsilon
                         * pmp->Y_m[jb + index1];
              }
           }
           lgGam += SumSIT;
              // Calculation of SIT sums - old variant
/*           SumSIT = 0.;
           if( pmp->EZ[j] > 0 )
           {       // this is a cation
              for( ia=0; ia<pmp->sitNan; ia++ )
                 SumSIT += (double)pmp->sitE[ icat*pmp->sitNan + ia ]
                        * I * pmp->Y_m[pmp->sitXan[ia]];
              lgGam += SumSIT;
           }
           else {  // this is an anion
              for( ic=0; ic<pmp->sitNcat; ic++ )
                 SumSIT += (double)pmp->sitE[ ic*pmp->sitNan + ian ]
                        * I * pmp->Y_m[pmp->sitXcat[ic]];
              lgGam += SumSIT;
           }
*/
      }
      else { // Neutral species
         if( pmp->DCC[j] != DC_AQ_SOLVENT ) // common salting-out coefficient ??
               lgGam = bgi * I;
            else // water-solvent - a0 - osmotic coefficient
               lgGam = 0.;
      }
      pmp->lnGam[j] = lgGam * lg_to_ln;
    } // j
//    if( ++icat != pmp->sitNcat || ++ian != pmp->sitNan )
//       Error( "SITgamma",
//          "Inconsistent numbers of cations and anions in gamma calculation" );
}

//----------------------------------------------------------------------------
// Extended Debye-Hueckel (EDH) model with a common ion-size parameter
//
void
TMulti::DebyeHueckel3Hel( int jb, int je, int jpb, int, int k )
{
    int j;
    double T, A, B, a0, a0c, I, sqI, bg, bgi, Z2, lgGam; //  molt;
    double Xw, Xaq, Nw, molT, Lgam, lnwxWat, lnGam;
    double Lam, SigTerm, Phi, lnActWat;
    float nPolicy;

    I= pmp->IC;
    if( I < TProfil::pm->pa.p.ICmin )
        return;
    T = pmp->Tc;
    A = (double)pmp->PMc[jpb+0];
    B = (double)pmp->PMc[jpb+1];
    bg = pmp->FitVar[0];   // Changed 07.06.05 for T,P-dep. b_gamma in DHH
//    bg = pmp->PMc[jpb+5];     // May be inappropriate for GEMIPM2K !
    a0c = (double)pmp->PMc[jpb+6];
    nPolicy = pmp->PMc[jpb+7];

    Xaq = pmp->XF[k]; // Mole amount of the whole aqueous phase
    Xw = pmp->XFA[k]; // Mole amount of water-solvent
    Nw = 1000./18.01528;
    molT = (Xaq-Xw)/(Nw/Xw);
//    molt = ( pmp->XF[0]-pmp->XFA[0] )*1000./18.01528/pmp->XFA[0]; // tot.molality
    lnwxWat = log(Xw/Xaq);
    sqI = sqrt( I );

#ifndef IPMGEMPLUGIN
    if( fabs(A) < 1e-9 )
    {
        A = 1.82483e6 * sqrt( (double)(tpp->RoW) ) /
            pow( T*(double)(tpp->EpsW), 1.5 );
//        pmp->PMc[jpb+0] = A;
    }
    if( fabs(B) < 1e-9 )
    {
        B = 50.2916 * sqrt( (double)(tpp->RoW) ) /
           sqrt( T*(double)(tpp->EpsW) );
//        pmp->PMc[jpb+1] = B;
    }
#else
    if( fabs(A) < 1e-9 )
        A = 1.82483e6 * sqrt( (double)(RoW_) ) /
           pow( T*(double)(EpsW_), 1.5 );
    if( fabs(B) < 1e-9 )
        B = 50.2916 * sqrt( (double)(RoW_) ) /
           sqrt( T*(double)(EpsW_) );
#endif
    ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "DebyeHueckel3Hel",
        "Error: A,B were not calculated - no values of RoW and EpsW !" );
    // Calculation of the EDH equation
    bgi = bg; // Common third parameter
    a0 = a0c; // Common ion-size parameter
    for( j=jb; j<je; j++ )
    {
//        bgi = bg; // Common third parameter
//        a0 = a0c; // Common ion-size parameter
        if( pmp->EZ[j] )
        { // Charged species
            Z2 = pmp->EZ[j]*pmp->EZ[j];
            lgGam = ( -A * sqI * Z2 ) / ( 1. + B * a0 * sqI ) + bgi * I ;
            pmp->lnGam[j] = lgGam * lg_to_ln;
            continue;
        }
        // Neutral species
        // <= -2 - gamma of both neutral aqueous species and water = 1.0
        // -1 - gamma of neutral species = 1.0, H2O activity coef. calculated
        // 0 - gamma of neutral species from salting out coefficient,
        //       H2O activity coefficient calculated
        //  >= 1 gamma of neutral species from salting out coefficient,
        //       H2O activity coefficient set to 1.0
        lgGam = 0.0;
        if( pmp->DCC[j] != DC_AQ_SOLVENT )
        {  // Calculation of gamma for neutral species except water
           if( nPolicy > -0.000001 )
  	       lgGam = bgi * I;
           pmp->lnGam[j] = lgGam * lg_to_ln;
           continue;
        }
        // Water-solvent
        lnGam = 0.0;  Lgam = 0.;
        if( nPolicy > -1.000001 && nPolicy < 0.000001 )
        {
       	   // Calculate activity coefficient of water solvent
	   // Lgam = -log10(1.+0.0180153*molT); // large gamma added (Helgeson 1981)
	   Lam = 1. + a0*B*sqI;
	   SigTerm = 3./(pow(a0,3.)*pow(B,3.)*pow(I,(3./2.)))*(Lam-1./Lam-2*log(Lam));
	   Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*I) - bgi*I/2.);
	   lnActWat = -Phi*molT/Nw;         
           lnGam = lnActWat - lnwxWat;      
//         lgGam = a0 * molt; // corrected: instead of I - tot.molality
        }
        pmp->lnGam[j] = lnGam;
    } // j
    return;
}

// Extended Debye-Hueckel eq. with Kielland ion-size and optional salting-out
// correction for neutral species
//
void
TMulti::DebyeHueckel2Kjel( int jb, int je, int jpb, int jdb, int k )
{
    int j;
    double T, A, B, a0, I, sqI, bg, bgi, Z2, lgGam, molt;
    float nPolicy;

    I= pmp->IC;
    if( I < TProfil::pm->pa.p.ICmin )
        return;
    T = pmp->Tc;
    A = (double)(pmp->PMc[jpb+0]);
    B = (double)(pmp->PMc[jpb+1]);
    bg = pmp->FitVar[0];   // Changed 07.06.05 for T,P-dep. b_gamma in DHH
//    bg = pmp->PMc[jpb+5];
//    a0c = pmp->PMc[jpb+6];
    nPolicy = (pmp->PMc[jpb+7]);

    molt = ( pmp->XF[0]-pmp->XFA[0] )*1000./18.01528/pmp->XFA[0]; // tot.molality
    sqI = sqrt( I );

#ifndef IPMGEMPLUGIN
    if( fabs(A) < 1e-9 )
    {
        A = 1.82483e6 * sqrt( (double)(tpp->RoW) ) /
           pow( T*(double)(tpp->EpsW), 1.5 );
//        pmp->PMc[jpb+0] = A;
    }
    if( fabs(B) < 1e-9 )
    {
        B = 50.2916 * sqrt( (double)(tpp->RoW) ) /
          sqrt( T*(double)(tpp->EpsW) );
//        pmp->PMc[jpb+1] = B;
    }
#else
    if( fabs(A) < 1e-9 )
        A = 1.82483e6 * sqrt( (double)(RoW_) ) /
          pow( T*(double)(EpsW_), 1.5 );
    if( fabs(B) < 1e-9 )
        B = 50.2916 * sqrt( (double)(RoW_) ) /
         sqrt( T*(double)(EpsW_) );
#endif
    ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "DebyeHueckel2Kjel",
        "Error: A,B were not calculated - no values of RoW and EpsW !" );
    // Calculation of EDH equation
    bgi = bg;
    for( j=jb; j<je; j++ )
    {
        a0 = (double)(pmp->DMc[jdb+j*pmp->LsMdc[k]]);
        if( pmp->EZ[j] && ( a0 < -1.000001 || a0 > -0.999999 ) ) // bugfix 26.07.07 DK
        { // Charged species
            Z2 = pmp->EZ[j]*pmp->EZ[j];
            lgGam = ( -A * sqI * Z2 ) / ( 1. + B * a0 * sqI ); // + bgi * I ;
        }
        else
        { // Neutral species
            if( nPolicy >= 0.0 )
            {
               if( a0 > 0.0 )
               {
                  if( pmp->DCC[j] != DC_AQ_SOLVENT ) // salting-out coefficient
                     lgGam = a0 * I;
                  else // water-solvent - a0 - rational osmotic coefficient
                     lgGam = a0 * molt; // corrected: instead of I - sum.molality
               }
               else {
                  if( a0 < -0.99 )
                      lgGam = 0.;
                  else if( fabs( a0 ) < 1e-9 )
                      lgGam = bgi * I;  // Average salting-out coeff.
                  else lgGam = a0 * I;  // Busenberg & Plummer
               }
            }
            else { // nPolicy < 0 - all gamma = 1 for neutral species
               lgGam = 0.;
            }
        }
        pmp->lnGam[j] = lgGam * lg_to_ln;
    } // j
}

// Debye-Hueckel limiting law
//
void
TMulti::DebyeHueckel1LL( int jb, int je, int )
{
    int j;
    double T, A, I, sqI, Z2, lgGam;
//    float nPolicy;

    I= pmp->IC;
    if( I < TProfil::pm->pa.p.ICmin )
        return;
    T = pmp->Tc;
//    A = pmp->PMc[jpb+0];
    sqI = sqrt( I );

#ifndef IPMGEMPLUGIN
//    if( fabs(A) < 1e-9 )
        A = 1.82483e6 * sqrt( (double)(tpp->RoW) ) /
          pow( T*(double)(tpp->EpsW), 1.5 );
#else
//    if( fabs(A) < 1e-9 )
        A = 1.82483e6 * sqrt( (double)(RoW_) ) /
          pow( T*(double)(EpsW_), 1.5 );
#endif
    ErrorIf( fabs(A) < 1e-9, "DebyeHueckel1LL",
        "Error: A was not calculated - no values of RoW and EpsW !" );
    // Calculation of DHLL equation
    for( j=jb; j<je; j++ )
    {
        if( pmp->EZ[j] )
        { // Charged species
            Z2 = pmp->EZ[j]*pmp->EZ[j];
            lgGam = ( -A * sqI * Z2 ); // / ( 1 + B * a0 * sqI ) + bgi * I ;
        }
        else  { // Neutral species
            lgGam = 0.;
        }
        pmp->lnGam[j] = lgGam * lg_to_ln;
    } // j
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Aqueous electrolyte  (Karpov's variant)
// Calculation of individual activity coefficients using extended
// Debye-Hueckel equation with common 3rd parameter (HKF81)
//  and individual (Kielland) ion size parameters
//
void TMulti::DebyeHueckel3Karp( int jb, int je, int jpb, int jdb, int k )
{
    int j;
    double T, A, B, a0, I, sqI, bg, bgi, Z2, lgGam, molt;
    float nPolicy;

    I= pmp->IC;
    if( I < TProfil::pm->pa.p.ICmin )
        return;
    T = pmp->Tc;
    A = (double)(pmp->PMc[jpb+0]);
    B = (double)(pmp->PMc[jpb+1]);
    bg = pmp->FitVar[0];   // Changed 07.06.05 for T,P-dep. b_gamma in DHH
//    bg = pmp->PMc[jpb+5];
//    a0c = pmp->PMc[jpb+6];
    nPolicy = (pmp->PMc[jpb+7]);

    molt = ( pmp->XF[0]-pmp->XFA[0] )*1000./18.01528/pmp->XFA[0]; /* tot.molality */
    sqI = sqrt( I );

#ifndef IPMGEMPLUGIN
    if( fabs(A) < 1e-9 )
    {
       A = 1.82483e6 * sqrt( (double)(tpp->RoW) ) /
         pow( T*(double)(tpp->EpsW), 1.5 );
//       pmp->PMc[jpb+0] = A;
    }
    if( fabs(B) < 1e-9 )
    {
       B = 50.2916 * sqrt( (double)(tpp->RoW) ) /
         sqrt( T*(double)(tpp->EpsW) );
//       pmp->PMc[jpb+1] = B;
    }
#else
    if( fabs(A) < 1e-9 )
        A = 1.82483e6 * sqrt( (double)(RoW_) ) /
          pow( T*(double)(EpsW_), 1.5 );
    if( fabs(B) < 1e-9 )
        B = 50.2916 * sqrt( (double)(RoW_) ) /
          sqrt( T*(double)(EpsW_) );
#endif
    ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "DebyeHueckel3Karp",
        "Error: A,B were not calculated - no values of RoW and EpsW !" );
    // Calculation of the EDH equation
//  bgi = bg;
    for( j=jb; j<je; j++ )
    {
        bgi = bg;
        a0 = (double)(pmp->DMc[jdb+j*pmp->LsMdc[k]]);
//        if( pmp->LsMdc[k] > 1 )
//        { // Individual bg coeff Truesdell-Jones (Parkhurst,1990)
//            bgi = pmp->DMc[jdb+j*pmp->LsMdc[k]+1];
//            if( !bgi )
//                bgi = bg;
//        }
        if( pmp->EZ[j] && ( a0 < -1.000001 || a0 > -0.999999 ) )
        { // Charged species
            Z2 = pmp->EZ[j]*pmp->EZ[j];
            lgGam = ( -A * sqI * Z2 ) / ( 1. + B * a0 * sqI ) + bgi * I ;
        }
        else
        { // Neutral species
            if( nPolicy >= 0.0 )
            {
               if( a0 > 0.0 )
               {
                  if( pmp->DCC[j] != DC_AQ_SOLVENT ) // Setchenow coefficient
                     lgGam = a0 * I;
                  else // water-solvent - a0 - rational osmotic coefficient
                     lgGam = a0 * molt; // corrected: instead of I - sum.molality
               }
               else {
                  if( a0 < -0.99 )
                      lgGam = 0.;
                  else if( fabs( a0 ) < 1e-9 )
                      lgGam = bgi * I;  // Average Setchenow coeff.
                  else lgGam = a0 * I;  // Busenberg & Plummer
               }
            }
            else { // nPolicy < 0 - all gamma = 1 for neutral species
               lgGam = 0.;
            }
        }
        pmp->lnGam[j] = lgGam * lg_to_ln;
    } // j
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Aqueous electrolyte
// Calculation of individual activity coefficients
// using the Davies equation with common 0.3 parameter and
// temperature-dependent A parameter
//
void TMulti::Davies03temp( int jb, int je, int )
{
    int j;
    double T, A, I, sqI, Z2, lgGam;

    I= pmp->IC;
    if( I < TProfil::pm->pa.p.ICmin )
        return;  // too low ionic strength

    T=pmp->Tc;
    sqI = sqrt( I );
//    if( fabs(A) < 1e-9 )
#ifndef IPMGEMPLUGIN
    A = 1.82483e6 * sqrt( (double)(tpp->RoW) ) /
      pow( T*(double)(tpp->EpsW), 1.5 );
#else
    A = 1.82483e6 * sqrt( (double)(RoW_) ) /
      pow( T*(double)(EpsW_), 1.5 );
#endif
//  at 25 C 1 bar: A = 0.5114
    ErrorIf( fabs(A) < 1e-9, "Davies03temp",
       "Error: A is not calculated - check values of RoW and EpsW !" );
    // Calculation of Davies equation: Langmuir 1997 p. 133
    for( j=jb; j<je; j++ )
    {
        if( pmp->EZ[j] )
        {   // Charged species
            Z2 = pmp->EZ[j]*pmp->EZ[j];
            lgGam = ( -A * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * I );
        }
        else
        { // Neutral species
          lgGam = 0;
        }
        pmp->lnGam[j] = lgGam * lg_to_ln;
    } // j
}

// non-ideal fluid mixtures
// ---------------------------------------------------------------------
// Churakov-Gottschalk (2004) calculation of pure gas/fluid component fugacity
// Added by D.Kulik on 15.02.2007
//
void
TMulti::CGofPureGases( int jb, int je, int, int jdb, int )
{
    double T, P, Fugacity = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
    float *Coeff, Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
                  Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;
    double X[1]={1.};
    int jdc, j, retCode = 0;

    TCGFcalc aCGF;
    P = pmp->Pc;
    T = pmp->Tc;

//    Coeff = pmp->DMc+jdb;
//    FugPure = (double*)malloc( NComp*sizeof(double) );
//    FugPure = pmp->Pparc+jb;

    for( jdc=0, j=jb; j<je; jdc++,j++)
    {
        Coeff = pmp->DMc+jdb+jdc*24;

        // Calling CG EoS for pure fugacity
        if( T >= 273.15 && T < 1e4 && P >= 1e-6 && P < 1e5 )
            retCode = aCGF.CGFugacityPT( Coeff, Eos4parPT, Fugacity, Volume, P, T );
        else {
            Fugacity = pmp->Pc;
            Volume = 8.31451*pmp->Tc/pmp->Pc;
            Coeff[15] = Coeff[0];
            if( Coeff[15] < 1. || Coeff[15] > 10. )
                Coeff[15] = 1.;                 // foolproof temporary
            Coeff[16] = Coeff[1];
            Coeff[17] = Coeff[2];
            Coeff[18] = Coeff[3];
            return;
        }

//        pmp->lnGmM[j] = log( Fugacity / pmp->Pc ); // Constant correction to G0 here!
pmp->GEX[j] = log( Fugacity / pmp->Pc );   // now here (since 26.02.2008)  DK
        pmp->Pparc[j] = Fugacity;          // Necessary only for performance
        pmp->Vol[j] = Volume * 10.;       // molar volume of pure fluid component, J/bar to cm3

        // passing corrected EoS coeffs to calculation of fluid mixtures
        Coeff[15] = Eos4parPT[0];
        if( Coeff[15] < 1. || Coeff[15] > 10. )
            Coeff[15] = 1.;                            // foolproof temporary
        Coeff[16] = Eos4parPT[1];
        Coeff[17] = Eos4parPT[2];
        Coeff[18] = Eos4parPT[3];
        
        aCGF.CGFugacityPT( Coeff, Eos4parPT1, Fugacity, Volume, P, T+T*aCGF.GetDELTA() );
         
        // passing corrected EoS coeffs for T+T*DELTA
        Coeff[19] = Eos4parPT1[0];
        if( Coeff[15] < 1. || Coeff[15] > 10. )
            Coeff[15] = 1.;                            // foolproof temporary
        Coeff[20] = Eos4parPT1[1];
        Coeff[21] = Eos4parPT1[2];
        Coeff[22] = Eos4parPT1[3];
        
        // Calculation of residual H and S
        aCGF.CGEntalpyRhoT( X, Eos4parPT, Eos4parPT1, 1, P, T, DeltaH, DeltaS);
        
    } // jdc, j

    if ( retCode )
    {
      char buf[150];
      sprintf(buf, "CG2004Fluid(): bad calculation of pure fugacities");
      Error( "E71IPM IPMgamma: ",  buf );
    }
}

// Churakov-Gottschalk (2004) multicomponent fluid mixing model
//
void
TMulti::ChurakovFluid( int jb, int je, int, int jdb, int k )
{
    double *FugCoefs;
    float *EoSparam, *EoSparam1, *Coeffs;
    int i, j, jj;
    double T, P, ro, DeltaH, DeltaS;
    TCGFcalc aCGF;
    P = pmp->Pc;
    T = pmp->Tc;

//    FugCoefs = (double*)malloc( pmp->L1[k]*sizeof(double) );
//    EoSparam = (float*)malloc( pmp->L1[k]*sizeof(double)*4 );
    FugCoefs =  new double[ pmp->L1[k] ];
    EoSparam =  new float[ pmp->L1[k]*4 ];
    EoSparam1 =  new float[ pmp->L1[k]*4 ];        
    Coeffs = pmp->DMc+jdb;

    // Copying T,P corrected coefficients
    for( j=0; j<pmp->L1[k]; j++)
    {   
    	for( i=0; i<4; i++)
          EoSparam[j*4+i] = Coeffs[j*24+i+15];
    	for( i=0; i<4; i++)
    	  EoSparam1[j*4+i] = Coeffs[j*24+i+19];
    }
//    EoSparam = pmp->DMc+jdb;
    if( T >= 273.15 && T < 1e4 && P >= 1e-6 && P < 1e5 )
    {
        ro = aCGF.CGActivCoefPT( pmp->X+jb, EoSparam, FugCoefs, pmp->L1[k],
             pmp->Pc, pmp->Tc );
        if (ro <= 0. )
        {
            delete[] FugCoefs;	
            delete[] EoSparam;
            delete[] EoSparam1;	
//           free( FugCoefs );
           char buf[150];
           sprintf(buf, "CGFluid(): bad calculation of density ro= %lg", ro);
           Error( "E71IPM IPMgamma: ",  buf );
        }
        // Phase volume of the fluid in cm3
        pmp->FVOL[k] = pmp->FWGT[k] / ro;
        // Get back residual H and S 
        aCGF.CGEntalpyRhoT( pmp->X+jb, EoSparam, EoSparam1, pmp->L1[k], pmp->Pc, pmp->Tc,
        		DeltaH, DeltaS );
        // Utilize residual enthalpy DeltaH, entropy DeltaS 
        // . . . . . . . . . . .
        // 
    }
    else // Setting Fugcoefs to 0 outside TP interval
      for( j=0; j<pmp->L1[k]; j++ )
        FugCoefs[ j ] = 0.0;

    for( jj=0, j=jb; j<je; j++, jj++ )
    {
        if( FugCoefs[jj] > 1e-23 /* && pmp->Pparc[j] > 1e-23 */ )
             pmp->lnGam[j] = log(FugCoefs[ jj ]/pmp->Pparc[j]);
        else
             pmp->lnGam[j] = 0;
    } // j
    delete[] EoSparam;
    delete[] EoSparam1;
    delete[] FugCoefs;	
//    free( FugCoefs );

}

#define MAXPRDCPAR 10
// ---------------------------------------------------------------------
// Entry to Peng-Robinson model for calculating pure gas fugacities
// Added by D.Kulik on 15.02.2007
//
void
TMulti::PRSVofPureGases( int jb, int je, int, int jdb, int, int )
{
    double /* *FugPure, */ Fugcoeff, Volume, DeltaH, DeltaS;
    float *Coeff; //  *BinPar;
    double Eos2parPT[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 } ;
    int j, jdc, NComp, retCode = 0;

    NComp = 1; // = pmp->L1[k];

    TPRSVcalc aPRSV( NComp, pmp->Pc, pmp->Tc );

//    Coeff = pmp->DMc+jdb;
//    FugPure = (double*)malloc( NComp*sizeof(double) );
//    FugPure = pmp->Pparc+jb;

    for( jdc=0, j=jb; j<je; jdc++,j++)
    {
         Coeff = pmp->DMc+jdb+jdc*12;	// increased from 10 to 12 (31.05.2008 TW)
         // Calling PRSV EoS for pure fugacity
         retCode = aPRSV.PRFugacityPT( pmp->Pc, pmp->Tc, Coeff,
                 Eos2parPT, Fugcoeff, Volume, DeltaH, DeltaS );

//    aW.twp->H +=  DeltaH;   // in J/mol - to be completed
//    aW.twp->S +=  DeltaS;   // to be completed
//         pmp->lnGmM[j] = log( Fugcoeff );    // Constant correction to G0 here!
pmp->GEX[j] = log( Fugcoeff );    // now here (since 26.02.2008) DK
         pmp->Pparc[j] = Fugcoeff * pmp->Pc; // Necessary only for performance
         pmp->Vol[j] = Volume * 10.;  // molar volume of pure fluid component, J/bar to cm3

//  passing corrected EoS coeffs to calculation of fluid mixtures (3 added, 31.05.2008 TW)
         Coeff[6] = Eos2parPT[0];	// a
         Coeff[7] = Eos2parPT[1];	// b
         Coeff[8] = Eos2parPT[2];	// sqrAL
         Coeff[9] = Eos2parPT[3];	// ac
         Coeff[10] = Eos2parPT[4];	// dALdT

    } // jdc, j

    if ( retCode )
    {
      char buf[150];
      sprintf(buf, "PRSVFluid(): bad calculation of pure fugacities");
      Error( "E71IPM IPMgamma: ",  buf );
    }

//    free( FugPure );
//    return retCode;
}

// ---------------------------------------------------------------------
// Entry to Peng-Robinson model for activity coefficients
// Added by Th.Wagner and D.Kulik on 19.07.2006, changed by DK on 15.02.2007
//
void
TMulti::PRSVFluid( int jb, int je, int jpb, int jdb, int k, int ipb )
{
    double *ActCoefs, PhVol, *FugPure;
    float *EoSparam, *BinPar;
    int j, jj, iRet, NComp, NPar, NPcoef, MaxOrd;
    short *aIPx;

    NComp = pmp->L1[k];

    TPRSVcalc aPRSV( NComp, pmp->Pc, pmp->Tc );

//    ActCoefs = (double*)malloc( NComp*sizeof(double) );
    ActCoefs =  new double[NComp];
    EoSparam = pmp->DMc+jdb;
//    FugPure = (double*)malloc( NComp*sizeof(double) );
    FugPure = pmp->Pparc + jb;
    BinPar = pmp->PMc+jpb;
    NPar = pmp->LsMod[k*3];      // Number of non-zero interaction parameters
    NPcoef = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    MaxOrd =  pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)
    aIPx = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters

    for( j=jb; j<je; j++)
       pmp->Wx[j] = pmp->X[j]/pmp->XF[k];

//    for( jj=0; jj<NComp; jj++)
//    {
//        FugPure[jj] = pmp->Pparc[jb+jj]; // left for debugging
//        FugPure[jj] = exp( pmp->lnGmM[jb+jj] ) * pmp->Pc ;
//    }

    iRet = aPRSV.PRActivCoefPT( NComp, pmp->Pc, pmp->Tc, pmp->Wx+jb, FugPure,
        BinPar, EoSparam, ActCoefs, PhVol, NPar, NPcoef, MaxOrd, aIPx );

    if ( iRet )
    {
      delete[] ActCoefs;	
//      free( ActCoefs );
//      free( FugPure );
      char buf[150];
      sprintf(buf, "PRSVFluid(): bad calculation");
      Error( "E71IPM IPMgamma: ",  buf );
    }
    // Phase volume of the fluid in cm3
    pmp->FVOL[k] = PhVol * 10.;

    for( jj=0, j=jb; j<je; j++, jj++ )
    {
        if( ActCoefs[jj] > 1e-23 )
             pmp->lnGam[j] = log( ActCoefs[ jj ]);
        else
             pmp->lnGam[j] = 0;
    } /* j */
   delete[] ActCoefs;	
//    free( ActCoefs );
//    free( FugPure );
}

// ------------------ condensed mixtures --------------------------
// Binary Redlich-Kister model - parameters (dimensionless)
// in ph_cf Phase opt.array 2x3, see also Phase module
// Implemented by KD on 31 July 2003
//
void
TMulti::RedlichKister( int jb, int, int jpb, int, int k )
{
  double a0, a1, a2, lnGam1, lnGam2, X1, X2;
  double gEX, vEX, hEX, sEX, cpEX, uEX;
  
  // load parameters
  a0 = (double)pmp->PMc[jpb+0];
  a1 = (double)pmp->PMc[jpb+1];  // in regular model should be 0
  a2 = (double)pmp->PMc[jpb+2];  // in regular model should be 0
  
  // load mole fractions
  X1 = pmp->X[jb] / pmp->XF[k];
  X2 = pmp->X[jb+1] / pmp->XF[k];
  
  // activity coeffs
  lnGam1 = X2*X2*( a0 + a1*(3.*X1-X2) + a2*(X1-X2)*(5.*X1-X2) );
  lnGam2 = X1*X1*( a0 - a1*(3.*X2-X1) + a2*(X2-X1)*(5.*X2-X1) );
  
  // bulk phase excess properties (added by TW, 29.05.2008)
  gEX = (X1*X2*( a0 + a1*(X1-X2) + a2*pow((X1-X2),2.) ))* pmp->RT;
  vEX = 0.0;
  uEX = (X1*X2*( a0 + a1*(X1-X2) + a2*pow((X1-X2),2.) ))* pmp->RT;
  sEX = 0.0;
  cpEX = 0.0;
  hEX = uEX;
  
  // assignments
  pmp->lnGam[jb] = lnGam1;
  pmp->lnGam[jb+1] = lnGam2;
}

// Binary Margules model - parameters (in J/mol) in ph_cf Phase opt.array 2x3
// See also Phase module
// Implemented by KD on 31 July 2003
//
void
TMulti::MargulesBinary( int jb, int, int jpb, int, int k )
{
  double T, P, WU1, WS1, WV1, WU2, WS2, WV2, WG1, WG2,
         a1, a2, lnGam1, lnGam2, X1, X2;
  double gEX, vEX, hEX, sEX, cpEX, uEX;

  // load parameters
  T = pmp->Tc;
  P = pmp->Pc;
  WU1 = (double)pmp->PMc[jpb+0];
  WS1 = (double)pmp->PMc[jpb+1];  // in J/K/mol, if unknown should be 0
  WV1 = (double)pmp->PMc[jpb+2];  // in J/bar if unknown should be 0
  WU2 = (double)pmp->PMc[jpb+3];
  WS2 = (double)pmp->PMc[jpb+4];  // if unknown should be 0
  WV2 = (double)pmp->PMc[jpb+5];  // if unknown should be 0
  
  // calculate parameters at (T,P)
  WG1 = WU1 - T*WS1 + P*WV1;
  WG2 = WU2 - T*WS2 + P*WV2;
  a1 = WG1 / pmp->RT;
  a2 = WG2 / pmp->RT;
  
  // load mole fractions
  X1 = pmp->X[jb] / pmp->XF[k];
  X2 = pmp->X[jb+1] / pmp->XF[k];
  
  // activity coefficients
  lnGam1 = (2.*a2-a1)*X2*X2 + 2.*(a1-a2)*X2*X2*X2;
  lnGam2 = (2.*a1-a2)*X1*X1 + 2.*(a2-a1)*X1*X1*X1;

  // bulk phase excess properties (extended by TW, 29.05.2008)
  gEX = X1*X2*( X2*WG1 + X1*WG2 );
  vEX = X1*X2*( X2*WV1 + X1*WV2 );
  uEX = X1*X2*( X2*WU1 + X1*WU2 );
  sEX = X1*X2*( X2*WS1 + X1*WS2 );
  cpEX = 0.0;
  hEX = uEX+vEX*P;
  
  // assignments
  pmp->lnGam[jb] = lnGam1;
  pmp->lnGam[jb+1] = lnGam2;
  pmp->FVOL[k] += vEX*10.;  // make consistent with TSolMod
  
 }

// Ternary regular Margules model - parameters (in J/mol)
// in ph_cf Phase opt.array 4x3; see also Phase module
// Implemented by KD on 31 July 2003
//
void
TMulti::MargulesTernary( int jb, int, int jpb, int, int k )
{
  double T, P, WU12, WS12, WV12, WU13, WS13, WV13, WU23, WS23, WV23,
         WU123, WS123, WV123, WG12, WG13, WG23, WG123,
         a12, a13, a23, a123, lnGam1, lnGam2, lnGam3, X1, X2, X3;
  double gEX, vEX, hEX, sEX, cpEX, uEX;

  // load parameters
  T = pmp->Tc;
  P = pmp->Pc;
  WU12 = (double)pmp->PMc[jpb+0];
  WS12 = (double)pmp->PMc[jpb+1];  // if unknown should be 0
  WV12 = (double)pmp->PMc[jpb+2];  // if unknown should be 0
  WU13 = (double)pmp->PMc[jpb+3];
  WS13 = (double)pmp->PMc[jpb+4];  // if unknown should be 0
  WV13 = (double)pmp->PMc[jpb+5];  // if unknown should be 0
  WU23 = (double)pmp->PMc[jpb+6];
  WS23 = (double)pmp->PMc[jpb+7];  // if unknown should be 0
  WV23 = (double)pmp->PMc[jpb+8];  // if unknown should be 0
  WU123 = (double)pmp->PMc[jpb+9];
  WS123 = (double)pmp->PMc[jpb+10];  // if unknown should be 0
  WV123 = (double)pmp->PMc[jpb+11];  // if unknown should be 0
  
  // calculate parameters at (T,P)
  WG12 = WU12 - T*WS12 + P*WV12;
  WG13 = WU13 - T*WS13 + P*WV13;
  WG23 = WU23 - T*WS23 + P*WV23;
  WG123 = WU123 - T*WS123 + P*WV123;
  a12 = WG12 / pmp->RT;
  a13 = WG13 / pmp->RT;
  a23 = WG23 / pmp->RT;
  a123 = WG123 / pmp->RT;
  
  // load mole fractions
  X1 = pmp->X[jb] / pmp->XF[k];
  X2 = pmp->X[jb+1] / pmp->XF[k];
  X3 = pmp->X[jb+2] / pmp->XF[k];

  // activity coefficients
  lnGam1 = a12*X2*(1.-X1) + a13*X3*(1.-X1) - a23*X2*X3
           + a123*X2*X3*(1.-2.*X1);
  lnGam2 = a23*X3*(1.-X2) + a12*X1*(1.-X2) - a13*X1*X3
           + a123*X1*X3*(1.-2.*X2);
  lnGam3 = a13*X1*(1.-X3) + a23*X2*(1.-X3) - a12*X1*X2
           + a123*X1*X2*(1.-2.*X3);
  
  // bulk phase excess properties (added by TW, 29.05.2008)
  gEX = X1*X2*WG12 + X1*X3*WG13 + X2*X3*WG23 + X1*X2*X3*WG123;
  vEX = X1*X2*WV12 + X1*X3*WV13 + X2*X3*WV23 + X1*X2*X3*WV123;
  uEX = X1*X2*WU12 + X1*X3*WU13 + X2*X3*WU23 + X1*X2*X3*WU123;
  sEX = X1*X2*WS12 + X1*X3*WS13 + X2*X3*WS23 + X1*X2*X3*WS123;
  cpEX = 0.0;
  hEX = uEX+vEX*P;
  
  // assignments
  pmp->lnGam[jb] = lnGam1;
  pmp->lnGam[jb+1] = lnGam2;
  pmp->lnGam[jb+2] = lnGam3;
  pmp->FVOL[k] += vEX*10.;  // make consistent with TSolMod

}


// ------------------------------------------------------------------------
// Wrapper calls for generic multi-component mixing models (see s_fgl.h and s_fgl2.cpp)
// Uses the TSolMod class by Th.Wagner and D.Kulik

void
TMulti::SolModParPT( int, int, int jpb, int jdb, int k, int ipb, char ModCode )
{
    int NComp, NPar, NPcoef, MaxOrd, NP_DC;
    float *aIPc, *aDCc;
    short *aIPx;
    double RhoW, EpsW;

    NComp = pmp->L1[k];          // Number of components in the phase
    NPar = pmp->LsMod[k*3];      // Number of interaction parameters
    NPcoef = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    MaxOrd =  pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)
    NP_DC = pmp->LsMdc[k]; // Number of non-ideality coeffs per one DC in multicomponent phase
    aIPx = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
                              // -> NPar x MaxOrd   added 07.12.2006   KD
    aIPc = pmp->PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    aDCc = pmp->DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    RhoW = pmp->denW;		// added 04.06.2008 (TW)
    EpsW = pmp->epsW;

    TSolMod aSM( NComp, NPar, NPcoef, MaxOrd, NP_DC, pmp->Tc, pmp->Pc, ModCode,
       aIPx, aIPc, aDCc, NULL, NULL, RhoW, EpsW, pmp->IC );

   // calculate P-T dependence of interaction parameters
    switch( ModCode )
    {
        case SM_VANLAAR:
             aSM.VanLaarPT();
             break;
        case SM_REGULAR:
             aSM.RegularPT();
             break;
        case SM_GUGGENM:
        	 aSM.RedlichKisterPT();
        	 break;
        case SM_NRTLLIQ:
        	 aSM.NRTL_PT();
        	 break;
        case SM_WILSLIQ:
        	 aSM.Wilson_PT();
        	 break;
        default:
             break;
    }
}

void
TMulti::SolModActCoeff( int jb, int, int jpb, int jdb, int k, int ipb,
            char ModCode )
{
    int NComp, NPar, NPcoef, MaxOrd, NP_DC;
    float *aIPc, *aDCc;
    double *aWx, *alnGam;
    short *aIPx;
    double Gex=0.0, Vex=0.0, Hex=0.0, Sex=0.0, CPex=0.0;
    double RhoW, EpsW, IonStr;

    NComp = pmp->L1[k];          // Number of components in the phase
    NPar = pmp->LsMod[k*3];      // Number of interaction parameters
    NPcoef = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    MaxOrd =  pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)
    NP_DC = pmp->LsMdc[k]; // Number of non-ideality coeffs per one DC in multicomponent phase[FIs]

    // These pointers provide direct access to parts of MULTI arrays related to this phase!
    aIPx = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
                              // -> NPar x MaxOrd   added 07.12.2006   KD
    aIPc = pmp->PMc+jpb;    // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    aDCc = pmp->DMc+jdb;    // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    aWx = pmp->Wx+jb;       // End member mole fractions
    alnGam = pmp->lnGam+jb; // End member ln activity coeffs
    RhoW = pmp->denW;		// added 04.06.2008 (TW)
    EpsW = pmp->epsW;
    IonStr = pmp->IC;

    TSolMod aSM( NComp, NPar, NPcoef, MaxOrd, NP_DC, pmp->Tc, pmp->Pc, ModCode,
       aIPx, aIPc, aDCc, aWx, alnGam, RhoW, EpsW, IonStr );
    // Extended constructor to connect to params, coeffs, and mole fractions

    switch( ModCode )
    {
        case SM_VANLAAR:
             aSM.VanLaarMixMod( Gex, Vex, Hex, Sex, CPex );
             break;
        case SM_REGULAR:
             aSM.RegularMixMod( Gex, Vex, Hex, Sex, CPex );
             break;
        case SM_GUGGENM:
        	 aSM.RedlichKisterMixMod( Gex, Vex, Hex, Sex, CPex );
        	 break;
        case SM_NRTLLIQ:
        	 aSM.NRTL_MixMod( Gex, Vex, Hex, Sex, CPex );
        	 break;
        case SM_WILSLIQ:
        	 aSM.Wilson_MixMod( Gex, Vex, Hex, Sex, CPex );
        	 break;
        default: // catch error here
              break;
    }
    // To add handling of excess properties for the phase
    // Gex, Vex, Hex, Sex, CPex

}



//--------------------- End of ipm_chemical3.cpp ---------------------------
