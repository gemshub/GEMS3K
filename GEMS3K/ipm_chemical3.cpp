//-------------------------------------------------------------------
// $Id: ipm_chemical3.cpp 690 2006-03-29 07:10:23Z gems $
//
// Copyright (C) 1992-2011  D.Kulik, T.Wagner, S.Dmitrieva, K.Chudnenko
//
// Implementation of chemistry-specific functions (concentrations,
// activity coefficients, adsorption models etc.)
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806); Kulik (2000), Geoch.Cosmoch.Acta 64,3161-3179
//
// This file is part of a GEM-Selektor (GEMS) v.3.x program
// environment for thermodynamic modeling in geochemistry and of the
// standalone GEMS3K code (define IPMGEMPLUGIN).
//
// This file may be distributed under the terms of the GEM-Selektor
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include <math.h>
#include "m_param.h"
#include "s_fgl.h"


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Linking DOD for executing Phase mixing model scripts
#ifdef IPMGEMPLUGIN
   void TMulti::pm_GC_ods_link( long int /*k*/, long int /*jb*/, long int /*jpb*/, long int /*jdb*/, long int /*ipb*/ )
   {
#else
	void TMulti::pm_GC_ods_link( long int k, long int jb, long int jpb, long int jdb, long int ipb )
	{

        ErrorIf( k < 0 || k >= pmp->FIs , "CalculateActivityCoefficients():", "Invalid link: k=0||>FIs" );
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
    aObj[ o_nppar].SetPtr(pmp->G0+jb );  // changed 10.12.2008 by DK
    aObj[o_nppar].SetN(  pmp->L1[k]);
//    aObj[ o_ngtn].SetPtr( pmp->G0+jb );
    aObj[ o_ngtn].SetPtr( pmp->fDQF+jb );     // changed 05.12.2006 by DK
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
    //
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


// Returns current value of smoothing factor for chemical potentials of highly non-ideal DCs
// added 18.06.2008 DK
double TMulti::SmoothingFactor( )
{
   if( pmp->FitVar[4] < 0 )
   {  // To start SIA mode (smart initial approximation)
      return 1.0;
   }
   if( pmp->FitVar[3] > 0 )
	   return pmp->FitVar[3];
   else
	   return pmp->FitVar[4];
}


// New correction of smoothing factor for highly non-ideal systems
// re-written 18.04.2009 DK+TW
// Smoothing function choice: AG >= 0.0001 and DGC > -0.0001: old f(IT)
//                            AG >= 0.0001 and DGC <= -0.0001: new f(1/IT)
//                            AG <= -0.0001 and DGC <= -0.0001: new f(1/CD)
// mode: 0 - taking single log(CD) value for calculation of smoothing factor SF;
//       1, 2, ...  taking log(CD) average from the moving window of length mode
//       (up to 5 consecutive values)
//
void TMulti::SetSmoothingFactor( long int mode )
{
    double TF=1., al, ag, dg, iim, irf; // rg=0.0;
    long int ir; //, Level, itqF, itq;

    ir = pmp->IT;
    irf = (double)ir;
    ag = prof->pa.p.AG; // pmp->FitVar[4];
    dg = prof->pa.p.DGC;
    iim = (double)prof->pa.p.IIM;

    if( dg > -0.0001 && ag >= 0.0001 ) // Smoothing used in the IPM-2 algorithm
    {					// with some improvements
        if(ag>1) ag=1;
        if(ag<0.1) ag=0.1;
        if(dg>0.15) dg=0.15;
        // if( irf > 1000. )
        //	irf = 1000;
        if( dg <= 0.0 )
          TF = ag;
        else
          TF = ag * ( 1 - pow(1-exp(-dg*irf),60.));
        if(TF < 1e-6 )
          TF = 1e-6;
    }
    else if( dg <= -0.0001 && ag >= 0.0001 )
    {
       // New sigmoid smoothing function of 1/IT
    	double logr, inv_r = 1., logr_m;
    	dg = fabs( dg );
    	if( pmp->IT )
    	  inv_r = 1./(double)pmp->IT;
        logr = log( inv_r );
        logr_m = log( 1./iim );
        al = dg + ( ag - dg ) / ( 1. + exp( logr_m - logr ) / dg );
        al += exp( log( 1. - ag ) + logr );
        if( al > 1. )
      	    al = 1.;
        TF = al;
    }
    else if( dg <= -0.0001 && ag <= -0.0001 )
    {
    	double dk, cd;   long int i;
    	dg = fabs( dg );
    	ag = fabs( ag );
    	dk = log( pmp->DXM );
    	// Checking the mode where it is called
    	switch( mode )
    	{
    	  default:
          case 0: // MassBalanceRefinement() after SolveSimplex()
    	  	     cd = log( pmp->PCI );
    	   	     break;
    	  case 1:
    	  case 2:
    	  case 3:
    	  case 4:
    	  case 5: // Getting average (log geometric mean) from sampled CD values
    	  	     cd = 0.0;
    	   	     for(i=0; i < mode; i++ )
    	   	    	 cd += pmp->logCDvalues[i];
                     cd /= (double)mode; // 5. - bugfix
    	   	     break;
    	}
        al = dg + ( ag - dg ) / ( 1. + exp( dk - cd ) / dg );
        al += exp( log( 1. - ag ) + cd );
        if( al > 1. )
      	    al = 1.;
        TF = al;
    }
//    if( pmp->IT )
      pmp->FitVar[3] = TF;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//   Function for converting internal lnGam[j] value into an external (phase-scale-specific)
//      Gamma[j] if DirFlag = 0 or external into internal value if DirFlag = 1
//  Returns the respectively corrected external gamma activity coefficient or internal lnGam
//  Returns trivial values (lnGam = 0 or Gamma = 1) when the respective component
//    amount is zero (X[j] == 0) (is this a correct policy for zeroed-off components?)
//
double
TMulti::PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag )
{
    double NonLogTerm = 0., NonLogTermW = 0., NonLogTermS = 0., MMC = 0.;
//    SPP_SETTING *pa = &prof->pa;

    switch( pmp->PHC[k] )
    {
      case PH_AQUEL:
           if( pmp->XF[k] && pmp->XFA[k] )
           {
                NonLogTerm = 1. - pmp->XFA[k]/pmp->XF[k];
                NonLogTermW = 2. - pmp->XFA[k]/pmp->XF[k] - pmp->XF[k]/pmp->XFA[k];
           }
           break;
      case PH_GASMIX:  case PH_FLUID:   case PH_PLASMA:   case PH_SIMELT:
      case PH_HCARBL:  case PH_SINCOND:  case PH_SINDIS:  case PH_LIQUID:
           break;
      case PH_POLYEL:
      case PH_SORPTION: // only sorbent end-members!
           if( pmp->XF[k] && pmp->XFA[k] )
           {
              for( long int jj=jb; jj<je; jj++ )
              {
                if( pmp->DCC[jj] == DC_SUR_CARRIER ||
                    pmp->DCC[jj] == DC_SUR_MINAL || pmp->DCC[jj] == DC_PEL_CARRIER )
                    MMC += pmp->MM[jj]*pmp->X[jj]/pmp->XFA[k];
                    // Weighted-average sorbent mole mass
              }
              NonLogTerm = 1. - pmp->XFA[k]/pmp->XF[k];  // Also for sorption phases
              NonLogTermS = 2. - pmp->XFA[k]/pmp->XF[k] - pmp->XF[k]/pmp->XFA[k];
           }
           break;
       default:
          break; // Phase class code error should be generated here!
    }

        if( DirFlag == 0 )
        {	 // Converting lnGam[j] into Gamma[j]
                if( !pmp->X[j] )
                        return 1.;
                double Gamma = 1.;
                double lnGamS = pmp->lnGam[j];

            switch( pmp->DCC[j] )
            { // Aqueous electrolyte
              case DC_AQ_PROTON: case DC_AQ_ELECTRON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                lnGamS += NonLogTerm;    // Correction by asymmetry term
                break;
                // calculate molar mass of solvent
            case DC_AQ_SOLVCOM:	    case DC_AQ_SOLVENT:
                lnGamS += NonLogTermW;
                break;
            case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:
            case DC_GAS_H2: case DC_GAS_N2:
                break;
            case DC_SOL_IDEAL:  case DC_SOL_MINOR:  case DC_SOL_MAJOR:
                break;
                // non-electrolyte condensed mixtures
            case DC_SCP_CONDEN: case DC_SUR_MINAL:
                break;
            case DC_SUR_CARRIER: case DC_PEL_CARRIER:
                lnGamS += NonLogTermS;
                break;
                // Sorption phases
            case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
            case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
            case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:  case DC_IESC_A:
            case DC_IEWC_B:
                lnGamS += NonLogTerm;
                break;
            default:
                break;
            }
        Gamma = exp( lnGamS );
            return Gamma;
        }
        else { // Converting Gamma[j] into lnGam[j]
                if( !pmp->X[j] )
                        return 0.;
                double Gamma = pmp->Gamma[j];
                double lnGam = 0.0;  // Cleanup by DK 5.12.2009
                if( Gamma != 1.0 && Gamma > pmp->lowPosNum )
                    lnGam = log( Gamma );
                switch( pmp->DCC[j] )
        { // Aqueous electrolyte
                   case DC_AQ_PROTON: case DC_AQ_ELECTRON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                        lnGam -= NonLogTerm;  // Correction by asymmetry term
                        break;
                   case DC_AQ_SOLVCOM:	    case DC_AQ_SOLVENT:
                        lnGam -= NonLogTermW;
                        break;
               case DC_GAS_COMP: case DC_GAS_H2O: case DC_GAS_CO2: case DC_GAS_H2: case DC_GAS_N2:
                                break;
                   case DC_SOL_IDEAL:  case DC_SOL_MINOR:  case DC_SOL_MAJOR:
                        break;
               case DC_SCP_CONDEN: case DC_SUR_MINAL:
                            break;
               case DC_SUR_CARRIER: case DC_PEL_CARRIER:
                        lnGam -= NonLogTermS;
                            break;
                                // Sorption phases
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
                   case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
                   case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:  case DC_IESC_A:
                   case DC_IEWC_B:
                        lnGam -= NonLogTerm;
                        break;
                    default:
                        break;
                }
            return lnGam;
        }
}

//--------------------------------------------------------------------------------
static double ICold=0.;
// Main call point for calculation of DC activity coefficients (lnGam vector)
//    formerly GammaCalc()
// Controls various built-in models, as well as generic Phase script calculation
// LinkMode is a parameter indicating the status of Gamma calculations:
// LINK_TP_MODE - calculation of equations depending on TP only;
// LINK_UX_MODE - calculation of equations depending on current
//      IPM approximation of the equilibrium state;
// LINK_PP_MODE - calculation of integral phase properties after GEMIPM has converged
//		needs to be implemented
// Returns: status code (0 if o.k., non-zero values if there were problems
//     with surface complexation models)
long int
TMulti::CalculateActivityCoefficients( long int LinkMode  )
{
    long int k, j, jb, je=0, jpb, jpe=0, jdb, jde=0, ipb, ipe=0, ja;
    char *sMod;
    long int statusGam=0, statusGC=0, statusSACT=0, SmMode = 0;
    double LnGam, pmpXFk;
    SPP_SETTING *pa = &prof->pa;

    // calculating concentrations of species in multi-component phases
    switch( LinkMode )
    {
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
                if(pmp->XF[k] < min( pmp->DSM, pmp->PhMinM ) ) // pmp->lowPosNum )   // workaround 10.03.2008 DK
            		pmp->Wx[j] = nxk;  // need this eventually to avoid problems with zero mole fractions
                pmp->fDQF[j] =0.0;  // cleaning fDQF in TP mode!
            	pmp->lnGmo[j] = pmp->lnGam[j]; // saving activity coefficients in TP mode
       	    }
            if( sMod[SGM_MODE] != SM_STNGAM )
                continue;  // The switch below is for built-in functions only!

            // the following section probably needs to be re-written to allow more flexible
            // combinations of fluid models for pure gases with gE mixing models,
            // scheme should probably be the same as in LINK_UX_MODE, 03.06.2008 (TW)
            switch( pmp->PHC[k] )
            {
               case PH_AQUEL:
               case PH_LIQUID:
               case PH_SINCOND:
               case PH_SINDIS:
               case PH_HCARBL:
               case PH_SIMELT:
               case PH_GASMIX:
               case PH_PLASMA:
               case PH_FLUID:
            	    SolModCreate( jb, je, jpb, jdb, k, ipb, sMod[SPHAS_TYP], sMod[MIX_TYP] ); // new solution models (TW, DK 2007)
            	    SolModParPT( k, sMod[SPHAS_TYP] );
            	    break;
              default:
                    break;
            }
        } // k
        break;

      case LINK_PP_MODE: // Mode of calculation of integral solution phase properties
    	for( k=0; k<pmp->FIs; k++ )
        { // loop on solution phases
            jb = je;
            je += pmp->L1[k];
            sMod = pmp->sMod[k];
    		switch( pmp->PHC[k] )
            {
              case PH_AQUEL:
              case PH_LIQUID:
              case PH_SINCOND:
              case PH_SINDIS:
              case PH_HCARBL:
              case PH_SIMELT:
              case PH_GASMIX:
              case PH_PLASMA:
              case PH_FLUID:
           	       SolModExcessProp( k, sMod[SPHAS_TYP] ); // extracting integral phase properties
           	       SolModIdealProp( jb, k, sMod[SPHAS_TYP] );
           	       SolModStandProp( jb, k, sMod[SPHAS_TYP] );
           	       SolModDarkenProp( jb, k, sMod[SPHAS_TYP] );
           	       break;
              default:
                       break;
            }
       } // k
       break;

    case LINK_UX_MODE:
    	// Getting actual smoothing parameter
    	SetSmoothingFactor( SmMode );
    	// calculating DC concentrations after this IPM iteration
        CalculateConcentrations( pmp->X, pmp->XF, pmp->XFA );
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
               {  long int ist;
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
        Error("CalculateActivityCoefficients()","Invalid LinkMode for a built-in solution model");
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
        ipb = ipe;
        ipe += pmp->LsMod[k*3]*pmp->LsMod[k*3+1];
        jpb = jpe;
        jpe += pmp->LsMod[k*3]*pmp->LsMod[k*3+2];
        jdb = jde;
        jde += pmp->LsMdc[k]*pmp->L1[k];

   if( LinkMode == LINK_UX_MODE && sMod[SGM_MODE] == SM_STNGAM )
   {    // check that SGM_MODE for adsorption or multi-site ideal SS is not SM_IDEAL in Phase records!
        switch( pmp->PHC[k] )
        {  // calculating activity coefficients using built-in functions
          case PH_AQUEL:   // DH III variant consistent with HKF
             if( pmpXFk > pmp->DSM && pmp->X[pmp->LO] > pmp->XwMinM && pmp->IC > pa->p.ICmin )
             {
                switch( sMod[SPHAS_TYP] )
                {
					case SM_AQDH3:
					case SM_AQDH2:
					case SM_AQDH1:
					case SM_AQDHS:
					case SM_AQDHH:
					case SM_AQDAV:
					case SM_AQSIT:
					case SM_AQPITZ:
					case SM_AQEXUQ:
						SolModActCoeff( k, sMod[SPHAS_TYP] );
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
                 if( sMod[SPHAS_TYP] == SM_CGFLUID )  // CG EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_PRFLUID )  // PRSV EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_SRFLUID )  // SRK EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_PR78FL )  // PR78 EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_CORKFL )  // CORK EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
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
					case SM_REDKIS:  // Redlich-Kister model (binary)
					case SM_MARGB:   // Subregular Margules model (binary)
					case SM_MARGT:   // Regular Margules model (ternary)
					case SM_GUGGENM: // Redlich-Kister model (multicomponent), 2007 (TW)
					case SM_VANLAAR: // VanLaar model (multicomponent), 2007 (TW)
					case SM_REGULAR: // Regular model (multicomponent), 2007 (TW)
					case SM_NRTLLIQ: // NRTL model (multicomponent), 03.06.2007 (TW)
					case SM_WILSLIQ: // Wilson model (multicomponent), 09.06.2007 (TW)
                       SolModActCoeff( k, sMod[SPHAS_TYP] );
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
                && (pmp->XFA[0] > pmp->XwMinM && pmp->XF[0] > pmp->DSM ))
                {
                    if( pmp->E )
                    {
                       statusGC = GouyChapman( jb, je, k );
                    // PoissonBoltzmann( q, jb, je, k )
                    }
                    // Calculating surface activity coefficient terms
                    statusSACT = SurfaceActivityCoeff(  jb, je, jpb, jdb, k );
                }
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

        case LINK_PP_MODE: // Mode of calculation of integral solution phase properties
        		switch( pmp->PHC[k] )
                {
                  case PH_AQUEL:
                  case PH_LIQUID:
                  case PH_SINCOND:
                  case PH_SINDIS:
                  case PH_HCARBL:
                  case PH_SIMELT:
                  case PH_GASMIX:
                  case PH_PLASMA:
                  case PH_FLUID:  // How to pull this stuff out of the script?
                	  // SolModExcessProp( k, sMod[SPHAS_TYP] ); // extracting integral phase properties
                	  // SolModIdealProp( jb, k, sMod[SPHAS_TYP] );
                	  // SolModStandProp( jb, k, sMod[SPHAS_TYP] );
                	  // SolModDarkenProp( jb, k, sMod[SPHAS_TYP] );
               	       break;
                  default:
                       break;
                }
                break;

        case LINK_UX_MODE:  // the model is dependent on current concentrations on IPM iteration
            switch( pmp->PHC[k] )
            {  //
              case PH_AQUEL:
            	  if(!(pmpXFk > pmp->DSM && pmp->X[pmp->LO] > pmp->XwMinM && pmp->IC > pa->p.ICmin ))
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
                      && (pmp->XFA[0] > pmp->XwMinM && pmp->XF[0] > pmp->DSM )))
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
            Error("CalculateActivityCoefficients()","Invalid LinkMode for a scripted solution model");
        } // end switch
#endif

END_LOOP:
        if( LinkMode == LINK_TP_MODE )  // TP mode - added 04.03.2008 by DK
        {
        	for( j=jb; j<je; j++ )
        	{
        	   if( pmp->XF[k] < pmp->DSM )   // workaround 10.03.2008 DK
        			pmp->Wx[j] = 0.0;               //
        	   LnGam = pmp->lnGmo[j];
                   pmp->lnGam[j] = LnGam;
                   if(  fabs( LnGam ) < 84. )
                       pmp->Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
                   else pmp->Gamma[j] = 1.0;
        	}
        }
        else if(LinkMode == LINK_UX_MODE )  // Bugfix! DK 06.04.11
        { // Real mode for activity coefficients
           double lnGamG;
           for( j=jb; j<je; j++ )
           {
             if( pmp->DCC[j] == DC_AQ_SURCOMP )  // Workaround for aqueous surface complexes DK 22.07.09
                pmp->lnGam[j] = 0.0;
             lnGamG = PhaseSpecificGamma( j, jb, je, k, 1 );
             LnGam = pmp->lnGam[j];
             if( fabs( lnGamG ) > 1e-9 )
            	LnGam += lnGamG;
             pmp->lnGmo[j] = LnGam;
             if( fabs( LnGam ) < 84. )   // before 26.02.08: < 42.
          	    pmp->Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
             else pmp->Gamma[j] = 1.0;

             pmp->F0[j] = DC_PrimalChemicalPotentialUpdate( j, k );
             pmp->G[j] = pmp->G0[j] + pmp->fDQF[j] + pmp->F0[j];
           }
        }
    }  // k - end loop over phases

    if( statusGC )
        return statusGC;
    if( statusSACT )
    	return statusSACT;
    return statusGam;
}


//--------------------------------------------------------------------------------
// Wrapper calls for multi-component mixing models
// using  TSolMod class
//
void TMulti::SolModCreate( long int jb, long int, long int jpb, long int jdb, long int k,
		long int ipb, char ModCode, char MixCode )
{
    long int NComp, NPar, NPcoef, MaxOrd, NP_DC;
    double *aIPc, *aDCc, *aWx, *alnGam, *aphVOL, *aZ, *aM;
    long int *aIPx;
    char *DCCp;

    NComp = pmp->L1[k];          // Number of components in the phase
    NPar = pmp->LsMod[k*3];      // Number of interaction parameters
    NPcoef = pmp->LsMod[k*3+2];  // and number of coefs per parameter in PMc table
    MaxOrd =  pmp->LsMod[k*3+1];  // max. parameter order (cols in IPx)
    NP_DC = pmp->LsMdc[k];       // Number of non-ideality coeffs per one DC in multicomponent phase

    if( phSolMod[k])
    	if(	phSolMod[k]->testSizes( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode ) )
    	{
    		phSolMod[k]->UpdatePT( pmp->Tc, pmp->Pc );
    		return; // using old allocation
    	}

    aIPx = pmp->IPx+ipb;   // Pointer to list of indexes of non-zero interaction parameters for non-ideal solutions
                              // -> NPar x MaxOrd   added 07.12.2006   KD
    aIPc = pmp->PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    aDCc = pmp->DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    aWx = pmp->Wx+jb;       // End member mole fractions
    alnGam = pmp->lnGam+jb; // End member ln activity coeffs
    aM = pmp->Y_m+jb;
    aZ = pmp->EZ+jb;
    aphVOL = pmp->FVOL+k;
    DCCp = pmp->DCC+jb;  // pointer to Dcomp class codes (added 02.05.2010 TW)

    TSolMod* mySM = 0;

   // creating instances of subclasses of TSolMod base class
    switch( ModCode )
    {

        case SM_OTHER:  // Hard-coded solid solution models (selected by phase name)
        {
                TModOther* myPT = new TModOther( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc, pmp->denW, pmp->epsW );
                myPT->GetPhaseName( pmp->SF[k] );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_VANLAAR:  // Van Laar solid solution model (multicomponent)
        {
                TVanLaar* myPT = new TVanLaar( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }
             // break;

        case SM_REGULAR:  // Regular solid solution model (multicomponent)
        {
                TRegular* myPT = new TRegular( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_GUGGENM:  // Redlich-Kister solid solution model (multicomponent)
        {
                TRedlichKister* myPT = new TRedlichKister( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_NRTLLIQ:  // NRTL liquid solution model (multicomponent)
        {
                TNRTL* myPT = new TNRTL( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_WILSLIQ:  // Wilson liquid solution model (multicomponent)
        {
                TWilson* myPT = new TWilson( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_MARGT:  // Margules ternary (regular) solid solution model
        {
                TMargules* myPT = new TMargules( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_MARGB:  // Margules binary (subregular) solid solution model
        {
                TSubregular* myPT = new TSubregular( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_REDKIS:  // Gugenheim binary (REdlich-Kister) solid solution
        {
                TGuggenheim* myPT = new TGuggenheim( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQPITZ:  // Pitzer aqueous electrolyte model (multicomponent)
        {
                TPitzer* myPT = new TPitzer( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc, pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQSIT:  // SIT aqueous electrolyte model (multicomponent)
        {
                TSIT* myPT = new TSIT( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc, pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQEXUQ:  // EUNIQUAC aqueous electrolyte model (multicomponent)
        {
                TEUNIQUAC* myPT = new TEUNIQUAC( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc, pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQDH3:  // extended Debye-Hueckel aqueous electrolyte model (Karpov version)
        {
                TKarpov* myPT = new TKarpov( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDH2:   // Debye-Hueckel aqueous electrolyte model
        {
                TDebyeHueckel* myPT = new TDebyeHueckel( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDH1:   // Debye-Hueckel limiting law aqueous electrolyte model
        {
                TLimitingLaw* myPT = new TLimitingLaw( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDHS:  // extended Debye-Hueckel aqueous electrolyte model (Shvarov version)
        {
                TShvarov* myPT = new TShvarov( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDHH:  // extended Debye-Hueckel aqueous electrolyte model (Helgeson version)
        {
                THelgeson* myPT = new THelgeson( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDAV:  // Davies aqueous electrolyte model (in NEA TDB version)
        {
                TDavies* myPT = new TDavies( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, aM, aZ, pmp->Tc, pmp->Pc,  pmp->denW, pmp->epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_PRFLUID:  // PRSV fluid mixture (multicomponent)
        {
                TPRSVcalc* myPT = new TPRSVcalc( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Pparc+jb,
                    pmp->fDQF+jb, pmp->Vol+jb, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_CGFLUID:  // CG fluid mixture (multicomponent)
        {
                TCGFcalc* myPT = new TCGFcalc( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL,
                    pmp->Pparc+jb, pmp->FWGT+k, pmp->X+jb,
                    pmp->fDQF+jb, pmp->Vol+jb, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_SRFLUID:  // SRK fluid mixture (multicomponent)
        {
                TSRKcalc* myPT = new TSRKcalc( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Pparc+jb,
                    pmp->fDQF+jb, pmp->Vol+jb, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_PR78FL:  // PR78 fluid mixture (multicomponent)
        {
                TPR78calc* myPT = new TPR78calc( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Pparc+jb,
                    pmp->fDQF+jb, pmp->Vol+jb, pmp->Tc, pmp->Pc );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_CORKFL:  // CORK fluid mixture (multicomponent)
        {
                TCORKcalc* myPT = new TCORKcalc( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
                    aIPx, aIPc, aDCc, aWx, alnGam, aphVOL, pmp->Pparc+jb,
                    pmp->fDQF+jb, pmp->Vol+jb, pmp->Tc, pmp->Pc, DCCp );
                mySM = (TSolMod*)myPT;
                break;
        }

        default:
                // mySM = new TSolMod( NComp, NPar, NPcoef, MaxOrd, NP_DC, ModCode, MixCode,
        	// aIPx, aIPc, aDCc,  aWx, alnGam, pmp->Tc, pmp->Pc );
        	break;
    }

  	if(phSolMod[k])
   	  delete phSolMod[k];
        phSolMod[k] = mySM; // set up new pointer for the solution model
}


// Wrapper call for calculation of temperature and pressure correction
// uses TSolMod class
void
TMulti::SolModParPT( long int k, char ModCode )
{
    // Extended constructor to connect to params, coeffs, and mole fractions
    switch( ModCode )
    {
		// solid and liquid solutions
        case SM_VANLAAR:
        case SM_REGULAR:
        case SM_GUGGENM:
        case SM_NRTLLIQ:
        case SM_WILSLIQ:

        // old solid solutions
        case SM_MARGT:
        case SM_MARGB:
        case SM_REDKIS:

        // aqueous DH models
        case SM_AQDH3:
        case SM_AQDH2:
        case SM_AQDH1:
        case SM_AQDHH:
        case SM_AQDHS:
        case SM_AQDAV:

        // aqueous models
        case SM_AQPITZ:
        case SM_AQSIT:
        case SM_AQEXUQ:

        // fluid (gas) models
        case SM_PRFLUID:
        case SM_CGFLUID:
        case SM_SRFLUID:
        case SM_PR78FL:
        case SM_CORKFL:

        {    ErrorIf( !phSolMod[k], "","Invalid index of phase");
              TSolMod* mySM = phSolMod[k];
              mySM->PTparam();
             break;
        }
        default:
              break;
    }
}


// Wrapper call for calculation of activity coefficients
// uses TSolMod class
void
TMulti::SolModActCoeff( long int k, char ModCode )
{
    switch( ModCode )
    {
		// solid and liquid solutions
        case SM_VANLAAR:
        case SM_REGULAR:
        case SM_GUGGENM:
        case SM_NRTLLIQ:
        case SM_WILSLIQ:

        // old solid solutions
        case SM_MARGT:
        case SM_MARGB:
        case SM_REDKIS:

        // aqueous DH models
        case SM_AQDH3:
        case SM_AQDH2:
        case SM_AQDH1:
        case SM_AQDHH:
        case SM_AQDHS:
        case SM_AQDAV:

        // aqueous models
        case SM_AQPITZ:
        case SM_AQSIT:
        case SM_AQEXUQ:

        // fluid (gas) models
        case SM_PRFLUID:
        case SM_CGFLUID:
        case SM_SRFLUID:
        case SM_PR78FL:
        case SM_CORKFL:

        {    ErrorIf( !phSolMod[k], "","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->MixMod();
             break;
        }
        default:
              break;
    }
}


// Wrapper call for calculation of bulk phase excess properties
// uses TSolMod class
void
TMulti::SolModExcessProp( long int k, char ModCode )
{
	// order of phase properties: G, H, S, CP, V, A, U
	long int j;
	double Gex, Hex, Sex, CPex, Vex, Aex, Uex;
	double zex[7];

	for (j =0; j<7; j++)
	{
		zex[j] = 0.0;
	}

	// insert cases for old solution and activity models
    switch( ModCode )
    {
        // solid and liquid solutions
        case SM_VANLAAR:
        case SM_REGULAR:
        case SM_GUGGENM:
        case SM_NRTLLIQ:
        case SM_WILSLIQ:

        // old solid solutions
        case SM_MARGT:
        case SM_MARGB:
        case SM_REDKIS:

        // aqueous DH models
        case SM_AQDH3:
        case SM_AQDH2:
        case SM_AQDH1:
        case SM_AQDHH:
        case SM_AQDHS:
        case SM_AQDAV:

        // aqueous models
        case SM_AQPITZ:
        case SM_AQSIT:
        case SM_AQEXUQ:

        // fluid (gas) models
        case SM_PRFLUID:
        case SM_CGFLUID:
        case SM_SRFLUID:
        case SM_PR78FL:
        case SM_CORKFL:

        {    ErrorIf( !phSolMod[k], "","Invalid index of phase");
              TSolMod* mySM = phSolMod[k];
              mySM->ExcessProp( zex );
              break;
        }
        default:
              break;
    }

    // assignments
    Gex = zex[0];
    Hex = zex[1];
    Sex = zex[2];
    CPex = zex[3];
    Vex = zex[4];
    Aex = zex[5];
    Uex = zex[6];
    pmp->GPh[k][2] = Gex;
    pmp->HPh[k][2] = Hex;
    pmp->SPh[k][2] = Sex;
    pmp->CPh[k][2] = CPex;
    pmp->VPh[k][2] = Vex;
    pmp->APh[k][2] = Aex;
    pmp->UPh[k][2] = Uex;

}


// Wrapper call for calculation of bulk phase ideal mixing properties
void
TMulti::SolModIdealProp( long int jb, long int k, char ModCode )
{
	// order of phase properties: G, H, S, CP, V, A, U
	long int j;
	double Gid, Hid, Sid, CPid, Vid, Aid, Uid;
	double zid[7];

	for (j=0; j<7; j++)
	{
		zid[j] = 0.0;
	}

	// needs to check for DC class state to catch ideal gas phase case
    switch( ModCode )
    {
		// check what solution phase (ideal gas?)
		// solid and liquid solutions
		case SM_VANLAAR:
		case SM_GUGGENM:
		case SM_REGULAR:
		case SM_NRTLLIQ:
		case SM_WILSLIQ:

		// old solid solutions
        case SM_MARGT:
        case SM_MARGB:
        case SM_REDKIS:

        // aqueous DH models
		case SM_AQDH3:
        case SM_AQDH2:
        case SM_AQDH1:
        case SM_AQDHH:
        case SM_AQDHS:
        case SM_AQDAV:

        // aqueous models
        case SM_AQPITZ:
        case SM_AQSIT:
        case SM_AQEXUQ:

        // fluid (gas) models
        case SM_PRFLUID:
        case SM_CGFLUID:
        case SM_SRFLUID:
        case SM_PR78FL:
        case SM_CORKFL:

        {    ErrorIf( !phSolMod[k], "","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->IdealProp( zid );
              break;
        }
		case SM_IDEAL:
		case SM_USERDEF:
		{
			IdealOneSite( jb, k, zid );
			break;
		}
		default:
			break;
    }

    // assignments
    Gid = zid[0];
    Hid = zid[1];
    Sid = zid[2];
    CPid = zid[3];
    Vid = zid[4];
    Aid = zid[5];
    Uid = zid[6];
    pmp->GPh[k][1] = Gid;
    pmp->HPh[k][1] = Hid;
    pmp->SPh[k][1] = Sid;
    pmp->CPh[k][1] = CPid;
    pmp->VPh[k][1] = Vid;
    pmp->APh[k][1] = Aid;
    pmp->UPh[k][1] = Uid;

}


// Wrapper call for retrieving bulk phase Darken quadratic terms
void
TMulti::SolModDarkenProp( long int jb, long int k, char ModCode )
{
	// order of phase properties: G, H, S, CP, V, A, U
	long int j;
	double Gdq, Hdq, Sdq, CPdq, Vdq, Adq, Udq;
	double zdq[7];

	for (j=0; j<7; j++)
	{
		zdq[j] = 0.0;
	}

	// data object for derivative properties needs to be added in Multi and DODs in scripts

	// assignments
    Gdq = zdq[0];
    Hdq = zdq[1];
    Sdq = zdq[2];
    CPdq = zdq[3];
    Vdq = zdq[4];
    Adq = zdq[5];
    Udq = zdq[6];
    pmp->GPh[k][3] = Gdq;
    pmp->HPh[k][3] = Hdq;
    pmp->SPh[k][3] = Sdq;
    pmp->CPh[k][3] = CPdq;
    pmp->VPh[k][3] = Vdq;
    pmp->APh[k][3] = Adq;
    pmp->UPh[k][3] = Udq;

}


// Wrapper call for retrieving bulk phase standard state terms
void
TMulti::SolModStandProp ( long int jb, long int k, char ModCode )
{
	// order of phase properties: G, H, S, CP, V, A, U
        double Gst=0., Hst=0., Sst=0., CPst=0., Vst=0., Ast=0., Ust=0.;

	// add if statement that checks DC class code (aqueous or not)

	// assignments
    pmp->GPh[k][0] = Gst;
    pmp->HPh[k][0] = Hst;
    pmp->SPh[k][0] = Sst;
    pmp->CPh[k][0] = CPst;
    pmp->VPh[k][0] = Vst;
    pmp->APh[k][0] = Ast;
    pmp->UPh[k][0] = Ust;

}


//-------------------------------------------------------------------------
// Internal memory allocation for TSolMod performance optimization
// (since version 2.3.0)
// added 26.11.2008 (SD)
void TMulti::Alloc_TSolMod( long int newFIs )
{
  if(  phSolMod && ( newFIs == sizeFIs) )
    return;

  Free_TSolMod();
  // alloc memory for all multicomponents phases
  phSolMod = new  TSolMod *[newFIs];
  sizeFIs = newFIs;
 for( long int ii=0; ii<newFIs; ii++ )
    	  phSolMod[ii] = 0;
}


void TMulti::Free_TSolMod()
{
  long int kk;

  if( phSolMod )
  {  for(  kk=0; kk<sizeFIs; kk++ )
      if( phSolMod[kk] )
           delete phSolMod[kk];

      delete[]  phSolMod;
  }
  phSolMod = 0;
  sizeFIs = 0;
}


// calculates ideal mixing properties of gas phases
void TMulti::IdealGas( long int jb, long int k, double *Zid )
{
	long int Nc, j;
	double T, P, R_T, R, s, sc, sp;
	double Gid, Hid, Sid, CPid, Vid, Aid, Uid;
	double *X;

	T = pmp->Tc;
	P = pmp->Pc;
	R_T = pmp->RT;
	R = R_T / T;
	Nc = pmp->L1[k];
	X = pmp->Wx+jb;

	s = 0.0;
	for (j=0; j<Nc; j++)
	{
		if ( X[j] > 1.0e-32 )
			s += X[j]*log(X[j]);
	}
	sc = (-1.)*R*s;
	sp = (-1.)*R*log(P);
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = sc + sp;

	// assignments (ideal mixing properties)
	Gid = Hid - Sid*T;
	Aid = Gid - Vid*P;
	Uid = Hid - Vid*P;
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
}


// calculates ideal mixing properties of condensed phases (one-site mixing)
void TMulti::IdealOneSite( long int jb, long int k, double *Zid )
{
	long int Nc, j;
	double T, P, R_T, R, s;
	double Gid, Hid, Sid, CPid, Vid, Aid, Uid;
	double *X;

	T = pmp->Tc;
	P = pmp->Pc;
	R_T = pmp->RT;
	R = R_T / T;
	Nc = pmp->L1[k];
	X = pmp->Wx+jb;

	s = 0.0;
	for (j=0; j<Nc; j++)
	{
		if ( X[j] > 1.0e-32 )
			s += X[j]*log(X[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R*s;
	Gid = Hid - T * Sid;
	Aid = Gid - Vid*P;
	Uid = Hid - Vid*P;

	// assignments
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
}


// calculates ideal mixing properties of condensed phases (multi-site mixing)
void TMulti::IdealMultiSite( long int jb, long int k, double *Zid )
{
	// reserved for future multi-site configurational property models
	//double Gid, Hid, Sid, CPid, Vid, Aid, Uid;


	// add assignments to Zid here
}


//--------------------- End of ipm_chemical3.cpp ---------------------------


