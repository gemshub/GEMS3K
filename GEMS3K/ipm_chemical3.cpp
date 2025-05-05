//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_chemical3.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, chemical potentials, etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, T.Wagner, S.Dmitrieva, K.Chudnenko
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

#include "ms_multi.h"
#include "v_detail.h"
#include "v_service.h"

/// Returns current value of smoothing factor for chemical potentials of highly non-ideal DCs
// added 18.06.2008 DK
double TMultiBase::SmoothingFactor( )
{
   if( pm.FitVar[4] < 0 )
   {  // To start SIA mode (smart initial approximation)
      return 1.0;
   }
   if( pm.FitVar[3] > 0 )
           return pm.FitVar[3];
   else
           return pm.FitVar[4];
}


/// New correction of smoothing factor for highly non-ideal systems.
// re-written 18.04.2009 DK+TW
/// Smoothing function choice: AG >= 0.0001 and DGC > -0.0001: old f(IT)
///                            AG >= 0.0001 and DGC <= -0.0001: new f(1/IT)
///                            AG <= -0.0001 and DGC <= -0.0001: new f(1/CD)
/// \param mode 0 - taking single log(CD) value for calculation of smoothing factor SF;
///       1, 2, ...  taking log(CD) average from the moving window of length mode
///       (up to 5 consecutive values)
///
void TMultiBase::SetSmoothingFactor( long int mode )
{
    double TF=1., al, ag, dg, iim, irf; // rg=0.0;
    long int ir; //, Level, itqF, itq;

    ir = pm.IT;
    irf = static_cast<double>(ir);
    ag = base_param()->AG; // pm.FitVar[4];
    dg = base_param()->DGC;
    iim = base_param()->IIM;

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
        if( pm.IT )
          inv_r = 1./static_cast<double>(pm.IT);
        logr = log( inv_r );
        logr_m = log( 1./iim );
        al = dg + ( ag - dg ) / ( 1. + exp( logr_m - logr ) / dg );
        if( !essentiallyEqual(ag, 1.0) ) {
           al += exp( log( 1. - ag ) + logr );
        }
        if( al > 1. )
      	    al = 1.;
        TF = al;
    }
    else if( dg <= -0.0001 && ag <= -0.0001 )
    {
    	double dk, cd;   long int i;
    	dg = fabs( dg );
    	ag = fabs( ag );
        dk = log( pm.DXM );
    	// Checking the mode where it is called
    	switch( mode )
    	{
    	  default:
          case 0: // MassBalanceRefinement() after SolveSimplex()
                     cd = log( pm.PCI );
    	   	     break;
    	  case 1:
    	  case 2:
    	  case 3:
    	  case 4:
    	  case 5: // Getting average (log geometric mean) from sampled CD values
    	  	     cd = 0.0;
    	   	     for(i=0; i < mode; i++ )
                         cd += pm.logCDvalues[i];
                     cd /= (static_cast<double>(mode)); // 5. - bugfix
    	   	     break;
    	}
        al = dg + ( ag - dg ) / ( 1. + exp( dk - cd ) / dg );
        if( !essentiallyEqual(ag, 1.0) ) {
           al += exp( log( 1. - ag ) + cd );
        }
        if( al > 1. )
      	    al = 1.;
        TF = al;
    }
//    if( pm.IT )
      pm.FitVar[3] = TF;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///  Function for converting internal lnGam[j] value into an external (phase-scale-specific)
///      Gamma[j] if DirFlag = 0 or external into internal value if DirFlag = 1.
///  Returns the respectively corrected external gamma activity coefficient or internal lnGam
///  Returns trivial values (lnGam = 0 or Gamma = 1) when the respective component
///    amount is zero (X[j] == 0) (is this a correct policy for zeroed-off components?)
//
double
TMultiBase::PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag )
{
    double NonLogTerm = 0., NonLogTermW = 0., NonLogTermS = 0.;//, MMC = 0.;

    if( k>=pm.FIs || pm.sMod[k][SPHAS_TYP] != SM_AQPITZ)
    {
        switch( pm.PHC[k] )
        {
        case PH_AQUEL:
            if( noZero( pm.XF[k] ) && noZero( pm.XFA[k] ) )
            {
                NonLogTerm = 1. - pm.XFA[k]/pm.XF[k];
                NonLogTermW = 2. - pm.XFA[k]/pm.XF[k] - pm.XF[k]/pm.XFA[k];
            }
            break;
        case PH_GASMIX:  case PH_FLUID:   case PH_PLASMA:   case PH_SIMELT:
        case PH_HCARBL:  case PH_SINCOND:  case PH_SINDIS:  case PH_LIQUID:
        case PH_IONEX:
            break;
        case PH_POLYEL:
        case PH_SORPTION: // only sorbent end-members!
            if( noZero( pm.XF[k] ) && noZero( pm.XFA[k] ) )
            {
//                for( long int jj=jb; jj<je; jj++ )
//                {
//                    if( pm.DCC[jj] == DC_SUR_CARRIER ||
//                            pm.DCC[jj] == DC_SUR_MINAL || pm.DCC[jj] == DC_PEL_CARRIER )
//                        MMC += pm.MM[jj]*pm.X[jj]/pm.XFA[k];
//                    // Weighted-average sorbent mole mass
//                }
                NonLogTerm = 1. - pm.XFA[k]/pm.XF[k];  // Also for sorption phases
                NonLogTermS = 2. - pm.XFA[k]/pm.XF[k] - pm.XF[k]/pm.XFA[k];
            }
            break;
        case PH_ADSORPT: // 'z' phase - adsorption on site balances
            break;
        default:
            break; // Phase class code error should be generated here!
        }
    }
#ifdef NOMUPNONLOGTERM
    NonLogTerm = 0.0;
    NonLogTermS = 0.0;
#endif
    if( DirFlag == 0 )
    {	 // Converting lnGam[j] into Gamma[j]
        if( approximatelyZero(pm.X[j]) && approximatelyZero(pm.XF[k]) )   // && !pm->XF[k]  added by DK 13.04.2012
            return 1.;
        double Gamma = 1.;
        double lnGamS = pm.lnGam[j];

        switch( pm.DCC[j] )
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
        case DC_SOL_MINDEP: case DC_SOL_MAJDEP:               case DC_SCM_SPECIES:
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
        if( fabs(lnGamS) < 706.8936 )
            Gamma = exp( lnGamS );
        else
            Gamma = 1.0; // Workaround - warning is needed
        return Gamma;
    }
    else { // Converting Gamma[j] into lnGam[j]
        if( approximatelyZero(pm.X[j]) && approximatelyZero(pm.XF[k]) )   // && !pm->XF[k]  added by DK 13.04.2012
            return 0.;
        double Gamma = pm.Gamma[j];
        double lnGam = 0.0;  // Cleanup by DK 5.12.2009
        if( !essentiallyEqual(Gamma, 1.0) && Gamma > pm.lowPosNum )
            lnGam = log( Gamma );
        switch( pm.DCC[j] )
        { // Aqueous electrolyte
        case DC_AQ_PROTON: case DC_AQ_ELECTRON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
            lnGam -= NonLogTerm;  // Correction by asymmetry term
            break;
        case DC_AQ_SOLVCOM:	    case DC_AQ_SOLVENT:
            lnGam -= NonLogTermW;
            break;
        case DC_GAS_COMP: case DC_GAS_H2O: case DC_GAS_CO2: case DC_GAS_H2: case DC_GAS_N2:
            break;
        case DC_SOL_IDEAL:  case DC_SOL_MINOR:  case DC_SOL_MAJOR: case DC_SOL_MINDEP:
        case DC_SOL_MAJDEP:                   case DC_SCM_SPECIES:
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
/// Main call point for calculation of DC activity coefficients (lnGam vector)
///    formerly GammaCalc().
/// Controls various built-in models, as well as generic Phase script calculation
/// LinkMode is a parameter indicating the status of Gamma calculations:
/// LINK_TP_MODE - calculation of equations depending on TP only;
/// LINK_UX_MODE - calculation of equations depending on current
///      IPM approximation of the equilibrium state;
/// LINK_PP_MODE - calculation of integral phase properties after GEMIPM has converged
///		needs to be implemented
/// \return status code (0 if o.k., non-zero values if there were problems
///     with surface complexation models)
long int
TMultiBase::CalculateActivityCoefficients( long int LinkMode  )
{
    long int k, j, jb, je=0, jpb, jdb, ipb,  jpe=0, jde=0, ipe=0;
    //long int  jmb, jme=0, jsb, jse=0;
    //long int jphl=0, jlphc=0, jdqfc=0,  jrcpc=0;
    char *sMod;
    long int statusGam=0, statusGC=0, statusSACT=0, SmMode = 0;
    double LnGam, pmpXFk;
    const BASE_PARAM *pa_p = base_param();

   ipm_logger->trace("CalculateActivityCoefficients {}", LinkMode);
    // calculating concentrations of species in multi-component phases
    switch( LinkMode )
    {
      case LINK_TP_MODE:  // Built-in functions depending on T,P only
      {
        long int  jdqfc=0,  jrcpc=0; // jphl=0, jlphc=0,
        long int  jmb, jme=0, jsb, jse=0;

        pm.FitVar[3] = 1.0;  // resetting the IPM smoothing factor

         for( k=0; k<pm.FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += pm.L1[k];
            if( pm.L1[k] == 1 && !( pm.PHC[k] == PH_GASMIX ||
                                    pm.PHC[k] == PH_PLASMA ||
                                    pm.PHC[k] == PH_FLUID ))  // SD 13/12/19
                 continue;
            // Indexes for extracting data from IPx, PMc and DMc arrays
            ipb = ipe;
            ipe += pm.LsMod[k*3]*pm.LsMod[k*3+1];
            jpb = jpe;
            jpe += pm.LsMod[k*3]*pm.LsMod[k*3+2];
            jdb = jde;
            jde += pm.LsMdc[k*3]*pm.L1[k];
            sMod = pm.sMod[k];

            jmb = jme;
            jme += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2]*pm.L1[k];
            jsb = jse;
            jse += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2];

                    double nxk = 1./pm.L1[k];
            for( j= jb; j<je; j++ )
    		{
                if(pm.XF[k] < std::min( pm.DSM, pm.PhMinM ) ) // pm.lowPosNum )   // workaround 10.03.2008 DK
                        pm.Wx[j] = nxk;  // need this eventually to avoid problems with zero mole fractions
                pm.fDQF[j] =0.0;  // cleaning fDQF in TP mode!
                pm.lnGmo[j] = pm.lnGam[j]; // saving activity coefficients in TP mode
       	    }
                // if( sMod[SGM_MODE] != SM_STNGAM ) This should not be the case anymore DK 24.11.2010
                // continue;  // The switch below is for built-in functions only!

            // the following section probably needs to be re-written to allow more flexible
            // combinations of fluid models for pure gases with gE mixing models,
            // scheme should probably be the same as in LINK_UX_MODE, 03.06.2008 (TW)
            switch( pm.PHC[k] )
            {
                case PH_AQUEL: case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
                case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID: case PH_ADSORPT:
                case PH_IONEX:
                    SolModCreate( jb, jmb, jsb, jpb, jdb, k, ipb,
                        sMod[SPHAS_TYP], sMod[MIX_TYP], /* jphl, jlphc,*/ jdqfc/*,  jrcpc*/  );
                    // new solution models (TW, DK 2007)
            	    SolModParPT( k, sMod[SPHAS_TYP] );
            	    break;
              default:
                    break;
            }

            // move pointers
 //           jphl  +=  (pm.LsPhl[k*2]*2);
 //           jlphc += (pm.LsPhl[k*2 ]*pm.LsPhl[k*2+1]);
            jdqfc += (pm.LsMdc2[k*3]*pm.L1[k]);
            jrcpc += (pm.LsMdc2[k*3+1]*pm.L1[k]);

          } // k
        }
        ipm_logger->trace("CalculateActivityCoefficients - LINK_TP_MODE");
        break;

      case LINK_PP_MODE: // Mode of calculation of integral solution phase properties
        for( k=0; k<pm.FIs; k++ )
        { // loop on solution phases
            jb = je;
            je += pm.L1[k];
            sMod = pm.sMod[k];
            switch( pm.PHC[k] )
            {
              case PH_AQUEL: case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
              case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:  case PH_ADSORPT:
              case PH_IONEX:
                  ipm_logger->trace("CalculateActivityCoefficients - LINK_PP_MODE");
                   SolModExcessProp( k, sMod[SPHAS_TYP] ); // extracting integral phase properties
                   SolModIdealProp(  /*jb,*/ k, sMod[SPHAS_TYP] );
                   SolModStandProp(  /*jb,*/ k, sMod[SPHAS_TYP] );
                   SolModDarkenProp( /*jb,*/ k/*, sMod[SPHAS_TYP]*/ );
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
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
        // cleaning activity coefficients
        for( j=0; j<pm.L; j++ )
        {
            pm.lnGam[j] = 0.;
            pm.Gamma[j] = 1.;
        }
        if( pm.E && pm.LO ) // checking electrostatics
        {
          IS_EtaCalc();  //  calculating charges and charge densities
          if( pm.FIat > 0 )
             for( k=0; k<pm.FIs; k++ )
             {
               if( pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION )
               {  long int ist;
                  for( ist=0; ist<pm.FIat; ist++ ) // loop over surface types
                  {
                     pm.XpsiA[k][ist] = 0.0;        // cleaning Psi before GouyChapman()
                     pm.XpsiB[k][ist] = 0.0;
                     pm.XpsiD[k][ist] = 0.0;
                  }  // ist
                }
             }  // k
         } // pm.E
        break;
    default:
        Error("CalculateActivityCoefficients()","Invalid LinkMode for a built-in solution model");
    }

    jpe=0; jde=0; ipe=0;
    je=0;
    for( k=0; k<pm.FI; k++ )
    { // loop on phases
        jb = je;
        je += pm.L1[k];
        if( pm.L1[k] == 1 )
            goto END_LOOP;
        sMod = pm.sMod[k];
            // if( sMod[SGM_MODE] == SM_IDEAL )  Comm.out 29.11.2010 by DK to introduce multi-site ideal models
            // goto END_LOOP;
        pmpXFk = 0.;  // Added 07.01.05 by KD
        for( j = jb; j < je; j++ )
            pmpXFk += pm.X[j];
        if( pm.XF[k] < pm.DSM ) // Bugfix by KD 09.08.2005 (bug report Th.Matschei)
            pm.XF[k] = pmpXFk;

        // Indexes for extracting data from IPx, PMc and DMc arrays
        ipb = ipe;                  // added 07.12.2006 by KD
        ipe += pm.LsMod[k*3]*pm.LsMod[k*3+1];
        jpb = jpe;
        jpe += pm.LsMod[k*3]*pm.LsMod[k*3+2];  // Changed 07.12.2006  by KD
        jdb = jde;
        jde += pm.LsMdc[k*3]*pm.L1[k];
        //jmb = jme;
        //jme += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2]*pm.L1[k];
        //jsb = jse;
        //jse += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2];

   if( LinkMode == LINK_UX_MODE && sMod[SGM_MODE] == SM_STNGAM )
   {    // check that SGM_MODE for adsorption or multi-site ideal SS is not SM_IDEAL in Phase records!
        switch( pm.PHC[k] )
        {  // calculating activity coefficients using built-in functions
          case PH_AQUEL:   // DH III variant consistent with HKF
             if( pmpXFk > pm.DSM && pm.X[pm.LO] > pm.XwMinM && pm.IC > pa_p->ICmin )
             {
                switch( sMod[SPHAS_TYP] )
                {
                    case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHS: case SM_AQDHH:
                    case SM_AQDAV: case SM_AQSIT: case SM_AQPITZ: case SM_AQEXUQ: case SM_AQELVIS:
						SolModActCoeff( k, sMod[SPHAS_TYP] );
						break;
					default:
						break;
                }
             }
             goto END_LOOP;
             break;
          case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
             if( pmpXFk > pm.DSM && pm.XF[k] >pa_p->PhMin )
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
                 if ( sMod[SPHAS_TYP] == SM_STFLUID )  // STP EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
             }
             goto END_LOOP;
             break;
         case PH_LIQUID: case PH_SIMELT: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
         case PH_ADSORPT: case PH_IONEX:
             if( pmpXFk > pm.DSM )
             {     // solid and liquid mixtures
                switch( sMod[SPHAS_TYP] )
                {
                    case SM_IDEAL:   // Ideal (multi-site) model (DK 29.11.2010)
                    case SM_BERMAN:  // Non-ideal (multi-site) model (DK 07.12.2010)
                    case SM_CEF:     // multi-site non-ideal ss model (CALPHAD) DK 15.08.2014
                    case SM_MBW:     // Modified Bragg-Williams model, [Vinograd et al. 2018]
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
                    case SM_SURCOM:   // new SCMs with site balances
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
                if( pm.PHC[0] == PH_AQUEL && pmpXFk > pm.DSM
                && (pm.XFA[0] > pm.XwMinM && pm.XF[0] > pm.DSM ))
                {
                    if( pm.E )
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

   if( !calculateActivityCoefficients_scripts( LinkMode, k, jb, jpb, jdb, ipb, pmpXFk ) )
       goto END_LOOP;

END_LOOP:
        if( LinkMode == LINK_TP_MODE )  // TP mode - added 04.03.2008 by DK
        {
        	for( j=jb; j<je; j++ )
        	{
                   if( pm.XF[k] < pm.DSM )   // workaround 10.03.2008 DK
                                pm.Wx[j] = 0.0;               //
                   LnGam = pm.lnGmo[j];
                   pm.lnGam[j] = LnGam;
                   if(  fabs( LnGam ) < 84. )
                       pm.Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
                   else pm.Gamma[j] = 1.0;
        	}
        }
        else if(LinkMode == LINK_UX_MODE )  // Bugfix! DK 06.04.11
        { // Real mode for activity coefficients
           double lnGamG;
           for( j=jb; j<je; j++ )
           {
             if( pm.DCC[j] == DC_AQ_SURCOMP )  // Workaround for aqueous surface complexes DK 22.07.09
             {   pm.lnGam[j] = 0.0; lnGamG = 0.0; }
             else
                 lnGamG = PhaseSpecificGamma( j, jb, je, k, 1 );
             LnGam = pm.lnGam[j];
             if( fabs( lnGamG ) > 1e-9 )
            	LnGam += lnGamG;
             pm.lnGmo[j] = LnGam;
             if( fabs( LnGam ) < 84. )   // before 26.02.08: < 42.
                    pm.Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
             else pm.Gamma[j] = 1.0;
             if( pm.DCC[j] == DC_AQ_SURCOMP )  // bugfix 18.11.2017 by DK for K aq species
                  pm.Gamma[j] = 1.0;
             pm.F0[j] = DC_PrimalChemicalPotentialUpdate( j, k );
             pm.G[j] = pm.G0[j] + pm.fDQF[j] + pm.F0[j];
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
/// Wrapper calls for creating multi-component mixing models for phases
/// using  TSolMod class. Now including multi-site ideal and scripted models
//
void TMultiBase::SolModCreate( long int jb, long int jmb, long int jsb, long int jpb, long int jdb,
                           long int k, long int ipb, char ModCode, char MixCode,
                           /* long int jphl, long int jlphc, */ long int jdqfc/*, long int  jrcpc*/)
{
    double *aZ, *aM;//, *aVol;
    //long int *aIPx;
    //char *DCCp;
    SolutionData sd;

    sd.NSpecies = pm.L1[k];          // Number of components (end members) in the phase
    sd.NParams = pm.LsMod[k*3];      // Number of interaction parameters
    sd.NPcoefs = pm.LsMod[k*3+2];    // and number of coefs per parameter in PMc table
    sd.MaxOrder =  pm.LsMod[k*3+1];  // max. parameter order (cols in IPx)
    sd.NPperDC = pm.LsMdc[k*3];      // Number of non-ideality coeffs per one DC in multicomponent phase
    sd.NSublat = pm.LsMdc[k*3+1];    // Number of site types (sublattices) for multi-site SS model
    sd.NMoiet = pm.LsMdc[k*3+2];     // Number of moieties for multi-site SS model
    sd.Mod_Code = ModCode;
    sd.Mix_Code = MixCode;

    //new objects to Phase 06/06/12
//    sd.NlPhs = pm.LsPhl[k*2];
//    sd.NlPhC = pm.LsPhl[k*2+1];
    sd.NDQFpDC = pm.LsMdc2[k*3];
//    sd.NrcPpDC = pm.LsMdc2[k*3+1];

    if( phSolMod[k])
        if(  phSolMod[k]->testSizes( &sd ) )
    	{
                phSolMod[k]->UpdatePT( pm.Tc, pm.P );
                return; // using old allocation and setup of the solution model
    	}

    // properties generic to all models
    sd.arIPx = pm.IPx+ipb;   // Pointer to list of indexes for non-ideal solutions -> NPar x MaxOrd
    sd.arIPc = pm.PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    sd.arDCc = pm.DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    sd.arWx = pm.Wx+jb;       // End member mole fractions
    sd.arlnGam = pm.lnGam+jb; // End member ln activity coeffs

    sd.arlnDQFt = pm.lnDQFt+jb; // End member ln activity coeffs
    sd.arlnRcpt = pm.lnRcpt+jb; // End member ln activity coeffs
    sd.arlnExet = pm.lnExet+jb; // End member ln activity coeffs
    sd.arlnCnft = pm.lnCnft+jb; // End member ln activity coeffs
    sd.arCTermt = pm.CTerms+jb; // End member coulombic terms

    sd.aphVOL = pm.FVOL+k;
    sd.DC_Codes = pm.DCC+jb;  // pointer to Dcomp class codes (added 02.05.2010 TW)
    sd.arMoiSN = pm.MoiSN+jmb;  // Pointer to sublattice-moiety multiplicity array
    sd.arSitFr = pm.SitFr+jsb;  // Pointer to sublattice-moiety multiplicity array
    sd.arGEX = pm.fDQF+jb;      // DQF parameters or pure-gas fugacities
    sd.arPparc = pm.Pparc+jb;
    sd.TP_Code = &pm.dcMod[jb];
    sd.T_k = pm.Tc;
    sd.P_bar = pm.P;

    //new objects to Phase 06/06/12
//    sd.arPhLin = pm.PhLin+jphl;
//    sd.lPhc = pm.lPhc+ jlphc;
    sd.arDQFc = pm.DQFc+ jdqfc;
//    sd.rcpc = pm.rcpc+ jrcpc;

    // specific properties
    aM = pm.Y_m+jb;
    aZ = pm.EZ+jb;
    sd.arVol = pm.Vol+jb;
    sd.arSM = pm.SM+jb;
    auto phase_name = char_array_to_string(pm.SF[k]+MAXSYMB, MAXPHNAME);
    strip(phase_name);
    sd.phaseName = phase_name;

    TSolMod* mySM = nullptr;

   // creating instances of subclasses of TSolMod base class
    switch( ModCode )
    {

        case SM_OTHER:  // Hard-coded solid solution models (selected by phase name)
        {
                TModOther* myPT = new TModOther( &sd, pm.denW, pm.epsW );
                mySM = myPT;
                break;
        }

        case SM_VANLAAR:  // Van Laar solid solution model (multicomponent)
        {
                TVanLaar* myPT = new TVanLaar( &sd );
                mySM = myPT;
                break;
        }
             // break;

        case SM_REGULAR:  // Regular solid solution model (multicomponent)
        {
                TRegular* myPT = new TRegular( &sd );
                mySM = myPT;
                break;
        }

        case SM_GUGGENM:  // Redlich-Kister solid solution model (multicomponent)
        {
                TRedlichKister* myPT = new TRedlichKister( &sd );
                mySM = myPT;
                break;
        }

        case SM_NRTLLIQ:  // NRTL liquid solution model (multicomponent)
        {
                TNRTL* myPT = new TNRTL( &sd );
                mySM = myPT;
                break;
        }

        case SM_WILSLIQ:  // Wilson liquid solution model (multicomponent)
        {
                TWilson* myPT = new TWilson( &sd );
                mySM = myPT;
                break;
        }

        case SM_MARGT:  // Margules ternary (regular) solid solution model
        {
                TMargules* myPT = new TMargules( &sd );
                mySM = myPT;
                break;
        }

        case SM_MARGB:  // Margules binary (subregular) solid solution model
        {
                TSubregular* myPT = new TSubregular( &sd );
                mySM = myPT;
                break;
        }

        case SM_REDKIS:  // Gugenheim binary (REdlich-Kister) solid solution
        {
                TGuggenheim* myPT = new TGuggenheim( &sd );
                mySM = myPT;
                break;
        }

        case SM_AQPITZ:  // Pitzer aqueous electrolyte model (multicomponent)
        {
                TPitzer* myPT = new TPitzer( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
                break;
        }

        case SM_AQSIT:  // SIT aqueous electrolyte model (multicomponent)
        {
                TSIT* myPT = new TSIT( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
                break;
        }

        case SM_AQEXUQ:  // EUNIQUAC aqueous electrolyte model (multicomponent)
        {
                TEUNIQUAC* myPT = new TEUNIQUAC( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
                break;
        }
/*
        case SM_AQELVIS:  // ELVIS aqueous electrolyte model (multicomponent)
        {
                TELVIS* myPT = new TELVIS( &sd, aM, aZ, pm.denW, pm.epsW );
		mySM = (TSolMod*)myPT;
                break;
        }
*/
        case SM_AQDH3:  // extended Debye-Hueckel aqueous electrolyte model (Karpov version)
        {
                TKarpov* myPT = new TKarpov( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_AQDH2:   // Debye-Hueckel aqueous electrolyte model
        {
                TDebyeHueckel* myPT = new TDebyeHueckel( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_AQDH1:   // Debye-Hueckel limiting law aqueous electrolyte model
        {
                TLimitingLaw* myPT = new TLimitingLaw( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_AQDHS:  // extended Debye-Hueckel aqueous electrolyte model (Shvarov version)
        {
                TShvarov* myPT = new TShvarov( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_AQDHH:  // extended Debye-Hueckel aqueous electrolyte model (Helgeson version)
        {
                THelgeson* myPT = new THelgeson( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_AQDAV:  // Davies aqueous electrolyte model (in NEA TDB version)
        {
                TDavies* myPT = new TDavies( &sd, aM, aZ, pm.denW, pm.epsW );
                mySM = myPT;
        	break;
        }

        case SM_PRFLUID:  // PRSV fluid mixture (multicomponent)
        {
                TPRSVcalc* myPT = new TPRSVcalc( &sd );
                mySM = myPT;
                break;
        }

        case SM_CGFLUID:  // CG fluid mixture (multicomponent)
        {
                TCGFcalc* myPT = new TCGFcalc( &sd, pm.FWGT+k, pm.X+jb  );
                mySM = myPT;
                break;
        }

        case SM_SRFLUID:  // SRK fluid mixture (multicomponent)
        {
                TSRKcalc* myPT = new TSRKcalc( &sd );
                mySM = myPT;
                break;
        }

        case SM_PR78FL:  // PR78 fluid mixture (multicomponent)
        {
                TPR78calc* myPT = new TPR78calc( &sd );
                mySM = myPT;
                break;
        }

        case SM_CORKFL:  // CORK fluid mixture (multicomponent)
        {
                TCORKcalc* myPT = new TCORKcalc( &sd );
                mySM = myPT;
                break;
        }

        case SM_STFLUID:  // STP fluid mixture (H2O-CO2)
        {
                TSTPcalc* myPT = new TSTPcalc ( &sd );
                mySM = myPT;
                break;
        }

        case SM_BERMAN:  // Non-ideal (multi-site) model
        {
                TBerman* myPT = new TBerman( &sd, pm.G0+jb );
                mySM = myPT;
                break;
        }
        case SM_CEF:  // Non-ideal (multi-site) model (CALPHAD)
        {
                TCEFmod* myPT = new TCEFmod( &sd, pm.G0+jb );
                mySM = myPT;
                break;
        }

        case SM_MBW:
        {
                TMBWmod* myPT = new TMBWmod( &sd, pm.G0+jb );
                mySM = myPT;
                break;
        }

        case SM_IDEAL:
        {
                TIdeal* myPT = new TIdeal( &sd );
                mySM = myPT;
                break;
        }
        // case SM_USERDEF:
        case SM_SURCOM:
        {
            TSCM_NEM* myPT = new TSCM_NEM( &sd );
            mySM = myPT;
            break;
        }
        default:
        	break;
    }
  	if(phSolMod[k])
            delete phSolMod[k];

     phSolMod[k] = mySM; // set up new pointer for the solution model
    if(TSolMod::solmod_logger->should_log(spdlog::level::debug)) {
     phSolMod[k]->to_json_file(std::string("solmod_")+std::to_string(k)+".json");
    }
}

/// Wrapper call for calculation of temperature and pressure correction
/// uses TSolMod class
void TMultiBase::SolModParPT( long int k, char ModCode )
{
    // Extended constructor to connect to params, coeffs, and mole fractions
    switch( ModCode )
    {  // solid and liquid solutions
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
    case SM_BERMAN:  case SM_CEF:  /* new ss models */   case SM_MBW: /* new mbw model */        case SM_SURCOM:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {
              ErrorIf( !phSolMod[k], "SolModParPT: ","Invalid index of phase");
              TSolMod* mySM = phSolMod[k];
              mySM->PTparam();
             break;
        }
        default:
              break;
    }
}

/// Wrapper call for calculation of activity coefficients
/// uses TSolMod class
void TMultiBase::SolModActCoeff( long int k, char ModCode )
{
    switch( ModCode )
    {   // solid and liquid solutions
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
    case SM_BERMAN:  case SM_CEF:  case SM_MBW:   case SM_SURCOM:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {
             ErrorIf( !phSolMod[k], "SolModActCoeff: ","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->MixMod();
             if(TSolMod::solmod_logger->should_log(spdlog::level::debug)) {
                 phSolMod[k]->to_text_file(std::string("solmod_act_coef_")+std::to_string(k)+".txt", true);
             }
             break;
        }
        default:
              break;
    }
}

/// Wrapper call for calculation of bulk phase excess properties
/// uses TSolMod class
void TMultiBase::SolModExcessProp( long int k, char ModCode )
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
        // case SM_USERDEF:  case SM_IDEAL:
        case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:   case SM_MBW:  case SM_SURCOM:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModExcessProp: ","Invalid index of phase");
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
    ipm_logger->debug("Assignment of calculated excess properties of mixing [2] Gex={} Hex={} Sex={} CPex={}", Gex, Hex, Sex, CPex);
    pm.GPh[k][2] = Gex;
    pm.HPh[k][2] = Hex;
    pm.SPh[k][2] = Sex;
    pm.CPh[k][2] = CPex;
    pm.VPh[k][2] = Vex;
    pm.APh[k][2] = Aex;
    pm.UPh[k][2] = Uex;
}


/// Wrapper call for calculation of bulk phase ideal mixing properties
void TMultiBase::SolModIdealProp( /*long int jb,*/ long int k, char ModCode )
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
    {   // check what solution phase (ideal gas?)
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:    case SM_MBW:   case SM_SURCOM:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModIdealProp: ","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->IdealProp( zid );
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
    ipm_logger->debug("Assignment of calculated excess properties of mixing [1] Gid={} Hid={} Sid={} CPid={}", Gid, Hid, Sid, CPid);
    pm.GPh[k][1] = Gid;
    pm.HPh[k][1] = Hid;
    pm.SPh[k][1] = Sid;
    pm.CPh[k][1] = CPid;
    pm.VPh[k][1] = Vid;
    pm.APh[k][1] = Aid;
    pm.UPh[k][1] = Uid;
}

/// Wrapper call for retrieving bulk phase Darken quadratic terms
void TMultiBase::SolModDarkenProp( /*long int jb,*/ long int k/*, char ModCode*/ )
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
    ipm_logger->debug("Assignment of calculated Darken terms of mixing [2] Gdq={} Hdq={} Sdq={} CPdq={}", Gdq, Hdq, Sdq, CPdq);
    pm.GPh[k][3] = Gdq;
    pm.HPh[k][3] = Hdq;
    pm.SPh[k][3] = Sdq;
    pm.CPh[k][3] = CPdq;
    pm.VPh[k][3] = Vdq;
    pm.APh[k][3] = Adq;
    pm.UPh[k][3] = Udq;

}

/// Wrapper call for retrieving bulk phase standard state terms
void TMultiBase::SolModStandProp ( /*long int jb,*/ long int k, char ModCode )
{
    // order of phase properties: G, H, S, CP, V, A, U
    long int j;
    double Gst=0., Hst=0., Sst=0., CPst=0., Vst=0., Ast=0., Ust=0.;
    double zst[7];

    // add if statement that checks DC class code (aqueous or not)
    for (j=0; j<7; j++)
    {
        zst[j] = 0.0;
    }
    // needs to check for DC class state to catch ideal gas phase case
    switch( ModCode )
    {   // check what solution phase (ideal gas?)
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:    case SM_MBW:   case SM_SURCOM:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModIdealProp: ","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->StandardProp( zst );
             break;
        }
        default:
            break;
    }
    // assignments
    Gst = zst[0];
    Hst = zst[1];
    Sst = zst[2];
    CPst = zst[3];
    Vst = zst[4];
    Ast = zst[5];
    Ust = zst[6];
    ipm_logger->debug("Assignment of calculated standard (reference) properties of mixed phase [0] Gst={} Hst={} Sst={} CPst={}", Gst, Hst, Sst, CPst);
    pm.GPh[k][0] = Gst;
    pm.HPh[k][0] = Hst;
    pm.SPh[k][0] = Sst;
    pm.CPh[k][0] = CPst;
    pm.VPh[k][0] = Vst;
    pm.APh[k][0] = Ast;
    pm.UPh[k][0] = Ust;

}

//-------------------------------------------------------------------------
/// Internal memory allocation for TSolMod performance optimization
/// (since version 2.3.0)
void TMultiBase::Alloc_TSolMod( long int newFIs )
{
  if(  phSolMod && ( newFIs == sizeFIs) )
    return;

  Free_TSolMod();
  // alloc memory for all multicomponents phases
  phSolMod = new  TSolMod *[newFIs];
  sizeFIs = newFIs;
 for( long int ii=0; ii<newFIs; ii++ )
          phSolMod[ii] = nullptr;
}

void TMultiBase::Free_TSolMod()
{
  long int kk;

  if( phSolMod )
  {
      for(  kk=0; kk<sizeFIs; kk++ )
        if( phSolMod[kk] )
           delete phSolMod[kk];
      delete[]  phSolMod;
  }
  phSolMod = nullptr;
  sizeFIs = 0;
}

/// Internal memory allocation for TSorpMod performance optimization
/// (since version 3.4.0)
void TMultiBase::Alloc_TSorpMod( long int newFIs )
{
  if(  phSorpMod && ( newFIs == sizeFIa) )
    return;

  Free_TSorpMod();
  // alloc memory for all multicomponents phases
  phSorpMod = new  TSorpMod *[newFIs];
  sizeFIa = newFIs;
  for( long int ii=0; ii<newFIs; ii++ )
          phSorpMod[ii] = nullptr;
}

void TMultiBase::Free_TSorpMod()
{
  long int kk;

  if( phSorpMod )
  {
      for(  kk=0; kk<sizeFIa; kk++ )
         if( phSorpMod[kk] )
           delete phSorpMod[kk];
      delete[]  phSorpMod;
  }
  phSorpMod = nullptr;
  sizeFIa = 0;
}

//--------------------- End of ipm_chemical3.cpp ---------------------------


