//-------------------------------------------------------------------
// $Id: s_fglbred.cpp 1140 2008-12-08 19:07:05Z wagner $
//
// Copyright (C) 2008  Th.Wagner, D.Kulik, S.Dmitrieva
//
// Implementation of  class
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "verror.h"
#include "s_fgl.h"


    //=============================================================================================
    // Van Laar model for solid solutions (c) TW March 2007
    // References:  Holland & Powell (2003)
    //=============================================================================================


    // Generic constructor for the TVanLaar class
    TModOther::TModOther( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
            long int NPperDC, char Mod_Code,
            long int *arIPx, double *arIPc, double *arDCc,
            double *arWx, double *arlnGam, double *aphVOL,
            double T_k, double P_bar, double dW, double eW ):
            	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
            			 Mod_Code, arIPx, arIPc, arDCc, arWx,
            			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
    {
      alloc_internal();

    }


    TModOther::~TModOther()
    {
      free_internal();
    }


    void TModOther::alloc_internal()
    {
    	   Wu = new double [NPar];
    	   Ws = new double [NPar];
    	   Wv = new double [NPar];
    	   Wpt = new double [NPar];
    	   Phi = new double [NComp];
    	   PsVol = new double [NComp];
    }


    void TModOther::free_internal()
    {
    	 if(Wu)  delete[]Wu;
    	 if(Ws)  delete[]Ws;
    	 if(Wv)  delete[]Wv;
    	 if(Wpt)  delete[]Wpt;
    	 if(Phi)  delete[]Phi;
    	 if(PsVol)  delete[]PsVol;
    }


    // Calculates T,P corrected binary interaction parameters
    long int TModOther::PTparam()
    {



    	/*
    	    if( !strncmp( PhaseName, "FELDSPAR", 8 ))
    	    	{
    	    		TModOther::PT_Feldspar();
    	    	} else if{ !strncmp( PhaseName, "Garnet", 6 ))
    	    	     TModOther::PT_Garnet();
    	    	} esle if ....
    	*/


    	long int ip;

        if ( NPcoef < 3 || NPar < 1 )
           return 1;

    // read P-T corrected interaction parameters
    	   for (ip=0; ip<NPar; ip++)
    	   {
    	     Wu[ip] = aIPc[NPcoef*ip];
    		 Ws[ip] = aIPc[NPcoef*ip+1];
    		 Wv[ip] = aIPc[NPcoef*ip+2];
    		 Wpt[ip] = Wu[ip]+ Ws[ip]*Tk + Wv[ip]*Pbar;
    	     aIP[ip] = Wpt[ip];
    		 // aIPc[NPcoef*ip+3] = Wpt[ip]; // obsolete
    	   }
    	   return 0;
    }


    // Calculates activity coefficients and excess functions
    long int TModOther::MixMod()
    {


       	/*
        	    if( !strncmp( PhaseName, "FELDSPAR", 8 ))
        	    	{
        	    		TModOther::MixMod_Feldspar();
        	    	} else if{ !strncmp( PhaseName, "Garnet", 6 ))
        	    	     TModOther::MixMod_Garnet();
        	    	} else if ....
        	*/


    	long int ip, j, i1, i2;
       double dj, dk;
       double sumPhi; // Sum of Phi terms
       double gE, vE, hE, sE, cpE, uE;

       if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
               return 1;

       // calculating Phi values
       sumPhi = 0.;
       for (j=0; j<NComp; j++)
       {
           PsVol[j] = aDCc[NP_DC*j];  // reading pseudo-volumes
           sumPhi +=  x[j]*PsVol[j];
       }

       if( fabs(sumPhi) < 1e-30 )
           return 2;    // to prevent zerodivide!

       for (j=0; j<NComp; j++)
           Phi[j] = x[j]*PsVol[j]/sumPhi;

       // calculate activity coefficients
       for (j=0; j<NComp; j++)      // index end members with j
       {
    	lnGamRT = 0.;
    	for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
    	{
            i1 = aIPx[MaxOrd*ip];
    	    i2 = aIPx[MaxOrd*ip+1];

       	    if( j == i1 )
    		dj = 1.;
    	    else
    		dj = 0.;
    	    if( j == i2 )
    		dk = 1.;
    	    else
    		dk = 0.;
    	    lnGamRT -= (dj-Phi[i1])*(dk-Phi[i2])*Wpt[ip]
                             *2.*PsVol[j]/(PsVol[i1]+PsVol[i2]);
    	}

        lnGam = lnGamRT/(R_CONST*Tk);
    	lnGamma[j] = lnGam;
    	}

       // calculate bulk phase excess properties
       gE = 0.0;
       vE = 0.0;
       hE = 0.0;
       sE = 0.0;
       cpE = 0.0;
       uE = 0.0;

       for (ip=0; ip<NPar; ip++)
       {
          i1 = aIPx[MaxOrd*ip];
          i2 = aIPx[MaxOrd*ip+1];
          gE += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wpt[ip];
          vE += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wv[ip];
          uE += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wu[ip];
          sE -= Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Ws[ip];
       }

       hE = uE+vE*Pbar;

       return 0;
    }






//--------------------- End of s_fgl3.cpp ---------------------------

