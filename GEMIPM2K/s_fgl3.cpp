//-------------------------------------------------------------------
// $Id: s_fglbred.cpp 1140 2008-12-08 19:07:05Z wagner $
//
// Copyright (C) 2008-2009  T.Wagner, D.Kulik, S.Dmitrieva
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
// Customized solid-solution and fluid models (c) TW March 2009
// References:
//=============================================================================================


// Generic constructor for the TVanLaar class
TModOther::TModOther( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
		long int NPperDC, char Mod_Code,
		long int *arIPx, double *arIPc, double *arDCc,
		double *arWx, double *arlnGam, double *aphVOL,
		double T_k, double P_bar, double *dW, double *eW ):
            	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
            			Mod_Code, arIPx, arIPc, arDCc, arWx,
            			arlnGam, aphVOL, T_k, P_bar )
{
	alloc_internal();
}


TModOther::~TModOther()
{
	free_internal();
}


void TModOther::alloc_internal()
{
    Gdqf = new double [NComp];
    Hdqf = new double [NComp];
    Sdqf = new double [NComp];
    CPdqf = new double [NComp];
    Vdqf = new double [NComp];
}


void TModOther::free_internal()
{
	delete[]Gdqf;
	delete[]Hdqf;
	delete[]Sdqf;
	delete[]CPdqf;
	delete[]Vdqf;
}


// calculates pure species properties (pure fugacities, DQF corrections)
long int TModOther::PureSpecies()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


// calculates T,P corrected binary interaction parameters
long int TModOther::PTparam()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


// calculates activity coefficients
long int TModOther::MixMod()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


// calculates excess properties
long int TModOther::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}




//--------------------- End of s_fgl3.cpp ---------------------------

