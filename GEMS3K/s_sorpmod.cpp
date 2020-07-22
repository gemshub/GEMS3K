//-------------------------------------------------------------------
// $Id: s_sorption.cpp 725 2012-10-02 15:43:37Z kulik $
//
/// \file s_sorpmod.cpp
/// Implementation of TSolMod derived classes
/// for activity models of mixing in site-balance-based surface complexation phases
///  (TSCM_NEM, TSCM_CCM, TSCM_DLM, TSCM_BSM, TSCM_TLM, TSCM_CDM, TSCM_ETLM, TSCM_BET)
//
// Copyright (C) 2017  D.Kulik
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
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "verror.h"
#include "s_solmod.h"

//=============================================================================================
// Implementation of ideal mixing model for multicomponent solid solutions
// References: Price 1989
// also used with scripted models to provide ideal mixing term in the multi-site case
// (c) DK/TW November 2010
//=============================================================================================


// Generic constructor for the TIdeal class
TSCM_NEM::TSCM_NEM( SolutionData *sd ):
                TSolMod( sd )
{
}

TSCM_NEM::~TSCM_NEM()
{
}

long int TSCM_NEM::PTparam()
{
   return 0;
}


/// Calculates ideal configurational terms in case of multi-site mixing
/// to preserve values computed in Phase scripts.
/// Only increments lnGamma[j] - may need to be cleaned before running MixMod
long int TSCM_NEM::MixMod()
{
   long int retCode, j;

   retCode = EDLmod();

   if(!retCode)
   {
      for(j=0; j<NComp; j++)
          lnGamma[j] += CTerm[j];
   }
   return 0;
}


/// calculates adsorption monolayer phase excess properties
long int TSCM_NEM::ExcessProp( double *Zex )
{

        // assignments (excess properties)
        Zex[0] = 0.;
        Zex[1] = 0.;
        Zex[2] = 0.;
        Zex[3] = 0.;
        Zex[4] = 0.;
        Zex[5] = 0.;
        Zex[6] = 0.;

        return 0;
}


/// calculates ideal mixing properties
long int TSCM_NEM::IdealProp( double *Zid )
{
     Hid = 0.0;
     CPid = 0.0;
     Vid = 0.0;
     Sid = SCM_conf_entropy();
     Gid = Hid - Sid*Tk;
     Aid = Gid - Vid*Pbar;
     Uid = Hid - Vid*Pbar;

     // assignments (ideal mixing properties)
     Zid[0] = Gid;
     Zid[1] = Hid;
     Zid[2] = Sid;
     Zid[3] = CPid;
     Zid[4] = Vid;
     Zid[5] = Aid;
     Zid[6] = Uid;

     return 0;
}

long int TSCM_NEM::EDLmod()
{
    long int j, retCode=0;

    if(MixCode != MR_UNDEF_ )
       return 1; // This should be 'N' code for NEM

    for(j=0; j<NComp; j++)
        CTerm[j] = 0.0;     // For NEM, Coulombic terms equal 0

    return retCode;
}

double TSCM_NEM::SCM_conf_entropy()
{
    long int j;
    double si = 0.0;

    for (j=0; j<NComp; j++)
    {
        if ( x[j] > 1.0e-32 )
            si += x[j]*log(x[j]);
    }
    auto Sid1 = (-1.)*R_CONST*si;
    return Sid1;
}

//--------------------- End of s_sorption.cpp ---------------------------

