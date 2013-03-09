//-------------------------------------------------------------------
// $Id$
//
// Stub implementation of TKinMet class for further development
//
// Copyright (C) 2012  D.Kulik, B.Thien, S.Dmitrieva
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
#include "verror.h"
//#include "s_fgl.h"
#include "s_kinmet.h"


// Generic constructor
TKinMet::TKinMet( const KinMetData *kmd )
{
  ;
}

// Destructor
TKinMet::~TKinMet()
{
  ;
}

// sets the specific surface area of the phase and 'parallel reactions' area fractions
long int TKinMet::UpdateFSA( const double *fSAf_p, const double As )
{

}

// returns modified specific surface area of the phase and 'parallel reactions' area fractions
double TKinMet::ModifiedFSA ( double *fSAf_p )
{

}

// sets new system TP state
long int TKinMet::UpdatePT ( const double T_k, const double P_bar )
{

}

// sets new time and time step
bool TKinMet::UpdateTime( const double Tau, const double dTau )
{

}

// Checks dimensions in order to re-allocate class instance, if necessary
bool TKinMet::testSizes( const KinMetData *kmd )
{

}













//--------------------- End of s_kinmet.cpp ---------------------------

