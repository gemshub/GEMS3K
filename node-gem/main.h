//--------------------------------------------------------------------
// $Id: main.cpp 517 2010-12-16 08:06:21Z gems $
//
// Demo test of usage of the TNode class for implementing a simple
// direct coupling scheme between FMT and GEM in a single-GEM-call
// fashion, assuming that the chemical speciation and all dynamic
// parameter data are kept in the FMT part, which calls GEMIPM
// calculation once per node.
// TNode class implements a  simple C/C++ interface between GEMIPM
// and FMT codes. Works with DATACH and work DATABR structures
//
// Copyright (c) 2006-2012 S.Dmytriyeva, D.Kulik, G.Kosakowski
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
#ifndef MAIN_H
#define MAIN_H

#include <time.h>
#include <math.h>
#include <string>
#include <iomanip> // setprecision
using namespace std;

#include "node.h"

class TMyTransport
{
 public:
    long int nNodes,  // Number of mass transport nodes
            nTimes,   // Number of time steps
            nIC,      // Number of chemical independent components
            nDC,      // Number of chemical dependent components
            nPH,      // Number of chemical phases
            nPS,      // Number of chemical phases-solutions
            nRecipes; // Number of different input node recipes to set boundary conditions

    long int *aNodeHandle,     // Node identification handles
             *aNodeStatusCH,   // Node status codes (changed after GEM calculation)
             *aIterDone;       // Number of GEM IPM iterations performed for each node
                               //   at the last time step
    double tm,      // time, s
           dt;      // time increment, s

    double *aT,     // Array of node temperatures T, Kelvin
           *aP,     // Array of node pressures P, Pa
           *aVs,    // Array of node volume V of reactive subsystem, m3
           *aMs,    // Array of node mass of reactive subsystem, kg
           *aGs,    // Array of node total Gibbs energy of reactive subsystems, J
           *aHs,    // Array of node total enthalpy of reactive subsystems, J (reserved)
           *aIC,    // Array of node effective aqueous ionic strengths, molal
           *apH,    // Array of node pH of aqueous solutions
           *ape,    // Array of node pe of aqueous solutions
           *aEh;    // Array of node Eh of aqueous solution, V

    double **axDC,  // Array of node mole amounts of dependent components (speciation)
           **agam,  // Array of node activity coefficients of dependent components
           **axPH,  // Array of node total mole amounts of all reactive phases
           **aaPH,  // Array of node specific surface areas of phases, m2/kg
           **avPS,  // Array of node total volumes of multicomponent phases, m3
           **amPS,  // Array of node total masses of multicomponent phases,kg
           **abPS,  // Array of node bulk compositions of multicomponent phases, moles
           **axPA,  // Array of node amount of carrier in asymmetric phases, moles
           **aaPh,  // Array of node surface areas of phases, m2
           **adul,  // Array of node upper restrictions to amounts of dependent components
           **adll,  // Array of node lower restrictions to amounts of dependent components
           **abIC,  // Array of node bulk mole amounts of independent components
           **arMB,  // Array of node mole balance residuals for independent components
           **auIC,  // Array of node chemical potentials of independent components (norm.)
           **abSP,  // Array for bulk composition of solid part of equilibrated sub-system
           **aomPH, // Array of stability indices of phases,log10 scale
         **amru,  // Array of upper metastability restrictions for amounts of phases  // needed for MARKS
         **amrl;  // Array of upper metastability restrictions for amounts of phases  // needed for MARKS
                  // MARK: Mineral-Aqueous Reaction Kinetic Simulations

        TMyTransport(   // Constructor (dynamic memory allocation)
           long int p_nNod,    // Number of nodes
           long int p_nTim,    // Number of time steps
           double p_Tim,       // Total maximum time
           double p_dTim,      // Time step duration
           long int p_nIC,     // Number of chemical independent components
           long int p_nDC,     // Number of chemical dependent components
           long int p_nPH,     // Number of chemical phases
           long int p_nPS,     // Number of chemical phases - solutions
           long int p_nRcps    // Number of different input node recipes to set boundary conditions
                );

        ~TMyTransport();  // Destructor of dynamic memory

         double OneTimeStepRun(   // Placeholder function for one transport time step, returns dt
            long int *ICndx,    // Indexes of mobile independent components
            long int nICndx     // Number of mobile independent components
                 );
         double OneTimeStepRun_CN(   // Placeholder function for one transport time step, returns dt
            long int *ICndx,    // Indexes of mobile independent components
            long int nICndx     // Number of mobile independent components
                 );
};

#endif // MAIN_H
