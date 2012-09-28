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
// Copyright (C) 2006,2012 S.Dmytriyeva, D.Kulik, G.Kosakowski
//
// This file is part of the GEMIPM2K code for thermodynamic modelling
// by Gibbs energy minimization
//
// See also http://gems.web.psi.ch/GEMS3K/
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------
#ifndef MAIN_H
#define MAIN_H

#include <time.h>
#include <math.h>
#include <string.h>

#include "node.h"

class TMyTransport
{
 public:
    long int nNodes,  // Number of mass transport nodes
            nTimes,   // Number of time steps
            nIC,      // Number of chemical independent components
            nDC,      // Number of chemical dependent components
            nPH,      // Number of chemical phases
            nPS;      // Number of chemical phases-solutions

    long int *aNodeHandle,     // Node identification handles
             *aNodeStatusCH,   // Node status codes (changed after GEM calculation)
             *aIterDone;       // Number of GEM IPM iterations performed for each node
                               //   at the last time step

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
           **abSP;  // Array for bulk composition of solid part of equilibrated sub-system

        TMyTransport(   // Constructor (dynamic memory allocation)
           long int p_nNod,    // Number of nodes
           long int p_nTim,    // Number of time steps
           long int p_nIC,     // Number of chemical independent components
           long int p_nDC,     // Number of chemical dependent components
           long int p_nPH,     // Number of chemical phases
           long int p_nPS      // Number of chemical phases - solutions
                );

        ~TMyTransport();  // Destructor of dynamic memory

         void OneTimeStepRun(   // Placeholder function for one transport time step
            double *stoich,     // Stoichiometry coefficients
            long int *ICndx,    // Indexes of mobile independent components
            long int nICndx     // Number of mobile independent components
                 );
};

#endif // MAIN_H
