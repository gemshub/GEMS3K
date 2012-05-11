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
// Copyright (C) 2006,2010 S.Dmytriyeva, D.Kulik, G.Kosakowski
//
// This file is part of the GEMIPM2K code for thermodynamic modelling
// by Gibbs energy minimization
//
// See also http://gems.web.psi.ch/
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

    long int *aNodeHandle,
             *aNodeStatusCH,
             *aIterDone;

    double *aT,
           *aP,
           *aVs,
           *aMs,
           *aGs,
           *aHs,
           *aIC,
           *apH,
           *ape,
           *aEh;

    double **axDC,
           **agam,
           **axPH,
           **aaPH,
           **avPS,
           **amPS,
           **abPS,
           **axPA,
           **aaPh,
           **adul,
           **adll,
           **abIC,
           **arMB,
           **auIC;

        TMyTransport( long int p_nNod, long int p_nTim, long int p_nIC, long int p_nDC,
                      long int p_nPH, long int p_nPS );

        ~TMyTransport();

        void OneTimeStepRun( double *stoich, long int *ICndx, long int nICndx );
};

#endif // MAIN_H
