//--------------------------------------------------------------------
// $Id$
//
/// \file node_copy.cpp
/// Copy constructor for TNode class
/// DATACH and DATABR structures allocations
//
// Copyright (c) 2017 S.Dmytriyeva, D.Kulik
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

#include <iomanip>
#include  <iostream>
//#include "io_arrays.h"
#include "node.h"
#include "gdatastream.h"
#include "num_methods.h"

#ifdef IPMGEMPLUGIN

TNode::TNode( const TNode& otherNode )
{
  CSD = 0;
  CNode = 0;

  allocMemory();
  dbr_file_name = otherNode.dbr_file_name;
  ipmlog_file_name = otherNode.ipmlog_file_name;

  // copy data from otherNode
  datach_copy( otherNode.CSD );
  databr_copy( otherNode.CNode );

  multi->copyMULTI( *otherNode.multi );

  // copy intervals for minimizatiom
   pmm->Pai[0] = CSD->Pval[0]/bar_to_Pa;
   pmm->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
   pmm->Pai[2] = getStep( pmm->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
   pmm->Pai[3] = CSD->Ptol/bar_to_Pa;

   pmm->Tai[0] = CSD->TKval[0]-C_to_K;
   pmm->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
   pmm->Tai[2] = getStep( pmm->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
   pmm->Tai[3] = CSD->Ttol;

  pmm->Fdev1[0] = 0.;
  pmm->Fdev1[1] = 1e-6;   // 24/05/2010 must be copy from GEMS3 structure
  pmm->Fdev2[0] = 0.;
  pmm->Fdev2[1] = 1e-6;

  cout << "copy constructor..." << endl;
}

// Copy CSD (DATACH structure) data from other structure.
void TNode::datach_copy( DATACH* otherCSD )
{
    // const data
    copyValues( &CSD->nIC,  &otherCSD->nIC, 14 );
    copyValues( &CSD->Ttol, &otherCSD->Ttol,4 );

    datach_realloc();
    databr_realloc();

    //dynamic data
    copyValues( CSD->nDCinPH, otherCSD->nDCinPH, CSD->nPH );
    //   if( CSD->nICb >0 )
     copyValues( CSD->xic, otherCSD->xic, CSD->nICb );
    copyValues( CSD->xdc, otherCSD->xdc, CSD->nDCb );
    copyValues( CSD->xph, otherCSD->xph, CSD->nPHb );

    copyValues( CSD->A, otherCSD->A, CSD->nIC*CSD->nDC );
    copyValues( CSD->ICmm, otherCSD->ICmm, CSD->nIC );
    copyValues( CSD->DCmm, otherCSD->DCmm, CSD->nDC );

    copyValues( CSD->TKval,  otherCSD->TKval,  CSD->nTp );
    copyValues( CSD->Psat,  otherCSD->Psat,CSD->nTp );
    copyValues( CSD->Pval,  otherCSD->Pval,  CSD->nPp );

    copyValues( CSD->ccIC, otherCSD->ccIC, CSD->nIC );
    copyValues( CSD->ccDC, otherCSD->ccDC, CSD->nDC );
    copyValues( CSD->ccPH, otherCSD->ccPH, CSD->nPH );

    if( CSD->ccPH[0] == PH_AQUEL )
       {
         copyValues( CSD->denW,  otherCSD->denW,  5*gridTP() );
         copyValues( CSD->denWg,  otherCSD->denWg,  5*gridTP() );
         copyValues( CSD->epsW, otherCSD->epsW, 5*gridTP() );
         copyValues( CSD->epsWg, otherCSD->epsWg, 5*gridTP() );
       }
    copyValues( CSD->G0,  otherCSD->G0,  CSD->nDC*gridTP() );
    copyValues( CSD->V0,  otherCSD->V0,  CSD->nDC*gridTP() );
    copyValues( CSD->H0,  otherCSD->H0,  CSD->nDC*gridTP() );
    copyValues( CSD->S0, otherCSD->S0, CSD->nDC*gridTP() );
    copyValues( CSD->Cp0, otherCSD->Cp0, CSD->nDC*gridTP() );
    copyValues( CSD->A0, otherCSD->A0, CSD->nDC*gridTP() );
    copyValues( CSD->U0, otherCSD->U0, CSD->nDC*gridTP() );
    if(  CSD->iGrd  )
         copyValues( CSD->DD, otherCSD->DD, CSD->nDCs*gridTP() );

    copyValues( (char *)CSD->ICNL, (char *)otherCSD->ICNL, MaxICN*CSD->nIC );
    copyValues( (char *)CSD->DCNL, (char *)otherCSD->DCNL, MaxDCN*CSD->nDC );
    copyValues( (char *)CSD->PHNL, (char *)otherCSD->PHNL, MaxPHN*CSD->nPH );
}

// Copy node (work DATABR structure) data from other DBR.
void TNode::databr_copy( DATABR* otherCNode )
{
    // const data
    copyValues( &CNode->NodeHandle, &otherCNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
    if( CNode->NodeStatusFMT != No_nodearray )
         copyValues( &CNode->TK, &otherCNode->TK, 32 );
    else
         copyValues( &CNode->TK, &otherCNode->TK, 15 );
#else
    fstream f_log(ipmLogFile().c_str(), ios::out|ios::app );
    ErrorIf(CNode->NodeStatusFMT != No_nodearray, ipmLogFile().c_str(),
         "Error reading work dataBR structure from binary file (No_nodearray)");
     copyValues( &CNode->TK, &otherCNode->TK, 15 );
#endif
    //dynamic data
    copyValues( CNode->bIC, otherCNode->bIC, CSD->nICb );
    copyValues( CNode->rMB, otherCNode->rMB, CSD->nICb );
    copyValues( CNode->uIC, otherCNode->uIC, CSD->nICb );
    copyValues( CNode->bSP, otherCNode->bSP, CSD->nICb );

    copyValues( CNode->xDC, otherCNode->xDC, CSD->nDCb );
    copyValues( CNode->gam, otherCNode->gam, CSD->nDCb );
    copyValues( CNode->dul, otherCNode->dul, CSD->nDCb );
    copyValues( CNode->dll, otherCNode->dll, CSD->nDCb );

    if( CSD->nAalp >0 )
        copyValues( CNode->aPH, otherCNode->aPH, CSD->nPHb );
    copyValues( CNode->xPH, otherCNode->xPH, CSD->nPHb );
    copyValues( CNode->vPS, otherCNode->vPS, CSD->nPSb );
    copyValues( CNode->mPS, otherCNode->mPS, CSD->nPSb );
    copyValues( CNode->bPS, otherCNode->bPS, CSD->nPSb*CSD->nICb );
    copyValues( CNode->xPA, otherCNode->xPA, CSD->nPSb );
    copyValues( CNode->amru, otherCNode->amru, CSD->nPSb );
    copyValues( CNode->amrl, otherCNode->amrl, CSD->nPSb );
    copyValues( CNode->omPH, otherCNode->omPH, CSD->nPHb );
}

#endif

// allocating DataCH structure
void TNode::datach_realloc()
{
  if( CSD->mLook == 1 &&  (CSD->nPp != CSD->nTp) )
     Error( "No-interpolation mode",
           "Different number of points for temperature and pressure ");

 CSD->nDCinPH = new long int[CSD->nPH];

 if( CSD->nICb >0 )
   CSD->xic = new long int[CSD->nICb];
 else  CSD->xic = 0;
 if( CSD->nDCb >0 )
   CSD->xdc = new long int[CSD->nDCb];
 else  CSD->xdc = 0;
 if( CSD->nPHb >0 )
   CSD->xph = new long int[CSD->nPHb];
 else  CSD->xph = 0;

  CSD->A = new double[CSD->nIC*CSD->nDC];
  CSD->ICmm = new double[CSD->nIC];
  CSD->DCmm = new double[CSD->nDC];
CSD->DCmm[0] = 0.0;   // Added by DK on 03.03.2007

  CSD->TKval = new double[CSD->nTp];
  CSD->Psat = new double[CSD->nTp];
  CSD->Pval = new double[CSD->nPp];

  CSD->denW = new double[ 5*gridTP()];
  CSD->denWg = new double[ 5*gridTP()];
  CSD->epsW = new double[ 5*gridTP()];
  CSD->epsWg = new double[ 5*gridTP()];

  CSD->G0 = new double[CSD->nDC*gridTP()];
  CSD->V0 = new double[CSD->nDC*gridTP()];
  CSD->H0 = new double[CSD->nDC*gridTP()];
  CSD->S0 = new double[CSD->nDC*gridTP()];
  CSD->Cp0 = new double[CSD->nDC*gridTP()];
  CSD->A0 = new double[CSD->nDC*gridTP()];
  CSD->U0 = new double[CSD->nDC*gridTP()];

  if(  CSD->iGrd  )
       CSD->DD = new double[CSD->nDCs*gridTP()];
  else
       CSD->DD = 0;
  CSD->ICNL = new char[CSD->nIC][MaxICN];
  CSD->DCNL = new char[CSD->nDC][MaxDCN];
  CSD->PHNL = new char[CSD->nPH][MaxPHN];

  CSD->ccIC = new char[CSD->nIC];
  CSD->ccDC = new char[CSD->nDC];
  CSD->ccPH = new char[CSD->nPH];
}

// free dynamic memory
void TNode::datach_free()
{
 if( CSD->nDCinPH )
  { delete[] CSD->nDCinPH;
    CSD->nDCinPH = 0;
  }
 if( CSD->xic )
  { delete[] CSD->xic;
    CSD->xic = 0;
  }
 if( CSD->xdc )
  { delete[] CSD->xdc;
    CSD->xdc = 0;
  }
 if( CSD->xph )
  { delete[] CSD->xph;
    CSD->xph = 0;
  }
 if( CSD->A )
  { delete[] CSD->A;
    CSD->A = 0;
  }
 if( CSD->ICmm )
  { delete[] CSD->ICmm;
    CSD->ICmm = 0;
  }
 if( CSD->DCmm )
  { delete[] CSD->DCmm;
    CSD->DCmm = 0;
  }

 if( CSD->TKval )
  { delete[] CSD->TKval;
    CSD->TKval = 0;
  }
 if( CSD->Psat )
  { delete[] CSD->Psat;
    CSD->Psat = 0;
  }
 if( CSD->Pval )
  { delete[] CSD->Pval;
    CSD->Pval = 0;
  }

 if( CSD->denW )
  { delete[] CSD->denW;
    CSD->denW = 0;
  }
 if( CSD->denWg )
  { delete[] CSD->denWg;
    CSD->denWg = 0;
  }
 if( CSD->epsW )
  { delete[] CSD->epsW;
    CSD->epsW = 0;
  }
 if( CSD->epsWg )
  { delete[] CSD->epsWg;
    CSD->epsWg = 0;
  }
 if( CSD->G0 )
  { delete[] CSD->G0;
    CSD->G0 = 0;
  }
 if( CSD->V0 )
  { delete[] CSD->V0;
    CSD->V0 = 0;
  }
 if( CSD->H0 )
  { delete[] CSD->H0;
    CSD->H0 = 0;
  }
 if( CSD->Cp0 )
  { delete[] CSD->Cp0;
    CSD->Cp0 = 0;
  }
  if( CSD->S0 )
  { delete[] CSD->S0;
     CSD->S0 = 0;
  }
  if( CSD->A0 )
  { delete[] CSD->A0;
     CSD->A0 = 0;
  }
  if( CSD->U0 )
  { delete[] CSD->U0;
     CSD->U0 = 0;
  }
  if( CSD->DD )
  { delete[] CSD->DD;
     CSD->DD = 0;
  }

 if( CSD->ICNL )
  { delete[] CSD->ICNL;
    CSD->ICNL = 0;
  }
 if( CSD->DCNL )
  { delete[] CSD->DCNL;
    CSD->DCNL = 0;
  }
 if( CSD->PHNL )
  { delete[] CSD->PHNL;
    CSD->PHNL = 0;
  }

 if( CSD->ccIC )
  { delete[] CSD->ccIC;
    CSD->ccIC = 0;
  }
 if( CSD->ccDC )
  { delete[] CSD->ccDC;
    CSD->ccDC = 0;
  }
 if( CSD->ccPH )
  { delete[] CSD->ccPH;
    CSD->ccPH = 0;
  }
 // delete[] CSD;
}

// Allocates DataBR structure
void TNode::databr_realloc( DATABR * CNode_)
{
  long int j,k;
  CNode_->bIC = new double[CSD->nICb];
  CNode_->rMB = new double[CSD->nICb];
  CNode_->uIC = new double[CSD->nICb];
  CNode_->bSP = new double[CSD->nICb];

  for(  j=0; j<CSD->nICb; j++ )
  {
      CNode_->rMB[j] = 0.;
      CNode_->uIC[j] = 0.;
      CNode_->bSP[j] = 0.;
   }

  CNode_->xDC = new double[CSD->nDCb];
  CNode_->gam = new double[CSD->nDCb];

  for(  j=0; j<CSD->nDCb; j++ )
  {
    CNode_->xDC[j] = 0.;
    CNode_->gam[j] = 1.;
  }

  //  default assignment
 CNode_->dul = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode_->dul[j] = 1.0e6;            // default assignment
 CNode_->dll = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode_->dll[j] = 0.0;              // default assignment

 if( CSD->nAalp >0 )
 {
    CNode_->aPH = new double[CSD->nPHb];
    for(  k=0; k<CSD->nPHb; k++ )
      CNode_->aPH[k] = 0.0;       // default assignment
 }
 else
    CNode_->aPH = 0;

 CNode_->xPH = new double[CSD->nPHb];
 CNode_->omPH = new double[CSD->nPHb];

 for(  k=0; k<CSD->nPHb; k++ )
 {
     CNode_->xPH[k] = 0.0;       // default assignment
     CNode_->omPH[k] = 0.0;
 }

 CNode_->vPS = new double[CSD->nPSb];
 CNode_->mPS = new double[CSD->nPSb];
 CNode_->bPS = new double[CSD->nPSb*CSD->nICb];
 CNode_->xPA = new double[CSD->nPSb];
 CNode_->amru = new double[CSD->nPSb];
 CNode_->amrl = new double[CSD->nPSb];

 for(  k=0; k<CSD->nPSb; k++ )
 {
     CNode_->vPS[k] = 0.0;
     CNode_->mPS[k] = 0.0;
     CNode_->xPA[k] = 0.0;
     for(  j=0; j<CSD->nICb; j++ )
        CNode_->bPS[k*CSD->nICb+j] = 0.0;
     CNode_->amru[k] = 1.0e6;
     CNode_->amrl[k] = 0.0;
 }
}

// free dynamic memory
DATABR * TNode::databr_free( DATABR *CNode_ )
{
  if( CNode_ == 0)
    CNode_ = CNode;

 if( CNode_->bIC )
 { delete[] CNode_->bIC;
   CNode_->bIC = 0;
 }
 if( CNode_->rMB )
 { delete[] CNode_->rMB;
   CNode_->rMB = 0;
 }
 if( CNode_->uIC )
 { delete[] CNode_->uIC;
   CNode_->uIC = 0;
 }

 if( CNode_->xDC )
  { delete[] CNode_->xDC;
    CNode_->xDC = 0;
  }
 if( CNode_->gam )
  { delete[] CNode_->gam;
    CNode_->gam = 0;
  }
 if( CNode_->dul )
   { delete[] CNode_->dul;
     CNode_->dul = 0;
   }
 if( CNode_->dll )
   { delete[] CNode_->dll;
     CNode_->dll = 0;
   }

 if( CNode_->aPH )
 { delete[] CNode_->aPH;
   CNode_->aPH = 0;
 }
 if( CNode_->xPH )
  { delete[] CNode_->xPH;
    CNode_->xPH = 0;
  }
 if( CNode_->omPH )
  { delete[] CNode_->omPH;
    CNode_->omPH = 0;
  }
 if( CNode_->vPS )
  { delete[] CNode_->vPS;
    CNode_->vPS = 0;
  }
 if( CNode_->mPS )
  { delete[] CNode_->mPS;
    CNode_->mPS = 0;
  }
 if( CNode_->bPS )
  { delete[] CNode_->bPS;
    CNode_->bPS = 0;
  }
 if( CNode_->xPA )
  { delete[] CNode_->xPA;
    CNode_->xPA = 0;
  }

 if( CNode_->bSP )
 { delete[] CNode_->bSP;
   CNode_->bSP = 0;
 }
 if( CNode_->amru )
  { delete[] CNode_->amru;
    CNode_->amru = 0;
  }
 if( CNode_->amrl )
  { delete[] CNode_->amrl;
    CNode_->amrl = 0;
  }
 delete CNode_;
 return NULL;
}

// set default values(zeros) for DATABR structure
void TNode::databr_reset( DATABR *CNode, long int level )
{
    //  FMT variables (units or dimensionsless) - to be used for storing them
    //  at the nodearray level = 0.; normally not used in the single-node FMT-GEM coupling
        CNode->Tm = 0.;
        CNode->dt = 0.;
#ifdef NODEARRAYLEVEL
        CNode->Dif = 0.;
        CNode->Vt = 0.;
        CNode->vp = 0.;
        CNode->eps = 0.;
        CNode->Km = 0.;
        CNode->Kf = 0.;
        CNode->S = 0.;
        CNode->Tr = 0.;
        CNode->h = 0.;
        CNode->rho = 0.;
        CNode->al = 0.;
        CNode->at = 0.;
        CNode->av = 0.;
        CNode->hDl = 0.;
        CNode->hDt = 0.;
        CNode->hDv = 0.;
        CNode->nto = 0.; //19
#endif
        if(level <1 )
          return;

   CNode->NodeHandle = 0;
   CNode->NodeTypeHY = normal;
   CNode->NodeTypeMT = normal;
   CNode->NodeStatusFMT = Initial_RUN;
   CNode->NodeStatusCH = NEED_GEM_AIA;
   CNode->IterDone = 0;      //6

// Chemical scalar variables
    CNode->TK = 0.;
    CNode->P = 0.;
    CNode->Vs = 0.;
    CNode->Vi = 0.;
    CNode->Ms = 0.;
    CNode->Mi = 0.;
    CNode->Gs = 0.;
    CNode->Hs = 0.;
    CNode->Hi = 0.;
    CNode->IC = 0.;
    CNode->pH = 0.;
    CNode->pe = 0.;
    CNode->Eh = 0.; //13

    if( level < 2 )
       return;

// Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
    CNode->bIC = 0;
    CNode->rMB = 0;
    CNode->uIC = 0;
    CNode->xDC = 0;
    CNode->gam = 0;
   CNode->dul = 0;
   CNode->dll = 0;
   CNode->aPH = 0;
   CNode->xPH = 0;
   CNode->vPS = 0;
   CNode->mPS = 0;
   CNode->bPS = 0;
   CNode->xPA = 0;
   CNode->bSP = 0;
   CNode->amru = 0;
   CNode->amrl = 0;
   CNode->omPH = 0;
}

// set default values(zeros) for DATACH structure
void TNode::datach_reset()
{
    CSD->nIC = 0;
    CSD->nDC = 0;
    CSD->nPH = 0;
    CSD->nPS = 0;
    CSD->nDCs = 0;
    CSD->nTp = 0;
    CSD->nPp = 0;
    CSD->iGrd = 0;
    CSD->nAalp = 0;
    CSD->nICb = 0;
    CSD->nDCb = 0;
    CSD->nPHb = 0;
    CSD->nPSb = 0;
    CSD->mLook = 0;
// Lists = 0; vectors and matrices
    CSD->nDCinPH = 0;
    CSD->xic = 0;
    CSD->xdc = 0;
    CSD->xph = 0;  //18

    CSD->TKval = 0;
    CSD->Psat = 0;
    CSD->Pval = 0;
    CSD->A = 0;
    CSD->Ttol = 0.;
    CSD->Ptol = 0.;
    CSD->dRes1 = 0.;
    CSD->dRes2 = 0.;
    CSD->ICmm = 0;
    CSD->DCmm = 0;
    CSD->DD = 0;
    CSD->denW = 0;
    CSD->epsW = 0;
    CSD->denWg = 0;
    CSD->epsWg = 0;
    CSD->G0 = 0;
    CSD->V0 = 0;
    CSD->S0 = 0;
    CSD->H0 = 0;
    CSD->Cp0 = 0;
    CSD->A0 = 0;
    CSD->U0 = 0;
    CSD->ICNL = 0;
    CSD->DCNL = 0;
    CSD->PHNL = 0;
    CSD->ccIC = 0;
    CSD->ccDC = 0;
    CSD->ccPH = 0;
}

