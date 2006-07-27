//-------------------------------------------------------------------
// $Id: ms_param.cpp 783 2006-07-27 11:18:26Z gems $
//
// Copyright  (C) 1992-2000 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation  of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifdef __unix__
#include <unistd.h>
#endif

#include <math.h>
#include "m_param.h"
#include "num_methods.h"
#include "gdatastream.h"
#include "node.h"

TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481;

enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

SPP_SETTING pa_ = {
    "GEM-Selektor v2-PSI: Controls and defaults for numeric modules",
    {
        1,  /* PC */  3,     /* PD */   3,   /* PRD */
        1,  /* PSM  */ 144,  /* DP */   15,   /* DW */
        0, /* was -2 DT */  0,     /* PLLG */   1,   /* PE */
        500,   /* IIM */
        1e-30, /* DG */   1e-8,  /* DHB */  1e-12,  /* DS */
        1e-5,  /* DK */  0.01,  /* DF */  1e-9,  /* DFM */
        1e-6,  /* DFYw */  1e-6,  /* DFYaq */    1e-6,  /* DFYid */
        1e-6,  /* DFYr,*/  1e-6,  /* DFYh,*/   1e-6,  /* DFYc,*/
        1e-12, /* DFYs, */  1e-17,  /* DB */   0.7,   /* AG */
        0.07,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
        0.01, /* GAS */  12.05,  /* DNS */  1e-8,  /* XwMin, */
        1e-8,  /* ScMin, */  1e-20, /* DcMin, */   1e-10, /* PhMin, */
        3e-5,  /* ICmin */   1e-7,  /* EPS */   1e-3,  /* IEPS */
        1e-4,  /* DKIN  */ 0,  /* tprn */
    },
}; /* SPP_SETTING */


void BASE_PARAM::write(ostream& oss)
{
    oss.write( (char*)&PC, 10*sizeof(short) );
    oss.write( (char*)&DG, 28*sizeof(double) );
    oss.write( (char*)&tprn, sizeof(char*) );
}

void SPP_SETTING::write(ostream& oss)
{
    oss.write( ver, TDBVERSION );
    p.write( oss );
}

TProfil::TProfil( TMulti* amulti )
{
    pa= pa_;
    multi = amulti;
    pmp = multi->GetPM();
}

// GEM IPM calculation of equilibrium state in MULTI
void TProfil::calcMulti()
{
    multi->MultiCalcInit( 0 );
    if( multi->AutoInitialApprox() == false )
    {
       //       TNode::na->GEM_printf( "calc_multi_test.ipm", 0, 0 );
        multi->MultiCalcIterations();
    }
}

void TProfil::outMulti( GemDataStream& ff, gstring& path  )
{
    ff.writeArray( &pa.p.PC, 10 );
    ff.writeArray( &pa.p.DG, 28 );
    multi->to_file( ff, path );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( GemDataStream& ff )
{
      ff.readArray( &pa.p.PC, 10 );
      ff.readArray( &pa.p.DG, 28 );
      multi->from_file( ff );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( const char* path )
{
      multi->from_text_file_gemipm( path);
}

// Load Thermodynamic Data from MTPARM to MULTI using LagranInterp
void TMulti::CompG0Load()
{
  int j, jj, k, jb, je=0;
  double Gg, Vv;
  float TC, P;

  DATACH  *dCH = TNode::na->pCSD();
//  DATABR  *dBR = TNodeArray::na->pCNode();

  TC = TNode::na->cT()-C_to_K;
  P = TNode::na->cP();

// if( dCH->nTp <=1 && dCH->nPp <=1 )
if( dCH->nTp <1 && dCH->nPp <1 )
   return;

  for( jj=0; jj<dCH->nTp; jj++)
    if( fabs( TC - dCH->Tval[jj] ) < dCH->Ttol )
    {
        TC = dCH->Tval[jj];
        break;
    }
  for( jj=0; jj<dCH->nPp; jj++)
   if( fabs( P - dCH->Pval[jj] ) < dCH->Ptol )
   {
        P = dCH->Pval[jj];
        break;
   }

//Test outpur ***********************************
//  fstream f_log("CompG0Load.txt", ios::out|ios::app );

//  f_log << "TC = " <<  TC << "  P =  " << P << endl;
//Test outpur ***********************************

 if( fabs( pmp->TC - TC ) < 1.e-10 &&
            fabs( pmp->P - P ) < 1.e-10 )
   return;    //T, P not changed


 pmp->T = pmp->Tc = TC + C_to_K;
 pmp->TC = pmp->TCc = TC;
 pmp->P = pmp->Pc = P;
 if( dCH->ccPH[0] == PH_AQUEL )
 {
   pmp->denW = LagranInterp( dCH->Pval, dCH->Tval, dCH->roW,
                          P, TC, dCH->nTp, dCH->nPp,1 );
   //       pmp->denWg = tpp->RoV;
   pmp->epsW = LagranInterp( dCH->Pval, dCH->Tval, dCH->epsW,
                          P, TC, dCH->nTp, dCH->nPp,1 );
   //       pmp->epsWg = tpp->EpsV;
 }
 else
 {
   pmp->denW = 1.;
   pmp->epsW = 78.;
 }

//Test outpur ***********************************
//  f_log << "roW = " <<  pmp->denW << "  epsW =  " << pmp->epsW << endl;
//Test outpur ***********************************

 pmp->RT = R_CONSTANT * pmp->Tc;
 pmp->FRT = F_CONSTANT/pmp->RT;
 pmp->lnP = 0.;
 if( P != 1. ) /* ??????? */
   pmp->lnP = log( P );

 for( k=0; k<pmp->FI; k++ )
 {
   jb = je;
   je += pmp->L1[k];
   /*load t/d data from DC */
    for( j=jb; j<je; j++ )
    {
      jj =  j * dCH->nPp * dCH->nTp;
      Gg = LagranInterp( dCH->Pval, dCH->Tval, dCH->G0+jj,
                          P, TC, dCH->nTp, dCH->nPp,1 );
      pmp->G0[j] = Cj_init_calc( Gg, j, k );
//Test outpur ***********************************
//  f_log << j  << " Gg = " <<  Gg  << "  GO =  " << pmp->G0[j] << endl;
//Test outpur ***********************************
    }
 }
 for( j=0; j<pmp->L; j++ )
 {
    jj =  j * dCH->nPp * dCH->nTp;
    Vv = LagranInterp( dCH->Pval, dCH->Tval, dCH->V0+jj,
                          P, TC, dCH->nTp, dCH->nPp, 1 );
    switch( pmp->PV )
    { /* make mol volumes of components */
       case VOL_CONSTR:
                    pmp->A[j*pmp->N] = Vv;
       case VOL_CALC:
       case VOL_UNDEF:
                    pmp->Vol[j] = Vv  * 10.;
                    break;
    }
//Test outpur ***********************************
//  f_log << j  << " Vv = " <<  Vv  << "  VO =  " << pmp->Vol[j] << endl;
//Test outpur ***********************************
 }
}

// GEM IPM calculation of equilibrium state in MULTI
void TMulti::MultiCalcInit( const char *key )
{
  short j,k;
  SPP_SETTING *pa = &TProfil::pm->pa;

    pmp->Ec = pmp->K2 = pmp->MK = 0;
    pmp->PZ = pa->p.DW; // IPM-2 default
    pmp->W1 = 0;
    pmp->is = 0;
    pmp->js = 0;
    pmp->next  = 0;
    pmp->ln5551 = 4.0165339;
    pmp->lowPosNum = pa->p.DcMin;
    pmp->logXw = -16.;
    pmp->logYFk = -9.;
    pmp->FitVar[4] = pa->p.AG;
    pmp->pRR1 = 0;      //IPM smoothing factor and level
    pmp->DX = pa->p.DK;

    pmp->T0 = 273.15;    // not used anywhere
    pmp->FX = 7777777.;
    pmp->YMET = 0;
    pmp->PCI = 0.0;

//    if( pmp->pESU  && pmp->pNP )     // problematic statement !!!!!!!!!
    {
//      unpackData(); // loading data from EqstatUnpack( key );
        pmp->IC = 0.;
        for( j=0; j< pmp->L; j++ )
            pmp->X[j] = pmp->Y[j];
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
    }
//    else
//        for( j=0; j<pmp->L; j++ )
//            pmp->Y[j] = 0.0;

    CompG0Load();
    for( j=0; j< pmp->L; j++ )
        pmp->G[j] = pmp->G0[j];
    // test phases - solutions and load models
    if( pmp->FIs )
    {
        for( j=0; j< pmp->Ls; j++ )
        {
            pmp->lnGmo[j] = pmp->lnGam[j];
            pmp->Gamma[j] = 1.0;
        }
        pmp->PD = pa->p.PD;
//        SolModLoad();
        GammaCalc( LINK_TP_MODE);
    }
    // recalc restrictions for DC quantities
    if( pmp->pULR && pmp->PLIM )
         Set_DC_limits(  DC_LIM_INIT );

   // dynamic arrays - begin load
    for( k=0; k<pmp->FI; k++ )
    {
        pmp->XFs[k] = pmp->YF[k];
        pmp->Falps[k] = pmp->Falp[k];
    }
}

//-------------------------------------------------------------------------
// internal functions

// read string as: "<characters>",
istream& f_getline(istream& is, gstring& str, char delim)
{
    char ch;
    is.get(ch);
    str="";

    while( is.good() && ( ch==' ' || ch=='\n' || ch== '\t') )
        is.get(ch);
    if(ch == '\"')
        is.get(ch);
    while( is.good() &&  ch!=delim && ch!= '\"' )
    {
        str += ch;
        is.get(ch);
    }
    while( is.good() &&  ch!=delim )
            is.get(ch);

   return is;
}

gstring u_makepath(const gstring& dir,
           const gstring& name, const gstring& ext)
{
    gstring Path(dir);
    if( dir != "")
      Path += "/";
    Path += name;
    Path += ".";
    Path += ext;

    return Path;
}

void u_splitpath(const gstring& Path, gstring& dir,
            gstring& name, gstring& ext)
{
    size_t pos = Path.rfind("/");
    if( pos != npos )
        dir = Path.substr(0, pos), pos++;
    else
        dir = "",    pos = 0;

    size_t pose = Path.rfind(".");
    if( pose != npos )
    {
        ext = Path.substr( pose+1, npos );
        name = Path.substr(pos, pose-pos);
    }
    else
    {
        ext = "";
        name = Path.substr(pos, npos);
    }
}

// ------------------ End of ms_param.cpp -----------------------




