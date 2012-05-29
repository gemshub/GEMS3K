//-------------------------------------------------------------------
// $Id$
//
// Copyright  (C) 1992,2012 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation  of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.3.x.x program
// environment for thermodynamic modeling in geochemistry and
// of the GEMS3K standalone code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
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

//TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481,
                             H2O_mol_to_kg = 55.50837344,
                               Min_phys_amount = 1.66e-24;

enum volume_code {  // Codes of volume parameter ???
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

SPP_SETTING pa_ = {
    " GEMS-GUI v.3.1 r.2147 (rc) " " GEMS3K v.3.1 r.681 (rc) ",
    {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     13000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
    },
}; // SPP_SETTING



void BASE_PARAM::write(fstream& oss)
{
  short arr[10];

  arr[0] = PC;
  arr[1] = PD;
  arr[2] = PRD;
  arr[3] = PSM;
  arr[4] = DP;
  arr[5] = DW;
  arr[6] = DT;
  arr[7] = PLLG;
  arr[8] = PE;
  arr[9] = IIM;

  oss.write( (char*)arr, 10*sizeof(short) );
  oss.write( (char*)&DG, 28*sizeof(double) );
  oss.write( (char*)&tprn, sizeof(char*) );
}

void SPP_SETTING::write(fstream& oss)
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

/// GEM IPM calculation of equilibrium state in MULTI
double TProfil::ComputeEquilibriumState( long int& RefinLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
 return multi->CalculateEquilibriumState( 0, NumIterFIA_, NumIterIPM_ );
}

/// Writing structure MULTI (GEM IPM work structure) to binary file
void TProfil::outMulti( GemDataStream& ff, gstring& /*path*/  )
{
	 short arr[10];

	  arr[0] = pa.p.PC;
	  arr[1] = pa.p.PD;
	  arr[2] = pa.p.PRD;
	  arr[3] = pa.p.PSM;
	  arr[4] = pa.p.DP;
	  arr[5] = pa.p.DW;
	  arr[6] = pa.p.DT;
	  arr[7] = pa.p.PLLG;
	  arr[8] = pa.p.PE;
	  arr[9] = pa.p.IIM;

	ff.writeArray( arr, 10 );
    ff.writeArray( &pa.p.DG, 28 );
    multi->to_file( ff );
}

/// Writing structure MULTI ( free format text file  )
void TProfil::outMultiTxt( const char *path, bool append  )
{
    multi->to_text_file( path, append );
}

/// Reading structure MULTI (GEM IPM work structure) from binary file
void TProfil::readMulti( GemDataStream& ff )
{
    //DATACH  *dCH = TNode::na->pCSD();
    DATACH  *dCH =  multi->node->pCSD();
    short arr[10];

	 ff.readArray( arr, 10 );
	  pa.p.PC = arr[0];
	  pa.p.PD = arr[1];
	  pa.p.PRD = arr[2];
	  pa.p.PSM = arr[3];
	  pa.p.DP = arr[4];
	  pa.p.DW = arr[5];
	  pa.p.DT = arr[6];
	  pa.p.PLLG = arr[7];
	  pa.p.PE = arr[8];
	  pa.p.IIM = arr[9];

      ff.readArray( &pa.p.DG, 28 );
      multi->from_file( ff );

      // copy intervals for minimizatiom
      if(  dCH->nPp > 1  )
      {
         pmp->Pai[0] = dCH->Pval[0];
         pmp->Pai[1] = dCH->Pval[dCH->nPp-1];
         pmp->Pai[2] = (pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
      }
      pmp->Pai[3] = dCH->Ptol;
      if(  dCH->nTp > 1  )
      {
         pmp->Tai[0] = dCH->TKval[0];
         pmp->Tai[1] = dCH->TKval[dCH->nTp-1];
         pmp->Tai[2] = (pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
      }
      pmp->Tai[3] = dCH->Ttol;

  }

/// Reading structure MULTI (GEM IPM work structure) from text file
void TProfil::readMulti( const char* path, DATACH  *dCH )
{
      multi->from_text_file_gemipm( path, dCH);
}


// Load Thermodynamic data from Database
// moved to multi from Project

/// Test and load thermodynamic data
void TMulti::CheckMtparam()
{
  double TK, P, PPa;

  //DATACH  *dCH = TNode::na->pCSD();
  //TK = TNode::na->cTK();
  //PPa = TNode::na->cP();

  DATACH  *dCH = node->pCSD();
  TK = node->cTK();
  PPa = node->cP();

  P = PPa/bar_to_Pa;

 //pmp->pTPD = 2;

 if( !load || fabs( pm.Tc - TK ) > dCH->Ttol
           || fabs( pm.P - P )  > dCH->Ptol  )
     pm.pTPD = 0;      //T, P is changed - problematic for UnSpace!

  load = true;
}

//-------------------------------------------------------------------------
// internal functions

void strip(string& str)
{
  string::size_type pos1 = str.find_first_not_of(' ');
  string::size_type pos2 = str.find_last_not_of(' ');
  str = str.substr(pos1 == string::npos ? 0 : pos1,
    pos2 == string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

// read string as: "<characters>"
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

const long int bGRAN = 20;

// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path,
		long int& nElem, char delim ))[fileNameLength]
{
  long int ii, bSize = bGRAN;
  char  (*filesList)[fileNameLength];
  char  (*filesListNew)[fileNameLength];
  filesList = new char[bSize][fileNameLength];
  gstring name;

// Get path
   gstring path_;
   gstring flst_name = f_name;
   unsigned long int pos = flst_name.rfind("/");
   path_ = "";
   if( pos < npos )
      path_ = flst_name.substr(0, pos+1);
   strncpy( Path, path_.c_str(), 256-fileNameLength);
   Path[255] = '\0';

//  open file stream for the file names list file
   fstream f_lst( f_name/*flst_name.c_str()*/, ios::in );
   ErrorIf( !f_lst.good(), f_name, "Fileopen error");

// Reading list of names from file
  nElem = 0;
  while( !f_lst.eof() )
  {
	f_getline( f_lst, name, delim);
    if( nElem >= bSize )
    {    bSize = bSize+bGRAN;
         filesListNew = new char[bSize][fileNameLength];
         for( ii=0; ii<nElem-1; ii++ )
		   strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	     delete[] filesList;
		 filesList =  filesListNew;
	}
    strncpy( filesList[nElem], name.c_str(), fileNameLength);
	filesList[nElem][fileNameLength-1] = '\0';
    nElem++;
  }

  // Realloc memory for reading size
  if( nElem != bSize )
  {
    filesListNew = new char[nElem][fileNameLength];
    for(  ii=0; ii<nElem; ii++ )
	  strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	delete[] filesList;
	filesList =  filesListNew;
  }

  return filesList;
}

//====================================================================================================
// old functions


//#include "io_arrays.h"
/*
void TProfil::test_G0_V0_H0_Cp0_DD_arrays( long int nT, long int nP )
{
  long int kk, jj, ii, l1, l2, xTP, lev = 11;
  double cT, cP, dT, dP;
  double *G0, *V0, *H0,* *Cp0, *S0, *A0, *U0,* *denW, *epsW, *denWg, *epsWg;
  DATACH  *CSD = TNode::na->pCSD();
  fstream ff("lagrange_T_11_3.out.txt", ios::out );
  ErrorIf( !ff.good() , "lagrange.out", "Fileopen error");

  G0 =  new double[CSD->nDC*nT*nP];
//   V0 =  new double[CSD->nDC*nT*nP];
//   H0 =  new double[TProfil::pm->mup->L*nT*nP];
  Cp0 = new double[CSD->nDC*nT*nP];
  S0 = new double[CSD->nDC*nT*nP];
//   A0 = new double[TProfil::pm->mup->L*nT*nP];
//   U0 = new double[TProfil::pm->mup->L*nT*nP];
  denW =  new double[5*nT*nP];
  epsW =  new double[5*nT*nP];

  cT = CSD->TCval[0];
  cP = CSD->Pval[0];
  if(nT > 1)
    dT = (CSD->TCval[CSD->nTp-1]-CSD->TCval[0])/(nT-1);
  else
       dT = 0;
  if(nP > 1 )
    dP = (CSD->Pval[CSD->nPp-1]-CSD->Pval[0])/(nP-1);
  else
         dP = 0;
  for(  ii=0; ii<nT; ii++)
  {
//     cT += dT;
    for(  jj=0; jj<nP; jj++)
    {
//       cP += dP;
      xTP = TNode::na->check_grid_TP( cT, cP ); // changed to Kelvin
      for( kk=0; kk<5; kk++)
      {
        l1 = ( kk * nP + jj) * nT + ii;
        l2  =  kk * CSD->nPp * CSD->nTp;

        if( xTP >= 0 )
         {
          denW[l1] = CSD->denW[l2+xTP];
          epsW[l1] = CSD->epsW[l2+xTP];
         }
         else
         {
           denW[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->denW+l2,
                           cP, cT, CSD->nTp, CSD->nPp,lev );
           epsW[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->epsW+l2,
                           cP, cT, CSD->nTp, CSD->nPp,lev );
         }
      }
      // copy to arrays
      for( kk=0; kk<CSD->nDC; kk++)
      {
          l1 = ( kk * nP + jj) * nT + ii;
          l2  =  kk * CSD->nPp * CSD->nTp;
         if( xTP >= 0 )
          {
                G0[l1] = CSD->G0[ l2+xTP];
                S0[l1] = CSD->S0[ l2+xTP];
                Cp0[l1] = CSD->Cp0[ l2+xTP];
           }
           else
           {
             G0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->G0+l2,
                             cP, cT, CSD->nTp, CSD->nPp,lev );
             S0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->S0+l2,
                             cP, cT, CSD->nTp, CSD->nPp,lev );
             Cp0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->Cp0+l2,
                             cP, cT, CSD->nTp, CSD->nPp,lev );
           }
       }
      cP += dP;
    }
    cT += dT;
  }

//  fstream ff("lagrange_T_5_3.out.txt", ios::out );
//  ErrorIf( !ff.good() , "lagrange.out", "Fileopen error");
 TPrintArrays  prar(0, 0, ff);
 prar.writeArray(  "denW", denW, 5**nP*nT, nP*nT );
 prar.writeArray(  "epsW", epsW, 5**nP*nT, nP*nT );
 prar.writeArray(  "S0_Ca+2", S0, nP*nT,  nP*nT ); //Ca+2
 prar.writeArray(  "S0_H2O", S0+13*nP*nT, nP*nT, nP*nT); //H2O@
 prar.writeArray(  "S0_Brc", S0+18*nP*nT, nP*nT, nP*nT ); //Brc
 prar.writeArray(  "G0_Ca+2", G0, nP*nT, nP*nT ); //Ca+2
 prar.writeArray(  "G0_H2O", G0+13*nP*nT, nP*nT, nP*nT ); //H2O@
 prar.writeArray(  "G0_Brc", G0+18*nP*nT, nP*nT, nP*nT ); //Brc
 prar.writeArray(  "Cp0_Ca+2", Cp0, nP*nT, nP*nT ); //Ca+2
 prar.writeArray(  "Cp0_H2O", Cp0+13*nP*nT, nP*nT, nP*nT ); //H2O@
 prar.writeArray(  "Cp0_Brc", Cp0+18*nP*nT, nP*nT, nP*nT ); //Brc
//  prar.writeArray(  "V0", V0, CSD->nDC*nP*nT, nP*nT );
//  prar.writeArray(  "G0", G0, CSD->nDC*nP*nT, nP*nT );

  delete[] G0;
//   delete[] V0;
//   delete[] H0;
  delete[] Cp0;
  delete[] S0;
//   delete[] A0;
//   delete[] U0;
  delete[] denW;
  delete[] epsW;
}
*/

//=============================================================================================
/* GEM IPM calculation of equilibrium state in MULTI
void TMulti::MultiCalcInit( const char* *key* )
{
  long int  j,k;

  SPP_SETTING *pa = &TProfil::pm->pa;
   // SD 09/03/09 strncpy( pmp->stkey, key, EQ_RKLEN );  // needed for the ipmlog error diagnostics
    pmp->Ec = pmp->K2 = pmp->MK = 0;
    pmp->PZ = pa->p.DW; // IPM-2 default
    pmp->W1 = 0;
    pmp->is = 0;
    pmp->js = 0;
    pmp->next  = 0;
    pmp->ln5551 = log( H2O_mol_to_kg ); //  4.016533882  ln(55.50837344)  4.0165339;
//    pmp->lowPosNum = pa->p.DcMin;  obsolete
    pmp->lowPosNum = Min_phys_amount;
    pmp->logXw = -16.;
    pmp->logYFk = -9.;
    pmp->FitVar[0] = 0.0640000030398369; // pa->aqPar[0]; setting T,P dependent b_gamma parameters
    pmp->DXM = pa->p.DK;
    RescaleToSize( true );  // Added to set default cutoffs/inserts 30.08.2009 DK
    pmp->T0 = 273.15;    // not used anywhere
    pmp->FX = 7777777.;
    pmp->YMET = 0;
//    pmp->PCI = 0.0;
    pmp->PCI = 1.0;

    // calculating mass of the system
    pmp->MBX = 0.0;
    for(long int i=0; i<pmp->N; i++ )
    {
        pmp->MBX += pmp->B[i] * pmp->Awt[i];
    }
    pmp->MBX /= 1000.;

    // optimization 08/02/2007 - allocation of A matrix index lists and IPM work arrays
    Alloc_internal(); // performance optimization 08/02/2007

    if(  pmp->pNP )     // Checking if this is SIA or AIA mode
    {   //  Smart Initial Approximation mode
        for( j=0; j< pmp->L; j++ )
          pmp->X[j] = pmp->Y[j];
 //       pmp->IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
        ConCalc( pmp->X, pmp->XF, pmp->XFA);
 //       pmp->IT = pmp->ITaia;
       if( pmp->FIs )
       {
           // Load activity coeffs for phases-solutions
          long int jb, je=0;
             for( k=0; k<pmp->FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += pmp->L1[k];
            // Load activity coeffs for phases-solutions
            for( j=jb; j< je; j++ )
            {
               pmp->lnGmo[j] = pmp->lnGam[j];
               if( fabs( pmp->lnGam[j] ) <= 84. )
      //                pmp->Gamma[j] = exp( pmp->lnGam[j] );
                      pmp->Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
               else pmp->Gamma[j] = 1.0;
             } // j
          }  // k
       }
    }
    else // Simplex initial approximation to be done (AIA mode)
    {
        for( j=0; j<pmp->L; j++ )
        {                           // cleaning work vectors
                pmp->X[j] = pmp->Y[j] = pmp->lnGam[j] = pmp->lnGmo[j] = 0.0;
                pmp->Gamma[j] = 1.0;
                pmp->MU[j] = 0.;// SD 06/12/2009
                pmp->XU[j] = 0.;// SD 06/12/2009
                pmp->EMU[j] = 0.;// SD 06/12/2009
                pmp->NMU[j] = 0.;// SD 06/12/2009
                pmp->W[j] = 0.;// SD 06/12/2009
                pmp->F[j] = 0.;// SD 06/12/2009
                pmp->F0[j] = 0.;// SD 06/12/2009
        }
    }
    CompG0Load(); // Loading thermodynamic data into MULTI structure


     pmp->PD = abs(pa->p.PD);
         //           SolModLoad() and Phase scripts cannot be used here!

    // recalculating kinetic restrictions for DC amounts
     if( pmp->pULR && pmp->PLIM )
          Set_DC_limits(  DC_LIM_INIT );

     bool AllPhasesPure = true;
     // checking if all phases are pure
     for( k=0; k < pmp->FI; k++ )
         if( pmp->L1[k] > 1 )
             AllPhasesPure = false;
     if( !AllPhasesPure )     // bugfix DK 09.03.2010   was if(!pmp->FIs)
    {
      long int jb, je=0;
     // Calc Eh, pe, pH,and other stuff
    if( pmp->E && pmp->LO && pmp->pNP )
    {
        ConCalc( pmp->X, pmp->XF, pmp->XFA);
        IS_EtaCalc();
        if( pmp->Lads )  // Calling this only when sorption models are present

        {
           for( k=0; k<pmp->FIs; k++ )
           { // loop on solution phases
              jb = je;
              je += pmp->L1[k];
              if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
              {
                     if( pmp->PHC[0] == PH_AQUEL && pmp->XF[k] > pmp->DSM
                       && (pmp->XFA[0] > pmp->ScMinM && pmp->XF[0] > pmp->XwMinM ))  // fixed 30.08.2009 DK
                       GouyChapman( jb, je, k );  // getting PSIs - electrical potentials on surface planes
              }
           }
        }
    }
    GammaCalc( LINK_TP_MODE);   // Computing DQF, FugPure and G wherever necessary
                                      // Activity coeffs are restored from lnGmo
 }
   else {  // no multi-component phases
        pmp->PD = 0;  pmp->pNP = 0;
  //      pmp->pIPN = 1;
    }

}
*/


// ------------------ End of ms_param.cpp -----------------------




