//-------------------------------------------------------------------
// $Id: ms_param.cpp 1186 2009-01-23 14:19:06Z gems $
//
// Copyright  (C) 1992,2007 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation  of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry and
// of the GEMIPM2K standalone code
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

enum volume_code {  // Codes of volume parameter ???
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

SPP_SETTING pa_ = {
  "GEMS-PSI v2.3.0: Controls and defaults for numeric modules",
  {
        1,  /* PC */  3,     /* PD */   3,   /* PRD */
        1,  /* PSM  */ 150,  /* DP */   15,   /* DW */
        0, /* DT */     0,   /* PLLG */   1,  /* PE */
        500,   /* IIM */
        1e-6, /* DG */   1e-8,  /* DHB */  1e-12,  /* DS */
        1e-4,  /* DK */  0.01,  /* DF */  0.1,  /* DFM */
        1e-6,  /* DFYw */  1e-6,  /* DFYaq */    1e-6,  /* DFYid */
        1e-6,  /* DFYr,*/  1e-6,  /* DFYh,*/   1e-6,  /* DFYc,*/
        1e-7, /* DFYs, */  1e-17,  /* DB */   0.7,   /* AG */
        0.07,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
        0.001, /* GAS */   12.05,  /* DNS */   1e-5,  /* XwMin, */
        1e-7,  /* ScMin, */  1e-19, /* DcMin, */   1e-10, /* PhMin, */
        1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
        1e-10,  /* DKIN  */ 0,  /* tprn */
  },
}; // SPP_SETTING


void BASE_PARAM::write(ostream& oss)
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
// Modified on 10.09.2007 to return elapsed GEMIPM2 runtime in seconds
// Modified on 15.11.2007 to return more detailed info on FIA and IPM iterations
// and precision refinement loops
// Modified on 22.01.2008 to implement "smart PIA" mode 
//
double TProfil::calcMulti( long int& RefinLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
    pmp = multi->GetPM();
    pmp->t_start = clock();
    pmp->t_end = pmp->t_start;
//    multi->MultiCalcInit( 0 );
    pmp->ITF = pmp->ITG = 0;  
FORCED_AIA:
    multi->MultiCalcInit( "GEMIPM2K" );
	if( pmp->pNP )
    { 
	   if( pmp->ITaia <=30 )       // Foolproof     
		  pmp->IT = 30;
	   else 
		  pmp->IT = pmp->ITaia;     // Setting number of iterations for the smoothing parameter  
    }
	if( multi->AutoInitialApprox( ) == false )
    {
    	multi->MultiCalcIterations( -1 );
    }
    if( pmp->MK == 2 )
    {
    	pmp->pNP = 0; 
    	pmp->MK = 0;
    	goto FORCED_AIA;  // Trying again with AIA set after bad PIA 
    }    
    RefinLoops_ = pmp->W1 + pmp->K2 - 1; // Prec.ref. + Selekt2() loops
    NumIterFIA_ = pmp->ITF;
    NumIterIPM_ = pmp->ITG;
    pmp->IT = pmp->ITG;   // This is to provide correct number of IPM iterations to upper levels

    if( pa.p.PRD < 0 && pa.p.PRD > -50 ) // max 50 loops
    {  // Test refinement loops for highly non-ideal systems Added here by KD on 15.11.2007
       long int pp, pNPo = pmp->pNP,  TotW1 = pmp->W1+pmp->K2-1,  
ITstart = 10,        TotIT = pmp->IT;  //  ITold = pmp->IT,
       pmp->pNP = 1;
       for( pp=0; pp < abs(pa.p.PRD); pp++ )
       {
    	  pmp->IT = ITstart; // Value is important for refinement in highly non-ideal systems!
          if( multi->AutoInitialApprox( ) == false )
          {
             multi->MultiCalcIterations( pp );
          }
          TotIT += pmp->IT - ITstart;
          TotW1 += pmp->W1+pmp->K2-1; 
       }  // end pp loop
       
       pmp->pNP = pNPo;
       pmp->IT = TotIT; // ITold;         

       RefinLoops_ = TotW1; 
       NumIterFIA_ = pmp->ITF;  //   TotITF;
       NumIterIPM_ = pmp->ITG;  //   TotITG;
    }       
pmp->t_end = clock();
pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);
return pmp->t_elap_sec;
}

void TProfil::outMulti( GemDataStream& ff, gstring& path  )
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
    multi->to_file( ff, path );
}

void TProfil::outMultiTxt( const char *path  )
{
    multi->to_text_file( path );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( GemDataStream& ff )
{
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
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( const char* path )
{
      multi->from_text_file_gemipm( path);
}

 bool load = false;

// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
// (only used in standalone GEMIPM2K version)
 //
void TMulti::CompG0Load()
{
  long int j, jj, k, xTP, jb, je=0;
  double Go, Gg, Vv;
  double TC, P;

  DATACH  *dCH = TNode::na->pCSD();
//  DATABR  *dBR = TNodeArray::na->pCNode();

  TC = TNode::na->cTC();
  P = TNode::na->cP();

// if( dCH->nTp <=1 && dCH->nPp <=1 )
  if( dCH->nTp <1 || dCH->nPp <1 || TNode::na->check_TP( TC, P ) == false )
  {
	  char buff[256]; 
	  sprintf( buff, " Temperature %g or pressure %g out of range, or no T/D data are provided\n", 
			  TC, P );
	  Error( "ECompG0Load: " , buff );
      return;
  }
   xTP = TNode::na->check_grid_TP( TC, P );

 if( load && fabs( pmp->TC - TC ) < 1.e-10 &&
            fabs( pmp->P - P ) < 1.e-10 )
   return;    //T, P not changed - problematic for UnSpace! 

 pmp->T = pmp->Tc = TC + C_to_K;
 pmp->TC = pmp->TCc = TC;
 pmp->P = pmp->Pc = P;
 if( dCH->ccPH[0] == PH_AQUEL )
 {
   if( xTP >= 0 )
   {
      pmp->denW = dCH->roW[xTP];
      pmp->epsW = dCH->epsW[xTP];
   }
   else
   {
       pmp->denW = LagranInterp( dCH->Pval, dCH->TCval, dCH->roW,
                          P, TC, dCH->nTp, dCH->nPp,1 );
        //       pmp->denWg = tpp->RoV;
       pmp->epsW = LagranInterp( dCH->Pval, dCH->TCval, dCH->epsW,
                          P, TC, dCH->nTp, dCH->nPp,1 );
       //       pmp->epsWg = tpp->EpsV;
   }
 }
 else
 {
   pmp->denW = 1.;
   pmp->epsW = 78.;
 }

 pmp->RT = R_CONSTANT * pmp->Tc;
 pmp->FRT = F_CONSTANT/pmp->RT;
 pmp->lnP = 0.;
 if( P != 1. ) // ???????
   pmp->lnP = log( P );

 for( k=0; k<pmp->FI; k++ )
 {
   jb = je;
   je += pmp->L1[k];
   // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
   // depending on the presence of these arrays in DATACH and Multi structures
    for( j=jb; j<je; j++ )
    {
      jj =  j * dCH->nPp * dCH->nTp;
      if( xTP >= 0 )
      {
        Go = dCH->G0[ jj+xTP];
        Vv = dCH->V0[ jj+xTP];
      }
     else
     {
       Go = LagranInterp( dCH->Pval, dCH->TCval, dCH->G0+jj,
                          P, TC, dCH->nTp, dCH->nPp,1 );
       Vv = LagranInterp( dCH->Pval, dCH->TCval, dCH->V0+jj,
                          P, TC, dCH->nTp, dCH->nPp, 1 );
     }
     if( pmp->tpp_G )
    	  pmp->tpp_G[j] = Go;
     if( pmp->Guns )
           Gg = pmp->Guns[j];
     else 
           Gg = 0.;	   
     pmp->G0[j] = Cj_init_calc( Go+Gg, j, k ); // Inside this function, pmp->YOF[k] can be added!
     switch( pmp->PV )
     { // put mol volumes of components into A matrix or into the vector of molar volumes
       // to be checked! 
       case VOL_CONSTR:
           if( pmp->Vuns )
              Vv += pmp->Vuns[jj];
           // pmp->A[j*pmp->N+xVol] = tpp->Vm[jj]+Vv;
             pmp->A[j*pmp->N] = Vv; // !!  error
       case VOL_CALC:
       case VOL_UNDEF:
    	     if( pmp->tpp_Vm )
    	    	  pmp->tpp_Vm[j] = Vv;
              if( pmp->Vuns )
                     Vv += pmp->Vuns[j];
 	          pmp->Vol[j] = Vv  * 10.;
              break;
     }
   }
 }
 load = true;
}

// GEM IPM calculation of equilibrium state in MULTI
void TMulti::MultiCalcInit( const char* key )
{
  long int  j,k, jb, je=0;;
  SPP_SETTING *pa = &TProfil::pm->pa;
    strncpy( pmp->stkey, key, EQ_RKLEN );  // needed for the ipmlog error diagnostics
    pmp->Ec = pmp->K2 = pmp->MK = 0;
    pmp->PZ = pa->p.DW; // IPM-2 default
    pmp->W1 = 0;
    pmp->is = 0;
    pmp->js = 0;
    pmp->next  = 0;
    pmp->ln5551 = 4.016533882;  //  ln(55.50837344)  4.0165339; 
    pmp->lowPosNum = pa->p.DcMin;
    pmp->logXw = -16.;
    pmp->logYFk = -9.;
//    pmp->FitVar[4] = pa->p.AG;
    pmp->FitVar[0] = 0.0640000030398369; // pa->aqPar[0]; setting T,P dependent b_gamma parameters 

    pmp->DX = pa->p.DK;

    pmp->T0 = 273.15;    // not used anywhere
    pmp->FX = 7777777.;
    pmp->YMET = 0;
    pmp->PCI = 0.0;

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
    {
        for( j=0; j< pmp->L; j++ )
          pmp->X[j] = pmp->Y[j];
 //       pmp->IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
        ConCalc( pmp->X, pmp->XF, pmp->XFA);
 //       pmp->IT = pmp->ITaia; 
    }
    else // Simplex initial approximation to be done (AIA mode)
    {
    	for( j=0; j<pmp->L; j++ )
    	{                           // cleaning work vectors
    		pmp->X[j] = pmp->Y[j] = pmp->lnGam[j] = pmp->lnGmo[j] = 0.0;
    		pmp->Gamma[j] = 1.0;
    	}
//    	pmp->FitVar[4] = pa->p.AG;
//        pmp->IT = 0;     // needed here to clean LINK_TP_MODE
    }

    CompG0Load(); // Loading thermodynamic data into MULTI structure

    // multicomponent phases and mixing models
 if( pmp->FIs )
 {
     long int k, jb, je=0; 
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
     pmp->PD = pa->p.PD;
	 //           SolModLoad();   Scripts cannot be used here!
    // Calculate Eh, pe, pH,and other stuff
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
    		       && (pmp->XFA[0] > pmp->lowPosNum && pmp->XF[0] > pa->p.XwMin ))
    		       GouyChapman( jb, je, k );  // getting PSIs - elecrtic potentials on surface planes               
    	      }                   
    	   }
        }  
    }
    //   double FitVar3 = pmp->FitVar[3];  // Reset the smoothing factor
    //   pmp->FitVar[3] = 1.0;
        GammaCalc( LINK_TP_MODE);   // Computing DQF, FugPure and G wherever necessary
                                       // Activity coeffs are restored from lnGmo 
    //   pmp->FitVar[3]=FitVar3;
}
 
    // recalculate kinetic restrictions for DC quantities
    if( pmp->pULR && pmp->PLIM )
         Set_DC_limits(  DC_LIM_INIT );

// dynamic demo arrays - do we need it here at all? - DK  14.03.2008
//    for( k=0; k<pmp->FI; k++ )
//    {
//       pmp->XFs[k] = pmp->YF[k];
//        pmp->Falps[k] = pmp->Falp[k];
//    }
}

//-------------------------------------------------------------------------
// internal functions

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
   long int pos = flst_name.rfind("/");
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


// ------------------ End of ms_param.cpp -----------------------




