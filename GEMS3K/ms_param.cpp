//-------------------------------------------------------------------
// $Id: ms_param.cpp 1392 2009-08-10 13:39:26Z gems $
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

// TProfil* TProfil::pm;

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





// test result GEM IPM calculation of equilibrium state in MULTI
long int TMulti::testMulti(  )
{
  if( pm.MK || pm.PZ )
  {
	if( pa.p.PSM == 2 )
	{
      fstream f_log("ipmlog.txt", ios::out|ios::app );
      f_log << "Warning " << pm.stkey << ": " <<  pm.errorCode << ":" << endl;
      f_log << pm.errorBuf << endl;
	}
   return 1L;
  }
  return 0L	;
}

// GEM IPM calculation of equilibrium state in MULTI
double TMulti::ComputeEquilibriumState( long int& RefinLoops_, long int& NumIterFIA_, long int& NumIterIPM_ , TNode* mynode)
{
 return CalculateEquilibriumState( 0, NumIterFIA_, NumIterIPM_ , mynode);
}

void TMulti::outMulti( GemDataStream& ff, gstring& /*path*/  )
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
    to_file( ff );
}

void TMulti::outMultiTxt( const char *path, bool append  )
{
    to_text_file( path, append );
}

// Reading structure MULTI (GEM IPM work structure)
void TMulti::readMulti( GemDataStream& ff , TNode* mynode)
{
    DATACH  *dCH = mynode->pCSD();
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
      from_file( ff );

      // copy intervals for minimizatiom
      if(  dCH->nPp > 1  )
      {
         pm.Pai[0] = dCH->Pval[0];
         pm.Pai[1] = dCH->Pval[dCH->nPp-1];
         pm.Pai[2] = (pm.Pai[1]-pm.Pai[0])/(double)dCH->nPp;
      }
      pm.Pai[3] = dCH->Ptol;
      if(  dCH->nTp > 1  )
      {
         pm.Tai[0] = dCH->TKval[0];
         pm.Tai[1] = dCH->TKval[dCH->nTp-1];
         pm.Tai[2] = (pm.Tai[1]-pm.Tai[0])/(double)dCH->nTp;
      }
      pm.Tai[3] = dCH->Ttol;

  }

// Reading structure MULTI (GEM IPM work structure)
void TMulti::readMulti( TNode *na, const char* path )
{
      from_text_file_gemipm( na, path);
}


// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
// (only used in standalone GEMIPM2K version)
 //
void TMulti::DC_LoadThermodynamicData(TNode* mynode)
{
  long int j, jj, k, xTP, jb, je=0;
  double Go, Gg, Vv, h0=0., S0 = 0., Cp0= 0., a0 = 0., u0 = 0.;
  double TK, P, PPa;

  DATACH  *dCH = mynode->pCSD();
//  DATABR  *dBR = TNodeArray::na->pCNode();

  TK = mynode->cTK();
  PPa = mynode->cP();
  P = PPa/bar_to_Pa;
// if( dCH->nTp <=1 && dCH->nPp <=1 )
  if( dCH->nTp <1 || dCH->nPp <1 || mynode->check_TP( TK, PPa ) == false )
  {
	  char buff[256];
	  sprintf( buff, " Temperature %g or pressure %g out of range, or no T/D data are provided\n",
			  TK, PPa );
	  Error( "ECompG0Load: " , buff );
      return;
  }
  xTP = mynode->check_grid_TP( TK, PPa );

 if( load && fabs( pm.Tc - TK ) < dCH->Ttol /*1.e-10*/ &&
            fabs( pm.P - P ) < dCH->Ptol /*1.e-10*/ )
   return;    //T, P not changed - problematic for UnSpace!

 pm.T = pm.Tc = TK;
 pm.TC = pm.TCc = TK-C_to_K;
 pm.P = pm.Pc = P;

// if( dCH->ccPH[0] == PH_AQUEL )
// {
   for( k=0; k<5; k++ )
   {
     jj =  k * dCH->nPp * dCH->nTp;
     if( xTP >= 0 )
      {
       pm.denW[k] = dCH->denW[jj+xTP]/1e3;
       pm.epsW[k] = dCH->epsW[jj+xTP];
       pm.denWg[k] = dCH->denWg[jj+xTP]/1e3;
       pm.epsWg[k] = dCH->epsWg[jj+xTP];
      }
     else
     {
       pm.denW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denW+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,6 )/1e3;// from test denW enough
       pm.epsW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsW+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );// from test epsW enough
       pm.denWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denWg+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 )/1e3;
       pm.epsWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsWg+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
     }
  }
// }
// else
// {
//   pm.denW[0] = 1.;
//   pm.epsW[0] = 78.;
// }

 pm.RT = R_CONSTANT * pm.Tc;
 pm.FRT = F_CONSTANT/pm.RT;
 pm.lnP = 0.;
 if( P != 1. ) // ???????
   pm.lnP = log( P );

 for( k=0; k<pm.FI; k++ )
 {
   jb = je;
   je += pm.L1[k];
   // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
   // depending on the presence of these arrays in DATACH and Multi structures
    for( j=jb; j<je; j++ )
    {
      jj =  j * dCH->nPp * dCH->nTp;
      if( xTP >= 0 )
      {
        Go = dCH->G0[ jj+xTP];
        Vv = dCH->V0[ jj+xTP]*1e5;
        if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
        if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
        if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
        if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
        if( dCH->U0 ) h0 = dCH->U0[ jj+xTP];
      }
     else
     {
       Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
       Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 5 )*1e5;
       if( dCH->S0 ) S0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->S0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 4 ); // from test S0[Ca+2] enough
       if( dCH->H0 ) h0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->H0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
       if( dCH->Cp0 ) Cp0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->Cp0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 3 ); // from test Cp0[Ca+2] not more
       if( dCH->A0 ) a0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->A0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
       if( dCH->U0 ) u0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->U0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
     }
     if( pm.tpp_G )
    	  pm.tpp_G[j] = Go;
     if( pm.Guns )
           Gg = pm.Guns[j];
     else
           Gg = 0.;
     pm.G0[j] = ConvertGj_toUniformStandardState( Go+Gg, j, k ); // Inside this function, pm.YOF[k] can be added!
     switch( pm.PV )
     { // put mol volumes of components into A matrix or into the vector of molar volumes
       // to be checked!
       case VOL_CONSTR:
           if( pm.Vuns )
              Vv += pm.Vuns[jj];
           // pm.A[j*pm.N+xVol] = tpp->Vm[jj]+Vv;
             pm.A[j*pm.N] = Vv; // !!  error
       case VOL_CALC:
       case VOL_UNDEF:
    	     if( pm.tpp_Vm )
    	    	  pm.tpp_Vm[j] = Vv;
              if( pm.Vuns )
                     Vv += pm.Vuns[j];
 	          pm.Vol[j] = Vv  * 10.;
              break;
     }
     if( pm.S0 ) pm.S0[j] = S0;
     if( pm.H0 ) pm.H0[j] = h0;
     if( pm.Cp0 ) pm.Cp0[j] = Cp0;
     if( pm.A0 ) pm.A0[j] = a0;
     if( pm.U0 ) pm.U0[j] = u0;
   }
 }
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

