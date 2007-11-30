//-------------------------------------------------------------------
// $Id: ms_multi_format.cpp 774 2006-07-26 08:45:45Z gems $
//
// Implementation of text writing/reading IPM, DCH and DBR files
//
// Copyright (C) 2006-2007 S.Dmytriyeva
//
// This file is part of the GEM-Vizor library and GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------
//#include  <iostream>

#include "io_arrays.h"
#include "m_param.h"
#include "node.h"
#include <iomanip>

#ifdef IPMGEMPLUGIN
  istream& f_getline(istream& is, gstring& str, char delim);
#endif

bool _comment = true;

//===================================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here.
//
outField MULTI_static_fields[8] =  {
  { "pa_PE" , 0 , 0 },
  { "PV" , 0 , 0 },
  { "PSOL" , 0 , 0 },
  { "PAalp" , 0 , 0 },
  { "PSigm" , 0 , 0 },
  { "Lads" , 0 , 0 },
  { "FIa" , 0 , 0 },
  { "FIat" , 0 , 0 }
};

outField MULTI_dynamic_fields[66] =  {
//read dynamic (array) data from the txt input file
   {  "sMod", 1 , 0 },
   {  "LsMod", 1 , 0 },
   {  "LsMdc", 1 , 0 },
   {  "B", 1 , 0 },
   {  "DCCW", 0 , 0 },  // placeholder - something else can be used here
   {  "Pparc", 0 , 0 },
   {  "GEX", 0 , 0 },
   {  "lnGmf", 0 , 0 },
   {  "RLC", 0 , 0 },
   {  "RSC", 0 , 0 },
   {  "DLL", 0 , 0 },
   {  "DUL", 0 , 0 },
   {  "Aalp", 0 , 0 },
   {  "Sigw", 0 , 0 },
   {  "Sigg", 0 , 0 },
   {  "YOF", 0 , 0 },
   {  "Nfsp", 1 , 0 },
   {  "MASDT", 1 , 0 },
   {  "C1", 1 , 0 },
   {  "C2", 1 , 0 },
   {  "C3", 1 , 0 },
   {  "pCh", 1 , 0 },
   {  "SATX", 1 , 0 },
   {  "MASDJ", 1 , 0 },
   {  "SCM", 1 , 0 },
   {  "SACT", 1 , 0 },
   {  "DCads", 1 , 0 },
   // static
   { "pa_DB" , 0, 0 },
   { "pa_DHB", 0 , 0 },
   { "pa_EPS" , 0 , 0 },
   { "pa_DK", 0 , 0 },
   { "pa_DF" , 0 , 0 },
   { "pa_DP", 0 , 0 },
   { "pa_IIM", 0 , 0 },
   { "pa_PD" , 0 , 0 },
   { "pa_PRD" , 0 , 0 },
   { "pa_AG" , 0 , 0 },
   { "pa_DGC" , 0 , 0 },
   { "pa_PSM" , 0 , 0 },
   { "pa_GAR" , 0 , 0 },
   { "pa_GAH" , 0 , 0 },
   { "pa_DS", 0 , 0 },
   { "pa_XwMin" , 0 , 0 },
   { "pa_ScMin", 0 , 0 },
   { "pa_DcMin" , 0 , 0 },
   { "pa_PhMin" , 0 , 0 },
   { "pa_ICmin" , 0 , 0 },
   { "pa_PC" , 0 , 0 },
   { "pa_DFM" , 0 , 0 },
   { "pa_DFYw" , 0 , 0 },
   { "pa_DFYaq" , 0 , 0 },
   { "pa_DFYid" , 0 , 0 },
   { "pa_DFYr" , 0 , 0 },
   { "pa_DFYh" , 0 , 0 },
   { "pa_DFYc" , 0 , 0 },
   { "pa_DFYs", 0 , 0 },
   { "pa_DW", 0 , 0 },
   { "pa_DT", 0 , 0 },
   { "pa_GAS", 0 , 0 },
   { "pa_DNS" , 0 , 0 },
   { "pa_IEPS" , 0 , 0 },
   { "pKin" , 0 , 0 },
   { "pa_DKIN" , 0 , 0 },
   { "mui" , 0 , 0 },
   { "muk" , 0 , 0 },
   { "muj" , 0 , 0 }
};


//===================================================================

void TMulti::to_text_file_gemipm( const char *path, bool addMui, bool with_comments )
{
  SPP_SETTING *pa = &TProfil::pm->pa;
   _comment = with_comments;
   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;

#ifndef IPMGEMPLUGIN
   PAalp = syp->PAalp;
   PSigm = syp->PSigm;
   EpsW = TProfil::pm->tpp->EpsW;
   RoW = TProfil::pm->tpp->RoW;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
   EpsW = EpsW_;
   RoW = RoW_;
#endif
  fstream ff( path, ios::out );
  ErrorIf( !ff.good() , path, "Fileopen error");
  TPrintArrays  prar(ff);

if( _comment )
{   ff << "# GEMIPM2K v. 2.2.0" << endl;
   ff << "# Prototype 28.11.2007" << endl;
   ff << "# Comments can be marked with # $ ; as the first character in the line" << endl << endl;
   ff << "# Template for the ipm-dat text input file for the internal MULTI data" << endl;
   ff << "# (should be read after the DATACH file and before DATABR files)" << endl << endl;
   ff << "# ID key of the initial chemical system definition" << endl;
}
  ff << "\"" << pmp->stkey << "\"" << endl << endl;

if( _comment )
{  ff << "## (1) Important flags that affect memory allocation" << endl;
   ff << "# PE: Flag for using electroneutrality condition in GEM IPM calculations " << endl;
}
   ff << left << setw(12) << "<pa_PE> " <<  right << setw(8) << pa->p.PE << endl;
//   ff << "# 'E'                1" << endl;
//   ff << left << setw(12) << "<E> " <<  right << setw(8) << pmp->E << endl;
   if( _comment )
     ff << "\n# PV: Flag for the volume balance constraint (on Vol IC)" << endl;
   ff << left << setw(12) << "<PV> " <<  right << setw(8) << pmp->PV << endl;
//   ff << "# These dimensions can be calculated from the DATACH information" << endl;
//   ff << "# 'Ls'              23" << endl;
//   ff << "# 'LO'              18" << endl;
//   ff << "# 'PG'               4" << endl;
//   ff << left << setw(12) << "<Ls> " <<  right << setw(8) << pmp->Ls << endl;
//   ff << left << setw(12) << "<LO> " <<  right << setw(8) << pmp->LO << endl;
//   ff << left << setw(12) << "<PG> " <<  right << setw(8) << pmp->PG << endl;
   if( _comment )
       ff << "\n# PSOL: Total number of DCs in liquid hydrocarbon phases" << endl;
   ff << left << setw(12) << "<PSOL> " <<  right << setw(8) << pmp->PSOL << endl;
//   ff << "# Do not know if this stuff is really necessary" << endl;
//   ff << "# 'GWAT'         55.51" << endl;
//   ff << "# 'EpsW'       78.2451" << endl;
//   ff << "# 'RoW'       0.997061" << endl << endl;
//   ff << left << setw(12) << "<GWAT> " <<  right << setw(8) << pmp->GWAT << endl;
//   ff << left << setw(12) << "<EpsW> " <<  right << setw(8) << EpsW << endl;
//   ff << left << setw(12) << "<RoW>  " <<  right << setw(8) << RoW << endl;
   if( _comment )
     ff << "\n# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases " << endl;
   ff << left << setw(12) << "<PAalp> " <<  right << setw(6) <<
      "\'" << PAalp << "\'" << endl;
   if( _comment )
    ff << "\n# PSigm: Flag for using (+) or ignoring (-) specific surface free energies of phase interfaces " << endl;
   ff << left << setw(12) << "<PSigm> " <<  right << setw(6) <<
      "\'" << PSigm << "\'" << endl;
if( _comment )
{  ff << "\n## (2) Important dimensionalities that affect memory allocation" << endl;
   ff << "# Lads: Total number of dependent components in sorption phases included in this system" << endl;
}
ff << left << setw(12) << "<Lads> " <<  right << setw(8) << pmp->Lads << endl;
if( _comment )
  ff << "# FIa: Number of sorption phases included in this system" << endl;
ff << left << setw(12) << "<FIa> " <<  right << setw(8) << pmp->FIa << endl;
if( _comment )
  ff << "# FIat: Allowed number of surface types per adsorption phase (default: 6 if FIa > 0)" << endl;
ff << left << setw(12) << "<FIat> " <<  right << setw(8) << pmp->FIat << endl << endl;
//   ff << left << setw(12) << "<FIat> " <<  right << setw(8) << pmp->FIat << endl;
//   ff << left << setw(12) << "<sitNc> " <<  right << setw(8) << pmp->sitNcat << endl;
//   ff << left << setw(12) << "<sitNa> " <<  right << setw(8) << pmp->sitNan << endl;
ff << "\n<END_DIM>\n\n";

// static data not affected by dimensionalities
  if( _comment )
  { ff << "## (3) Controls of the numerical behavior of the GEM IPM algorithm" << endl;
    ff << "#      - Need to be changed only in rare special cases" << endl;
    ff << "# DB - Minimum amount of independent component in bulk composition (except charge Zz), moles" << endl;
  }
   ff << left << setw(12) << "<pa_DB> " <<  right << setw(8) << pa->p.DB << endl;
   if( _comment )
     ff << "\n# DHB: Maximum allowed mass balance residual (moles) for major independent components" << endl;
   ff << left << setw(12) << "<pa_DHB> " << right << setw(8) <<  pa->p.DHB << endl;
   if( _comment )
     ff << "\n# EPS: Precision criterion of the simplex() procedure for automatic initial approximation" << endl;
   ff << left << setw(12) << "<pa_EPS> " <<  right << setw(8) << pa->p.EPS << endl;
   if( _comment )
     ff << "\n# DK: IPM convergence threshold for the Dikin criterion (1e-6 < DK < 1e-4)" << endl;
   ff << left << setw(12) << "<pa_DK> " <<  right << setw(8) << pa->p.DK << endl;
   if( _comment )
     ff << "\n# DB: Cutoff (threshold) for the amount of phase for phase elimination " << endl;
   ff << left << setw(12) << "<pa_DS> " << right << setw(8) <<  pa->p.DS << endl;
   if( _comment )
     ff << "\n# DF: Threshold for Karpov's criterion (Fa > DF) for a lost stable phase to be inserted" << endl;
   ff << left << setw(12) << "<pa_DF> " <<  right << setw(8) << pa->p.DF << endl;
   if( _comment )
     ff << "# DFM: Threshold for Karpov's criterion (Fa < -DFM) for a present unstable phase to be eliminated" << endl;
   ff << left << setw(12) << "<pa_DFM> " <<  right << setw(8) << pa->p.DFM << endl;
   if( _comment )
     ff << "\n# DP: Maximum allowed number of iterations in the EnterFeasibleDomain() procedure" << endl;
   ff << left << setw(12) << "<pa_DP> " << right << setw(8) << pa->p.DP << endl;
   if( _comment )
     ff << "\n# IIM: Maximum allowed number of iterations in the MainIPM_Descent() procedure" << endl;
   ff << left << setw(12) << "<pa_IIM> " << right << setw(8) <<  pa->p.IIM << endl;
   if( _comment )
     ff << "\n# PD: Control on calling built-in Debye-Hueckel() and other models for aqueous activity coefficients" << endl;
   ff << left << setw(12) << "<pa_PD> " <<  right << setw(8) << pa->p.PD << endl;
   if( _comment )
   { ff << "\n# PRD: Positive: control of calling IPM_gamma() on iterations of GEM EFD and IPM (default 3)" << endl;
     ff <<   "#      Negative: the number of additional EFD-IPM loops to improve the GEM final solution" << endl;
   }
   ff << left << setw(12) << "<pa_PRD> " <<  right << setw(8) << pa->p.PRD << endl;
   if( _comment )
     ff << "\n# AG,DGC: Smoothing parameters controlling convergence in highly non-ideal systems " << endl;
   ff << left << setw(12) << "<pa_AG> " <<  right << setw(8) << pa->p.AG << endl;
   ff << left << setw(12) << "<pa_DGC> " <<  right << setw(8) << pa->p.DGC << endl;
   if( _comment )
     ff << "\n# PSM: Flag for using initial activity coefficients in simplex() approximation (1-enable, 0-disable)" << endl;
   ff << left << setw(12) << "<pa_PSM> " <<  right << setw(8) << pa->p.PSM << endl;
   if( _comment )
     ff << "# Activity coefficient values for simplex(): GAR for major and GAH for minor components" << endl;
   ff << left << setw(12) << "<pa_GAR> " <<  right << setw(8) << pa->p.GAR << endl;
   ff << left << setw(12) << "<pa_GAH> " <<  right << setw(8) << pa->p.GAH << endl;
   if( _comment )
   {  ff << "\n# _Min: Cutoffs for elimination of: Xw water solvent; Sc - solid sorbent; Dc - solution or surface species; " << endl;
      ff << "# Ph - non-electrolyte solution phases" << endl;
   }
   ff << left << setw(12) << "<pa_XwMin> " <<  right << setw(8) << pa->p.XwMin << endl;
   ff << left << setw(12) << "<pa_ScMin> " <<  right << setw(8) << pa->p.ScMin << endl;
   ff << left << setw(12) << "<pa_DcMin> " <<  right << setw(8) << pa->p.DcMin << endl;
   ff << left << setw(12) << "<pa_PhMin> " <<  right << setw(8) << pa->p.PhMin << endl;
   if( _comment )
     ff << "\n# ICmin: Minimal ionic strength, below which the aqueous activity coefficients are set to 1" << endl;
   ff << left << setw(12) << "<pa_ICmin> " <<  right << setw(8) << pa->p.ICmin << endl;
   if( _comment )
     ff << "\n# PC: Mode of Selekt2() procedure operation" << endl;
   ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   if( _comment )
   {  ff << "# DFY: Insertion amounts used after the simplex() initial approximation and in Selekt2() algorithm" << endl;
      ff << "# DFYw - water solvent;" << endl;
   }
   ff << left << setw(12) << "<pa_DFYw> " <<  right << setw(8) << pa->p.DFYw << endl;
   if( _comment )
     ff << "# DFYaq - aqueous species;" << endl;
   ff << left << setw(12) << "<pa_DFYaq> " <<  right << setw(8) << pa->p.DFYaq << endl;
   if( _comment )
     ff << "\n# DFYid - ideal solution components;" << endl;
//   ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   ff << left << setw(12) << "<pa_DFYid> " <<  right << setw(8) << pa->p.DFYid << endl;
   if( _comment )
     ff << "# DFYr - major solution components;" << endl;
   ff << left << setw(12) << "<pa_DFYr> " <<  right << setw(8) << pa->p.DFYr << endl;
   if( _comment )
     ff << "# DFYh - minor solution components;" << endl;
   ff << left << setw(12) << "<pa_DFYh> " <<  right << setw(8) << pa->p.DFYh << endl;
   if( _comment )
     ff << "# DFYc - single-component phase; " << endl;
   ff << left << setw(12) << "<pa_DFYc> " <<  right << setw(8) << pa->p.DFYc << endl;
   if( _comment )
     ff << "# Insertion amount of single-component phase (Selekt2() algorithm only)" << endl;
   ff << left << setw(12) << "<pa_DFYs> " << right << setw(8) <<  pa->p.DFYs << endl;
   if( _comment )
   {  ff << "\n# Parameters for high-accuracy IPM-2 algorithm " << endl;
      ff << "# DW: Number of the IPM-2 enhancement loops for high-accuracy mass balance (from 0 to 14)" << endl;
   }
   ff << left << setw(12) << "<pa_DW> " << right << setw(8) << pa->p.DW  << endl;
   if( _comment )
     ff << "# DT: Exponent for dual-thermodynamic restoring of low amounts of solutes (+2 to -5), relative to DHBM" << endl;
   ff << left << setw(12) << "<pa_DT> " << right << setw(8) << pa->p.DT  << endl;
   if( _comment )
     ff << "\n# GAS: IPM2 balance accuracy control factor DHBM[i]/b[i]" << endl;
   ff << left << setw(12) << "<pa_GAS> " << right << setw(8) <<  pa->p.GAS << endl;
//   ff << "# 'pa_DG'        1e-30  now internal in LU decomposition procedure from JAMA-TNT" << endl << endl;
//  ff << left << setw(12) << "<pa_DG> " <<  right << setw(8) << pa->p.DG << endl;
   if( _comment )
     ff << "# DNS: Standard surface density (nm-2) for calculating activity of surface species" << endl;
   ff << left << setw(12) << "<pa_DNS> " <<  right << setw(8) << pa->p.DNS << endl;
   if( _comment )
     ff << "# IEPS: Control parameter of SACT calculation in sorption/surface complexation models" << endl;
   ff << left << setw(12) << "<pa_IEPS> " <<  right << setw(8) << pa->p.IEPS << endl;
   if( _comment )
     ff << "\n# pKin: Flag for using metastability/kinetic constraints on calculated DC amounts (see dll, dul arrays)" << endl;
   ff << left << setw(12) << "<pKin> " <<  right << setw(8) << pmp->PLIM << endl;
   if( _comment )
     ff << "# DKIN: Tolerance on amount of DC with two-side metastability constraints set in dll, dul (moles) " << endl;
   ff << left << setw(12) << "<pa_DKIN> " <<  right << setw(8) << pa->p.DKIN << endl;
//   ff << "# 'pa_PLLG'          0 used only in GEMS-PSI shell" << endl;
//   ff << left << setw(12) << "<pa_PLLG> " <<  right << setw(8) << pa->p.PLLG << endl;

//dynamic arrays
if( pm.FIs > 0 && pm.Ls > 0 )
{
  if( _comment )
  {   ff << "\n## (4) Initial data for multicomponent phases (see DATACH file for dimension nPHs)" << endl;
      ff << "# sMod: Codes for mixing models of multicomponent phases";
  }
  prar.writeArray(  "sMod", pmp->sMod[0], pmp->FIs, 6 );

int LsModSum;
int LsIPxSum;
int LsMdcSum;
getLsModsum( LsModSum, LsIPxSum );
getLsMdcsum( LsMdcSum );

  if( _comment )
  {  ff << "\n\n# LsMod: Dimensions for parameters of non-ideal mixing models for each multicomponent phase" << endl;
     ff << "# Number of parameters per phase";
  }
  prar.writeArray(  "LsMod", pmp->LsMod, pmp->FIs*3, 3);

if(LsIPxSum )
 {
   if( _comment )
      ff << "\n\n# IPx: List of indexes of interaction parameters for non-ideal solutions ";
  prar.writeArray(  "IPxPH", pmp->IPx,  LsIPxSum);
}
  if(LsModSum )
   {
     if( _comment )
        ff << "\n\n# PMc: Parameters of non-ideal mixing models for multicomponent phases ";
    prar.writeArray(  "PMc", pmp->PMc,  LsModSum);
   }
   if( _comment )
     ff << "\n\n# LsMdc: Number of parameters per component of the phase";
   prar.writeArray(  "LsMdc", pmp->LsMdc, pmp->FIs);
   if(LsMdcSum )
   {   if( _comment )
          ff << "\n\n# DMc: Parameters of non-ideal mixing models for components in the phase ";
    prar.writeArray(  "DMc", pmp->DMc,  LsMdcSum);
   }
  }
  if( _comment )
  {  ff << "\n\n## (5) Some data arrays which are not provided in DATACH and DATABR files" << endl;
     ff << "# B: Full total bulk composition of the initial system (vector b) - see DATACH file for dimension nIC";
  }
  prar.writeArray(  "B", pmp->B,  pmp->N);
  if( _comment )
   {  ff << "\n\n# Initial data for DCs - see DATACH file for dimensions nDC, nDCs" << endl;
      ff << "# Pparc: Partial pressures (fugacities) of dependent components";
   }
  prar.writeArray(  "Pparc", pmp->Pparc,  pmp->L);
 //  ff << "\n\n# This is not necessary - can be calculated from G0 ???????????";
 // prar.writeArray(  "G0", pmp->G0,  pmp->L);
   if( _comment )
      ff << "\n\n# GEX: for dependent components, G0 increments";
  prar.writeArray(  "GEX", pmp->GEX,  pmp->L);
   if( _comment )
      ff << "\n\n# lnGmf:  Fixed (initial) activity coefficients";
  prar.writeArray(  "lnGmf", pmp->lnGmf,  pmp->L);
//   if( pmp->E )
//   {
//     if( _comment )
//        ff << "\n\n# DC Unit formula charges - can be extracted from the stoich. matrix ????";
 //   prar.writeArray(  "EZ", pmp->EZ,  pmp->L);
//   }
   if( _comment )
   {  ff << "\n\n# (6) Section for metastability/ kinetic constraints" << endl;
      ff << "# RLC: Codes of metastability/kinetic constraints for DCs";
   }
  prar.writeArray(  "RLC", pmp->RLC, pmp->L, 1 );
   if( _comment )
     ff << "\n\n# RSC: Units of metastability/kinetic constraints for DCs (see vectors dul, dll)";
  prar.writeArray(  "RSC", pmp->RSC, pmp->L, 1 );
   if( _comment )
     ff << "\n\n# DLL: Vector of lower metastability constraints on DC amounts in the system";
  prar.writeArray(  "DLL", pmp->DLL,  pmp->L);
   if( _comment )
     ff << "\n\n# DUL: Vector of upper metastability constraints on DC amounts in the system";
  prar.writeArray(  "DUL", pmp->DUL,  pmp->L);
   if( _comment )
   {  ff << "\n\n# (7) Initial data for phases" << endl;
      ff << "\n# Aalp: Specific surface areas of phases (whole list), set 0 if unknown";
   }
  prar.writeArray(  "Aalp", pmp->Aalp,  pmp->FI);
   if( PSigm != S_OFF )
   {
      if( _comment )
         ff << "\n\n# Sigw Specific surface free energy for phase-water interface (J/m2; 0 if unknown)";
     prar.writeArray(  "Sigw", pmp->Sigw,  pmp->FI);
      if( _comment )
         ff << "\n\n# Sigg Specific surface free energy for phase-gas interface (for future versions)";
     prar.writeArray(  "Sigg", pmp->Sigg,  pmp->FI);
   }
   if( _comment )
      ff << "\n\n# YOF: Surface energy or metastability parameters for phases (J/g)";
  prar.writeArray(  "YOF", pmp->YOF,  pmp->FI);
   if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
      if( _comment )
      {  ff << "\n\n# (8) Initial data for sorption" << endl;
         ff << "\n# Nfsp: Function of sorbent surface allocated to surface types";
      }
     prar.writeArray(  "Nfsp", &pmp->Nfsp[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# MASDT: Total maximum site  density per surface type, mkmol/g";
     prar.writeArray(  "MASDT", &pmp->MASDT[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# C1: Inner capacitance density parameter C1 (TLM, BSM, CCM), F/m2";
     prar.writeArray(  "C1", &pmp->XcapA[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# C2: Outer capacitance density parameter C2 (TLM, 3LM), F/m2";
     prar.writeArray(  "C2", &pmp->XcapB[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# C3: Third capacitance density parameter C3 (reserved)";
     prar.writeArray(  "C3", &pmp->XcapF[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# pCh: Density of permanent surface type charge (mkeq/m2)";
     prar.writeArray(  "pCh", &pmp->Xetaf[0][0], pmp->FIs*pmp->FIat, pmp->FIat);

     if( _comment )
     {   ff << "\n# SATX: Setup of surface sites and specres: link to";
         ff << "\n# [0] surface type; [1] sorbent emd member;";
         ff << "\n# [2] surface site in surf. type; [3] surface EDL plane";
     }
     prar.writeArray(  "SATX", &pmp->SATX[0][0], pmp->Lads*4, 4);
      if( _comment )
      {  ff << "\n# MASDJ: Parameters of surface binding model:";
         ff << "\n# [0] max site density mmol/g; [1] charge allocated to 0 plane;";
         ff << "\n# [2] charge allocated to beta -or third plane; [3] Frumkin interaction parameter;";
         ff << "\n# [4] dentateness or CN; [5] reserved isoterm parameter.";
      }
     prar.writeArray(  "MASDJ", &pmp->MASDJ[0][0], pmp->Lads*DFCN, DFCN);
      if( _comment )
         ff << "\n# SCM: Classifier of EDL models applied to surface types.";
     prar.writeArray(  "SCM", pmp->SCM[0], pmp->FIs, pmp->FIat );
      if( _comment )
         ff << "\n# SACT: Classifier of applied SACT terms.";
     prar.writeArray(  "SACT", pmp->SATT, pmp->Lads, 1 );
      if( _comment )
         ff << "\n# DCads: Classifier of species in sorption phase.";
     prar.writeArray(  "DCads", pmp->DCC3, pmp->Lads, 1 );
    }
//outArray( ff, "Vol", pmp->Vol,  pmp->L);
//outArray( ff, "G0", pmp->G0,  pmp->L);
//outArray( ff, "PUL", pmp->PUL,  pmp->L);
//outArray( ff, "PLL", pmp->PLL,  pmp->L);
//outArray( ff, "lnGam", pmp->lnGam,  pmp->L);
//outArray( ff, "F0", pmp->F0,  pmp->L);

/*
   if( pm.sitNcat*pm.sitNcat )
    prar.writeArray(  "sitE", pmp->sitE, pmp->sitNcat*pmp->sitNan );
   if( pm.sitNcat )
    prar.writeArray(  "sitXc", pmp->sitXcat, pmp->sitNcat );
   if( pm.sitNan )
     prar.writeArray(  "sitXa", pmp->sitXan, pmp->sitNan );
*/

if( addMui )
{
  if( _comment )
    ff << "\n\n# mui: IC indices in RMULTS IC list";
  prar.writeArray(  "mui", pmp->mui,  pmp->N);
  if( _comment )
    ff << "\n\n# muk: Phase indices in RMULTS phase list";
  prar.writeArray(  "muk", pmp->muk,  pmp->FI);
  if( _comment )
    ff << "\n\n# muj: DC indices in RMULTS DC list";
  prar.writeArray(  "muj", pmp->muj,  pmp->L);
}

 if( _comment )
   ff << "\n\n# End of file" << endl;

}

void TMulti::from_text_file_gemipm( const char *path )
{
  SPP_SETTING *pa = &TProfil::pm->pa;
  DATACH  *dCH = TNode::na->pCSD();
  int ii, nfild;

   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;

  memset( &pm.N, 0, 38*sizeof(short));
  memset( &pm.TC, 0, 55*sizeof(double));
  // get sizes from DATACH
  pmp->TC = pmp->TCc = 25.;
  pmp->T = pmp->Tc =298.15;
  pmp->P = pmp->Pc = 1.;
  pmp->N = pmp->NR = dCH->nIC;
  pmp->L = dCH->nDC;
  pmp->FI = dCH->nPH;
  pmp->FIs = dCH->nPS;
  pmp->Ls = 0; //dCH->nDCs;
  for( ii=0; ii<dCH->nPS; ii++)
  {
    pmp->Ls += dCH->nDCinPH[ii];
    if( dCH->ccPH[ii] == 'a' )
     pmp->LO = pmp->Ls-1;
    if( dCH->ccPH[ii] == 'g' || dCH->ccPH[ii] == 'p' || dCH->ccPH[ii] == 'f')
      pmp->PG = dCH->nDCinPH[ii];
  }
  // setup default constants
  pa->p.PE =  pmp->E = 1;
  pmp->PV = 0;
  pmp->PSOL = 0;
  PAalp = '+';
  PSigm = '+';
  pmp->Lads = 0;
  pmp->FIa = 0;
  pmp->FIat = 0; //6
  pmp->PLIM  = 1;

  // reads sizes and constants from txt file
  fstream ff( path, ios::in );
  ErrorIf( !ff.good() , path, "Fileopen error");

// static data
   TReadArrays  rdar( 8, MULTI_static_fields, ff);
   gstring str;
   rdar.skipSpace();
   f_getline( ff, str, '\n');
   memcpy( pmp->stkey, str.c_str(), EQ_RKLEN );

   nfild = rdar.findNext();
   while( nfild >=0 )
   {
     switch( nfild )
     {
       case 0: rdar.readArray("pa_PE" , &pa->p.PE, 1);
                 pmp->E = pa->p.PE;
              break;
       case 1: rdar.readArray("PV" , &pmp->PV, 1);
              break;
       case 2: rdar.readArray("PSOL" , &pmp->PSOL, 1);
              break;
       case 3: rdar.readArray("PAalp" , &PAalp, 1, 1);
              break;
       case 4: rdar.readArray("PSigm" , &PSigm, 1, 1);
              break;
       case 5: rdar.readArray("Lads" , &pmp->Lads, 1);
              break;
       case 6: rdar.readArray("FIa" , &pmp->FIa, 1);
              break;
       case 7: rdar.readArray("FIat" , &pmp->FIat, 1);
              break;
    }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from the MULTY structure";
    Error( "Error", ret);
  }

//   if( dCH->ccPH[0] == PH_AQUEL )
//   {
//     RoW = dCH->roW[0];
//     EpsW = dCH->epsW[0];
//   }
//   else
//  {
    RoW = (float)0.99706137180;
    EpsW = (float)78.245147705;
//  }

#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
   EpsW_ = EpsW;
   RoW_ =  RoW;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#endif

// get dynamic data from DATACH file
  for( ii=0; ii<dCH->nPH; ii++)
    pmp->L1[ii] = dCH->nDCinPH[ii];

  for( ii=0; ii<dCH->nIC*dCH->nDC; ii++)
    pmp->A[ii] = dCH->A[ii];

  if( pmp->EZ )
  { int iZ=-1;
    for(  ii=0; ii<dCH->nDC; ii++ )
     if( dCH->ccIC[ii] == IC_CHARGE )
         break;
    if( ii< dCH->nDC )
    { iZ = ii;
      for( ii=0; ii<dCH->nDC; ii++)
          pmp->EZ[ii] = pmp->A[pmp->N*ii+iZ];
    }
  }

  for( ii=0; ii< dCH->nIC; ii++ )
  { pmp->Awt[ii]  = dCH->ICmm[ii];
    memset(pmp->SB[ii], ' ', MaxICN*sizeof(char));
    memcpy( pmp->SB[ii], dCH->ICNL[ii], MaxICN*sizeof(char) );
    pmp->SB[ii][MaxICN] = dCH->ccIC[ii];
    pmp->ICC[ii] =  dCH->ccIC[ii];
  }

if( fabs(dCH->DCmm[0]) < 1e-32 )  // Restore DCmm if skipped from the DCH file
  for( int jj=0; jj< dCH->nDC; jj++ )  // Added by DK on 03.03.2007
  {
    dCH->DCmm[jj] = 0.0;
    for( ii=0; ii< dCH->nIC; ii++ )
       dCH->DCmm[jj] += dCH->ICmm[ii]*dCH->A[jj*dCH->nIC+ii];
  }

  for( ii=0; ii< dCH->nDC; ii++ )
  {
    pmp->MM[ii] = dCH->DCmm[ii];
    pmp->DCC[ii] = dCH->ccDC[ii];
    memcpy( pmp->SM[ii], dCH->DCNL[ii], MaxDCN*sizeof(char));
  }

  for( ii=0; ii< dCH->nPH; ii++ )
  {
     memset( pmp->SF[ii], ' ', MaxPHN*sizeof(char) );
     memcpy( pmp->SF[ii]+4, dCH->PHNL[ii], MaxPHN*sizeof(char));
     pmp->SF[ii][0] = dCH->ccPH[ii];
     pmp->PHC[ii] = dCH->ccPH[ii];
  }

// !!!!  memcpy( pmp->DCCW, dCH->ccDCW, dCH->nDC*sizeof(char));
  // set up DCCW
  ConvertDCC();

//reads dynamic values from txt file
   TReadArrays  rddar( 66, MULTI_dynamic_fields, ff);

// set up array flags for permanent fields

 if( !( pm.FIs > 0 && pm.Ls > 0 ) )
 {
   rddar.setNoAlws( short(0) /*"sMod"*/);
   rddar.setNoAlws( 1 /*"LsMod"*/);
   rddar.setNoAlws( 2 /*"LsMdc"*/);
 }
 if( PSigm == S_OFF )
 {
   rddar.setNoAlws( 13 /*"Sigw"*/);
   rddar.setNoAlws( 14 /*"Sigg"*/);
 }
 if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
 { /* ADSORPTION AND ION EXCHANGE */
   rddar.setNoAlws( 16 /*"Nfsp"*/);
   rddar.setNoAlws( 17 /*"MASDT"*/);
   rddar.setNoAlws( 18 /*"C1"*/);
   rddar.setNoAlws( 19 /*"C2"*/);
   rddar.setNoAlws( 20 /*"C3"*/);
   rddar.setNoAlws( 21 /*"pCh"*/);
   rddar.setNoAlws( 22 /*"SATX"*/);
   rddar.setNoAlws( 23 /*"MASDJ"*/);
   rddar.setNoAlws( 24 /*"SCM"*/);
   rddar.setNoAlws( 25 /*"SACT"*/);
   rddar.setNoAlws( 26 /*"DCads"*/);
 }

  // Read dynamic arrays
  nfild = rddar.findNext();
  while( nfild >=0 )
  {
    switch( nfild )
    { case 0: if( !pmp->sMod )
                Error( "Error", "Array sMod is not used in this problem");
              rddar.readArray( "sMod" , pmp->sMod[0], pmp->FIs, 6 );
              break;
      case 1:{ if( !pmp->LsMod )
                Error( "Error", "Array LsMod is not used in this problem");
              rddar.readArray( "LsMod" , pmp->LsMod, pmp->FIs*3) ;
              int LsModSum;
              int LsIPxSum;
              getLsModsum( LsModSum, LsIPxSum );
              if(LsIPxSum )
              { rddar.readNext( "IPxPH");
                if(!pmp->IPx )
                  pmp->IPx = new short[LsIPxSum];
                rddar.readArray( "IPxPH", pmp->IPx,  LsIPxSum);
              }
              if(LsModSum )
              { rddar.readNext( "PMc");
                if(!pmp->PMc )
                  pmp->PMc = new float[LsModSum];
                rddar.readArray( "PMc", pmp->PMc,  LsModSum);
              }
              break;
             }
      case 2: { if( !pmp->LsMdc )
                   Error( "Error", "Array LsMdc not used in this problem");
                rddar.readArray( "LsMdc" , pmp->LsMdc, pmp->FIs );
                int LsMdcSum;
                getLsMdcsum( LsMdcSum );
                if(LsMdcSum )
                { rddar.readNext( "DMc");
                  if(!pmp->DMc )
                     pmp->DMc = new float[LsMdcSum];
                  rddar.readArray( "DMc", pmp->DMc,  LsMdcSum);
                }
                break;
              }
      case 3: rddar.readArray( "B", pmp->B,  pmp->N);
              break;
      case 4: rddar.readArray( "DCCW", pmp->DCCW,  pmp->L, 1);
              break;
      case 5: rddar.readArray( "Pparc", pmp->Pparc,  pmp->L);
              break;
      case 6: rddar.readArray( "GEX", pmp->GEX,  pmp->L);
              break;
      case 7: rddar.readArray( "lnGmf", pmp->lnGmf,  pmp->L);
              break;
      case 8: rddar.readArray( "RLC", pmp->RLC, pmp->L, 1 );
              break;
      case 9: rddar.readArray( "RSC", pmp->RSC, pmp->L, 1 );
              break;
      case 10: rddar.readArray( "DLL", pmp->DLL,  pmp->L);
              break;
      case 11: rddar.readArray( "DUL", pmp->DUL,  pmp->L);
              break;
      case 12: rddar.readArray( "Aalp", pmp->Aalp,  pmp->FI);
              break;
      case 13: if( !pmp->Sigw )
                Error( "Error", "Array Sigw not used in this problem");
              rddar.readArray( "Sigw", pmp->Sigw,  pmp->FI);
              break;
      case 14: if( !pmp->Sigg )
                Error( "Error", "Array Sigg not used in this problem");
              rddar.readArray( "Sigg", pmp->Sigg,  pmp->FI);
              break;
      case 15: rddar.readArray( "YOF", pmp->YOF,  pmp->FI);
              break;
      case 16: if( !pmp->Nfsp )
                Error( "Error", "Array Nfsp not used in this problem");
              rddar.readArray( "Nfsp", &pmp->Nfsp[0][0], pmp->FIs*pmp->FIat);
              break;
      case 17: if( !pmp->MASDT )
                Error( "Error", "Array MASDT not used in this problem");
              rddar.readArray( "MASDT", &pmp->MASDT[0][0], pmp->FIs*pmp->FIat);
              break;
      case 18: if( !pmp->XcapA )
                Error( "Error", "Array XcapA not used in this problem");
              rddar.readArray( "C1", &pmp->XcapA[0][0], pmp->FIs*pmp->FIat);
              break;
      case 19: if( !pmp->XcapB )
                Error( "Error", "Array XcapB not used in this problem");
              rddar.readArray( "C2", &pmp->XcapB[0][0], pmp->FIs*pmp->FIat);
              break;
      case 20: if( !pmp->XcapF )
                Error( "Error", "Array XcapF not used in this problem");
              rddar.readArray( "C3", &pmp->XcapF[0][0], pmp->FIs*pmp->FIat);
              break;
      case 21: if( !pmp->Xetaf )
                Error( "Error", "Array Xetaf not used in this problem");
              rddar.readArray( "pCh", &pmp->Xetaf[0][0], pmp->FIs*pmp->FIat);
              break;
      case 22: if( !pmp->SATX )
                Error( "Error", "Array SATX not used in this problem");
              rddar.readArray( "SATX", &pmp->SATX[0][0], pmp->Lads*4);
              break;
      case 23: if( !pmp->MASDJ )
                Error( "Error", "Array MASDJ not used in this problem");
              rddar.readArray( "MASDJ", &pmp->MASDJ[0][0], pmp->Lads*DFCN);
              break;
      case 24: if( !pmp->SCM )
                Error( "Error", "Array SCM not used in this problem");
              rddar.readArray( "SCM", pmp->SCM[0], pmp->FIs, pmp->FIat );
              break;
      case 25: if( !pmp->SATT )
                Error( "Error", "Array SATT not used in this problem");
              rddar.readArray( "SACT", pmp->SATT, pmp->Lads, 1 );
              break;
      case 26: if( !pmp->DCC3 )
                Error( "Error", "Array DCC3 not used in this problem");
               rddar.readArray( "DCads", pmp->DCC3, pmp->Lads, 1 );
               break;
      case 27: rddar.readArray( "pa_DB" , &pa->p.DB, 1);
               break;
      case 28: rddar.readArray("pa_DHB", &pa->p.DHB, 1);
               break;
      case 29: rddar.readArray("pa_EPS" , &pa->p.EPS, 1);
               break;
      case 30: rddar.readArray("pa_DK" , &pa->p.DK, 1);
               break;
      case 31: rddar.readArray("pa_DF" , &pa->p.DF, 1);
               break;
      case 32: rddar.readArray("pa_DP", &pa->p.DP, 1);
               break;
      case 33: rddar.readArray("pa_IIM", &pa->p.IIM, 1);
               break;
      case 34: rddar.readArray("pa_PD" , &pa->p.PD, 1);
               break;
      case 35: rddar.readArray("pa_PRD" , &pa->p.PRD, 1);
               break;
      case 36: rddar.readArray("pa_AG" , &pa->p.AG, 1);
               break;
      case 37: rddar.readArray("pa_DGC" , &pa->p.DGC, 1);
               break;
      case 38: rddar.readArray("pa_PSM" , &pa->p.PSM, 1);
               break;
      case 39: rddar.readArray("pa_GAR" , &pa->p.GAR, 1);
               break;
      case 40: rddar.readArray("pa_GAH" , &pa->p.GAH, 1);
               break;
      case 41: rddar.readArray("pa_DS", &pa->p.DS, 1);
               break;
      case 42: rddar.readArray("pa_XwMin" , &pa->p.XwMin, 1);
               break;
      case 43: rddar.readArray("pa_ScMin" , &pa->p.ScMin, 1);
               break;
      case 44: rddar.readArray("pa_DcMin" , &pa->p.DcMin, 1);
               break;
      case 45: rddar.readArray("pa_PhMin" , &pa->p.PhMin, 1);
               break;
      case 46: rddar.readArray("pa_ICmin" , &pa->p.ICmin, 1);
               break;
      case 47: rddar.readArray("pa_PC" , &pa->p.PC, 1);
               break;
      case 48: rddar.readArray("pa_DFM" , &pa->p.DFM, 1);
               break;
      case 49: rddar.readArray("pa_DFYw" , &pa->p.DFYw, 1);
               break;
      case 50: rddar.readArray("pa_DFYaq" , &pa->p.DFYaq, 1);
               break;
      case 51: rddar.readArray("pa_DFYid" , &pa->p.DFYid, 1);
               break;
      case 52: rddar.readArray("pa_DFYr" , &pa->p.DFYr, 1);
               break;
      case 53: rddar.readArray("pa_DFYh" , &pa->p.DFYh, 1);
               break;
      case 54: rddar.readArray("pa_DFYc" , &pa->p.DFYc, 1);
               break;
      case 55: rddar.readArray("pa_DFYs", &pa->p.DFYs, 1);
               break;
      case 56: rddar.readArray("pa_DW", &pa->p.DW , 1);
               break;
      case 57: rddar.readArray("pa_DT", &pa->p.DT , 1);
               break;
      case 58: rddar.readArray("pa_GAS", &pa->p.GAS, 1);
               break;
      case 59: rddar.readArray("pa_DNS" , &pa->p.DNS, 1);
               break;
      case 60: rddar.readArray("pa_IEPS" , &pa->p.IEPS, 1);
               break;
      case 61: rddar.readArray("pKin" , &pmp->PLIM, 1);
               break;
      case 62: rddar.readArray("pa_DKIN" , &pa->p.DKIN, 1);
               break;
      case 63: rddar.readArray("mui" , pmp->mui, pmp->N);
               break;
      case 64: rddar.readArray("muk" , pmp->muk, pmp->FI);
               break;
      case 65: rddar.readArray("muj" , pmp->muj, pmp->L);
               break;
    }
    nfild = rddar.findNext();
  }
/*
   if( pm.sitNcat*pm.sitNcat )
     rddar.readArray( "sitE", pmp->sitE, pmp->sitNcat*pmp->sitNan );
   if( pm.sitNcat )
     rddar.readArray( "sitXc", pmp->sitXcat, pmp->sitNcat );
   if( pm.sitNan )
     rddar.readArray( "sitXa", pmp->sitXan, pmp->sitNan );
*/
 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from the MULTY structure";
    Error( "Error", ret);
  }
}

//=============================================================================
// ms_multi_format.cpp

