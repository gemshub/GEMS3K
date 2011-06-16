//-------------------------------------------------------------------
// $Id: ms_multi_format.cpp 774 2006-07-26 08:45:45Z gems $
//
// Implementation of text writing/reading IPM, DCH and DBR files
//
// Copyright (C) 2006,2009 S.Dmytriyeva,D.Kulik
//
// This file is part of the GEM-Vizor library and GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
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

outField MULTI_dynamic_fields[69] =  {
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
   { "pa_DG", 0 , 0 },
   { "pa_DNS" , 0 , 0 },
   { "pa_IEPS" , 0 , 0 },
   { "pKin" , 0 , 0 },
   { "pa_DKIN" , 0 , 0 },
   { "mui" , 0 , 0 },
   { "muk" , 0 , 0 },
   { "muj" , 0 , 0 },
   { "pa_PLLG" , 0 , 0 },
   { "tMin" , 0 , 0 }
};


//===================================================================

void TMulti::to_text_file_gemipm( const char *path, bool addMui,
		bool with_comments, bool brief_mode )
{
  SPP_SETTING *pa = &TProfil::pm->pa;
   _comment = with_comments;
   char PAalp;
   char PSigm;

#ifndef IPMGEMPLUGIN
   PAalp = syp->PAalp;
   PSigm = syp->PSigm;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
#endif
  fstream ff( path, ios::out );
  ErrorIf( !ff.good() , path, "Fileopen error");
  TPrintArrays  prar( 69, MULTI_dynamic_fields, ff);
// set up array flags for permanent fields
   if( !( pm.FIs > 0 && pm.Ls > 0 ) )
   {
     prar.setNoAlws( (long int)(0) /*"sMod"*/);
     prar.setNoAlws( 1 /*"LsMod"*/);
     prar.setNoAlws( 2 /*"LsMdc"*/);
   }
   if( PSigm == S_OFF )
   {
     prar.setNoAlws( 13 /*"Sigw"*/);
     prar.setNoAlws( 14 /*"Sigg"*/);
   }
   if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
   { /* ADSORPTION AND ION EXCHANGE */
     prar.setNoAlws( 16 /*"Nfsp"*/);
     prar.setNoAlws( 17 /*"MASDT"*/);
     prar.setNoAlws( 18 /*"C1"*/);
     prar.setNoAlws( 19 /*"C2"*/);
     prar.setNoAlws( 20 /*"C3"*/);
     prar.setNoAlws( 21 /*"pCh"*/);
     prar.setNoAlws( 22 /*"SATX"*/);
     prar.setNoAlws( 23 /*"MASDJ"*/);
     prar.setNoAlws( 24 /*"SCM"*/);
     prar.setNoAlws( 25 /*"SACT"*/);
     prar.setNoAlws( 26 /*"DCads"*/);
   }

if( _comment )
{   ff << "# GEMIPM2K v. 3.0 rev. 530(1848)" << endl;
   ff << "# Comments can be marked with # $ ; as the first character in the line" << endl << endl;
   ff << "# Template for the ipm-dat text input file for the internal MULTI data" << endl;
   ff << "# (should be read after the DATACH file and before DATABR files)" << endl << endl;
   ff << "# ID key of the initial chemical system definition" << endl;
}
  ff << "\"" << pmp->stkey << "\"" << endl;

 if( _comment )
     ff << "\n## (1) Important flags that affect memory allocation" << endl;

 if(!brief_mode || pa->p.PE != pa_.p.PE )
 { if( _comment )
      ff << "# PE: Flag for using electroneutrality condition in GEM IPM calculations " << endl;
   ff << left << setw(12) << "<pa_PE> " <<  right << setw(8) << pa->p.PE << endl;
 }

//   ff << "# 'E'                1" << endl;
//   ff << left << setw(12) << "<E> " <<  right << setw(8) << pmp->E << endl;
 if( _comment )
     ff << "\n# PV: Flag for the volume balance constraint (on Vol IC) - for indifferent equilibria at P_Sat" << endl;
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
//   ff << "# 'GWAT'         55.50837344" << endl;
//   ff << left << setw(12) << "<GWAT> " <<  right << setw(8) << pmp->GWAT << endl;
   if( _comment )
     ff << "\n# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases " << endl;
   ff << left << setw(12) << "<PAalp> " <<  right << setw(6) <<
      "\'" << PAalp << "\'" << endl;
   if( _comment )
    ff << "\n# PSigm: Flag for using (+) or ignoring (-) specific surface free energies  " << endl;
   ff << left << setw(12) << "<PSigm> " <<  right << setw(6) <<
      "\'" << PSigm << "\'" << endl;
  if( !brief_mode || pmp->FIat > 0 || pmp->Lads > 0 )
  { if( _comment )
    {  ff << "\n## (2) Important dimensionalities that affect memory allocation" << endl;
       ff << "# Lads: Total number of Dependent Components in sorption phases included into this system" << endl;
    }
   ff << left << setw(12) << "<Lads> " <<  right << setw(8) << pmp->Lads << endl;
   if( _comment )
     ff << "# FIa: Number of sorption phases included in this system (0 if no sorption phases are included)" << endl;
   ff << left << setw(12) << "<FIa> " <<  right << setw(8) << pmp->FIa << endl;
   if( _comment )
      ff << "# FIat: Maximum number of surface types per adsorption phase (if FIa > 0, FIat must be set to default value of 6)" << endl;
   ff << left << setw(12) << "<FIat> " <<  right << setw(8) << pmp->FIat << endl << endl;
//   ff << left << setw(12) << "<sitNc> " <<  right << setw(8) << pmp->sitNcat << endl;
//   ff << left << setw(12) << "<sitNa> " <<  right << setw(8) << pmp->sitNan << endl;
    } // brief_mode
   ff << "\n<END_DIM>\n";

// static data not affected by dimensionalities
   if( _comment )
   {  ff << "\n## (3) Tolerances and controls of the numerical behavior of GEM IPM-2 kernel" << endl;
      ff << "#      - Need to be changed only in rare special cases (see gems_ipm.html)" << endl;
   }
   if(!brief_mode || pa->p.DB != pa_.p.DB )
   { if( _comment )
       ff << "# DB: Minimum amount of Independent Component in the bulk system composition (except charge Zz) (moles)" << endl;
     ff << left << setw(12) << "<pa_DB> " <<  right << setw(8) << pa->p.DB << endl;
   }
   if(!brief_mode || pa->p.DHB != pa_.p.DHB )
   { if( _comment )
      ff << "\n# DHB: Maximum allowed relative mass balance residual for ICs ( 1e-9 to 1e-15 ) { 1e-12 } " << endl;
     ff << left << setw(12) << "<pa_DHB> " << right << setw(8) <<  pa->p.DHB << endl;
   }
   if(!brief_mode || pa->p.EPS != pa_.p.EPS )
   { if( _comment )
      ff << "\n# EPS: Tolerance of the SolveSimplex() convergence (1e-6 to 1e-14) { 1e-10 }" << endl;
    ff << left << setw(12) << "<pa_EPS> " <<  right << setw(8) << pa->p.EPS << endl;
   }
   if(!brief_mode || pa->p.DK != pa_.p.DK )
   { if( _comment )
      ff << "\n# DK: Tolerance threshold for the Dikin's criterion of IPM convergence (1e-6 to 1e-4) { 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_DK> " <<  right << setw(8) << pa->p.DK << endl;
   }
   if(!brief_mode || pa->p.DS != pa_.p.DS )
   { if( _comment )
      ff << "\n# DS: Cutoff min. amount of stable phase in GEM IPM primal solution (1e-8 to 1e-24) { 1e-20 }" << endl;
     ff << left << setw(12) << "<pa_DS> " << right << setw(8) <<  pa->p.DS << endl;
   }
   if(!brief_mode || pa->p.DF != pa_.p.DF )
   { if( _comment )
      ff << "\n# DF: Tolerance DF for Karpov's criterion (Fa > DF) for a lost stable phase to be inserted { 0.01 }" << endl;
     ff << left << setw(12) << "<pa_DF> " <<  right << setw(8) << pa->p.DF << endl;
   }
   if(!brief_mode || pa->p.DFM != pa_.p.DFM )
   { if( _comment )
       ff << "# DFM: Tolerance for Karpov's criterion (Fa < -DFM) for a present unstable phase to be eliminated { 0.01 } " << endl;
     ff << left << setw(12) << "<pa_DFM> " <<  right << setw(8) << pa->p.DFM << endl;
   }
   if(!brief_mode || pa->p.DP != pa_.p.DP )
   {  if( _comment )
       ff << "\n# DP: Maximal number of iterations in MassBalanceRefinement() procedure (20 to 130) { 90 }" << endl;
     ff << left << setw(12) << "<pa_DP> " << right << setw(8) << pa->p.DP << endl;
   }
   if(!brief_mode || pa->p.IIM != pa_.p.IIM )
   { if( _comment )
      ff << "\n# Maximum allowed number of iterations in one GEM IPM descent run (100 to 9999) { 7000 }" << endl;
     ff << left << setw(12) << "<pa_IIM> " << right << setw(8) <<  pa->p.IIM << endl;
   }
   if(!brief_mode || pa->p.PD != pa_.p.PD )
   { if( _comment )
       ff << "\n# PD: Mode of calculation of activity coefficients ( 1 -IPM, 2 +EFD, 3 IPM ) { 2 } " << endl;
     ff << left << setw(12) << "<pa_PD> " <<  right << setw(8) << pa->p.PD << endl;
   }
   if(!brief_mode || pa->p.PRD != pa_.p.PRD )
   { if( _comment )
       ff << "\n# PRD: Disable (0) or activate (-4 or less- max.dec.exp.for DC amount correction) SpeciationCleanup() { -4 }" << endl;
     ff << left << setw(12) << "<pa_PRD> " <<  right << setw(8) << pa->p.PRD << endl;
   }
   if(!brief_mode || pa->p.AG != pa_.p.AG )
   { if( _comment )
      ff << "\n# AG: Smoothing parameter 1 for non-ideal primal chemical potential increments -1 to +1 { 1.0 }" << endl;
     ff << left << setw(12) << "<pa_AG> " <<  right << setw(8) << pa->p.AG << endl;
   }
   if(!brief_mode || pa->p.DGC != pa_.p.DGC )
   { if( _comment )
      ff << "\n# DGC: Smoothing parameter 2 (exponent in smoothing function (-1 to +1) { -0.98 or 0.001 for adsorption }" << endl;
     ff << left << setw(12) << "<pa_DGC> " <<  right << setw(8) << pa->p.DGC << endl;
   }
   if(!brief_mode || pa->p.PSM != pa_.p.PSM )
   { if( _comment )
      ff << "\n# PSM: Level of diagnostic messages { 0- disabled (no ipmlog file); 1- normal; 2-including warnings }" << endl;
     ff << left << setw(12) << "<pa_PSM> " <<  right << setw(8) << pa->p.PSM << endl;
   }
   if(!brief_mode || pa->p.GAR != pa_.p.GAR )
   { if( _comment )
      ff << "# GAR: Initial activity coefficient value for major (M) species in a solution phase at LPP AIA { 1 }" << endl;
     ff << left << setw(12) << "<pa_GAR> " <<  right << setw(8) << pa->p.GAR << endl;
   }
   if(!brief_mode || pa->p.GAH != pa_.p.GAH )
   { if( _comment )
      ff << "# GAH: Initial activity coefficient value for minor (J) species in a solution phase at LPP AIA { 1000 }" << endl;
     ff << left << setw(12) << "<pa_GAH> " <<  right << setw(8) << pa->p.GAH << endl;
   }
   if(!brief_mode)
    if( _comment )
     {  ff << "\n# _Min: Cutoff amounts for elimination of: Xw - water-solvent { 1e-11 }; Sc - solid sorbent {1e-11}; " << endl;
        ff <<   "#       Dc - solution- or surface species { 1e-30 }; Ph - non-electrolyte solution phase with all its components { 1e-20 }" << endl;
     }
   if(!brief_mode || pa->p.XwMin != pa_.p.XwMin )
    ff << left << setw(12) << "<pa_XwMin> " <<  right << setw(8) << pa->p.XwMin << endl;
   if(!brief_mode || pa->p.ScMin != pa_.p.ScMin )
     ff << left << setw(12) << "<pa_ScMin> " <<  right << setw(8) << pa->p.ScMin << endl;
   if(!brief_mode || pa->p.DcMin != pa_.p.DcMin )
     ff << left << setw(12) << "<pa_DcMin> " <<  right << setw(8) << pa->p.DcMin << endl;
   if(!brief_mode || pa->p.PhMin != pa_.p.PhMin )
     ff << left << setw(12) << "<pa_PhMin> " <<  right << setw(8) << pa->p.PhMin << endl;

   if(!brief_mode || pa->p.ICmin != pa_.p.ICmin )
   { if( _comment )
      ff << "\n# ICmin: Cutoff value of effective molal ionic strength to disable aq-gamma calculation (1e-6 to 1e-3) { 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_ICmin> " <<  right << setw(8) << pa->p.ICmin << endl;
   }
   if(!brief_mode || pa->p.PC != pa_.p.PC )
   { if( _comment )
      ff << "\n# PC: Mode of Phase Selection: 1 old (Select-2), 2 new (PSSC)  { 2 }" << endl;
     ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   }
   if( _comment )
      ff << "# DFY: Insertion mole amounts used after the LPP AIA and in PhaseSelection() algorithm" << endl;
   if(!brief_mode || pa->p.DFYw != pa_.p.DFYw )
   { if( _comment )
      ff << "# DFYw: Insertion mole amount for water-solvent at Simplex()->MassBalanceRefinement() bridge { 1e-5 }" << endl;
      ff << left << setw(12) << "<pa_DFYw> " <<  right << setw(8) << pa->p.DFYw << endl;
   }
   if(!brief_mode || pa->p.DFYaq != pa_.p.DFYaq )
   { if( _comment )
      ff << "# DFYaq: Insertion mole amount for aqueous species at Simplex()->MBR() bridge { 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_DFYaq> " <<  right << setw(8) << pa->p.DFYaq << endl;
   }
   if(!brief_mode || pa->p.DFYid != pa_.p.DFYid )
   { if( _comment )
      ff << "\n# DFYid: Insertion mole amount for species of ideal solutions at Simplex()->MBR() bridge { 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_DFYid> " <<  right << setw(8) << pa->p.DFYid << endl;
   }
   if(!brief_mode || pa->p.DFYr != pa_.p.DFYr )
   { if( _comment )
      ff << "# DFYr: Insertion mole amount for a major species in a solution at Simplex()->MBR()bridge { 1e-5 }" << endl;
    ff << left << setw(12) << "<pa_DFYr> " <<  right << setw(8) << pa->p.DFYr << endl;
   }
   if(!brief_mode || pa->p.DFYh != pa_.p.DFYh )
   { if( _comment )
      ff << "# DFYh: Insertion mole amount for a junior species in a solution at Simplex()->MBR() bridge{ 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_DFYh> " <<  right << setw(8) << pa->p.DFYh << endl;
   }
   if(!brief_mode || pa->p.DFYc != pa_.p.DFYc )
   { if( _comment )
      ff << "# DFYc:  Insertion mole amount for a single-component phase at Simplex()->MBR() bridge { 1e-5 }" << endl;
     ff << left << setw(12) << "<pa_DFYc> " <<  right << setw(8) << pa->p.DFYc << endl;
   }
   if(!brief_mode || pa->p.DFYs != pa_.p.DFYs )
   { if( _comment )
      ff << "# DFYs: Insertion mole amount for a single-component phase in PSSC() { 1e-6 }" << endl;
     ff << left << setw(12) << "<pa_DFYs> " << right << setw(8) <<  pa->p.DFYs << endl;
   }
   if( _comment )
     ff << "\n# Parameters for high-accuracy IPM algorithm " << endl;
   if(!brief_mode || pa->p.DW != pa_.p.DW )
   {    if( _comment )
      ff << "# DW: Activate (1) or disable (0) the error condition on DP - maximum allowed number of MBR() iterations { 1 }" << endl;
   ff << left << setw(12) << "<pa_DW> " << right << setw(8) << pa->p.DW  << endl;
   }
   if(!brief_mode || pa->p.DT != pa_.p.DT )
   { if( _comment )
      ff << "# DT: DHB is rel.max.MB cutoff for all ICs (0) or for major ICs: dec.exponent (<-6) of abs.MB cutoff; (1) for DHB also as abs.cutoff { 1 }" << endl;
     ff << left << setw(12) << "<pa_DT> " << right << setw(8) << pa->p.DT  << endl;
   }
   if(!brief_mode || pa->p.GAS != pa_.p.GAS )
   { if( _comment )
       ff << "\n# GAS: Threshold for primal-dual norm.chem.pot.difference used in SpeciationCleanup() { 0.0001 }" << endl;
     ff << left << setw(12) << "<pa_GAS> " << right << setw(8) <<  pa->p.GAS << endl;
   }
   if(!brief_mode || pa->p.DG != pa_.p.DG )
   { if( _comment )
          ff << "# Total number of moles used in internal re-scaling of the system (disabled if < 1e-4) { 1e3 }" << endl;
     ff << left << setw(12) << "<pa_DG> " <<  right << setw(8) << pa->p.DG << endl;
   }
   if(!brief_mode || pa->p.DNS != pa_.p.DNS )
   { if( _comment )
       ff << "# DNS: Standard surface density (nm-2) for calculating activity of surface species { 12.05 nm-2 }" << endl;
     ff << left << setw(12) << "<pa_DNS> " <<  right << setw(8) << pa->p.DNS << endl;
   }
   if(!brief_mode || pa->p.IEPS != pa_.p.IEPS )
   { if( _comment )
       ff << "# IEPS: Tolerance for calculation of surface activity coefficient terms for surface species { 1e-3 }" << endl;
     ff << left << setw(12) << "<pa_IEPS> " <<  right << setw(8) << pa->p.IEPS << endl;
   }
   if(!brief_mode )
   { if( _comment )
       ff << "\n# pKin:Flag for using metastability constraints on calculated amounts of Dependent Components { 1 } " << endl;
     ff << left << setw(12) << "<pKin> " <<  right << setw(8) << pmp->PLIM << endl;
   }
   if(!brief_mode || pa->p.DKIN != pa_.p.DKIN )
   { if( _comment )
      ff << "# DKIN: Tolerance for non-trivial metastability restrictions on amounts of dependent components, moles { 1e-8 } " << endl;
     ff << left << setw(12) << "<pa_DKIN> " <<  right << setw(8) << pa->p.DKIN << endl;
   }
   if(!brief_mode || pa->p.PLLG != pa_.p.PLLG )
   { if( _comment )
       ff << "# pa_PLLG: Tolerance for checking changes in dual solution after PhaseSelect(), 1 to 100 { 10 }" << endl;
     ff << left << setw(12) << "<pa_PLLG> " <<  right << setw(8) << pa->p.PLLG << endl;
   }

   if(!brief_mode || pmp->tMin != G_TP_ )
   { if( _comment )
       ff << "# tMin: Type of thermodynamic potential to minimize" << endl;
     ff << left << setw(12) << "<tMin> " <<  right << setw(8) << pmp->tMin << endl;
   }

//dynamic arrays
if( pm.FIs > 0 && pm.Ls > 0 )
{
  if( _comment )
    ff << "\n## (4) Initial data for multicomponent phases (see DATACH file for dimension nPHs)" << endl;
  if( _comment )
     ff << "# sMod: Codes for built-in mixing models of multicomponent phases [nPS*6]";
    prar.writeArray(  "sMod", pmp->sMod[0], pmp->FIs, 6L );

long int LsModSum;
long int LsIPxSum;
long int LsMdcSum;
getLsModsum( LsModSum, LsIPxSum );
getLsMdcsum( LsMdcSum );

   if( _comment )
    {  ff << "\n\n# LsMod: Dimensions of <IPxPH> and <PMc> arrays [nPS*3]" << endl;

    }
    prar.writeArray(  "LsMod", pmp->LsMod, pmp->FIs*3, 3L);

  if(LsIPxSum )
  {
   if( _comment )
      ff << "\n\n# IPxPH:  Collected indexation table for interaction parameters of non-ideal solutions.";
   prar.writeArray(  "IPxPH", pmp->IPx,  LsIPxSum);
  }
  if(LsModSum )
   {
     if( _comment )
        ff << "\n\n# PMc: Collected interaction parameter coefficients for the (built-in) non-ideal mixing models";
    prar.writeArray(  "PMc", pmp->PMc,  LsModSum);
   }
   if( _comment )
     ff << "\n\n# LsMdc: Number of parameters per component of the phase for the non-ideal mixing models [nPS]";
   prar.writeArray(  "LsMdc", pmp->LsMdc, pmp->FIs);
   if(LsMdcSum )
   {   if( _comment )
          ff << "\n\n# DMc: Collected parameters per phase component for the non-ideal mixing models ";
    prar.writeArray(  "DMc", pmp->DMc,  LsMdcSum);
   }
} // sMod
  if( _comment )
  {  ff << "\n\n## (5) Some data arrays which are not provided in DATACH and DATABR files" << endl;
     ff << "# B: Full total bulk composition of the initial system (vector b)  (moles) [nIC]";
  }
  prar.writeArray(  "B", pmp->B,  pmp->N);

  if( _comment )
     ff << "\n\n# Initial data for DCs - see DATACH file for dimensions nDC, nDCs" << endl;

  if(!brief_mode || prar.getAlws("Pparc" ))
  {
	  if( _comment )
	    ff << "# Pparc: Partial pressures or fugacities of pure Dependent Components [nDC]";
      prar.writeArray(  "Pparc", pmp->Pparc,  pmp->L);
  }
 //  ff << "\n\n# This is not necessary - can be calculated from G0 ???????????";
 // prar.writeArray(  "G0", pmp->G0,  pmp->L);
  if(!brief_mode || prar.getAlws("GEX" ))
  {
   if( _comment )
      ff << "\n\n# fDQF: DQF parameters or pure gas fugacities in (J/mol/(RT) [nDC]";
   prar.writeArray(  "fDQF", pmp->fDQF,  pmp->L);
  }
  if(!brief_mode || prar.getAlws("lnGmf" ))
  { if( _comment )
      ff << "\n\n# lnGmf:  Natural logarithms of DC (activity coefficients) to be used for correcting g0(T,P) [nDC]";
    prar.writeArray(  "lnGmf", pmp->lnGmf,  pmp->L);
  }
  if( _comment )
     ff << "\n\n# (6) Section for metastability/ kinetic constraints" << endl;

   if(!brief_mode || prar.getAlws("RLC" ))
   {   if( _comment )
         ff << "# RLC: Code of metastability constraints for DCs [nDC]";
     prar.writeArray(  "RLC", pmp->RLC, pmp->L, 1L );
   }
   if(!brief_mode || prar.getAlws("RSC" ))
   { if( _comment )
      ff << "\n\n# RSC: Units of metastability/kinetic constraints for DCs [nDC]";
     prar.writeArray(  "RSC", pmp->RSC, pmp->L, 1L );
   }
   if(!brief_mode || prar.getAlws("DLL" ))
   { if( _comment )
        ff << "\n\n# DLL: Full vector of lower metastability constraints on DC amounts <xDC> in the system (moles) [nDC]";
     prar.writeArray(  "DLL", pmp->DLL,  pmp->L);
   }
   if(!brief_mode || prar.getAlws("DUL" ))
   {  if( _comment )
       ff << "\n\n# DUL:Full vector of upper metastability constraints on DC amounts <xDC> in the system (moles) [nDC]";
      prar.writeArray(  "DUL", pmp->DUL,  pmp->L);
   }
   if( _comment )
     ff << "\n\n# (7) Initial data for phases" << endl;
   if(!brief_mode || prar.getAlws("Aalp" ))
   { if( _comment )
      ff << "\n# Aalp: Full vector of specific surface areas of phases (m2/g) [nPH]";
     prar.writeArray(  "Aalp", pmp->Aalp,  pmp->FI);
   }
   if( PSigm != S_OFF )
   {
     if(!brief_mode || prar.getAlws("Sigw" ))
     { if( _comment )
         ff << "\n\n# Sigw: Specific surface free energy for phase-water interface (J/m2) [nPH]";
        prar.writeArray(  "Sigw", pmp->Sigw,  pmp->FI);
     }
     if(!brief_mode || prar.getAlws("Sigg" ))
     { if( _comment )
         ff << "\n\n# Sigg: Specific surface free energy for phase-gas interface (J/m2) (not yet used) [nPH]";
       prar.writeArray(  "Sigg", pmp->Sigg,  pmp->FI);
     }
   }
   if(!brief_mode || prar.getAlws("YOF" ))
   {  if( _comment )
        ff << "\n\n# YOF: Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition)  [nPH]";
      prar.writeArray(  "YOF", pmp->YOF,  pmp->FI);
   }
   if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
      if( _comment )
        ff << "\n\n# (8) Initial data for sorption" << endl;
      if(!brief_mode || prar.getAlws("Nfsp" ))
      { if( _comment )
         ff << "\n# Nfsp: Fractions of the sorbent specific surface area allocated to surface types [nPS*6]";
        prar.writeArray(  "Nfsp", &pmp->Nfsp[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("MASDT" ))
      { if( _comment )
         ff << "\n# MASDT: Total maximum site  density per surface type (mkmol/g) [nPS*6]";
        prar.writeArray(  "MASDT", &pmp->MASDT[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("C1" ))
      { if( _comment )
           ff << "\n# C1: Inner capacitance density parameter C1 (F/m2) [nPS*6]";
        prar.writeArray(  "C1", &pmp->XcapA[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("C2" ))
      { if( _comment )
          ff << "\n# C2: Outer capacitance density parameter C2 (F/m2) [nPS*6]";
        prar.writeArray(  "C2", &pmp->XcapB[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("C3" ))
      { if( _comment )
          ff << "\n# C3: Third capacitance density parameter C3  (F/m2) [nPS*6]";
        prar.writeArray(  "C3", &pmp->XcapF[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("pCh" ))
      { if( _comment )
          ff << "\n# pCh: Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [nPS*6]";
       prar.writeArray(  "pCh", &pmp->Xetaf[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      }
      if(!brief_mode || prar.getAlws("SATX" ))
      { if( _comment )
        {   ff << "\n# SATX: Setup of surface sites and species (will be applied separately within each sorption phase) [Lads*4]";
            ff << "\n# [0] surface type; [1] sorbent emd member;";
            ff << "\n# [2] surface site in surf. type; [3] surface EDL plane";
        }
       prar.writeArray(  "SATX", &pmp->SATX[0][0], pmp->Lads*4, 4L);
     }
     if(!brief_mode || prar.getAlws("MASDJ" ))
     { if( _comment )
       {  ff << "\n# MASDJ: Parameters of surface species in surface complexation models [Lads*6]";
          ff << "\n# [0] max site density mmol/g; [1] charge allocated to 0 plane;";
          ff << "\n# [2] charge allocated to beta -or third plane; [3] Frumkin interaction parameter;";
          ff << "\n# [4] dentateness or CN; [5] reserved isoterm parameter.";
       }
       prar.writeArray(  "MASDJ", &pmp->MASDJ[0][0], pmp->Lads*DFCN, (long int)DFCN);
     }
     if(!brief_mode || prar.getAlws("SCM" ))
     { if( _comment )
         ff << "\n# SCM: Classifier of built-in electrostatic models applied to surface types in sorption phases [nPS*6]";
       prar.writeArray(  "SCM", pmp->SCM[0], pmp->FIs, pmp->FIat );
     }
     if(!brief_mode || prar.getAlws("SACT" ))
     { if( _comment )
         ff << "\n# SACT: Classifier of applied SACT equations (isotherm corrections) [Lads].";
       prar.writeArray(  "SACT", pmp->SATT, pmp->Lads, 1L );
     }
     if(!brief_mode || prar.getAlws("DCads" ))
     { if( _comment )
         ff << "\n# DCads: Classifier of DCs involved in sorption phases [Lads]";
      prar.writeArray(  "DCads", pmp->DCC3, pmp->Lads, 1L );
     }
    }
/*
 * outArray( ff, "Vol", pmp->Vol,  pmp->L);
   outArray( ff, "G0", pmp->G0,  pmp->L);
   outArray( ff, "PUL", pmp->PUL,  pmp->L);
   outArray( ff, "PLL", pmp->PLL,  pmp->L);
   outArray( ff, "lnGam", pmp->lnGam,  pmp->L);
   outArray( ff, "F0", pmp->F0,  pmp->L);

   if( pm.sitNcat*pm.sitNcat )
    prar.writeArray(  "sitE", pmp->sitE, pmp->sitNcat*pmp->sitNan );
   if( pm.sitNcat )
    prar.writeArray(  "sitXc", pmp->sitXcat, pmp->sitNcat );
   if( pm.sitNan )
     prar.writeArray(  "sitXa", pmp->sitXan, pmp->sitNan );
*/

 if( addMui && !brief_mode )
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
 ff << endl;
 if( _comment )
   ff << "\n# End of file" << endl;

}

void TMulti::from_text_file_gemipm( const char *path )
{
  SPP_SETTING *pa = &TProfil::pm->pa;
  DATACH  *dCH = TNode::na->pCSD();
  long int ii, nfild;

   //static values
   char PAalp;
   char PSigm;

#ifdef IPMGEMPLUGIN
   set_def();
#endif
  //mem_set( &pm.N, 0, 38*sizeof(long int));
  //mem_set( &pm.TC, 0, 55*sizeof(double));
  // get sizes from DATACH
  pmp->TC = pmp->TCc = 25.;
  pmp->T = pmp->Tc =298.15;
  pmp->P = pmp->Pc = 1.;
  pmp->N = pmp->NR = dCH->nIC;
  pmp->L = dCH->nDC;
  pmp->FI = dCH->nPH;
  pmp->FIs = dCH->nPS;
  //
  pmp->Ls = 0; //dCH->nDCs;
  for( ii=0; ii<dCH->nPS; ii++)
  {
    pmp->Ls += dCH->nDCinPH[ii];
    if( dCH->ccPH[ii] == 'a' )
     pmp->LO = pmp->Ls-1;
    if( dCH->ccPH[ii] == 'g' || dCH->ccPH[ii] == 'p' || dCH->ccPH[ii] == 'f')
      pmp->PG = dCH->nDCinPH[ii];
  }

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
   copyValues( pmp->stkey, (char * )str.c_str(), EQ_RKLEN );

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

#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#else
   dyn_new();
#endif

// get dynamic data from DATACH file
  for( ii=0; ii<dCH->nPH; ii++)
    pmp->L1[ii] = dCH->nDCinPH[ii];

  for( ii=0; ii<dCH->nIC*dCH->nDC; ii++)
    pmp->A[ii] = dCH->A[ii];

  if( pmp->EZ )
  { long int iZ=-1;
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
  { pmp->Awt[ii]  = dCH->ICmm[ii]*1e3;
    fillValue(pmp->SB[ii], ' ', MaxICN );
    copyValues( pmp->SB[ii], dCH->ICNL[ii], min(MaxICN,(long int)MAXICNAME));
    pmp->SB[ii][MaxICN] = dCH->ccIC[ii];
    pmp->ICC[ii] =  dCH->ccIC[ii];
  }

if( fabs(dCH->DCmm[0]) < 1e-32 )  // Restore DCmm if skipped from the DCH file
  for( long int jj=0; jj< dCH->nDC; jj++ )  // Added by DK on 03.03.2007
  {
    dCH->DCmm[jj] = 0.0;
    for( ii=0; ii< dCH->nIC; ii++ )
       dCH->DCmm[jj] += dCH->ICmm[ii]*dCH->A[jj*dCH->nIC+ii];
  }

  for( ii=0; ii< dCH->nDC; ii++ )
  {
    pmp->MM[ii] = dCH->DCmm[ii]*1e3;
    pmp->DCC[ii] = dCH->ccDC[ii];
    copyValues( pmp->SM[ii], dCH->DCNL[ii], min(MaxDCN,(long int)MAXDCNAME) );
  }

  for( ii=0; ii< dCH->nPH; ii++ )
  {
	  fillValue( pmp->SF[ii], ' ', MAXPHNAME+MAXSYMB );
	  copyValues( pmp->SF[ii]+MAXSYMB, dCH->PHNL[ii], min(MaxPHN,(long int)MAXPHNAME) );
     pmp->SF[ii][0] = dCH->ccPH[ii];
     pmp->PHC[ii] = dCH->ccPH[ii];
  }

// !!!!  copyValues( pmp->DCCW, dCH->ccDCW, dCH->nDC);
  // set up DCCW
  ConvertDCC();

//reads dynamic values from txt file
   TReadArrays  rddar( 69, MULTI_dynamic_fields, ff);

// set up array flags for permanent fields

 if( !( pm.FIs > 0 && pm.Ls > 0 ) )
 {
   rddar.setNoAlws( (long int)(0) /*"sMod"*/);
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
              long int LsModSum;
              long int LsIPxSum;
              getLsModsum( LsModSum, LsIPxSum );
              if(LsIPxSum )
              { rddar.readNext( "IPxPH");
#ifdef IPMGEMPLUGIN
              if(!pmp->IPx )
                  pm.IPx = new long int[LsIPxSum];
#else
                 pm.IPx = (long int *)aObj[ o_wi_ipxpm ].Alloc(LsIPxSum, 1, L_);
#endif
                rddar.readArray( "IPxPH", pmp->IPx,  LsIPxSum);
              }
              if(LsModSum )
              { rddar.readNext( "PMc");
#ifdef IPMGEMPLUGIN
              if(!pmp->PMc )
                  pm.PMc = new double[LsModSum];
#else
               pm.PMc = (double *)aObj[ o_wi_pmc].Alloc( LsModSum, 1, D_);
#endif
                rddar.readArray( "PMc", pmp->PMc,  LsModSum);
              }
              break;
             }
      case 2: { if( !pmp->LsMdc )
                   Error( "Error", "Array LsMdc not used in this problem");
                rddar.readArray( "LsMdc" , pmp->LsMdc, pmp->FIs );
                long int LsMdcSum;
                getLsMdcsum( LsMdcSum );
                if(LsMdcSum )
                { rddar.readNext( "DMc");
#ifdef IPMGEMPLUGIN
                if(!pmp->DMc )
                     pm.DMc = new double[LsMdcSum];
#else
                pm.DMc = (double *)aObj[ o_wi_dmc].Alloc( LsMdcSum, 1, D_ );
#endif
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
      case 6: rddar.readArray( "fDQF", pmp->fDQF,  pmp->L);
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
      case 59: rddar.readArray("pa_DG" , &pa->p.DG, 1);
               break;
      case 60: rddar.readArray("pa_DNS" , &pa->p.DNS, 1);
               break;
      case 61: rddar.readArray("pa_IEPS" , &pa->p.IEPS, 1);
               break;
      case 62: rddar.readArray("pKin" , &pmp->PLIM, 1);
               break;
      case 63: rddar.readArray("pa_DKIN" , &pa->p.DKIN, 1);
               break;
      case 64: rddar.readArray("mui" , pmp->mui, pmp->N);
               break;
      case 65: rddar.readArray("muk" , pmp->muk, pmp->FI);
               break;
      case 66: rddar.readArray("muj" , pmp->muj, pmp->L);
               break;
      case 67: rddar.readArray("pa_PLLG" , &pa->p.PLLG, 1);
               break;
      case 68: rddar.readArray("tMin" , &pmp->tMin, 1);
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

