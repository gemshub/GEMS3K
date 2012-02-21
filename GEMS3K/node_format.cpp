//--------------------------------------------------------------------
// $Id: node_format.cpp 684 2005-11-23 13:17:15Z gems $
//
// C/C++ interface for writing/reading DBR and DCH files
// Works  with DATACH and DATABR structures
//
// Copyright (C) 2006,2011 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of the GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#include <iomanip>
#include  <iostream>
#include "io_arrays.h"
#include "node.h"
#include "gdatastream.h"

extern bool _comment;
extern const char* _GEMIPM_version_stamp;

//===============================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here.
//
outField DataBR_fields[51] =  {
  { "NodeHandle",  0, 0 },
  { "NodeTypeHY",  0, 0 },
  { "NodeTypeMT",  0, 0 },
  { "NodeStatusFMT",  1, 0 },
  { "NodeStatusCH",  1, 0 },
  { "IterDone",  0, 0 },
  { "TK",   1, 0 },
  { "P",  1, 0 },
  { "Vs",  0, 0 },
  { "Vi",   0, 0 },
  { "Ms",   0, 0 },
  { "Mi",   0, 0 },
  { "Hs",   0, 0 },
  { "Hi",   0, 0 },
  { "Gs",   0, 0 },
  { "IS",   0, 0 },
  { "pH",   0, 0 },
  { "pe",   0, 0 },
  { "Eh",   0, 0 },
  { "Tm",   0, 0 },
  { "dt",    0, 0 },
  { "Dif",   0, 0 },
  { "Vt",    0, 0 },
  { "vp",   0, 0 },
  { "eps",   0, 0 },
  { "Km",    0, 0 },
  { "Kf",   0, 0 },
  { "S",   0, 0 },
  { "Tr",   0, 0 },
  { "h",   0, 0 },
  { "rho",   0, 0 },
  { "al",   0, 0 },
  { "at",    0, 0 },
  { "av",   0, 0 },
  { "hDl",    0, 0 },
  { "hDt",   0, 0 },
  { "hDv",   0, 0 },
  { "nto",   0, 0 },
// dynamic arrays
  { "bIC",     1, 0 },
  { "rMB",     0, 0 },
  { "uIC",     0, 0 },
  { "xDC",     0, 0 },
  { "gam",    0, 0 },
  { "dll",    0, 0 },   // changed to non-obligatory by KD on 19.12.2006
  { "dul",    0, 0 },   // changed to non-obligatory by KD on 19.12.2006
  { "aPH",    0, 0 },   // changed to non-obligatory by KD on 19.12.2006
  { "xPH",    0, 0 },
  { "vPS",    0, 0 },
  { "mPS",   0, 0 },
  { "bPS",   0, 0 },
  { "xPA",   0, 0 }
};

outField DataCH_static_fields[13] =  {
  { "nIC",   1, 0 },
  { "nDC",   1, 0 },
  { "nPH",   1, 0 },
  { "nPS",    1, 0 },
  { "nDCs",   1, 0 },
  { "nICb",    1, 0 },
  { "nDCb",   1, 0 },
  { "nPHb",   1, 0 },
  { "nPSb",   1, 0 },
  { "nTp",    1, 0 },
  { "nPp",    1, 0 },
  { "iGrd",    1, 0 },
  { "fAalp",   1, 0 }
};

outField DataCH_dynamic_fields[29] =  { //+4
   { "xic",  1, 0 },
   { "xdc",  1, 0 },
   { "xph",  1, 0 },
   { "ICNL",  1, 0 },
   { "ccIC",  1, 0 },
   { "ICmm",  1, 0 },
   { "DCNL",  1, 0 },
   { "ccDC",  1, 0 },
   { "DCmm",  0, 0 },   // Changed to non-obligatory by DK on 3.05.2007
   { "PHNL",  1, 0 },
   { "ccPH",  1, 0 },
   { "nDCinPH",  1, 0 },
   { "A",  1, 0 },
   { "Ttol",  0, 0 },
   { "TKval",  1, 0 },
   { "Ptol",  0, 0 },
   { "Pval",  1, 0 },
   { "denW",  1, 0 },
   { "denWg",  1, 0 },
   { "epsW",  1, 0 },
   { "epsWg",  1, 0 },
//   { "visW",  1, 0 },
   { "V0",  1, 0 },
   { "G0",  1, 0 },
   { "H0", 0, 0 },
   { "S0",  0, 0 },
   { "Cp0",  0, 0 },
   { "A0",  0, 0 },
   { "U0",  0, 0 },
   { "DD",  0, 0 }    // Depending on iGrd flag
};


//===============================================================

void TNode::databr_to_text_file( fstream& ff, bool with_comments, bool brief_mode, const char* path )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");
  bool _comment = with_comments;

  TPrintArrays  prar(51, DataBR_fields, ff);

   if( _comment )
   {
      ff << "# " << _GEMIPM_version_stamp << endl << "# File: " << path << endl;
      ff << "# Comments can be marked with # $ ;" << endl << endl;
      ff << "# Template for the dbr-dat text input file for DATABR (node) data" << endl;
      ff << "# (should be read only after the DATACH and the IPM-DAT files)" << endl << endl;
      ff << "#Section (scalar-1): Controls of the GEM IPM operation and data exchange" << endl;
   }
   if(!brief_mode || prar.getAlws("NodeHandle" ))
   {  if( _comment )
        ff << "# NodeHandle: Node identification handle" << endl;
      ff << left << setw(17) << "<NodeHandle> " <<  CNode->NodeHandle << endl;
   }
if( CNode->NodeStatusFMT != No_transport )
{
  if(!brief_mode || prar.getAlws("NodeTypeHY" ))
  {  if( _comment )
       ff << "# NodeTypeHY:  Node type code (hydraulic), not used on TNode level ; see typedef NODETYPE" << endl;
     ff << left << setw(17) << "<NodeTypeHY> " <<  CNode->NodeTypeHY << endl;
  }
  if(!brief_mode || prar.getAlws("NodeTypeMT" ))
  { if( _comment )
      ff << "# NodeTypeMT:  Node type (mass transport), not used on TNode level; see typedef NODETYPE" << endl;
    ff << left << setw(17) << "<NodeTypeMT> " <<  CNode->NodeTypeMT << endl;
  }
}
  if(!brief_mode || prar.getAlws("NodeStatusFMT" ))
  { if( _comment )
      ff << "# NodeStatusFMT:  Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT" << endl;
    ff << left << setw(17) << "<NodeStatusFMT> " <<  CNode->NodeStatusFMT << endl;
  }
  if(!brief_mode || prar.getAlws("NodeStatusCH" ))
  {
    if( _comment )
      ff << "# NodeStatusCH: Node status code in GEM (input and output); see typedef NODECODECH" << endl;
    ff << left << setw(17) << "<NodeStatusCH> " <<  CNode->NodeStatusCH << endl;
  }
  if(!brief_mode || prar.getAlws("NodeStatusCH" ))
  {
   if( _comment )
      ff << "# IterDone:  Number of iterations performed by GEM IPM in the last run - GEM output" << endl;
   ff << left << setw(17) << "<IterDone> " <<  CNode->IterDone << endl;
  }
  if( _comment )
      ff << "\n## (2) Chemical scalar variables" << endl;
  if(!brief_mode || prar.getAlws("TK" ))
  { if( _comment )
         ff << "# TK: Node temperature T (Kelvin). This value must always be provided - GEM input." << endl;
    ff << left << setw(7) << "<TK> " <<  CNode->TK << endl;
  }
  if(!brief_mode || prar.getAlws("P" ))
  {  if( _comment )
         ff << "# P:  Node Pressure P (Pa)  GEM input" << endl;
     ff << left << setw(7) << "<P> " <<  CNode->P << endl;
  }
  if(!brief_mode || prar.getAlws("Vs" ))
  { if( _comment )
         ff << "# Vs: Volume V of reactive subsystem (m3)" << endl;
    ff << left << setw(7) << "<Vs> " << CNode->Vs << endl;
  }
if( CNode->NodeStatusFMT != No_transport )
{
  if(!brief_mode || prar.getAlws("Vi" ))
  { if( _comment )
         ff << "# Vi: Volume of inert subsystem (m3)" << endl;
    ff << left << setw(7) << "<Vi> " <<  CNode->Vi << endl;
  }
}
  if(!brief_mode || prar.getAlws("Ms" ))
  { if( _comment )
         ff << "# Ms:  Mass of reactive subsystem (kg) - GEM output" << endl;
    ff << left << setw(7) << "<Ms> " <<  CNode->Ms << endl;
  }
if( CNode->NodeStatusFMT != No_transport )
{
  if(!brief_mode || prar.getAlws("Mi" ))
  { if( _comment )
         ff << "# Mi: Mass of inert subsystem (kg)" << endl;
    ff << left << setw(7) << "<Mi> " <<  CNode->Mi << endl;
  }
}
  if(!brief_mode || prar.getAlws("Hs" ))
  { if( _comment )
         ff << "# Hs:  Total enthalpy of reactive subsystem (J) (reserved)" << endl;
    ff << left << setw(7) << "<Hs> " <<  CNode->Hs << endl;
  }
if( CNode->NodeStatusFMT != No_transport )
{
  if(!brief_mode || prar.getAlws("Hi" ))
  { if( _comment )
         ff << "# Hi:  Total enthalpy of inert subsystem (J) (reserved, can be used only in FMT part) " << endl;
    ff << left << setw(7) << "<Hi> " <<  CNode->Hi << endl;
  }
}
  if(!brief_mode || prar.getAlws("Gs" ))
  { if( _comment )
         ff << "# Gs: Total Gibbs energy of the reactive subsystem (J/RT) (normalized) - GEM output" << endl;
    ff << left << setw(7) << "<Gs> " <<  CNode->Gs << endl;
  }
if( CSD->ccPH[0] == PH_AQUEL )
{
  if(!brief_mode || prar.getAlws("IS" ))
  {	if( _comment )
         ff << "# IS: Effective aqueous ionic strength (molal)  - GEM output" << endl;
    ff << left << setw(7) << "<IS> " <<  CNode->IC << endl;
  }
  if(!brief_mode || prar.getAlws("pH" ))
  { if( _comment )
         ff << "# pH: pH of aqueous solution in the activity scale - GEM output" << endl;
    ff << left << setw(7) << "<pH> " <<  CNode->pH << endl;
  }
  if(!brief_mode || prar.getAlws("pe" ))
  {    if( _comment )
         ff << "# pe: pe of aqueous solution in the activity scale - GEM output " << endl;
     ff << left << setw(7) << "<pe> " <<  CNode->pe << endl;
  }
  if(!brief_mode || prar.getAlws("Eh" ))
  { if( _comment )
         ff << "# Eh: Eh of aqueous solution (V) - GEM output" << endl;
    ff << left << setw(7) << "<Eh> " <<  CNode->Eh << endl;
  }
}
if( CNode->NodeStatusFMT != No_transport )
{
	if( _comment )
	    ff << "\n## (3) FMT scalar variables (used only in NodeArray, not used in GEM)" << endl;
	if(!brief_mode || prar.getAlws("Tm" ))
    {  if( _comment )
          ff << "# Tm: Actual total simulation time (s)" << endl;
       ff << left << setw(7) << "<Tm> " <<  CNode->Tm << endl;
    }
	if(!brief_mode || prar.getAlws("dt" ))
    {  if( _comment )
        ff << "# dt:  Actual time step (s)" << endl;
      ff << left << setw(7) << "<dt> " <<  CNode->dt << endl;
    }
	if(!brief_mode || prar.getAlws("Dif" ))
    { if( _comment )
       ff << "# Dif: General diffusivity of disolved matter (m2/s)" << endl;
     ff << left << setw(7) << "<Dif> " <<  CNode->Dif << endl;
    }
	if(!brief_mode || prar.getAlws("Vt" ))
    { if( _comment )
       ff << "# Vt: Total volume of the node (m3)" << endl;
     ff << left << setw(7) << "<Vt> " <<  CNode->Vt << endl;
    }
	if(!brief_mode || prar.getAlws("vp" ))
    { if( _comment )
       ff << "# vp:  Advection velocity (in pores)  (m/s)" << endl;
      ff << left << setw(7) << "<vp> " <<  CNode->vp << endl;
    }
	if(!brief_mode || prar.getAlws("eps" ))
    { if( _comment )
       ff << "#  eps: Effective (actual) porosity normalized to 1" << endl;
      ff << left << setw(7) << "<eps> " <<  CNode->eps << endl;
    }
	if(!brief_mode || prar.getAlws("Km" ))
    { if( _comment )
       ff << "# Km: Actual permeability (m2)" << endl;
      ff << left << setw(7) << "<Km> " <<  CNode->Km << endl;
    }
	if(!brief_mode || prar.getAlws("Kf" ))
    { if( _comment )
       ff << "# Kf: Actual Darcy`s constant (m2/s)" << endl;
     ff << left << setw(7) << "<Kf> " <<  CNode->Kf << endl;
    }
	if(!brief_mode || prar.getAlws("S" ))
    { if( _comment )
       ff << "# S: Specific storage coefficient, dimensionless" << endl;
      ff << left << setw(7) << "<S> " <<  CNode->S << endl;
    }
	if(!brief_mode || prar.getAlws("Tr" ))
    { if( _comment )
       ff << "# Tr:  Transmissivity (m2/s)" << endl;
      ff << left << setw(7) << "<Tr> " <<  CNode->Tr << endl;
    }
	if(!brief_mode || prar.getAlws("h" ))
    { if( _comment )
       ff << "# h:  Actual hydraulic head (hydraulic potential) (m)" << endl;
      ff << left << setw(7) << "<h> " <<  CNode->h << endl;
    }
	if(!brief_mode || prar.getAlws("rho" ))
    { if( _comment )
       ff << "# rho:  Actual carrier density for density-driven flow (kg/m3)" << endl;
     ff << left << setw(7) << "<rho> " <<  CNode->rho << endl;
    }
	if(!brief_mode || prar.getAlws("al" ))
    { if( _comment )
       ff << "# al: Specific longitudinal dispersivity of porous media (m)" << endl;
      ff << left << setw(7) << "<al> " <<  CNode->al << endl;
    }
	if(!brief_mode || prar.getAlws("at" ))
    { if( _comment )
       ff << "# at:  Specific transversal dispersivity of porous media (m)." << endl;
     ff << left << setw(7) << "<at> " <<  CNode->at << endl;
    }
	if(!brief_mode || prar.getAlws("av" ))
    { if( _comment )
       ff << "# av:  Specific vertical dispersivity of porous media (m). " << endl;
     ff << left << setw(7) << "<av> " <<  CNode->av << endl;
    }
	if(!brief_mode || prar.getAlws("hDl" ))
    { if( _comment )
       ff << "# hDl: Hydraulic longitudinal dispersivity (m2/s)" << endl;
     ff << left << setw(7) << "<hDl> " <<  CNode->hDl << endl;
    }
	if(!brief_mode || prar.getAlws("hDt" ))
    { if( _comment )
       ff << "# hDt: Hydraulic transversal dispersivity (m2/s)" << endl;
     ff << left << setw(7) << "<hDt> " <<  CNode->hDt << endl;
    }
	if(!brief_mode || prar.getAlws("hDv" ))
    { if( _comment )
       ff << "# hDv: Hydraulic vertical dispersivity (m2/s)" << endl;
     ff << left << setw(7) << "<hDv> " <<  CNode->hDv << endl;
    }
	if(!brief_mode || prar.getAlws("nto" ))
    { if( _comment )
       ff << "# nto:  Tortuosity factor (dimensionless)" << endl;
      ff << left << setw(7) << "<nto> " <<  CNode->nto << endl;
    }
}
   if( _comment )
   {   ff << "\n### Arrays - for dimensions and index lists, see Section (2) of DATACH file" << endl << endl;
       ff << "## (4) IC data section";
       prar.writeArray(  NULL, CSD->ICNL[0], CSD->nIC, MaxICN );
       ff << endl;
   }
  if(!brief_mode || prar.getAlws("bIC" ))
  {  if( _comment )
       ff << "# bIC:  Bulk composition of (reactive part of) the system - main GEM input (amounts of IC in moles) [nICb]";
     prar.writeArray(  "bIC",  CNode->bIC, CSD->nICb );
  }
  if(!brief_mode || prar.getAlws("rMB" ))
  { if( _comment )
       ff << "\n\n# rMB: Mass balance residuals (moles) [nICb] - GEM output";
    prar.writeArray(  "rMB",  CNode->rMB, CSD->nICb );
  }
  if(!brief_mode || prar.getAlws("uIC" ))
  { if( _comment )
       ff << "\n\n# uIC: Chemical potentials of ICs (dual GEM solution) - GEM output, normalized scale [nICb]";
    prar.writeArray(  "uIC",  CNode->uIC, CSD->nICb );
  }

  if( _comment )
  {    ff << "\n\n## (5) DC data section";
       prar.writeArray(  NULL, CSD->DCNL[0], CSD->nDC, MaxDCN );
       ff << endl;
  }

  if(!brief_mode || prar.getAlws("xDC" ))
  { if( _comment )
      ff << "# xDC:  Speciation - amounts of DCs in equilibrium state - primal GEM solution (moles) [nDCb] - GEM output";
    prar.writeArray(  "xDC",  CNode->xDC, CSD->nDCb );
  }
  if(!brief_mode || prar.getAlws("gam" ))
  { if( _comment )
       ff << "\n\n# gam:  Activity coefficients of DCs in their respective phases [nDCb] - GEM output";
    prar.writeArray(  "gam",  CNode->gam, CSD->nDCb );
  }
  if(!brief_mode || prar.getAlws("dll" ))
  { if( _comment )
       ff << "\n\n# dll: Lower metastability constraints on amounts of DCs (moles) [nDCb] - GEM input";
    prar.writeArray(  "dll",  CNode->dll, CSD->nDCb );
  }
  if(!brief_mode || prar.getAlws("dul" ))
  { if( _comment )
       ff << "\n\n# dul:  Upper metastability constraints on amounts of DCs (moles) [nDCb] - GEM input";
    prar.writeArray(  "dul",  CNode->dul, CSD->nDCb );
  }
   if( _comment )
   {    ff << "\n\n## (6) Phase data section";
        prar.writeArray(  NULL, CSD->PHNL[0], CSD->nPH, MaxPHN );
        ff << endl;
   }
  if(!brief_mode || prar.getAlws("aPH" ))
  {
    if( _comment )
     ff << "# aPH: Specific surface areas of phases (m2/kg) [nPHb] - GEM input";
    prar.writeArray(  "aPH",  CNode->aPH, CSD->nPHb );
  }
  if(!brief_mode || prar.getAlws("xPH" ))
  { if( _comment )
        ff << "\n\n# xPH: Amounts of phases in equilibrium state (moles) [nPHb] - GEM output";
    prar.writeArray(  "xPH",  CNode->xPH, CSD->nPHb );
  }
  if(!brief_mode || prar.getAlws("vPS" ))
  { if( _comment )
        ff << "\n\n# vPS: Volumes of multicomponent phases (m3) [nPSb] - GEM output";
    prar.writeArray(  "vPS",  CNode->vPS, CSD->nPSb );
  }
  if(!brief_mode || prar.getAlws("mPS" ))
  { if( _comment )
        ff << "\n\n# mPS: Masses of multicomponent phases (kg) [nPSb] - GEM output";
    prar.writeArray(  "mPS",  CNode->mPS, CSD->nPSb );
  }
  if(!brief_mode || prar.getAlws("xPA" ))
  { if( _comment )
       ff << "\n\n# xPA: Amount of carrier (sorbent or solvent) in multicomponent phases [nPSb] - GEM output";
    prar.writeArray(  "xPA",  CNode->xPA, CSD->nPSb );
  }
  if(!brief_mode || prar.getAlws("bPS" ))
  {  if( _comment )
     {
	  ff << "\n\n# bPS: Bulk elemental compositions of multicomponent phases (moles) [nPSb*nICb]- GEM output";
	  prar.writeArray(  NULL, CSD->ICNL[0], CSD->nIC, MaxICN );
//	  ff << endl;
      }
     prar.writeArray(  "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb, CSD->nICb );
  }
 ff << endl;
  if( _comment )
   {     ff << "\n# reserved" << endl;
         ff << "\n# End of file"<< endl;
   }
}

// Reading work dataBR structure from text file
void TNode::databr_from_text_file( fstream& ff )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

 // mem_set( &CNode->Tm, 0, 19*sizeof(double));
 databr_reset( CNode );
 TReadArrays  rdar( 51, DataBR_fields, ff);
 long int nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
    case 0: rdar.readArray( "NodeHandle",  &CNode->NodeHandle, 1);
            break;
    case 1: rdar.readArray( "NodeTypeHY",  &CNode->NodeTypeHY, 1);
            break;
    case 2: rdar.readArray( "NodeTypeMT",  &CNode->NodeTypeMT, 1);
            break;
    case 3: rdar.readArray( "NodeStatusFMT",  &CNode->NodeStatusFMT, 1);
            break;
    case 4: rdar.readArray( "NodeStatusCH",  &CNode->NodeStatusCH, 1);
            break;
    case 5: rdar.readArray( "IterDone",  &CNode->IterDone, 1);
            break;
    case 6: rdar.readArray( "TK",  &CNode->TK, 1);
            break;
    case 7: rdar.readArray( "P",  &CNode->P, 1);
            break;
    case 8: rdar.readArray( "Vs", &CNode->Vs, 1);
            break;
    case 9: rdar.readArray( "Vi",  &CNode->Vi, 1);
            break;
    case 10: rdar.readArray( "Ms",  &CNode->Ms, 1);
            break;
    case 11: rdar.readArray( "Mi",  &CNode->Mi, 1);
            break;
    case 12: rdar.readArray( "Hs",  &CNode->Hs, 1);
            break;
    case 13: rdar.readArray( "Hi",  &CNode->Hi, 1);
            break;
    case 14: rdar.readArray( "Gs",  &CNode->Gs, 1);
             break;
    case 15: rdar.readArray( "IS",  &CNode->IC, 1);
            break;
    case 16: rdar.readArray( "pH",  &CNode->pH, 1);
            break;
    case 17: rdar.readArray( "pe",  &CNode->pe, 1);
            break;
    case 18: rdar.readArray( "Eh",  &CNode->Eh, 1);
            break;
    case 19: rdar.readArray( "Tm",  &CNode->Tm, 1);
            break;
    case 20: rdar.readArray( "dt",  &CNode->dt, 1);
            break;
    case 21: rdar.readArray( "Dif",  &CNode->Dif, 1);
            break;
    case 22: rdar.readArray( "Vt",  &CNode->Vt, 1);
            break;
    case 23: rdar.readArray( "vp",  &CNode->vp, 1);
            break;
    case 24: rdar.readArray( "eps",  &CNode->eps, 1);
            break;
    case 25: rdar.readArray( "Km",  &CNode->Km, 1);
            break;
    case 26: rdar.readArray( "Kf",  &CNode->Kf, 1);
            break;
    case 27: rdar.readArray( "S",  &CNode->S, 1);
            break;
    case 28: rdar.readArray( "Tr",  &CNode->Tr, 1);
            break;
    case 29: rdar.readArray( "h",  &CNode->h, 1);
            break;
    case 30: rdar.readArray( "rho",  &CNode->rho, 1);
            break;
    case 31: rdar.readArray( "al",  &CNode->al, 1);
            break;
    case 32: rdar.readArray( "at",  &CNode->at, 1);
            break;
    case 33: rdar.readArray( "av",  &CNode->av, 1);
            break;
    case 34: rdar.readArray( "hDl",  &CNode->hDl, 1);
            break;
    case 35: rdar.readArray( "hDt",  &CNode->hDt, 1);
            break;
    case 36: rdar.readArray( "hDv",  &CNode->hDv, 1);
            break;
    case 37: rdar.readArray( "nto",  &CNode->nto, 1);
            break;
    case 38: rdar.readArray( "bIC",  CNode->bIC, CSD->nICb );
            break;
    case 39: rdar.readArray( "rMB",  CNode->rMB, CSD->nICb );
            break;
    case 40: rdar.readArray( "uIC",  CNode->uIC, CSD->nICb );
            break;
    case 41: rdar.readArray( "xDC",  CNode->xDC, CSD->nDCb );
            break;
    case 42: rdar.readArray( "gam",  CNode->gam, CSD->nDCb );
            break;
    case 43: rdar.readArray( "dll",  CNode->dll, CSD->nDCb );
            break;
    case 44: rdar.readArray( "dul",  CNode->dul, CSD->nDCb );
            break;
    case 45: rdar.readArray( "aPH",  CNode->aPH, CSD->nPHb );
            break;
    case 46: rdar.readArray( "xPH",  CNode->xPH, CSD->nPHb );
            break;
    case 47: rdar.readArray( "vPS",  CNode->vPS, CSD->nPSb );
            break;
    case 48: rdar.readArray( "mPS",  CNode->mPS, CSD->nPSb );
            break;
    case 49: rdar.readArray( "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
            break;
    case 50: rdar.readArray( "xPA",  CNode->xPA, CSD->nPSb );
   }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataBR structure";
    Error( "Error", ret);
  }
}


void TNode::datach_to_text_file( fstream& ff, bool with_comments, bool brief_mode, const char* path )
{
// fstream ff("DataCH.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");
  bool _comment = with_comments;
  TPrintArrays  prar(29, DataCH_dynamic_fields, ff);
  if( CSD->nIC == CSD->nICb )
	  prar.setNoAlws( "xic");
  if(CSD->nDC == CSD->nDCb )
	  prar.setNoAlws( 1 /*"xdc"*/);
  if(CSD->nPH == CSD->nPHb )
	  prar.setNoAlws( 2 /*"xph"*/);

  if( _comment )
  {
     ff << "# " << _GEMIPM_version_stamp << endl << "# File: " << path << endl;
     ff << "# Comments are marked with # $ ;" << endl;
     ff << "\n# Template for the dch-dat text input file for DATACH data " << endl;
     ff << "# (should be read first, before the IPM-DAT file and DATABR files)" << endl;
     ff << "\n## (1) Dimensions for memory allocation" << endl;
     ff << "# nIC: Number of Independent Components (stoichiometry units, usually chemical elements and charge)"<< endl;
  }
  ff << left << setw(7) << "<nIC> " <<  CSD->nIC << endl;
  if( _comment )
     ff << "# nDC: Total number of Dependent Components (chemical species made of Independent Components)" << endl;
  ff << left << setw(7) << "<nDC> " <<  CSD->nDC << endl;
  if( _comment )
     ff << "# nPH: Number of phases (into which Dependent Components are grouped)" << endl;
  ff << left << setw(7) << "<nPH> " <<  CSD->nPH << endl;
  if( _comment )
     ff << "# nPS: Number of phases-solutions (multicomponent phases) in the chemical system definition (<= nPH)" << endl;
  ff << left << setw(7) << "<nPS> " <<  CSD->nPS << endl;
  if( _comment )
     ff << "# nDCs: Number of Dependent Components in phases-solutions (multicomponent phases)" << endl;
  ff << left << setw(7) << "<nDCs> " <<  CSD->nDCs << endl;

  if( _comment )
  {  ff << "\n## (2) Databridge configuration section (for memory allocation)" << endl;
     ff << "# nICb: Number of Independent Components kept in the DBR file and DATABR memory structure (<= nIC)" << endl;
  }
  ff << left << setw(7) << "<nICb> " <<  CSD->nICb << endl;
  if( _comment )
     ff << "# nDCb: Number of Dependent Components kept in the DBR file and DATABR memory structure (<=nDC)" << endl;
  ff << left << setw(7) << "<nDCb> " <<  CSD->nDCb << endl;
  if( _comment )
     ff << "# nPHb: Number of Phases to be kept in the DBR file and DATABR structure (<=nPH)" << endl;
  ff << left << setw(7) << "<nPHb> " <<  CSD->nPHb << endl;
  if( _comment )
     ff << "# nPSb: Number of Phases-solutions to be kept in the DBR file and DATABR memory structure (<=nPS)" << endl;
  ff << left << setw(7) << "<nPSb> " <<  CSD->nPSb << endl;

  if( _comment )
  {   ff << "\n## (3) Dimensions for thermodynamic data arrays" << endl;
      ff << "# nTp: Number of temperature grid points in interpolation lookup arrays, 1 or more" << endl;
  }
  ff << left << setw(7) << "<nTp> " <<  CSD->nTp << endl;
  if( _comment )
     ff << "# nPp: Number of pressure grid points in interpolation lookup arrays, 1 or more" << endl;
  ff << left << setw(7) << "<nPp> " <<  CSD->nPp << endl;
  if( _comment )
   {  ff << "# iGrd: Flag for selection Diffusition coefficients array provided in the DCH file." << endl;
   }
  ff << left << setw(7) << "<iGrd> " <<  CSD->iGrd << endl;
  if( _comment )
    ff << "# fAalp: Flag for keeping specific surface areas of phases in DATABR structure (1) or ignoring them (0)" << endl;
  ff << left << setw(7) << "<fAalp> " <<  CSD->nAalp << endl;

  ff<< "\n<END_DIM>\n";

// dynamic arrays - must follow static data
  if( _comment )
     ff << "\n## (4) Databridge configuration section (for memory allocation)";
  if(!brief_mode || prar.getAlws("xIC" ))
  { if( _comment )
      ff << "\n# xIC: DATACH access index list for Independent Components kept in the DATABR structure and in DBR files [nICb]";
    prar.writeArray(  "xic", CSD->xic, CSD->nICb);
  }
  if(!brief_mode || prar.getAlws("xDC" ))
  { if( _comment )
     ff << "\n# xDC: DATACH access index list of Dependent Components kept in the DATABR  structure and in DBR files [nDCb]";
    prar.writeArray(  "xdc", CSD->xdc, CSD->nDCb);
  }
  if(!brief_mode || prar.getAlws("xPH" ))
  { if( _comment )
      ff << "\n# xPH: DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]";
    prar.writeArray(  "xph", CSD->xph, CSD->nPHb);
  }

  if( _comment )
     ff << "\n\n## (5) Independent components section";
  if(!brief_mode || prar.getAlws("ICNL" ))
  {
     if( _comment )
         ff << "\n# ICNL: Name list of Independent Components (up to 4 characters per name) [nIC]";
      prar.writeArray(  "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
  }
  if(!brief_mode || prar.getAlws("ccIC" ))
  { if( _comment )
      ff << "\n# ccIC: Codes of Independent Components [nIC]";
    prar.writeArray(  "ccIC", CSD->ccIC, CSD->nIC, 1L );
  }
  if(!brief_mode || prar.getAlws("ICmm" ))
  { if( _comment )
      ff << "\n# ICmm: Atomic (molar) masses of Independent Components  (kg/mol) [nIC]";
    prar.writeArray(  "ICmm", CSD->ICmm, CSD->nIC);
  }

  if( _comment )
    ff << "\n\n## (6) Dependent Components section (codes and names)";
  if(!brief_mode || prar.getAlws("DCNL" ))
  {	  if( _comment )
       ff << "\n# DCNL: Name list of Dependent Components (up to 16 characters per name) [nDC]";
     prar.writeArray(  "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
  }
  if(!brief_mode || prar.getAlws("ccDC" ))
  { if( _comment )
     ff << "\n# ccDC: Type codes of Dependent Components [nDC]";
    prar.writeArray(  "ccDC", CSD->ccDC, CSD->nDC, 1L );
  }
  if(!brief_mode || prar.getAlws("DCmm" ))
  { if( _comment )
      ff << "\n\n# DCmm: Molar masses of Dependent Components (kg/mol) [nDC]";
     prar.writeArray(  "DCmm", CSD->DCmm, CSD->nDC);
  }
  if( _comment )
    ff << "\n\n## (7) Phases section" << endl;
  if(!brief_mode || prar.getAlws("PHNL" ))
  { if( _comment )
      ff << "# PHNL: List of Phase names (up to 16 characters per name) [nPH]";
    prar.writeArray(  "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
  }
  if(!brief_mode || prar.getAlws("ccPH" ))
  { if( _comment )
      ff << "\n# ccPH: Phase aggregate state codes [nPH]";
    prar.writeArray(  "ccPH", CSD->ccPH, CSD->nPH, 1L );
  }
  if(!brief_mode || prar.getAlws("nDCinPH" ))
  {  if( _comment )
       ff << "\n# nDCinPH: This vector tells how many Dependent Components is included in each phase [nPH]";
     prar.writeArray(  "nDCinPH", CSD->nDCinPH, CSD->nPH);
  }

  if( _comment )
    ff << "\n\n# (8) Data section for DCs";
  if(!brief_mode || prar.getAlws("A" ))
  { if( _comment )
     ff << "\n# A: Stoichiometry matrix A for Dependent Components, [nDC*nIC]";
    prar.writeArray(  "A", CSD->A, CSD->nDC*CSD->nIC, CSD->nIC );
  }
  ff << endl;
  if( _comment )
    ff << "\n## (9) Thermodynamic data section";
  if(!brief_mode || prar.getAlws("Ttol" ))
  { if( _comment )
     ff << "\n# Ttol: Tolerance for the temperature interpolation (K)" << endl;
    ff << left << setw(7) << "<Ttol> " <<  CSD->Ttol;
  }
  if(!brief_mode || prar.getAlws("TKval" ))
  { if( _comment )
      ff << "\n# Tval: Temperature values for the interpolation grid (K) for the lookup arrays of thermodynamic data [nTp]";
    prar.writeArray(  "TKval", CSD->TKval, CSD->nTp );
  }
  ff << endl;
  if(!brief_mode || prar.getAlws("Ptol" ))
  {  if( _comment )
      ff << "\n# Ptol: Tolerance for the pressure interpolation (Pa)" << endl;
    ff << left << setw(7) << "<Ptol> " <<  CSD->Ptol;
  }
  if(!brief_mode || prar.getAlws("Pval" ))
  { if( _comment )
      ff << "\n# Pval: Pressure values for the interpolation grid (Pa) for the lookup arrays of thermodynamic data [nPp]";
    prar.writeArray(  "Pval", CSD->Pval, CSD->nPp );
  }

  if( CSD->ccPH[0] == PH_AQUEL )
  {
    if(!brief_mode || prar.getAlws("denW" ))
    {  if( _comment )
         ff << "\n\n# denW: Lookup array for the density of water-solvent (kg/m3) [5*nPp*nTp]";
      prar.writeArray(  "denW", CSD->denW, 5*(CSD->nPp*CSD->nTp), CSD->nPp*CSD->nTp );
    }
    if(!brief_mode || prar.getAlws("denWg" ))
    {  if( _comment )
         ff << "\n\n# denWg: Optional lookup array for the density of water vapour (kg/m3) [5*nPp*nTp]";
      prar.writeArray(  "denWg", CSD->denWg, 5*(CSD->nPp*CSD->nTp), CSD->nPp*CSD->nTp );
    }
    if(!brief_mode || prar.getAlws("epsW" ))
    {  if( _comment )
        ff << "\n\n# epsW: Lookup array for the dielectric constant of water-solvent (dimensionless) [5*nPp*nTp]";
      prar.writeArray(  "epsW", CSD->epsW,  5*(CSD->nPp*CSD->nTp), CSD->nPp*CSD->nTp );
    }
    if(!brief_mode || prar.getAlws("epsWg" ))
    {  if( _comment )
        ff << "\n\n# epsWg: Optional lookup array for the dielectric constant of water vapour [5*nPp*nTp]";
      prar.writeArray(  "epsWg", CSD->epsWg, 5*(CSD->nPp*CSD->nTp),  CSD->nPp*CSD->nTp );
    }
  }
  if(!brief_mode || prar.getAlws("V0" ))
  { if( _comment )
      ff << "\n\n# V0: Obligatory lookup array for (standard) molar volumes of Dependent Components (J/Pa) [nDC*nPp*nTp]";
    prar.writeArray(  "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp, CSD->nPp*CSD->nTp );
  }
  if(!brief_mode || prar.getAlws("G0" ))
  { if( _comment )
     ff << "\n\n# G0: Obligatory lookup array for DC molar Gibbs energy function g(T,P) (J/mol) [nDC*nPp*nTp]";
    prar.writeArray(  "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp, CSD->nPp*CSD->nTp );
  }
 if(!brief_mode || prar.getAlws("H0" ))
    {  if( _comment )
        ff << "\n\n# H0: Optional lookup array for DC molar enthalpy h(T,P) (J/mol) [nDC*nPp*nTp]";
       prar.writeArray(  "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp,CSD->nPp*CSD->nTp );
    }
  if(!brief_mode || prar.getAlws("S0" ))
    { if( _comment )
       ff << "\n\n# S0: Optional lookup array for the DC absolute entropy function (J/K/mol) [nDC*nPp*nTp] ";
      prar.writeArray(  "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp,  CSD->nPp*CSD->nTp  );
    }
  if(!brief_mode || prar.getAlws("Cp0" ))
	 {  if( _comment )
            ff << "\n\n# Cp0: Optional lookup array for DC heat capacity function (J/K/mol) [nDC*nPp*nTp]";
        prar.writeArray(  "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp, CSD->nPp*CSD->nTp  );
	 }
 if(!brief_mode || prar.getAlws("A0" ))
	 {  if( _comment )
            ff << "\n\n# A0: Optional lookup array for Helmholtz energy of DC (J/mol) reserved [nDC*nPp*nTp]";
        prar.writeArray(  "A0", CSD->A0, CSD->nDC*CSD->nPp*CSD->nTp, CSD->nPp*CSD->nTp  );
	 }
 if(!brief_mode || prar.getAlws("U0" ))
	 {  if( _comment )
            ff << "\n\n# U0: Optional lookup array for Internal energy of DC (J/K/mol) [nDC*nPp*nTp]";
        prar.writeArray(  "U0", CSD->U0, CSD->nDC*CSD->nPp*CSD->nTp, CSD->nPp*CSD->nTp  );
	 }

  if( CSD->iGrd  )
  {
    if(!brief_mode || prar.getAlws("DD" ))
    { if( _comment )
        ff << "\n\n# DD: Lookup array for diffusion coefficients of DCs (reserved) [nDC*nPp*nTp]";
      prar.writeArray(  "DD", CSD->DD, CSD->nDCs*CSD->nPp*CSD->nTp,  CSD->nPp*CSD->nTp);
    }
  }
  ff << endl;
  if( _comment )
      ff << "\n# End of file";
}

// Reading dataCH structure from text file
void TNode::datach_from_text_file(fstream& ff)
{
  long int ii;
// fstream ff("DataCH.out", ios::in );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

// static arrays
 TReadArrays  rdar( 13, DataCH_static_fields, ff);
 long int nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
    case 0: rdar.readArray( "nIC", &CSD->nIC, 1);
            break;
    case 1: rdar.readArray( "nDC", &CSD->nDC, 1);
            break;
    case 2: rdar.readArray( "nPH", &CSD->nPH, 1);
            break;
    case 3: rdar.readArray( "nPS", &CSD->nPS, 1);
            break;
    case 4: rdar.readArray( "nDCs", &CSD->nDCs, 1);
            break;
    case 5: rdar.readArray( "nICb", &CSD->nICb, 1);
            break;
    case 6: rdar.readArray( "nDCb", &CSD->nDCb, 1);
            break;
    case 7: rdar.readArray( "nPHb", &CSD->nPHb, 1);
            break;
    case 8: rdar.readArray( "nPSb", &CSD->nPSb, 1);
            break;
    case 9: rdar.readArray( "nTp", &CSD->nTp, 1);
            break;
    case 10: rdar.readArray( "nPp", &CSD->nPp, 1);
            break;
    case 11: rdar.readArray( "iGrd", &CSD->iGrd, 1);
            break;
    case 12: rdar.readArray( "fAalp", &CSD->nAalp, 1);
  }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataCH structure";
    Error( "Error", ret);
  }

  datach_realloc();
  databr_realloc();

//dynamic data
 TReadArrays  rddar( 29, DataCH_dynamic_fields, ff);

   if( CSD->iGrd  )
      rddar.setNoAlws( 28 /*"DD"*/);
/*   if( CSD->iGrd <= 4 )
      rddar.setNoAlws( 27 *"U0"*);
   if( CSD->iGrd <= 3 )
      rddar.setNoAlws( 26 *"A0"*);
   if( CSD->iGrd <= 2 )
      rddar.setNoAlws( 25 *"Cp0"*);
   if( CSD->iGrd <= 0 )
      rddar.setNoAlws( 23 *"H0"*);
   if( CSD->iGrd <= 1 )
      rddar.setNoAlws( 24 *"S0"*);
*/
// default set up
  for( ii=0; ii< CSD->nDCs*CSD->nPp*CSD->nTp; ii++ )
  {
    if( CSD->DD ) CSD->DD[ii] = 0.;
    if( CSD->Cp0) CSD->Cp0[ii] = 0.;
    if( CSD->H0 ) CSD->H0[ii] = 0.;
    if( CSD->S0 ) CSD->S0[ii] = 0.;
    if( CSD->A0 ) CSD->A0[ii] = 0.;
    if( CSD->U0 ) CSD->U0[ii] = 0.;
  }
  CSD->Ttol = 0.1;
  CSD->Ptol = 10000;
  if( CSD->nIC == CSD->nICb )
  {
    rddar.setNoAlws( "xic" );
    for( ii=0; ii< CSD->nICb; ii++ )
      CSD->xic[ii] = ii;
  }
  if(CSD->nDC == CSD->nDCb )
  {
    rddar.setNoAlws( 1 /*"xdc"*/);
    for( ii=0; ii< CSD->nDCb; ii++ )
      CSD->xdc[ii] = ii;
  }
  if(CSD->nPH == CSD->nPHb )
  {
    rddar.setNoAlws( 2 /*"xph"*/);
    for( ii=0; ii< CSD->nPHb; ii++ )
      CSD->xph[ii] = ii;
  }

  nfild = rddar.findNext();
  while( nfild >=0 )
  {
   switch( nfild )
   {
    case 0: rddar.readArray( "xic", CSD->xic, CSD->nICb);
            break;
    case 1: rddar.readArray( "xdc", CSD->xdc, CSD->nDCb);
            break;
    case 2: rddar.readArray( "xph", CSD->xph, CSD->nPHb);
            break;
    case 3: rddar.readArray( "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
            break;
    case 4: rddar.readArray( "ccIC", CSD->ccIC, CSD->nIC, 1 );
            break;
    case 5: rddar.readArray( "ICmm", CSD->ICmm, CSD->nIC);
            break;
    case 6: rddar.readArray( "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
            break;
    case 7: rddar.readArray( "ccDC", CSD->ccDC, CSD->nDC, 1 );
            break;
    case 8: rddar.readArray( "DCmm", CSD->DCmm, CSD->nDC);
            break;
    case 9: rddar.readArray( "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
            break;
    case 10: rddar.readArray( "ccPH", CSD->ccPH, CSD->nPH, 1 );
            break;
    case 11: rddar.readArray( "nDCinPH", CSD->nDCinPH, CSD->nPH);
            break;
    case 12: rddar.readArray( "A", CSD->A, CSD->nDC*CSD->nIC );
            break;
    case 13: rddar.readArray( "Ttol", &CSD->Ttol, 1);
            break;
    case 14: rddar.readArray( "TKval", CSD->TKval, CSD->nTp );
            break;
    case 15: rddar.readArray( "Ptol", &CSD->Ptol, 1);
            break;
    case 16: rddar.readArray( "Pval", CSD->Pval, CSD->nPp );
              break;
    case 17: if( !CSD->denW )
                   Error( "Error", "Array denW is not allocated in DCH!");
             rddar.readArray( "denW", CSD->denW, 5*CSD->nPp*CSD->nTp );
              break;
    case 18: if( !CSD->denWg )
                   Error( "Error", "Array denWg is not allocated in DCH!");
             rddar.readArray( "denWg", CSD->denWg, 5*CSD->nPp*CSD->nTp );
              break;
    case 19: if( !CSD->epsW )
                   Error( "Error", "Array epsW is not allocated in DCH!");
             rddar.readArray( "epsW", CSD->epsW,  5*CSD->nPp*CSD->nTp );
            break;
    case 20: if( !CSD->epsWg )
                   Error( "Error", "Array epsWg is not allocated in DCH!");
             rddar.readArray( "epsWg", CSD->epsWg,  5*CSD->nPp*CSD->nTp );
            break;
    case 21: rddar.readArray( "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
            break;
    case 22: rddar.readArray( "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp );
              break;
    case 23: if( !CSD->H0 )
                   Error( "Error", "Array HO is not allocated in DCH!");
            rddar.readArray( "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp);
            break;
    case 24: if( !CSD->S0 )
                   Error( "Error", "Array S0 is not allocated in DCH!");
            rddar.readArray( "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp);
            break;
    case 25: if( !CSD->Cp0 )
                   Error( "Error", "Array CpO is not allocated in DCH!");
            rddar.readArray( "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp );
            break;
    case 26: if( !CSD->A0 )
                   Error( "Error", "Array AO is not allocated in DCH!");
            rddar.readArray( "A0", CSD->A0, CSD->nDC*CSD->nPp*CSD->nTp );
            break;
    case 27: if( !CSD->U0 )
                   Error( "Error", "Array UO is not allocated in DCH!");
            rddar.readArray( "U0", CSD->U0, CSD->nDC*CSD->nPp*CSD->nTp );
            break;
    case 28: if( !CSD->DD )
                    Error( "Error", "Array DD is not allocated in DCH!");
            rddar.readArray( "DD", CSD->DD, CSD->nDCs*CSD->nPp*CSD->nTp);
           break;
  }
     nfild = rddar.findNext();
 }

  // Set up flags
  if( CSD->ccPH[0] != PH_AQUEL )
  {
        rddar.setNoAlws( 17 /*"denW"*/);
        rddar.setNoAlws( 18 /*"denWg"*/);
        rddar.setNoAlws( 19 /*"epsW"*/);
        rddar.setNoAlws( 20 /*"epsWg"*/);
  }

 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataCH structure";
    Error( "Error", ret);
  }

}

//---------------------------------------------------------------
// new i/o structures

// Writing DataCH to binary file
void TNode::datach_to_file( GemDataStream& ff )
{
// const data
   ff.writeArray( &CSD->nIC, 14 );
   ff.writeArray( &CSD->Ttol, 4 );

//dynamic data
   ff.writeArray( CSD->nDCinPH, CSD->nPH );
//   if( CSD->nICb >0 )
   ff.writeArray( CSD->xic, CSD->nICb );
   ff.writeArray( CSD->xdc, CSD->nDCb );
   ff.writeArray( CSD->xph, CSD->nPHb );

   ff.writeArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.writeArray( CSD->ICmm, CSD->nIC );
   ff.writeArray( CSD->DCmm, CSD->nDC );

   ff.writeArray( CSD->TKval,  CSD->nTp );
   ff.writeArray( CSD->Pval,  CSD->nPp );

   ff.writeArray( CSD->ccIC, CSD->nIC );
   ff.writeArray( CSD->ccDC, CSD->nDC );
   ff.writeArray( CSD->ccPH, CSD->nPH );

   if( CSD->ccPH[0] == PH_AQUEL )
   { ff.writeArray( CSD->denW,  5*CSD->nPp*CSD->nTp );
     ff.writeArray( CSD->denWg,  5*CSD->nPp*CSD->nTp );
     ff.writeArray( CSD->epsW, 5*CSD->nPp*CSD->nTp );
     ff.writeArray( CSD->epsWg, 5*CSD->nPp*CSD->nTp );
   }
   ff.writeArray( CSD->G0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->S0, CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->Cp0, CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->A0, CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->U0, CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd  )
      ff.writeArray( CSD->DD, CSD->nDCs*CSD->nPp*CSD->nTp );

   ff.writeArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.writeArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.writeArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );
}

// Reading DataCH structure from binary file
void TNode::datach_from_file( GemDataStream& ff )
{
// const data
   ff.readArray( &CSD->nIC, 14 );
   ff.readArray( &CSD->Ttol, 4 );

  datach_realloc();
  databr_realloc();

//dynamic data
   ff.readArray( CSD->nDCinPH, CSD->nPH );
//   if( CSD->nICb >0 )
   ff.readArray( CSD->xic, CSD->nICb );
   ff.readArray( CSD->xdc, CSD->nDCb );
   ff.readArray( CSD->xph, CSD->nPHb );

   ff.readArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.readArray( CSD->ICmm, CSD->nIC );
   ff.readArray( CSD->DCmm, CSD->nDC );

   ff.readArray( CSD->TKval,  CSD->nTp );
   ff.readArray( CSD->Pval,  CSD->nPp );

   ff.readArray( CSD->ccIC, CSD->nIC );
   ff.readArray( CSD->ccDC, CSD->nDC );
   ff.readArray( CSD->ccPH, CSD->nPH );

   if( CSD->ccPH[0] == PH_AQUEL )
   {
	 ff.readArray( CSD->denW,  5*CSD->nPp*CSD->nTp );
	 ff.readArray( CSD->denWg,  5*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->epsW, 5*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->epsWg, 5*CSD->nPp*CSD->nTp );
   }
   ff.readArray( CSD->G0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.readArray( CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->S0, CSD->nDC*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->Cp0, CSD->nDC*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->A0, CSD->nDC*CSD->nPp*CSD->nTp );
     ff.readArray( CSD->U0, CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd  )
     ff.readArray( CSD->DD, CSD->nDCs*CSD->nPp*CSD->nTp );

   ff.readArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.readArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.readArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );


}

// allocating DataCH structure
void TNode::datach_realloc()
{
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
  CSD->Pval = new double[CSD->nPp];

  CSD->denW = new double[ 5*CSD->nPp*CSD->nTp];
  CSD->denWg = new double[ 5*CSD->nPp*CSD->nTp];
  CSD->epsW = new double[ 5*CSD->nPp*CSD->nTp];
  CSD->epsWg = new double[ 5*CSD->nPp*CSD->nTp];

  CSD->G0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->V0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->H0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->S0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->Cp0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->A0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->U0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];

  if(  CSD->iGrd  )
       CSD->DD = new double[CSD->nDCs*CSD->nPp*CSD->nTp];
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

// writing DataBR to binary file
void TNode::databr_to_file( GemDataStream& ff )
{
// const data
   ff.writeArray( &CNode->NodeHandle, 6 );
   ff.writeArray( &CNode->TK, 32 );

//dynamic data
   ff.writeArray( CNode->bIC, CSD->nICb );
   ff.writeArray( CNode->rMB, CSD->nICb );
   ff.writeArray( CNode->uIC, CSD->nICb );

   ff.writeArray( CNode->xDC, CSD->nDCb );
   ff.writeArray( CNode->gam, CSD->nDCb );
   ff.writeArray( CNode->dul, CSD->nDCb );
   ff.writeArray( CNode->dll, CSD->nDCb );

   if( CSD->nAalp >0 )
        ff.writeArray( CNode->aPH, CSD->nPHb );
   ff.writeArray( CNode->xPH, CSD->nPHb );
   ff.writeArray( CNode->vPS, CSD->nPSb );
   ff.writeArray( CNode->mPS, CSD->nPSb );
   ff.writeArray( CNode->bPS, CSD->nPSb*CSD->nICb );
   ff.writeArray( CNode->xPA, CSD->nPSb );

   CNode->dRes1 = 0;
//   datach_to_text_file();
//   databr_to_text_file();
}

// Reading work dataBR structure from binary file
void TNode::databr_from_file( GemDataStream& ff )
{
// const data
   ff.readArray( &CNode->NodeHandle, 6 );
   ff.readArray( &CNode->TK, 32 );

//dynamic data
   ff.readArray( CNode->bIC, CSD->nICb );
   ff.readArray( CNode->rMB, CSD->nICb );
   ff.readArray( CNode->uIC, CSD->nICb );

   ff.readArray( CNode->xDC, CSD->nDCb );
   ff.readArray( CNode->gam, CSD->nDCb );
   ff.readArray( CNode->dul, CSD->nDCb );
   ff.readArray( CNode->dll, CSD->nDCb );

   if( CSD->nAalp >0 )
        ff.readArray( CNode->aPH, CSD->nPHb );
   ff.readArray( CNode->xPH, CSD->nPHb );
   ff.readArray( CNode->vPS, CSD->nPSb );
   ff.readArray( CNode->mPS, CSD->nPSb );
   ff.readArray( CNode->bPS, CSD->nPSb*CSD->nICb );
   ff.readArray( CNode->xPA, CSD->nPSb );

   CNode->dRes1 = 0;

}

// Allocates DataBR structure
void TNode::databr_realloc()
{
  long int j;
  CNode->bIC = new double[CSD->nICb];
  CNode->rMB = new double[CSD->nICb];
  CNode->uIC = new double[CSD->nICb];

 CNode->xDC = new double[CSD->nDCb];
 CNode->gam = new double[CSD->nDCb];

for(  j=0; j<CSD->nDCb; j++ )
   CNode->gam[j] = 1.;               //  default assignment
 CNode->dul = new double[CSD->nDCb];
for(  j=0; j<CSD->nDCb; j++ )
   CNode->dul[j] = 1.0e6;            // default assignment
 CNode->dll = new double[CSD->nDCb];
for(  j=0; j<CSD->nDCb; j++ )
   CNode->dll[j] = 0.0;              // default assignment

 if( CSD->nAalp >0 )
 {
    CNode->aPH = new double[CSD->nPHb];
    for( long int k=0; k<CSD->nPHb; k++ )
      CNode->aPH[k] = 0.0;       // default assignment
 }
// else
//    CNode->aPH = 0;

 CNode->xPH = new double[CSD->nPHb];
 CNode->vPS = new double[CSD->nPSb];
 CNode->mPS = new double[CSD->nPSb];
 CNode->bPS = new double[CSD->nPSb*CSD->nICb];
 CNode->xPA = new double[CSD->nPSb];

 CNode->dRes1 = 0;
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
   CNode->dRes1 = 0;
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
	CSD->uRes1 = 0;
// Lists = 0; vectors and matrices
	CSD->nDCinPH = 0;
	CSD->xic = 0;
	CSD->xdc = 0;
	CSD->xph = 0;  //18

	CSD->TKval = 0;
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
//-----------------------End of node_format.cpp--------------------------
