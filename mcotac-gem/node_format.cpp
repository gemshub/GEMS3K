//--------------------------------------------------------------------
// $Id: node.cpp 684 2005-11-23 13:17:15Z gems $
//
// C/C++ interface between GEM IPM and FMT node array
// Working whith DATACH and DATABR structures
//
// Copyright (C) 2004-2005 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#include <iomanip>
#include  <iostream>
#include "io_arrays.h"
#include "node.h"
#include "gdatastream.h"

extern bool _comment;
//===============================================================

outField DataBR_fields[51] =  {
  { "NodeHandle",  0, 0 },
  { "NodeTypeHY",  0, 0 },
  { "NodeTypeMT",  0, 0 },
  { "NodeStatusFMT",  0, 0 },
  { "NodeStatusCH",  1, 0 },
  { "IterDone",  0, 0 },
  { "T",   1, 0 },
  { "P",  1, 0 },
  { "Vs",  1, 0 },
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
  { "dll",    1, 0 },
  { "dul",    1, 0 },
  { "aPH",    1, 0 },
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

outField DataCH_dynamic_fields[25] =  {
   { "xIC",  1, 0 },
   { "xDC",  1, 0 },
   { "xPH",  1, 0 },
   { "ICNL",  1, 0 },
   { "ccIC",  1, 0 },
   { "ICmm",  1, 0 },
   { "DCNL",  1, 0 },
   { "ccDC",  1, 0 },
   { "DCmm",  1, 0 },
   { "PHNL",  1, 0 },
   { "ccPH",  1, 0 },
   { "nDCinPH",  1, 0 },
   { "A",  1, 0 },
   { "Ttol",  0, 0 },
   { "Tval",  1, 0 },
   { "Ptol",  0, 0 },
   { "Pval",  1, 0 },
   { "roW",  1, 0 },
   { "epsW",  1, 0 },
   { "V0",  1, 0 },
   { "G0",  1, 0 },
   { "H0", 1, 0 },
   { "S0",  1, 0 },
   { "Cp0",  1, 0 },
   { "DD",  0, 0 }
};

//===============================================================

void TNode::databr_to_text_file( fstream& ff )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

  TPrintArrays  prar(ff);

   if( _comment )
   {  ff << "# GEMIPM2K v. 0.98" << endl;
      ff << "# Prototype 12.12.2006" << endl;
      ff << "# Comments can be marked with # $ ;" << endl << endl;
      ff << "# Template for the dbr-dat text input file for DATABR (node) data" << endl;
      ff << "# (should be read only after the DATACH and the IPM-DAT files)" << endl << endl;
      ff << "#Section (scalar-1): Controls of the GEM IPM operation and data exchange" << endl;
   }
   if( _comment )
      ff << "# Node identification handle" << endl;
   ff << left << setw(17) << "<NodeHandle> " <<  CNode->NodeHandle << endl;
   if( _comment )
      ff << "# Node type (hydraulic); see typedef NODETYPE" << endl;
   ff << left << setw(17) << "<NodeTypeHY> " <<  CNode->NodeTypeHY << endl;
   if( _comment )
      ff << "# Node type (mass transport); see typedef NODETYPE" << endl;
   ff << left << setw(17) << "<NodeTypeMT> " <<  CNode->NodeTypeMT << endl;
   if( _comment )
      ff << "# Node status code FMT; see typedef NODECODEFMT" << endl;
   ff << left << setw(17) << "<NodeStatusFMT> " <<  CNode->NodeStatusFMT << endl;
   if( _comment )
      ff << "# Node status code CH;  see typedef NODECODECH" << endl;
   ff << left << setw(17) << "<NodeStatusCH> " <<  CNode->NodeStatusCH << endl;
   if( _comment )
      ff << "# Number of iterations performed by IPM (output)" << endl;
   ff << left << setw(17) << "<IterDone> " <<  CNode->IterDone << endl;
   ff << endl;
   if( _comment )
      ff << "##Section (scalar-2): Chemical scalar variables" << endl;
   if( _comment )
         ff << "# Temperature T, K" << endl;
   ff << left << setw(7) << "<T> " <<  CNode->T << endl;
   if( _comment )
         ff << "# Pressure P, bar" << endl;
   ff << left << setw(7) << "<P> " <<  CNode->P << endl;
   if( _comment )
         ff << "# Volume V of reactive subsystem, m3 (GEM output)" << endl;
   ff << left << setw(7) << "<Vs> " << CNode->Vs << endl;
   if( _comment )
         ff << "# Volume Vi of inert subsystem, m3" << endl;
   ff << left << setw(7) << "<Vi> " <<  CNode->Vi << endl;
   if( _comment )
         ff << "# Mass Ms of reactive subsystem,  kg " << endl;
   ff << left << setw(7) << "<Ms> " <<  CNode->Ms << endl;
   if( _comment )
         ff << "# Mass Mi of inert subsystem, kg" << endl;
   ff << left << setw(7) << "<Mi> " <<  CNode->Mi << endl;
   if( _comment )
         ff << "# Enthalpy Hs of reactive subsystem, J " << endl;
   ff << left << setw(7) << "<Hs> " <<  CNode->Hs << endl;
   if( _comment )
         ff << "# Enthalpy Hi of inert subsystem, J " << endl;
   ff << left << setw(7) << "<Hi> " <<  CNode->Hi << endl;
   if( _comment )
         ff << "# Gibbs energy Gs of reactive subsystem, J" << endl;
   ff << left << setw(7) << "<Gs> " <<  CNode->Gs << endl;
   if( _comment )
         ff << "# Effective aqueous ionic strength IS, molal" << endl;
   ff << left << setw(7) << "<IS> " <<  CNode->IC << endl;
   if( _comment )
         ff << "# pH of aqueous solution " << endl;
   ff << left << setw(7) << "<pH> " <<  CNode->pH << endl;
   if( _comment )
         ff << "# pe of aqueous solution" << endl;
   ff << left << setw(7) << "<pe> " <<  CNode->pe << endl;
   if( _comment )
         ff << "# Eh of aqueous solution, V" << endl;
   ff << left << setw(7) << "<Eh> " <<  CNode->Eh << endl;
   ff << endl;
   if( _comment )
       ff << "## FMT scalar variables (used only on the level of NodeArray)" << endl;
   if( _comment )
       ff << "# actual total simulation time Tm, s" << endl;
   ff << left << setw(7) << "<Tm> " <<  CNode->Tm << endl;
   if( _comment )
       ff << "# actual time step dt, s" << endl;
   ff << left << setw(7) << "<dt> " <<  CNode->dt << endl;
   if( _comment )
       ff << "# General diffusivity Dif of disolved matter in the mode, m2/s" << endl;
   ff << left << setw(7) << "<Dif> " <<  CNode->Dif << endl;
   if( _comment )
       ff << "# total volume Vt of the node (voxel), m3" << endl;
   ff << left << setw(7) << "<Vt> " <<  CNode->Vt << endl;
   if( _comment )
       ff << "# advection velocity vp in this node, m/s" << endl;
   ff << left << setw(7) << "<vp> " <<  CNode->vp << endl;
   if( _comment )
       ff << "#  effective (actual) porosity eps, normalized to 1" << endl;
   ff << left << setw(7) << "<eps> " <<  CNode->eps << endl;
   if( _comment )
       ff << "# actual permeability Km, m2" << endl;
   ff << left << setw(7) << "<Km> " <<  CNode->Km << endl;
   if( _comment )
       ff << "# actual DARCY`s constant Kf, m2/s" << endl;
   ff << left << setw(7) << "<Kf> " <<  CNode->Kf << endl;
   if( _comment )
       ff << "# specific storage coefficient S, dimensionless" << endl;
   ff << left << setw(7) << "<S> " <<  CNode->S << endl;
   if( _comment )
       ff << "# transmissivity Tr, m2/s" << endl;
   ff << left << setw(7) << "<Tr> " <<  CNode->Tr << endl;
   if( _comment )
       ff << "# actual hydraulic head h (hydraulic potential), m" << endl;
   ff << left << setw(7) << "<h> " <<  CNode->h << endl;
   if( _comment )
       ff << "# actual carrier density rho for density-driven flow, kg/m3" << endl;
   ff << left << setw(7) << "<rho> " <<  CNode->rho << endl;
   if( _comment )
       ff << "# specific longitudinal dispersivity al of porous media, m" << endl;;
   ff << left << setw(7) << "<al> " <<  CNode->al << endl;
   if( _comment )
       ff << "# specific transversal dispersivity at of porous media, m" << endl;;
   ff << left << setw(7) << "<at> " <<  CNode->at << endl;
   if( _comment )
       ff << "# specific vertical dispersivity av of porous media, m" << endl;;
   ff << left << setw(7) << "<av> " <<  CNode->av << endl;
   if( _comment )
       ff << "# hydraulic longitudinal dispersivity hDl, m2/s" << endl;;
   ff << left << setw(7) << "<hDl> " <<  CNode->hDl << endl;
   if( _comment )
       ff << "# hydraulic transversal dispersivity hDt, m2/s" << endl;;
   ff << left << setw(7) << "<hDt> " <<  CNode->hDt << endl;
   if( _comment )
       ff << "# hydraulic vertical dispersivity hDv, m2/s" << endl;;
   ff << left << setw(7) << "<hDv> " <<  CNode->hDv << endl;
   if( _comment )
       ff << "# tortuosity factor nto, dimensionless" << endl;;
   ff << left << setw(7) << "<nto> " <<  CNode->nto << endl;
   ff << endl;

   if( _comment )
   {   ff << "### Arrays - for dimensions and index lists, see Section (2) of DATACH file" << endl << endl;
       ff << "## IC data section" << endl;
       ff << "# Bulk composition bIC of the reactive part of the node (GEM input, moles)";
   }
  prar.writeArray(  "bIC",  CNode->bIC, CSD->nICb );
   if( _comment )
       ff << "\n\n# Mass balance residuals rMB of GEM solution (GEM output, moles)";
  prar.writeArray(  "rMB",  CNode->rMB, CSD->nICb );
   if( _comment )
       ff << "\n\n# Dual chemical potentials uIC (GEM output, normalized)";
  prar.writeArray(  "uIC",  CNode->uIC, CSD->nICb );
   if( _comment )
   {    ff << "\n\n## DC data section" << endl;
        ff << "# Speciation xDC (amounts of DCs in equilibrium state) - GEM output, moles";
   }
  prar.writeArray(  "xDC",  CNode->xDC, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Activity coefficients gam of Dependent Components, GEM output";
  prar.writeArray(  "gam",  CNode->gam, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Lower metastability constraints dll on amounts in xDC (GEM input, moles)";
  prar.writeArray(  "dll",  CNode->dll, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Upper metastability constraints dul on amounts in xDC (GEM input, moles)";
  prar.writeArray(  "dul",  CNode->dul, CSD->nDCb );
   if( _comment )
   {    ff << "\n\n## Phase data section" << endl;
        ff << "# Specific surface areas of phases aPH (m2/g) - GEM input";
   }
  prar.writeArray(  "aPH",  CNode->aPH, CSD->nPHb );
   if( _comment )
        ff << "\n\n# Amounts of phases in equilibrium state xPH (GEM output, moles)";
  prar.writeArray(  "xPH",  CNode->xPH, CSD->nPHb );
   if( _comment )
        ff << "\n\n# Volumes of the multicomponent phases vPS (cm3), GEM output";
  prar.writeArray(  "vPS",  CNode->vPS, CSD->nPSb );
   if( _comment )
        ff << "\n\n# Masses of the multicomponent phases mPS (g), GEM output";
  prar.writeArray(  "mPS",  CNode->mPS, CSD->nPSb );
   if( _comment )
        ff << "\n\n# Bulk elemental compositions of multicomponent phases bPS (GEM output, moles)";
  prar.writeArray(  "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
   if( _comment )
        ff << "\n\n# Amounts of carrier xPA (sorbent or solvent) in multicomponent phases";
  prar.writeArray(  "xPA",  CNode->xPA, CSD->nPSb );
   if( _comment )
   {     ff << "\n\n# reserved" << endl;
         ff << "\n# End of file"<< endl;
   }
}

// Reading work dataBR structure from text file
void TNode::databr_from_text_file( fstream& ff )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

 memset( &CNode->Tm, 0, 19*sizeof(double));
 TReadArrays  rdar( 51, DataBR_fields, ff);
 short nfild = rdar.findNext();
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
    case 6: rdar.readArray( "T",  &CNode->T, 1);
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
  { ret += " - fields must be readed from DataBR structure";
    Error( "Error", ret);
  }
}

void TNode::datach_to_text_file( fstream& ff )
{
// fstream ff("DataCH.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

  TPrintArrays  prar(ff);

  if( _comment )
  {  ff << "# GEMIPM2K v. 0.98" << endl;
     ff << "# Prototype 12.12.2006" << endl;
     ff << "# Comments are marked with # $ ;" << endl;
     ff << "\n# Template for the dch-dat text input file for DATACH data " << endl;
     ff << "# (should be read first, before the IPM-DAT file and DATABR files)" << endl;
     ff << "\n## (1) Dimensions for memory allocation" << endl;
     ff << "# nIC: Number of Independent Components"<< endl;
  }
  ff << left << setw(7) << "<nIC> " <<  CSD->nIC << endl;
  if( _comment )
     ff << "# nDC: Number of Dependent Components" << endl;
  ff << left << setw(7) << "<nDC> " <<  CSD->nDC << endl;
  if( _comment )
     ff << "# nPH: Number of Phases" << endl;
  ff << left << setw(7) << "<nPH> " <<  CSD->nPH << endl;
  if( _comment )
     ff << "# nPS: Number of Phases-solutions (<= nPH)" << endl;
  ff << left << setw(7) << "<nPS> " <<  CSD->nPS << endl;
  if( _comment )
     ff << "# nDCs: Number of Dependent Components in Phases-solutions" << endl;
  ff << left << setw(7) << "<nDCs> " <<  CSD->nDCs << endl;

  if( _comment )
  {  ff << "\n## (2) Databridge configuration section (for memory allocation)" << endl;
     ff << "# nICb: number of ICs to be kept in DATABR structure (<= nIC)" << endl;
  }
  ff << left << setw(7) << "<nICb> " <<  CSD->nICb << endl;
  if( _comment )
     ff << "# nDCb: number of DCs to be kept in DATABR structure (<=nDC)" << endl;
  ff << left << setw(7) << "<nDCb> " <<  CSD->nDCb << endl;
  if( _comment )
     ff << "# nPHb: number of Phases to be kept in DATABR structure (<=nPH)" << endl;
  ff << left << setw(7) << "<nPHb> " <<  CSD->nPHb << endl;
  if( _comment )
     ff << "# nPSb: number of Phases-solutions to be kept in DATABR structure" << endl;
  ff << left << setw(7) << "<nPSb> " <<  CSD->nPSb << endl;

  if( _comment )
  {   ff << "\n## (3) Dimensions for thermodynamic data arrays" << endl;
      ff << "# nTp: Number of temperature points in the interpolation grid array" << endl;
  }
  ff << left << setw(7) << "<nTp> " <<  CSD->nTp << endl;
  if( _comment )
     ff << "# nPp: Number of pressure points in the interpolation grid array" << endl;
  ff << left << setw(7) << "<nPp> " <<  CSD->nPp << endl;
  if( _comment )
   {  ff << "# iGrd: flag for DC array setup: 0 - only V0 and G0; 1 - plus H0; 2 - plus S0; 3 - plus Cp0;" << endl;
      ff << "# 4 - plus A0 (Helmholtz)" << endl;
   }
  ff << left << setw(7) << "<iGrd> " <<  CSD->iGrd << endl;
  if( _comment )
    ff << "# fAalp: Flag for keeping specific surface areas in DATABR structures/files" << endl;
  ff << left << setw(7) << "<fAalp> " <<  CSD->nAalp << endl;

  ff<< "\n<END_DIM>\n";

// dynamic arrays - must follow after static data
  if( _comment )
  {   ff << "\n## (4) Databridge configuration section (for memory allocation)";
      ff << "\n# xIC: indexes of ICs to be kept in DATABR structure";
  }
 prar.writeArray(  "xIC", CSD->xIC, CSD->nICb);
  if( _comment )
    ff << "\n# xDC: indexes of DCs to be kept in DATABR structure";
 prar.writeArray(  "xDC", CSD->xDC, CSD->nDCb);
  if( _comment )
    ff << "\n# xPH: indexes of Phases to be kept in DATABR structure";
 prar.writeArray(  "xPH", CSD->xPH, CSD->nPHb);

  if( _comment )
     ff << "\n\n## (5) Independent components section";
         ff << "\n# ICNL: List of names of Independent Components";
 prar.writeArray(  "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
   if( _comment )
    ff << "\n# ccIC: List of class codes for Independent Components";
 prar.writeArray(  "ccIC", CSD->ccIC, CSD->nIC, 1 );
  if( _comment )
    ff << "\n# ICmm: Atomic (molar) masses of Independent Components, g/mol";
 prar.writeArray(  "ICmm", CSD->ICmm, CSD->nIC);

  if( _comment )
    ff << "\n\n## (6) Dependent components section (codes and names)";
         ff << "\n# DCNL: List of names of Dependent Components";
 prar.writeArray(  "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
  if( _comment )
    ff << "\n# ccDC: class codes of Dependent Components";
 prar.writeArray(  "ccDC", CSD->ccDC, CSD->nDC, 1 );
 if( _comment )
   ff << "\n\n# DCmm: Molar masses of DCs ";
 prar.writeArray(  "DCmm", CSD->DCmm, CSD->nDC);

  if( _comment )
  {  ff << "\n\n## (7) Phases section" << endl;
     ff << "# PHNL: Phase name list (without a g s ...)";
  }
 prar.writeArray(  "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
  if( _comment )
    ff << "\n# ccPH: Phase aggregate state code list";
prar.writeArray(  "ccPH", CSD->ccPH, CSD->nPH, 1 );
  if( _comment )
    ff << "\n# nDCinPH: Vector L1 telling how many DCs is included in each phase";
prar.writeArray(  "nDCinPH", CSD->nDCinPH, CSD->nPH);

  if( _comment )
  {  ff << "\n\n# (8) Data section for DCs";
     ff << "\n# A: Stoichiometry matrix for DCs - one column per IC, row per DC";
  }
 prar.writeArray(  "A", CSD->A, CSD->nDC*CSD->nIC, CSD->nIC );

  if( _comment )
  {  ff << "\n\n## (8) Thermodynamic data section";
     ff << "\n# Ttol: Tolerance for the interpolation over temperature (K)" << endl;
  }
  ff << left << setw(7) << "<Ttol> " <<  CSD->Ttol;
  if( _comment )
    ff << "\n# Tval: Grid temperatures for the interpolation";
 prar.writeArray(  "Tval", CSD->Tval, CSD->nTp );
  if( _comment )
    ff << "\n\n# Ptol: Tolerance for the interpolation over pressure (K)" << endl;
  ff << left << setw(7) << "<Ptol> " <<  CSD->Ptol;
  if( _comment )
      ff << "\n# Grid pressures for the interpolation";
 prar.writeArray(  "Pval", CSD->Pval, CSD->nPp );

  if( CSD->ccPH[0] == PH_AQUEL )
  { if( _comment )
      ff << "\n\n# roW: Grid array for the density of water-solvent (g/cm3)";
   prar.writeArray(  "roW", CSD->roW, CSD->nPp*CSD->nTp );
    if( _comment )
      ff << "\n\n# epsW: grid array for the dielectric constant of water-solvent";
   prar.writeArray(  "epsW", CSD->epsW,  CSD->nPp*CSD->nTp );
  }
  if( _comment )
    ff << "\n\n# V0: Grid array for the molar volumes of Dependent Components (J/bar)";
 prar.writeArray(  "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp,
                                       CSD->nPp*CSD->nTp );
  if( _comment )
     ff << "\n\n# G0: Grid array for DC molar Gibbs energy function (J/mol)";
 prar.writeArray(  "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp,
                                 CSD->nPp*CSD->nTp );

  if( CSD->iGrd > 0 )
  {
    if( _comment )
      ff << "\n\n# H0: Grid array for DC molar enthalpy function (J/mol)";
   prar.writeArray(  "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp );
  }
  if( CSD->iGrd > 1 )
  {
    if( _comment )
      ff << "\n\n# S0: Grid array for DC absolute entropy function (J/K/mol)";
   prar.writeArray(  "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp  );
  }
  if( CSD->iGrd > 2 )
  {
     if( _comment )
      ff << "\n\n# Cp0: Grid array for DC heat capacity function (J/K/mol)";
    prar.writeArray(  "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp  );
  }
  if( _comment )
    ff << "\n\n# DD: Diffusion coefficients for DCs (reserved)";
  prar.writeArray(  "DD", CSD->DD, CSD->nDCs);

  if( _comment )
      ff << "\n\n# End of file";
}

// Reading dataCH structure from text file
void TNode::datach_from_text_file(fstream& ff)
{
  int ii;
// fstream ff("DataCH.out", ios::in );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

// static arrays
 TReadArrays  rdar( 13, DataCH_static_fields, ff);
 short nfild = rdar.findNext();
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
  { ret += " - fields must be readed from DataCH structure";
    Error( "Error", ret);
  }

  datach_realloc();
  databr_realloc();


//dynamic data
 TReadArrays  rddar( 25, DataCH_dynamic_fields, ff);

// Set up flags
   if( CSD->ccPH[0] == PH_AQUEL )
   {
      rddar.setNoAlws( 17 /*"roW"*/);
      rddar.setNoAlws( 18 /*"epsW"*/);
   }
   if( CSD->iGrd <= 2 )
      rddar.setNoAlws( 23 /*"Cp0"*/);
   if( CSD->iGrd <= 0 )
      rddar.setNoAlws( 21 /*"H0"*/);
   if( CSD->iGrd <= 1 )
      rddar.setNoAlws( 22 /*"S0"*/);

// default set up
  for( ii=0; ii< CSD->nDCs; ii++ )
    CSD->DD[ii] = 0.;
  CSD->Ttol = 0.1;
  CSD->Ptol = 0.1;
  if( CSD->nIC == CSD->nICb )
  {
    rddar.setNoAlws( "xIC");
    for( ii=0; ii< CSD->nICb; ii++ )
      CSD->xIC[ii] = (short)ii;
  }
  if(CSD->nDC == CSD->nDCb )
  {
    rddar.setNoAlws( 1 /*"xDC"*/);
    for( ii=0; ii< CSD->nDCb; ii++ )
      CSD->xDC[ii] = (short)ii;
  }
  if(CSD->nPH == CSD->nPHb )
  {
    rddar.setNoAlws( 2 /*"xPH"*/);
    for( ii=0; ii< CSD->nPHb; ii++ )
      CSD->xPH[ii] = (short)ii;
  }

  nfild = rddar.findNext();
  while( nfild >=0 )
  {
   switch( nfild )
   {
    case 0: rddar.readArray( "xIC", CSD->xIC, CSD->nICb);
            break;
    case 1: rddar.readArray( "xDC", CSD->xDC, CSD->nDCb);
            break;
    case 2: rddar.readArray( "xPH", CSD->xPH, CSD->nPHb);
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
    case 14: rddar.readArray( "Tval", CSD->Tval, CSD->nTp );
            break;
    case 15: rddar.readArray( "Ptol", &CSD->Ptol, 1);
            break;
    case 16: rddar.readArray( "Pval", CSD->Pval, CSD->nPp );
              break;
    case 17: if( !CSD->roW )
                   Error( "Error", "Array roW not used in task");
             rddar.readArray( "roW", CSD->roW, CSD->nPp*CSD->nTp );
              break;
    case 18: if( !CSD->epsW )
                   Error( "Error", "Array epsW not used in task");
             rddar.readArray( "epsW", CSD->epsW,  CSD->nPp*CSD->nTp );
            break;
    case 19: rddar.readArray( "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
            break;
    case 20: rddar.readArray( "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp );
              break;
    case 21: if( !CSD->H0 )
                   Error( "Error", "Array HO not used in task");
            rddar.readArray( "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp);
            break;
    case 22: if( !CSD->S0 )
                   Error( "Error", "Array SO not used in task");
            rddar.readArray( "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp);
            break;
    case 23: if( !CSD->Cp0 )
                   Error( "Error", "Array CpO not used in task");
            rddar.readArray( "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp );
            break;
   case 24: rddar.readArray( "DD", CSD->DD, CSD->nDCs);
              break;
  }
     nfild = rddar.findNext();
 }

 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be readed from DataCH structure";
    Error( "Error", ret);
  }

}

//---------------------------------------------------------------

// new structures i/o

// Writting DataCH to binary file
void TNode::datach_to_file( GemDataStream& ff )
{
// const data
   ff.writeArray( &CSD->nIC, 14 );
   ff.writeArray( &CSD->Ttol, 4 );

//dynamic data
   ff.writeArray( CSD->nDCinPH, CSD->nPH );
//   if( CSD->nICb >0 )
   ff.writeArray( CSD->xIC, CSD->nICb );
   ff.writeArray( CSD->xDC, CSD->nDCb );
   ff.writeArray( CSD->xPH, CSD->nPHb );

   ff.writeArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.writeArray( CSD->ICmm, CSD->nIC );
   ff.writeArray( CSD->DCmm, CSD->nDC );
   ff.writeArray( CSD->DD, CSD->nDCs );

   ff.writeArray( CSD->Tval,  CSD->nTp );
   ff.writeArray( CSD->Pval,  CSD->nPp );

   ff.writeArray( CSD->ccIC, CSD->nIC*sizeof(char) );
   ff.writeArray( CSD->ccDC, CSD->nDC*sizeof(char) );
   ff.writeArray( CSD->ccPH, CSD->nPH*sizeof(char) );

   if( CSD->ccPH[0] == PH_AQUEL )
   { ff.writeArray( CSD->roW,  CSD->nPp*CSD->nTp );
     ff.writeArray( CSD->epsW, CSD->nPp*CSD->nTp );
   }
   ff.writeArray( CSD->G0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.writeArray( CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 0 )
      ff.writeArray( CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 1 )
      ff.writeArray( CSD->S0, CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 2 )
      ff.writeArray( CSD->Cp0, CSD->nDC*CSD->nPp*CSD->nTp );

   ff.writeArray( (char *)CSD->ICNL, MaxICN*CSD->nIC*sizeof(char) );
   ff.writeArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC*sizeof(char) );
   ff.writeArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH*sizeof(char) );
}

// Reading dataCH structure from binary file
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
   ff.readArray( CSD->xIC, CSD->nICb );
   ff.readArray( CSD->xDC, CSD->nDCb );
   ff.readArray( CSD->xPH, CSD->nPHb );

   ff.readArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.readArray( CSD->ICmm, CSD->nIC );
   ff.readArray( CSD->DCmm, CSD->nDC );
   ff.readArray( CSD->DD, CSD->nDCs );

   ff.readArray( CSD->Tval,  CSD->nTp );
   ff.readArray( CSD->Pval,  CSD->nPp );

   ff.readArray( CSD->ccIC, CSD->nIC*sizeof(char) );
   ff.readArray( CSD->ccDC, CSD->nDC*sizeof(char) );
   ff.readArray( CSD->ccPH, CSD->nPH*sizeof(char) );

   if( CSD->ccPH[0] == PH_AQUEL )
   {  ff.readArray( CSD->roW,  CSD->nPp*CSD->nTp );
      ff.readArray( CSD->epsW, CSD->nPp*CSD->nTp );
   }
   ff.readArray( CSD->G0,  CSD->nDC*CSD->nPp*CSD->nTp );
   ff.readArray( CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 0 )
     ff.readArray( CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 1 )
     ff.readArray( CSD->S0, CSD->nDC*CSD->nPp*CSD->nTp );
   if(  CSD->iGrd > 2 )
     ff.readArray( CSD->Cp0, CSD->nDC*CSD->nPp*CSD->nTp );

   ff.readArray( (char *)CSD->ICNL, MaxICN*CSD->nIC*sizeof(char) );
   ff.readArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC*sizeof(char) );
   ff.readArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH*sizeof(char) );


}

// allocate DataCH structure
void TNode::datach_realloc()
{
 CSD->nDCinPH = new short[CSD->nPH];

 if( CSD->nICb >0 )
   CSD->xIC = new short[CSD->nICb];
 else  CSD->xIC = 0;
 if( CSD->nDCb >0 )
   CSD->xDC = new short[CSD->nDCb];
 else  CSD->xDC = 0;
 if( CSD->nPHb >0 )
   CSD->xPH = new short[CSD->nPHb];
 else  CSD->xPH = 0;

  CSD->A = new float[CSD->nIC*CSD->nDC];
  CSD->ICmm = new double[CSD->nIC];
  CSD->DCmm = new double[CSD->nDC];
  CSD->DD = new double[CSD->nDCs];

  CSD->Tval = new float[CSD->nTp];
  CSD->Pval = new float[CSD->nPp];

  CSD->roW = new double[ CSD->nPp*CSD->nTp];
  CSD->epsW = new double[ CSD->nPp*CSD->nTp];

  CSD->G0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  CSD->V0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];

  if(  CSD->iGrd > 0 )
    CSD->H0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  else
    CSD->H0 = 0;
  if(  CSD->iGrd > 1 )
      CSD->S0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  else
      CSD->S0 = 0;
  if(  CSD->iGrd > 2 )
    CSD->Cp0 = new double[CSD->nDC*CSD->nPp*CSD->nTp];
  else
    CSD->Cp0 = 0;

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
 if( CSD->xIC )
  { delete[] CSD->xIC;
    CSD->xIC = 0;
  }
 if( CSD->xDC )
  { delete[] CSD->xDC;
    CSD->xDC = 0;
  }
 if( CSD->xPH )
  { delete[] CSD->xPH;
    CSD->xPH = 0;
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
 if( CSD->DD )
  { delete[] CSD->DD;
    CSD->DD = 0;
  }

 if( CSD->Tval )
  { delete[] CSD->Tval;
    CSD->Tval = 0;
  }
 if( CSD->Pval )
  { delete[] CSD->Pval;
    CSD->Pval = 0;
  }

 if( CSD->roW )
  { delete[] CSD->roW;
    CSD->roW = 0;
  }
 if( CSD->epsW )
  { delete[] CSD->epsW;
    CSD->epsW = 0;
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
   ff.writeArray( &CNode->T, 32 );

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
   ff.readArray( &CNode->T, 32 );

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

// allocate DataBR structure
void TNode::databr_realloc()
{
  CNode->bIC = new double[CSD->nICb];
  CNode->rMB = new double[CSD->nICb];
  CNode->uIC = new double[CSD->nICb];

 CNode->xDC = new double[CSD->nDCb];
 CNode->gam = new double[CSD->nDCb];
 for( int ii=0; ii<CSD->nDCb; ii++ )
   CNode->gam[ii] = 1.;
 CNode->dul = new double[CSD->nDCb];
 CNode->dll = new double[CSD->nDCb];

 if( CSD->nAalp >0 )
     CNode->aPH = new double[CSD->nPHb];
 else
    CNode->aPH = 0;

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
  memset( &CNode_->NodeHandle, 0, 6*sizeof(short));
  memset( &CNode_->T, 0, 32*sizeof(double));

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

 delete[] CNode_;
 return NULL;
}


//-----------------------End of node_format.cpp--------------------------
