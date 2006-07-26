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
#include "node.h"
#include "gdatastream.h"

extern bool _comment;

void TNode::databr_to_text_file( fstream& ff )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

   if( _comment )
   {  ff << "# GEMIPM2K v. 0.725" << endl;
      ff << "# Prototype 12.07.2006" << endl;
      ff << "# Comments can be marked with #" << endl << endl;
      ff << "# Template for the dbr-dat text input file for DATABR data" << endl;
      ff << "# (should be read only after the DATACH and the IPM-DAT files)" << endl << endl;
      ff << "#'sCon' Node handle and status code" << endl;
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
      ff << "# Number of iterations performed by IPM" << endl;
   ff << left << setw(17) << "<IterDone> " <<  CNode->IterDone << endl;
   ff << endl;
   if( _comment )
      ff << "##'dCon' Chemical scalar variables" << endl;
   if( _comment )
         ff << "# Temperature T, K" << endl;
   ff << left << setw(7) << "<T> " <<  CNode->T << endl;
   if( _comment )
         ff << "# Pressure P, bar" << endl;
   ff << left << setw(7) << "<P> " <<  CNode->P << endl;
   if( _comment )
         ff << "# Volume V of reactive subsystem, cm3" << endl;
   ff << left << setw(7) << "<Vs> " << CNode->Vs << endl;
   if( _comment )
         ff << "# Volume of inert subsystem, cm3" << endl;
   ff << left << setw(7) << "<Vi> " <<  CNode->Vi << endl;
   if( _comment )
         ff << "# Mass of reactive subsystem,  g " << endl;
   ff << left << setw(7) << "<Ms> " <<  CNode->Ms << endl;
   if( _comment )
         ff << "# Mass of inert subsystem, g" << endl;
   ff << left << setw(7) << "<Mi> " <<  CNode->Mi << endl;
   if( _comment )
         ff << "# Gibbs energy of reactive subsystem, J" << endl;
   ff << left << setw(7) << "<Gs> " <<  CNode->Gs << endl;
   if( _comment )
         ff << "# Enthalpy of reactive subsystem, J " << endl;
   ff << left << setw(7) << "<Hs> " <<  CNode->Hs << endl;
   if( _comment )
         ff << "# Enthalpy of inert subsystem, J  " << endl;
   ff << left << setw(7) << "<Hi> " <<  CNode->Hi << endl;
   if( _comment )
         ff << "# Effective aqueous ionic strength, molal" << endl;
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
       ff << "# actual total simulation time, s" << endl;
   ff << left << setw(7) << "<Tm> " <<  CNode->Tm << endl;
   if( _comment )
       ff << "# actual time step" << endl;
   ff << left << setw(7) << "<dt> " <<  CNode->dt << endl;
   if( _comment )
       ff << "# General diffusivity of disolved matter in the mode" << endl;
   ff << left << setw(7) << "<Dif> " <<  CNode->Dif << endl;
   if( _comment )
       ff << "# total volume of the node (voxel), m3" << endl;
   ff << left << setw(7) << "<Vt> " <<  CNode->Vt << endl;
   if( _comment )
       ff << "# advection velocity (in pores) in this node" << endl;
   ff << left << setw(7) << "<vp> " <<  CNode->vp << endl;
   if( _comment )
       ff << "#  effective (actual) porosity normalized to 1" << endl;
   ff << left << setw(7) << "<eps> " <<  CNode->eps << endl;
   if( _comment )
       ff << "# actual permeability, m2" << endl;
   ff << left << setw(7) << "<Km> " <<  CNode->Km << endl;
   if( _comment )
       ff << "# actual DARCY`s constant, m2/s" << endl;
   ff << left << setw(7) << "<Kf> " <<  CNode->Kf << endl;
   if( _comment )
       ff << "# specific storage coefficient, dimensionless" << endl;
   ff << left << setw(7) << "<S> " <<  CNode->S << endl;
   if( _comment )
       ff << "# transmissivity m2/s" << endl;
   ff << left << setw(7) << "<Tr> " <<  CNode->Tr << endl;
   if( _comment )
       ff << "# actual hydraulic head (hydraulic potential), m" << endl;
   ff << left << setw(7) << "<h> " <<  CNode->h << endl;
   if( _comment )
       ff << "# actual carrier density for density-driven flow, g/cm3" << endl;
   ff << left << setw(7) << "<rho> " <<  CNode->rho << endl;
   if( _comment )
       ff << "# specific longitudinal dispersivity of porous media, m" << endl;;
   ff << left << setw(7) << "<al> " <<  CNode->al << endl;
   if( _comment )
       ff << "# specific transversal dispersivity of porous media, m" << endl;;
   ff << left << setw(7) << "<at> " <<  CNode->at << endl;
   if( _comment )
       ff << "# specific vertical dispersivity of porous media, m" << endl;;
   ff << left << setw(7) << "<av> " <<  CNode->av << endl;
   if( _comment )
       ff << "# hydraulic longitudinal dispersivity, m2/s, diffusities from chemical database" << endl;;
   ff << left << setw(7) << "<hDl> " <<  CNode->hDl << endl;
   if( _comment )
       ff << "# hydraulic transversal dispersivity, m2/s" << endl;;
   ff << left << setw(7) << "<hDt> " <<  CNode->hDt << endl;
   if( _comment )
       ff << "# hydraulic vertical dispersivity, m2/s" << endl;;
   ff << left << setw(7) << "<hDv> " <<  CNode->hDv << endl;
   if( _comment )
       ff << "# tortuosity factor" << endl;;
   ff << left << setw(7) << "<nto> " <<  CNode->nto << endl;
   ff << endl;
   if( _comment )
   {   ff << "### Arrays - for dimensions and index lists, see Section (2) of DATACH file" << endl << endl;
       ff << "## IC data section" << endl;
       ff << "# Bulk composition (partial) of the system - GEM input (moles)";
   }
   outArray( ff, "bIC",  CNode->bIC, CSD->nICb );
   if( _comment )
       ff << "\n\n# Mass balance residuals - GEM output (moles)";
   outArray( ff, "rMB",  CNode->rMB, CSD->nICb );
   if( _comment )
       ff << "\n\n# Chemical potentials of ICs (dual GEM solution) - GEM output, normalized";
   outArray( ff, "uIC",  CNode->uIC, CSD->nICb );
   if( _comment )
   {    ff << "\n\n## DC data section" << endl;
        ff << "# Speciation - amount of DCs in equilibrium state - GEM output";
   }
   outArray( ff, "xDC",  CNode->xDC, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Activity coefficients of DCs - GEM output";
   outArray( ff, "gam",  CNode->gam, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Lower metastability constraints on amounts of DCs - GEM input, moles";
   outArray( ff, "dll",  CNode->dll, CSD->nDCb );
   if( _comment )
       ff << "\n\n# Upper metastability constraints on amounts of DCs - GEM input, moles";
   outArray( ff, "dul",  CNode->dul, CSD->nDCb );
   if( _comment )
   {    ff << "\n\n## Phase data section" << endl;
        ff << "# Specific surface areas of phases (m2/g) - GEM input";
   }
   outArray( ff, "aPH",  CNode->aPH, CSD->nPHb );
   if( _comment )
        ff << "\n\n# Amounts of phases in equilibrium state - GEM output (moles)";
   outArray( ff, "xPH",  CNode->xPH, CSD->nPHb );
   if( _comment )
        ff << "\n\n# Volumes of the multicomponent phases (cm3), GEM output";
   outArray( ff, "vPS",  CNode->vPS, CSD->nPSb );
   if( _comment )
        ff << "\n\n# masses of multicomponent phases (g), GEM output";
   outArray( ff, "mPS",  CNode->mPS, CSD->nPSb );
   if( _comment )
        ff << "\n\n# bulk elemental compositions of multicomponent phases - GEM output, moles";
   outArray( ff, "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
   if( _comment )
        ff << "\n\n# amounts of carrier (sorbent or solvent) in multicomponent phases";
   outArray( ff, "xPA",  CNode->xPA, CSD->nPSb );
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

  inArray( ff, "NodeHandle",  &CNode->NodeHandle, 1);
  inArray( ff, "NodeTypeHY",  &CNode->NodeTypeHY, 1);
  inArray( ff, "NodeTypeMT",  &CNode->NodeTypeMT, 1);
  inArray( ff, "NodeStatusFMT",  &CNode->NodeStatusFMT, 1);
  inArray( ff, "NodeStatusCH",  &CNode->NodeStatusCH, 1);
  inArray( ff, "IterDone",  &CNode->IterDone, 1);
  inArray( ff, "T",  &CNode->T, 1);
  inArray( ff, "P",  &CNode->P, 1);
  inArray( ff, "Vs", &CNode->Vs, 1);
  inArray( ff, "Vi",  &CNode->Vi, 1);
  inArray( ff, "Ms",  &CNode->Ms, 1);
  inArray( ff, "Mi",  &CNode->Mi, 1);
  inArray( ff, "Gs",  &CNode->Gs, 1);
  inArray( ff, "Hs",  &CNode->Hs, 1);
  inArray( ff, "Hi",  &CNode->Hi, 1);
  inArray( ff, "IS",  &CNode->IC, 1);
  inArray( ff, "pH",  &CNode->pH, 1);
  inArray( ff, "pe",  &CNode->pe, 1);
  inArray( ff, "Eh",  &CNode->Eh, 1);
  inArray( ff, "Tm",  &CNode->Tm, 1);
  inArray( ff, "dt",  &CNode->dt, 1);
  inArray( ff, "Dif",  &CNode->Dif, 1);
  inArray( ff, "Vt",  &CNode->Vt, 1);
  inArray( ff, "vp",  &CNode->vp, 1);
  inArray( ff, "eps",  &CNode->eps, 1);
  inArray( ff, "Km",  &CNode->Km, 1);
  inArray( ff, "Kf",  &CNode->Kf, 1);
  inArray( ff, "S",  &CNode->S, 1);
  inArray( ff, "Tr",  &CNode->Tr, 1);
  inArray( ff, "h",  &CNode->h, 1);
  inArray( ff, "rho",  &CNode->rho, 1);
  inArray( ff, "al",  &CNode->al, 1);
  inArray( ff, "at",  &CNode->at, 1);
  inArray( ff, "av",  &CNode->av, 1);
  inArray( ff, "hDl",  &CNode->hDl, 1);
  inArray( ff, "hDt",  &CNode->hDt, 1);
  inArray( ff, "hDv",  &CNode->hDv, 1);
  inArray( ff, "nto",  &CNode->nto, 1);

// dynamic arrays
  inArray( ff, "bIC",  CNode->bIC, CSD->nICb );
  inArray( ff, "rMB",  CNode->rMB, CSD->nICb );
  inArray( ff, "uIC",  CNode->uIC, CSD->nICb );

  inArray( ff, "xDC",  CNode->xDC, CSD->nDCb );
  inArray( ff, "gam",  CNode->gam, CSD->nDCb );
  inArray( ff, "dll",  CNode->dll, CSD->nDCb );
  inArray( ff, "dul",  CNode->dul, CSD->nDCb );

  inArray( ff, "aPH",  CNode->aPH, CSD->nPHb );
  inArray( ff, "xPH",  CNode->xPH, CSD->nPHb );
  inArray( ff, "vPS",  CNode->vPS, CSD->nPSb );
  inArray( ff, "mPS",  CNode->mPS, CSD->nPSb );
  inArray( ff, "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
  inArray( ff, "xPA",  CNode->xPA, CSD->nPSb );
}

void TNode::datach_to_text_file( fstream& ff )
{
// fstream ff("DataCH.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

  if( _comment )
  {  ff << "# GEMIPM2K v. 0.725" << endl;
     ff << "# Prototype 12.07.2006" << endl;
     ff << "# Comments marked with #" << endl;
     ff << "\n# Template for the dch-dat text input file for DATACH data " << endl;
     ff << "# (should be read first, before the IPM-DAT file and DATABR files)" << endl;
     ff << "\n## (1) Dimensions for memory allocation" << endl;
     ff << "# Number of ICs"<< endl;
  }
  ff << left << setw(7) << "<nIC> " <<  CSD->nIC << endl;
  if( _comment )
     ff << "# Number of DCs" << endl;
  ff << left << setw(7) << "<nDC> " <<  CSD->nDC << endl;
  if( _comment )
     ff << "# Number of phases" << endl;
  ff << left << setw(7) << "<nPH> " <<  CSD->nPH << endl;
  if( _comment )
     ff << "# Number of phases-solutions" << endl;
  ff << left << setw(7) << "<nPS> " <<  CSD->nPS << endl;
  if( _comment )
     ff << "# Number of DCs in phases-solutions" << endl;
  ff << left << setw(7) << "<nDCs> " <<  CSD->nDCs << endl;

  if( _comment )
  {  ff << "\n## (2) Databridge configuration section (for memory allocation)" << endl;
     ff << "# number of ICs to be kept in DATABR structure" << endl;
  }
  ff << left << setw(7) << "<nICb> " <<  CSD->nICb << endl;
  if( _comment )
     ff << "# number of DCs to be kept in DATABR structure" << endl;
  ff << left << setw(7) << "<nDCb> " <<  CSD->nDCb << endl;
  if( _comment )
     ff << "# number of Phases to be kept in DATABR structure" << endl;
  ff << left << setw(7) << "<nPHb> " <<  CSD->nPHb << endl;
  if( _comment )
     ff << "# number of Phases-solutions to be kept in DATABR structure" << endl;
  ff << left << setw(7) << "<nPSb> " <<  CSD->nPSb << endl;

  if( _comment )
  {   ff << "\n## (3) Thermodynamic data arrays dimensions" << endl;
      ff << "# Number of temperature points in the interpolation grid array" << endl;
  }
  ff << left << setw(7) << "<nTp> " <<  CSD->nTp << endl;
  if( _comment )
     ff << "# Number of pressure points in the interpolation grid array" << endl;
  ff << left << setw(7) << "<nPp> " <<  CSD->nPp << endl;
  if( _comment )
   {  ff << "# flag for DC array setup: 0 - only V0 and G0; 1 - plus H0; 2 - plus S0; 3 - plus Cp0;" << endl;
      ff << "# 4 - plus A0 (Helmholtz)" << endl;
   }
  ff << left << setw(7) << "<iGrd> " <<  CSD->iGrd << endl;
  if( _comment )
    ff << "# Flag for keeping specific surface areas in DATABR structures/files" << endl;
  ff << left << setw(7) << "<fAalp> " <<  CSD->nAalp << endl;

// dynamic arrays
  if( _comment )
  {   ff << "\n## (4) Databridge configuration section (for memory allocation)";
      ff << "\n# number of ICs to be kept in DATABR structure";
  }
  outArray( ff, "xIC", CSD->xIC, CSD->nICb);
  if( _comment )
    ff << "\n# number of DCs to be kept in DATABR structure";
  outArray( ff, "xDC", CSD->xDC, CSD->nDCb);
  if( _comment )
    ff << "\n# number of Phases to be kept in DATABR structure";
  outArray( ff, "xPH", CSD->xPH, CSD->nPHb);

  if( _comment )
     ff << "\n\n## (5) Independent components section";
  outArray( ff, "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
   if( _comment )
    ff << "\n# too much info";
  outArray( ff, "ccIC", CSD->ccIC, CSD->nIC, 1 );
  if( _comment )
    ff << "\n# Atomic (molar) masses of IC, g/mol";
  outArray( ff, "ICmm", CSD->ICmm, CSD->nIC);

  if( _comment )
    ff << "\n\n## (6) Dependent components section (codes and names)";
  outArray( ff, "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
  if( _comment )
    ff << "\n# DC class codes";
  outArray( ff, "ccDC", CSD->ccDC, CSD->nDC, 1 );

  if( _comment )
  {  ff << "\n\n## (7) Phases section" << endl;
     ff << "# Phase name list (without a g s ...)";
  }
  outArray( ff, "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
  if( _comment )
    ff << "\n# Phase aggregate state code list";
 outArray( ff, "ccPH", CSD->ccPH, CSD->nPH, 1 );
  if( _comment )
    ff << "\n# Vector L1 telling how many DCs is included in each phase";
 outArray( ff, "nDCinPH", CSD->nDCinPH, CSD->nPH);

  if( _comment )
  {  ff << "\n\n# (8) Data section for DCs";
     ff << "\n# Stoichiometry matrix for DCs - one column per IC, row per DC";
  }
  outArray( ff, "A", CSD->A, CSD->nDC*CSD->nIC, CSD->nIC );
  if( _comment )
    ff << "\n\n# Molar masses of DCs ";
  outArray( ff, "DCmm", CSD->DCmm, CSD->nDC);
  if( _comment )
    ff << "\n\n# Diffusion coefficients for DCs";
  outArray( ff, "DD", CSD->DD, CSD->nDCs);

  if( _comment )
  {  ff << "\n\n## (8) Thermodynamic data section";
     ff << "\n# Tolerance for the interpolation over temperature (K)" << endl;
  }
  ff << left << setw(7) << "<Ttol> " <<  CSD->Ttol;
  if( _comment )
    ff << "\n# Temperatures for the grid";
  outArray( ff, "Tval", CSD->Tval, CSD->nTp );
  if( _comment )
    ff << "\n\n# Tolerance for the interpolation over pressure (K)" << endl;
  ff << left << setw(7) << "<Ptol> " <<  CSD->Ptol;
  if( _comment )
      ff << "\n# Pressures for the grid";
  outArray( ff, "Pval", CSD->Pval, CSD->nPp );

  if( CSD->ccPH[0] == PH_AQUEL )
  { if( _comment )
      ff << "\n\n# Grid array for density of water-solvent";
    outArray( ff, "roW", CSD->roW, CSD->nPp*CSD->nTp );
    if( _comment )
      ff << "\n\n# Grid array for diel. const. of water-solvent";
    outArray( ff, "epsW", CSD->epsW,  CSD->nPp*CSD->nTp );
  }
  if( _comment )
    ff << "\n\n# Grid array for DC molar volumes (J/bar)";
  outArray( ff, "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp,
                                       CSD->nPp*CSD->nTp );
  if( _comment )
     ff << "\n\n# Grid array for DC molar Gibbs energy function (J/mol)";
  outArray( ff, "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp,
                                 CSD->nPp*CSD->nTp );

  if( CSD->iGrd > 0 )
  {
    if( _comment )
      ff << "\n\n# Grid array for DC molar enthalpy function (J/mol)";
    outArray( ff, "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp );
  }
  if( CSD->iGrd > 1 )
  {
    if( _comment )
      ff << "\n\n# Grid array for DC absolute entropy function (J/mol)";
    outArray( ff, "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp  );
  }
  if( CSD->iGrd > 1 )
  {
     if( _comment )
      ff << "\n\n# Grid array for DC heat capacity function (J/mol)";
     outArray( ff, "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp,
                                        CSD->nPp*CSD->nTp  );
  }
  if( _comment )
      ff << "\n\n# End of file";
}

// Reading dataCH structure from text file
void TNode::datach_from_text_file(fstream& ff)
{
// fstream ff("DataCH.out", ios::in );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

  inArray( ff, "nIC", &CSD->nIC, 1);
  inArray( ff, "nDC", &CSD->nDC, 1);
  inArray( ff, "nPH", &CSD->nPH, 1);
  inArray( ff, "nPS", &CSD->nPS, 1);
  inArray( ff, "nDCs", &CSD->nDCs, 1);
  inArray( ff, "nICb", &CSD->nICb, 1);
  inArray( ff, "nDCb", &CSD->nDCb, 1);
  inArray( ff, "nPHb", &CSD->nPHb, 1);
  inArray( ff, "nPSb", &CSD->nPSb, 1);
  inArray( ff, "nTp", &CSD->nTp, 1);
  inArray( ff, "nPp", &CSD->nPp, 1);
  inArray( ff, "iGrd", &CSD->iGrd, 1);
  inArray( ff, "fAalp", &CSD->nAalp, 1);

  datach_realloc();
  databr_realloc();

//dynamic data

   inArray( ff, "xIC", CSD->xIC, CSD->nICb);
   inArray( ff, "xDC", CSD->xDC, CSD->nDCb);
   inArray( ff, "xPH", CSD->xPH, CSD->nPHb);

   inArray( ff, "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
   inArray( ff, "ccIC", CSD->ccIC, CSD->nIC, 1 );
   inArray( ff, "ICmm", CSD->ICmm, CSD->nIC);

   inArray( ff, "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
   inArray( ff, "ccDC", CSD->ccDC, CSD->nDC, 1 );

   inArray( ff, "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
   inArray( ff, "ccPH", CSD->ccPH, CSD->nPH, 1 );
   inArray( ff, "nDCinPH", CSD->nDCinPH, CSD->nPH);

   inArray( ff, "A", CSD->A, CSD->nDC*CSD->nIC );
   inArray( ff, "DCmm", CSD->DCmm, CSD->nDC);
   inArray( ff, "DD", CSD->DD, CSD->nDCs);

   inArray( ff, "Ttol", &CSD->Ttol, 1);
   inArray( ff, "Tval", CSD->Tval, CSD->nTp );
   inArray( ff, "Ptol", &CSD->Ptol, 1);
   inArray( ff, "Pval", CSD->Pval, CSD->nPp );

   if( CSD->ccPH[0] == PH_AQUEL )
   {
     inArray( ff, "roW", CSD->roW, CSD->nPp*CSD->nTp );
     inArray( ff, "epsW", CSD->epsW,  CSD->nPp*CSD->nTp );
   }
   inArray( ff, "V0", CSD->V0,  CSD->nDC*CSD->nPp*CSD->nTp );
   inArray( ff, "G0", CSD->G0, CSD->nDC*CSD->nPp*CSD->nTp );
   if( CSD->iGrd > 0 )
     inArray( ff, "H0", CSD->H0,  CSD->nDC*CSD->nPp*CSD->nTp);
   if( CSD->iGrd > 1 )
     inArray( ff, "S0", CSD->S0,CSD->nDC*CSD->nPp*CSD->nTp);
   if( CSD->iGrd > 1 )
      inArray( ff, "Cp0", CSD->Cp0,CSD->nDC*CSD->nPp*CSD->nTp );
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

   ff.writeArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.writeArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.writeArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );

   ff.writeArray( CSD->ccIC, CSD->nIC );
   ff.writeArray( CSD->ccDC, CSD->nDC );
   ff.writeArray( CSD->ccPH, CSD->nPH );

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

   ff.readArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.readArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.readArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );

   ff.readArray( CSD->ccIC, CSD->nIC );
   ff.readArray( CSD->ccDC, CSD->nDC );
   ff.readArray( CSD->ccPH, CSD->nPH );

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
