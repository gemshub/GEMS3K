//--------------------------------------------------------------------
// $Id$
//
/// \file datach_formats.cpp
/// Interface for writing/reading DBR and DCH I/O files of GEMS3K
//
// Copyright (c) 2023 S.Dmytriyeva, D.Kulik
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

#include "datach_api.h"
#include "io_template.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "io_keyvalue.h"
#include "gdatastream.h"
#include "v_detail.h"
#include "m_const_base.h"

extern const char* _GEMIPM_version_stamp;

namespace  dbr_dch_api {

//===============================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here. The fourth field contains
// the text of the comment for this data object, optionally written into the
// text-format output DBR or DCH file.
//
std::vector<io_formats::outField> DataBR_fields =  {  // [f_lga+1/*60*/]
  { "NodeHandle",  0, 0, 1, "# NodeHandle: Node identification handle"},
  { "NodeTypeHY",  0, 0, 1, "# NodeTypeHY:  Node type code (hydraulic), not used on TNode level; see typedef NODETYPE" },
  { "NodeTypeMT",  0, 0, 1, "# NodeTypeMT:  Node type (mass transport), not used on TNode level; see typedef NODETYPE" },
  { "NodeStatusFMT",  1, 0, 1, "# NodeStatusFMT:  Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT" },
  { "NodeStatusCH",  1, 0, 1, "# NodeStatusCH: Node status code and control in GEM input and output; see typedef NODECODECH"},
  { "IterDone",  0, 0, 1, "# IterDone:  Number of iterations performed by GEM IPM in the last run (GEM output)" },
  { "TK",   1, 0, 1, "# TK: Node temperature T, Kelvin. This value must always be provided (GEM input)"  },
  { "P",    1, 0, 1, "# P:  Node Pressure P, Pa. This value must always be provided (GEM input)" },
  { "Vs",   0, 0, 1, "# Vs: Volume V of reactive subsystem, m3 (GEM output)" },
  { "Vi",   0, 0, 1, "# Vi: Volume of inert subsystem, m3 (mass transport)" },
  { "Ms",   0, 0, 1, "# Ms: Mass of reactive subsystem, kg (GEM output)" },
  { "Mi",   0, 0, 1, "# Mi: Mass of inert subsystem, kg (mass transport)" },
  { "Hs",   0, 0, 1, "# Hs: Total enthalpy of reactive subsystem, J (reserved)" },
  { "Hi",   0, 0, 1, "# Hi: Total enthalpy of inert subsystem, J (reserved, mass transport) " },
  { "Gs",   0, 0, 1, "# Gs: Total Gibbs energy of the reactive subsystem, J/(RT) (GEM output)" },
  { "IS",   0, 0, 1, "# IS: Effective aqueous ionic strength, molal (GEM output)" },
  { "pH",   0, 0, 1, "# pH: pH of aqueous solution in molal activity scale (GEM output)" },
  { "pe",   0, 0, 1, "# pe: pe of aqueous solution in molal activity scale (GEM output)" },
  { "Eh",   0, 0, 1, "# Eh: Eh of aqueous solution, V (GEM output)" },
  { "Tm",   0, 0, 1, "# Tm: Actual total simulation time, s (kinetics, metastability, transport)" },
  { "dt",   0, 0, 1, "# dt: Actual time step, s (kinetics, metastability, transport)" },
//#ifdef NODEARRAYLEVEL
// Scalar parameters below are only used at TNodeArray level
  { "Dif",  0, 0, 1, "# Dif: General diffusivity of disolved matter, m2/s (mass transport)" },
  { "Vt",   0, 0, 1, "# Vt: Total volume of the node, m3 (mass transport)" },
  { "vp",   0, 0, 1, "# vp: Advection velocity in pores, m/s (mass transport)" },
  { "eps",  0, 0, 1, "# eps: Effective (actual) porosity normalized to 1 (mass transport)" },
  { "Km",   0, 0, 1, "# Km: Actual permeability, m2 (mass transport)" },
  { "Kf",   0, 0, 1, "# Kf: Actual Darcy`s constant, (m2/s (mass transport)" },
  { "S",    0, 0, 1, "# S: Specific storage coefficient, dimensionless (mass transport)" },
  { "Tr",   0, 0, 1, "# Tr: Transmissivity, m2/s (mass transport)" },
  { "h",    0, 0, 1, "# h: Actual hydraulic head (hydraulic potential), m (mass transport)" },
  { "rho",  0, 0, 1, "# rho: Actual carrier density for density-driven flow, kg/m3 (mass transport)" },
  { "al",   0, 0, 1, "# al: Specific longitudinal dispersivity of porous media, m (mass transport)" },
  { "at",   0, 0, 1, "# at: Specific transversal dispersivity of porous media, m (mass transport)" },
  { "av",   0, 0, 1, "# av: Specific vertical dispersivity of porous media, m (mass transport)" },
  { "hDl",  0, 0, 1, "# hDl: Hydraulic longitudinal dispersivity, m2/s (mass transport)" },
  { "hDt",  0, 0, 1, "# hDt: Hydraulic transversal dispersivity, m2/s (mass transport)" },
  { "hDv",  0, 0, 1, "# hDv: Hydraulic vertical dispersivity, m2/s (mass transport)" },
  { "nto",  0, 0, 1, "# nto: Tortuosity factor, dimensionless (mass transport)" },
//#endif
 // dynamic arrays (53-38=15)
  { "bIC",  1, 0, nICbi, "# bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]" },
  { "rMB",  0, 0, nICbi, "\n# rMB: Mass balance residuals, moles (GEM output) [nICb]" },
  { "uIC",  0, 0, nICbi, "\n# uIC: Chemical potentials of ICs in equilibrium (dual solution), J/(RT) (GEM output) [nICb]" },
  { "xDC",  0, 0, nDCbi, "# xDC: Speciation - amounts of DCs in equilibrium (primal solution), moles (GEM output/input) [nDCb]" },
  { "gam",  0, 0, nDCbi, "\n# gam: Activity coefficients of DCs in their respective phases (GEM output/input) [nDCb]" },
  { "dll",  0, 0, nDCbi, "\n# dll: Lower metastability restrictions on amounts of DCs, moles (GEM input) [nDCb]" },
  { "dul",  0, 0, nDCbi, "\n# dul: Upper metastability constraints on amounts of DCs, moles (GEM input) [nDCb]" },
  { "aPH",  0, 0, nPHbi, "# aPH: Specific surface areas of phases, m2/kg (GEM input) [nPHb]" },
  { "xPH",  0, 0, nPHbi, "\n# xPH: Amounts of phases in equilibrium state, moles (GEM output) [nPHb]" },
  { "vPS",  0, 0, nPSbi, "\n# vPS: Volumes of multicomponent phases, m3 (GEM output) [nPSb]" },
  { "mPS",  0, 0, nPSbi, "\n# mPS: Masses of multicomponent phases, kg (GEM output) [nPSb]" },
  { "bPS",  0, 0, nPSbnICbi, "\n# bPS: Bulk elemental compositions of multicomponent phases, moles (GEM output) [nPSb*nICb]"},
  { "xPA",  0, 0, nPSbi, "\n# xPA: Amount of carrier (sorbent or solvent) in multicomponent phases, moles (GEM output) [nPSb]" },
  { "bSP",  0, 0, nICbi, "\n# bSP: Output bulk composition of the equilibrium solid part of the system, moles " },
    { "amru",  0, 0, nPSbi, "\n# amru: Upper AMRs on amounts of multi-component phases (mol) [nPSb]  " },
    { "amrl",  0, 0, nPSbi, "\n# amrl: Lower AMRs on amounts of multi-component phases (mol) [nPSb]" },
  { "omPH",  0, 0, nPHbi, "\n# omPH: stability (saturation) indices of phases in log10 scale, can change in GEM [nPHb] " },

// only for VTK format output
    { "mPH",  0, 0, nPHbi, "# mPH: Masses of phases in equilibrium, kg [nPHb]" },
    { "vPH",  0, 0, nPHbi, "# vPH: Volumes of phases in equilibrium, m3 [nPHb]" },
    { "m_t",  0, 0, nICbi, "# m_t: Total dissolved molality of independent components, m [nICb]" },
    { "con",  0, 0, nDCbi, "# con: DC concentrations in phases (molal, mole fraction) [nDCb]" },
    { "mju",  0, 0, nDCbi, "# mju: DC chemical potentials in equilibrium, J/mol [nDCb]" },
    { "lga",  0, 0, nDCbi, "# lga: DC activities in equilibrium, in log10 scale [nDCb]" }
};

std::vector<io_formats::outField> DataCH_static_fields =  {  // [14]
  { "nIC",   1, 0, 0, "# nIC: Number of Independent Components (usually chemical elements and charge)" },
  { "nDC",   1, 0, 0, "# nDC: Number of Dependent Components (chemical species made of Independent Components)" },
  { "nPH",   1, 0, 0, "# nPH: Number of phases (into which Dependent Components are grouped)" },
  { "nPS",   1, 0, 0, "# nPS: Number of phases-solutions (multicomponent phases) <= nPH" },
  { "nDCs",  1, 0, 0, "# nDCs: Number of Dependent Components in phases-solutions <= nDC" },
  { "nICb",  1, 0, 0, "# nICb: Number of ICs kept in the DBR file and DATABR memory structure (<= nIC)" },
  { "nDCb",  1, 0, 0, "# nDCb: Number of DCs kept in the DBR file and DATABR memory structure (<=nDC)"  },
  { "nPHb",  1, 0, 0, "# nPHb: Number of phases kept in the DBR file and DATABR structure (<=nPH)" },
  { "nPSb",  1, 0, 0, "# nPSb: Number of phases-solutions kept in the DBR file and DATABR structure (<=nPS)" },
  { "nTp",   1, 0, 0, "# nTp: Number of temperature grid points in lookup arrays for data interpolation, >=1" },
  { "nPp",   1, 0, 0, "# nPp: Number of pressure grid points in lookup arrays for data interpolation, >=1" },
  { "iGrd",  1, 0, 0, "# iGrd: Flag for allocation of array of diffusition coefficients in DATACH structure (DCH file)" },
  { "fAalp", 1, 0, 0, "# fAalp: Flag for keeping specific surface areas of phases in DATABR structure (1) or ignoring them (0)" },
  { "mLook", 1, 0, 0, "# mLook: Lookup mode: 0 interpolation over nTp*nPp grid; 1 data for T,P pairs, no interpolation"}
};

std::vector<io_formats::outField> DataCH_dynamic_fields =  { //  [30] +4
   { "xic",   1, 0, 0, "# xIC: DATACH access index list for ICs kept in the DATABR structure and in DBR files [nICb]" },
   { "xdc",   1, 0, 0, "# xDC: DATACH access index list of DCs kept in the DATABR  structure and in DBR files [nDCb]" },
   { "xph",   1, 0, 0, "# xPH: DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]" },
   { "ICNL",  1, 0, 0, "\n# ICNL: Name list of ICs (Independent Components, up to 5 characters per name) [nIC]" },
   { "ccIC",  1, 0, 0, "# ccIC: Class codes of ICs (Independent Components) [nIC]" },
   { "ICmm",  1, 0, 0, "# ICmm: Atomic (molar) masses of ICs,  kg/mol [nIC]" },
   { "DCNL",  1, 0, 0, "\n# DCNL: Name list of DCs (Dependent Components, up to 16 characters per name) [nDC]" },
   { "ccDC",  1, 0, 0, "# ccDC: Class codes of DCs (Dependent Components) [nDC]" },
   { "DCmm",  0, 0, 0, "\n# DCmm: Molar masses of DCs, kg/mol [nDC]" },
   { "PHNL",  1, 0, 0, "# PHNL: List of Phase names (up to 16 characters per name) [nPH]" },
   { "ccPH",  1, 0, 0, "# ccPH: Codes of phase aggregate state [nPH]" },
   { "nDCinPH",  1, 0, 0, "# nDCinPH: Number of DCs included in each phase [nPH]" },
   { "A",     1, 0, 0, "# A: Stoichiometry matrix A (expanded formulae) for DCs [nDC*nIC]"},
   { "Ttol",  0, 0, 0, "# Ttol: Tolerance for the temperature interpolation, K" },
   { "TKval", 1, 0, 0, "# TKval: Temperature values, K for lookup arrays of thermodynamic data [nTp]" },
   { "Ptol",  0, 0, 0, "# Ptol: Tolerance for the pressure interpolation, Pa" },
   { "Pval",  1, 0, 0, "# Pval: Pressure values, Pa for lookup arrays of thermodynamic data [nPp]" },
   { "denW",  1, 0, 0, "\n# denW: Look-up array for the density of water-solvent, kg/m3, and its derivatives [5*nPp*nTp]" },
   { "denWg", 1, 0, 0, "\n# denWg: Look-up array for the density of water vapour, kg/m3, and its derivatives [5*nPp*nTp]" },
   { "epsW",  1, 0, 0, "\n# epsW: Look-up array for the dielectric constant of water-solvent and its derivatives [5*nPp*nTp]" },
   { "epsWg", 1, 0, 0, "\n# epsWg: Look-up array for the dielectric constant of water vapour and its derivatives [5*nPp*nTp]" },
//   { "visW",  1, 0, 0 },
   { "V0",    1, 0, 0, "\n# V0: Look-up array for DC (standard) molar volumes, J/Pa [nDC*nPp*nTp]" },
   { "G0",    1, 0, 0, "\n# G0: Look-up array for DC molar Gibbs energy function g(T,P), J/mol [nDC*nPp*nTp]" },
   { "H0",    0, 0, 0, "\n# H0: Look-up array for DC molar enthalpy h(T,P), J/mol [nDC*nPp*nTp]" },
   { "S0",    0, 0, 0, "\n# S0: Look-up array for DC absolute entropy S(T,P), J/K/mol [nDC*nPp*nTp] " },
   { "Cp0",   0, 0, 0, "\n# Cp0: Look-up array for DC heat capacity Cp(T,P), J/K/mol [nDC*nPp*nTp]" },
   { "A0",    0, 0, 0, "\n# A0: reserved: Look-up array for DC Helmholtz energy function, J/mol [nDC*nPp*nTp]" },
   { "U0",    0, 0, 0, "\n# U0: reserved: Look-up array for DC internal energy function, J/mol [nDC*nPp*nTp]" },
   { "DD",    0, 0, 0, "\n# DD: reserved: Look-up array for DC diffusion coefficients [nDC*nPp*nTp]" },
   { "Psat",  0, 0, 0, "# Psat: Pressure Pa at saturated H2O vapour at given temperature [nTp]" }
};

//===============================================================

template<typename TIO>
void databr_to_text_file(const DATACH* CSD, const DATABR* CNode, TIO& out_format, bool with_comments, bool brief_mode)
{
    bool _comment = with_comments;

    out_format.put_head( GEMS3KGenerator::gen_dbr_name(out_format.set_name(), 0, CNode->NodeHandle), "dbr");
    io_formats::TPrintArrays<TIO>  prar( f_omph+1/*55*/, DataBR_fields, out_format );

    if( _comment )
    {
        prar.writeComment( _comment, std::string("# ") + _GEMIPM_version_stamp );
        //prar.writeComment( _comment, std::string("# File: ") + path );
        prar.writeComment( _comment, "# Comments can be marked with # $ ; as the first character in the line" );
        prar.writeComment( _comment, "# DBR text input file for node system recipe and speciation data" );
        prar.writeComment( _comment, "# (should be read only after the DCH and the IPM files)\n" );
        prar.writeComment( _comment, "# (1): Flags controlling GEM IPM-3 operation and data exchange");
    }

    prar.writeField(f_NodeHandle, CNode->NodeHandle, _comment, brief_mode  );
    prar.writeField(f_NodeTypeHY, CNode->NodeTypeHY, _comment, brief_mode  );
    prar.writeField(f_NodeTypeMT, CNode->NodeTypeMT, _comment, brief_mode  );
    prar.writeField(f_NodeStatusFMT, CNode->NodeStatusFMT, _comment, brief_mode  );
    prar.writeField(f_NodeStatusCH, CNode->NodeStatusCH, _comment, brief_mode  );
    prar.writeField(f_IterDone, CNode->IterDone, _comment, brief_mode  );

    if( _comment )
        prar.writeComment( _comment, "\n## (2) Chemical scalar properies of the node system");

    prar.writeField(f_TK, CNode->TK, _comment, brief_mode  );
    prar.writeField(f_P, CNode->P, _comment, brief_mode  );

    prar.writeField(f_Vs, CNode->Vs, _comment, brief_mode  );
    prar.writeField(f_Vi, CNode->Vi, _comment, brief_mode  );

    prar.writeField(f_Ms, CNode->Ms, _comment, brief_mode  );
    prar.writeField(f_Mi, CNode->Mi, _comment, brief_mode  );


    prar.writeField(f_Hs, CNode->Hs, _comment, brief_mode  );
    prar.writeField(f_Hi, CNode->Hi, _comment, brief_mode  );

    prar.writeField(f_Gs, CNode->Gs, _comment, brief_mode  );

    if( CSD->ccPH[0] == PH_AQUEL )
    {
        prar.writeField(f_IS, CNode->IC, _comment, brief_mode  );
        prar.writeField(f_pH, CNode->pH, _comment, brief_mode  );
        prar.writeField(f_pe, CNode->pe, _comment, brief_mode  );
        prar.writeField(f_Eh, CNode->Eh, _comment, brief_mode  );
    }

    prar.writeField(f_Tm, CNode->Tm, _comment, brief_mode  );
    prar.writeField(f_dt, CNode->dt, _comment, brief_mode  );

#ifdef NODEARRAYLEVEL
    if( CNode->NodeStatusFMT != No_nodearray /*TNodeArray::na->nNodes() > 1*/ )
    {
        if( _comment )
            prar.writeComment( _comment, "\n## (3) Scalar mass-trasport properties (used only at NodeArray level)");
        prar.writeField(f_Dif, CNode->Dif, _comment, brief_mode  );
        prar.writeField(f_Vt, CNode->Vt, _comment, brief_mode  );
        prar.writeField(f_vp, CNode->vp, _comment, brief_mode  );
        prar.writeField(f_eps, CNode->eps, _comment, brief_mode  );
        prar.writeField(f_Km, CNode->Km, _comment, brief_mode  );
        prar.writeField(f_Kf, CNode->Kf, _comment, brief_mode  );
        prar.writeField(f_S, CNode->S, _comment, brief_mode  );
        prar.writeField(f_Tr, CNode->Tr, _comment, brief_mode  );
        prar.writeField(f_h, CNode->h, _comment, brief_mode  );
        prar.writeField(f_rho, CNode->rho, _comment, brief_mode  );
        prar.writeField(f_al, CNode->al, _comment, brief_mode  );
        prar.writeField(f_at, CNode->at, _comment, brief_mode  );
        prar.writeField(f_av, CNode->av, _comment, brief_mode  );
        prar.writeField(f_hDl, CNode->hDl, _comment, brief_mode  );
        prar.writeField(f_hDt, CNode->hDt, _comment, brief_mode  );
        prar.writeField(f_hDv, CNode->hDv, _comment, brief_mode  );
        prar.writeField(f_nto, CNode->nto, _comment, brief_mode  );
    }
#endif

    if( _comment )
    {
        prar.writeComment( _comment, "\n### Arrays: for dimensions and index lists, see Section (2) of DCH file\n" );
        prar.writeComment( _comment, "## (4) Data for Independent Components");
        prar.writeArray( "", CSD->ICNL[0], CSD->nIC, MaxICN );
    }

    prar.writeArray(  f_bIC,  CNode->bIC, CSD->nICb, -1L,_comment, brief_mode );
    prar.writeArray(  f_rMB,  CNode->rMB, CSD->nICb, -1L,_comment, brief_mode );
    prar.writeArray(  f_uIC,  CNode->uIC, CSD->nICb, -1L,_comment, brief_mode );
    prar.writeArray(  f_bSP,  CNode->bSP, CSD->nICb, -1L,_comment, brief_mode );

    if( _comment )
    {
        prar.writeComment( _comment, "\n## (5) Data for Dependent Components");
        prar.writeArray(  "", CSD->DCNL[0], CSD->nDC, MaxDCN );
    }

    prar.writeArray(  f_xDC,  CNode->xDC, CSD->nDCb, -1L,_comment, brief_mode  );
    prar.writeArray(  f_gam,  CNode->gam, CSD->nDCb, -1L,_comment, brief_mode  );
    prar.writeArray(  f_dll,  CNode->dll, CSD->nDCb, -1L,_comment, brief_mode  );
    prar.writeArray(  f_dul,  CNode->dul, CSD->nDCb, -1L,_comment, brief_mode  );

    if( _comment )
    {
        prar.writeComment( _comment, "\n## (6) Data for Phases");
        prar.writeArray(  "", CSD->PHNL[0], CSD->nPH, MaxPHN );
    }

    prar.writeArray(  f_aPH,  CNode->aPH, CSD->nPHb, -1L,_comment, brief_mode );
    prar.writeArray(  f_xPH,  CNode->xPH, CSD->nPHb, -1L,_comment, brief_mode );
    prar.writeArray(  f_vPS,  CNode->vPS, CSD->nPSb, -1L,_comment, brief_mode );
    prar.writeArray(  f_mPS,  CNode->mPS, CSD->nPSb, -1L,_comment, brief_mode );
    prar.writeArray(  f_xPA,  CNode->xPA, CSD->nPSb, -1L,_comment, brief_mode );
    prar.writeArray(  f_amru,  CNode->amru, CSD->nPSb, -1L,_comment, brief_mode );
    prar.writeArray(  f_amrl,  CNode->amrl, CSD->nPSb, -1L,_comment, brief_mode );
    prar.writeArray(  f_omph,  CNode->omPH, CSD->nPHb, -1L,_comment, brief_mode );

    if(!brief_mode || prar.getAlws( f_bPS ))
    {
        if( _comment )
        {
            prar.writeComment( _comment,  DataBR_fields[f_bPS].comment.c_str() );
            prar.writeArray(  "", CSD->ICNL[0], CSD->nIC, MaxICN );
        }
        prar.writeArray(  f_bPS,  CNode->bPS, CSD->nPSb*CSD->nICb, CSD->nICb,false, brief_mode );
    }

    out_format.dump(  _comment );
}

// Reading work dataBR structure from text file
template<typename TIO>
void databr_from_text_file(const DATACH* CSD, DATABR* CNode, TIO& in_format )
{
#ifndef NODEARRAYLEVEL
    double tmpVal;
#endif

    databr_free_internal(CNode);
    databr_realloc(CSD, CNode);
    // mem_set( &CNode->Tm, 0, 19*sizeof(double));
    databr_reset(CNode, 0);

    io_formats::TReadArrays<TIO>  rdar(f_omph+1/*55*/, DataBR_fields, in_format);

    long int nfild = rdar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f_NodeHandle: rdar.readArray( "NodeHandle",  &CNode->NodeHandle, 1);
            break;
        case f_NodeTypeHY: rdar.readArray( "NodeTypeHY",  &CNode->NodeTypeHY, 1);
            break;
        case f_NodeTypeMT: rdar.readArray( "NodeTypeMT",  &CNode->NodeTypeMT, 1);
            break;
        case f_NodeStatusFMT: rdar.readArray( "NodeStatusFMT",  &CNode->NodeStatusFMT, 1);
            break;
        case f_NodeStatusCH: rdar.readArray( "NodeStatusCH",  &CNode->NodeStatusCH, 1);
            break;
        case f_IterDone: rdar.readArray( "IterDone",  &CNode->IterDone, 1);
            break;
        case f_TK: rdar.readArray( "TK",  &CNode->TK, 1);
            break;
        case f_P: rdar.readArray( "P",  &CNode->P, 1);
            break;
        case f_Vs: rdar.readArray( "Vs", &CNode->Vs, 1);
            break;
        case f_Vi: rdar.readArray( "Vi",  &CNode->Vi, 1);
            break;
        case f_Ms: rdar.readArray( "Ms",  &CNode->Ms, 1);
            break;
        case f_Mi: rdar.readArray( "Mi",  &CNode->Mi, 1);
            break;
        case f_Hs: rdar.readArray( "Hs",  &CNode->Hs, 1);
            break;
        case f_Hi: rdar.readArray( "Hi",  &CNode->Hi, 1);
            break;
        case f_Gs: rdar.readArray( "Gs",  &CNode->Gs, 1);
            break;
        case f_IS: rdar.readArray( "IS",  &CNode->IC, 1);
            break;
        case f_pH: rdar.readArray( "pH",  &CNode->pH, 1);
            break;
        case f_pe: rdar.readArray( "pe",  &CNode->pe, 1);
            break;
        case f_Eh: rdar.readArray( "Eh",  &CNode->Eh, 1);
            break;
        case f_Tm: rdar.readArray( "Tm",  &CNode->Tm, 1);
            break;
        case f_dt: rdar.readArray( "dt",  &CNode->dt, 1);
            break;
#ifdef NODEARRAYLEVEL
        case f_Dif: rdar.readArray( "Dif",  &CNode->Dif, 1);
            break;
        case f_Vt: rdar.readArray( "Vt",  &CNode->Vt, 1);
            break;
        case f_vp: rdar.readArray( "vp",  &CNode->vp, 1);
            break;
        case f_eps: rdar.readArray( "eps",  &CNode->eps, 1);
            break;
        case f_Km: rdar.readArray( "Km",  &CNode->Km, 1);
            break;
        case f_Kf: rdar.readArray( "Kf",  &CNode->Kf, 1);
            break;
        case f_S: rdar.readArray( "S",  &CNode->S, 1);
            break;
        case f_Tr: rdar.readArray( "Tr",  &CNode->Tr, 1);
            break;
        case f_h: rdar.readArray( "h",  &CNode->h, 1);
            break;
        case f_rho: rdar.readArray( "rho",  &CNode->rho, 1);
            break;
        case f_al: rdar.readArray( "al",  &CNode->al, 1);
            break;
        case f_at: rdar.readArray( "at",  &CNode->at, 1);
            break;
        case f_av: rdar.readArray( "av",  &CNode->av, 1);
            break;
        case f_hDl: rdar.readArray( "hDl",  &CNode->hDl, 1);
            break;
        case f_hDt: rdar.readArray( "hDt",  &CNode->hDt, 1);
            break;
        case f_hDv: rdar.readArray( "hDv",  &CNode->hDv, 1);
            break;
        case f_nto: rdar.readArray( "nto",  &CNode->nto, 1);
            break;
#else
        case f_Dif: rdar.readArray( "Dif",  &tmpVal, 1);
            break;
        case f_Vt: rdar.readArray( "Vt",  &tmpVal, 1);
            break;
        case f_vp: rdar.readArray( "vp",  &tmpVal, 1);
            break;
        case f_eps: rdar.readArray( "eps",  &tmpVal, 1);
            break;
        case f_Km: rdar.readArray( "Km",  &tmpVal, 1);
            break;
        case f_Kf: rdar.readArray( "Kf",  &tmpVal, 1);
            break;
        case f_S: rdar.readArray( "S",  &tmpVal, 1);
            break;
        case f_Tr: rdar.readArray( "Tr",  &tmpVal, 1);
            break;
        case f_h: rdar.readArray( "h",  &tmpVal, 1);
            break;
        case f_rho: rdar.readArray( "rho",  &tmpVal, 1);
            break;
        case f_al: rdar.readArray( "al",  &tmpVal, 1);
            break;
        case f_at: rdar.readArray( "at",  &tmpVal, 1);
            break;
        case f_av: rdar.readArray( "av",  &tmpVal, 1);
            break;
        case f_hDl: rdar.readArray( "hDl",  &tmpVal, 1);
            break;
        case f_hDt: rdar.readArray( "hDt",  &tmpVal, 1);
            break;
        case f_hDv: rdar.readArray( "hDv",  &tmpVal, 1);
            break;
        case f_nto: rdar.readArray( "nto",  &tmpVal, 1);
            break;
#endif
        case f_bIC: rdar.readArray( "bIC",  CNode->bIC, CSD->nICb );
            break;
        case f_rMB: rdar.readArray( "rMB",  CNode->rMB, CSD->nICb );
            break;
        case f_uIC: rdar.readArray( "uIC",  CNode->uIC, CSD->nICb );
            break;
        case f_xDC: rdar.readArray( "xDC",  CNode->xDC, CSD->nDCb );
            break;
        case f_gam: rdar.readArray( "gam",  CNode->gam, CSD->nDCb );
            break;
        case f_dll: rdar.readArray( "dll",  CNode->dll, CSD->nDCb );
            break;
        case f_dul: rdar.readArray( "dul",  CNode->dul, CSD->nDCb );
            break;
        case f_aPH: rdar.readArray( "aPH",  CNode->aPH, CSD->nPHb );
            break;
        case f_xPH: rdar.readArray( "xPH",  CNode->xPH, CSD->nPHb );
            break;
        case f_vPS: rdar.readArray( "vPS",  CNode->vPS, CSD->nPSb );
            break;
        case f_mPS: rdar.readArray( "mPS",  CNode->mPS, CSD->nPSb );
            break;
        case f_bPS: rdar.readArray( "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
            break;
        case f_xPA: rdar.readArray( "xPA",  CNode->xPA, CSD->nPSb );
            break;
        case f_bSP: rdar.readArray( "bSP",  CNode->bSP, CSD->nICb );
            break;
        case f_amru: rdar.readArray( "amru",  CNode->amru, CSD->nPSb );
            break;
        case f_amrl: rdar.readArray( "amrl",  CNode->amrl, CSD->nPSb );
            break;
        case f_omph: rdar.readArray( "omPH",  CNode->omPH, CSD->nPHb );
            break;
        }
        nfild = rdar.findNext();
    }

    // testing read
    std::string ret = rdar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from DataBR structure";
        Error( "Error", ret);
    }
}

//===============================================================

template<typename TIO>
void datach_to_text_file(const DATACH* CSD, TIO& out_format, bool use_thermofun, bool with_comments, bool brief_mode )
{
    bool _comment = with_comments;

    out_format.put_head( GEMS3KGenerator::gen_dch_name( out_format.set_name() ), "dch");
    io_formats::TPrintArrays<TIO>  prar1(14, DataCH_static_fields, out_format );
    io_formats::TPrintArrays<TIO>  prar( 30, DataCH_dynamic_fields, out_format );

    if( CSD->nIC == CSD->nICb )
        prar.setNoAlws( f_xic);
    if(CSD->nDC == CSD->nDCb )
        prar.setNoAlws( f_xdc);
    if(CSD->nPH == CSD->nPHb )
        prar.setNoAlws( f_xph );

    if( _comment )
    {
        prar.writeComment( _comment, std::string( "# ") + _GEMIPM_version_stamp );
        //prar.writeComment( _comment, std::string("# File: ")+ path  );
        prar.writeComment( _comment, "# Comments can be marked with # $ ; as the first character in the line");
        prar.writeComment( _comment, "# DCH text input file (should be read before IPM and DBR files)");
        prar.writeComment( _comment, "\n## (1) Dimensions for memory allocation");
    }
    prar1.writeField(f_nIC, CSD->nIC, _comment, brief_mode  );
    prar1.writeField(f_nDC, CSD->nDC, _comment, brief_mode  );
    prar1.writeField(f_nPH, CSD->nPH, _comment, brief_mode  );
    prar1.writeField(f_nPS, CSD->nPS, _comment, brief_mode  );
    prar1.writeField(f_nDCs, CSD->nDCs, _comment, brief_mode  );

    if( _comment )
        prar.writeComment( _comment, "\n## (2) Dimensions for DBR node recipe (memory allocation)");
    prar1.writeField(f_nICb, CSD->nICb, _comment, brief_mode  );
    prar1.writeField(f_nDCb, CSD->nDCb, _comment, brief_mode  );
    prar1.writeField(f_nPHb, CSD->nPHb, _comment, brief_mode  );
    prar1.writeField(f_nPSb, CSD->nPSb, _comment, brief_mode  );

    if( _comment )
        prar.writeComment( _comment, "\n## (3) Dimensions for thermodynamic data arrays");
    prar1.writeField(f_nTp, CSD->nTp, _comment, brief_mode  );
    prar1.writeField(f_nPp, CSD->nPp, _comment, brief_mode  );
    prar1.writeField(f_iGrd, CSD->iGrd, _comment, brief_mode  );
    prar1.writeField(f_fAalp, CSD->nAalp, _comment, brief_mode  );
    prar1.writeField(f_mLook, CSD->mLook, _comment, brief_mode  );

    prar.writeComment( true,  "\n<END_DIM>\n" );

    // dynamic arrays - must follow static data
    if( _comment )
        prar.writeComment( _comment, "## (4) DBR node recipe connection index lists");
    prar.writeArray(  f_xic, CSD->xic, CSD->nICb, -1L,_comment, brief_mode);
    prar.writeArray(  f_xdc, CSD->xdc, CSD->nDCb, -1L,_comment, brief_mode);
    prar.writeArray(  f_xph, CSD->xph, CSD->nPHb, -1L,_comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n## (5) Independent Components and their properties");
    if(!brief_mode || prar.getAlws( f_ICNL ))
    {
        if( _comment )
            prar.writeComment( _comment, "# ICNL: List of Independent Component names (<=4 characters per name) [nIC]");
        prar.writeArray(  "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
    }
    prar.writeArrayF(  f_ccIC, CSD->ccIC, CSD->nIC, 1L,_comment, brief_mode );
    prar.writeArray(  f_ICmm, CSD->ICmm, CSD->nIC, -1L,_comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n## (6) Dependent Components and their codes");
    if(!brief_mode || prar.getAlws( f_DCNL ))
    {
        if( _comment )
            prar.writeComment( _comment, "# DCNL: Name list of Dependent Components (<=16 characters per name) [nDC]");
        prar.writeArray(  "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
    }
    prar.writeArrayF(  f_ccDC, CSD->ccDC, CSD->nDC, 1L,_comment, brief_mode );
    prar.writeArray(  f_DCmm, CSD->DCmm, CSD->nDC, -1L,_comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n## (7) Phases and their codes");
    if(!brief_mode || prar.getAlws( f_PHNL ))
    {
        if( _comment )
            prar.writeComment( _comment, "# PHNL: List of Phase names (<=16 characters per name) [nPH]");
        prar.writeArray(  "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
    }
    prar.writeArrayF(  f_ccPH, CSD->ccPH, CSD->nPH, 1L,_comment, brief_mode );
    prar.writeArray(  f_nDCinPH, CSD->nDCinPH, CSD->nPH, -1L,_comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n# (8) Data for Dependent Components");
    prar.writeArray(  f_A, CSD->A, CSD->nDC*CSD->nIC, CSD->nIC, _comment, brief_mode );
    prar.writeComment( true,  "");

    if( _comment )
        prar.writeComment( _comment, "## (9) Thermodynamic data for Dependent Components");
    prar.writeField(  f_Ttol, CSD->Ttol, _comment, brief_mode  );
    prar.writeArray(  f_TKval, CSD->TKval, CSD->nTp, -1L,_comment, brief_mode );
    prar.writeArray(  f_Psat, CSD->Psat, CSD->nTp, -1L,_comment, brief_mode );
    prar.writeComment( true,  "");

    prar.writeField(  f_Ptol, CSD->Ptol, _comment, brief_mode  );
    prar.writeArray(  f_Pval, CSD->Pval, CSD->nPp,  -1L,_comment, brief_mode );

    if( CSD->ccPH[0] == PH_AQUEL )
    {
        prar.writeArray(  f_denW, CSD->denW, 5*(gridTP(CSD)), gridTP(CSD), _comment, brief_mode );
        prar.writeArray(  f_denWg, CSD->denWg, 5*(gridTP(CSD)), gridTP(CSD), _comment, brief_mode );
        prar.writeArray(  f_epsW, CSD->epsW,  5*(gridTP(CSD)), gridTP(CSD), _comment, brief_mode );
        prar.writeArray(  f_epsWg, CSD->epsWg, 5*(gridTP(CSD)),  gridTP(CSD), _comment, brief_mode );
    }
    //if( !isAllZero( CSD->V0,  CSD->nDC*gridTP(CSD) ))
    prar.writeArray(  f_V0, CSD->V0,  CSD->nDC*gridTP(CSD), gridTP(CSD), _comment, brief_mode );
    prar.writeArray(  f_G0, CSD->G0, CSD->nDC*gridTP(CSD), gridTP(CSD), _comment, brief_mode );
    if( !isAllZero( CSD->H0,  CSD->nDC*gridTP(CSD) ))
        prar.writeArray(  f_H0, CSD->H0,  CSD->nDC*gridTP(CSD),gridTP(CSD), _comment, brief_mode );
    if( !isAllZero( CSD->S0,  CSD->nDC*gridTP(CSD) ))
        prar.writeArray(  f_S0, CSD->S0,CSD->nDC*gridTP(CSD),  gridTP(CSD), _comment, brief_mode  );
    if( !isAllZero( CSD->Cp0,  CSD->nDC*gridTP(CSD) ))
        prar.writeArray(  f_Cp0, CSD->Cp0,CSD->nDC*gridTP(CSD), gridTP(CSD), _comment, brief_mode  );
    if( !isAllZero( CSD->A0,  CSD->nDC*gridTP(CSD) ))
        prar.writeArray(  f_A0, CSD->A0, CSD->nDC*gridTP(CSD), gridTP(CSD), _comment, brief_mode  );
    if( !isAllZero( CSD->U0,  CSD->nDC*gridTP(CSD) ))
        prar.writeArray(  f_U0, CSD->U0, CSD->nDC*gridTP(CSD), gridTP(CSD), _comment, brief_mode  );

    if( CSD->iGrd  )
        prar.writeArray(  f_DD, CSD->DD, CSD->nDCs*gridTP(CSD),  gridTP(CSD),  _comment, brief_mode);

    out_format.dump( _comment );
}

// Reading dataCH structure from text file
template<typename TIO>
void datach_from_text_file(DATACH* CSD, TIO& in_format, bool use_thermofun)
{
    long int ii;

    // static arrays
    io_formats::TReadArrays<TIO> rdar( 14, DataCH_static_fields, in_format);

    long int nfild = rdar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f_nIC: rdar.readArray( "nIC", &CSD->nIC, 1);
            break;
        case f_nDC: rdar.readArray( "nDC", &CSD->nDC, 1);
            break;
        case f_nPH: rdar.readArray( "nPH", &CSD->nPH, 1);
            break;
        case f_nPS: rdar.readArray( "nPS", &CSD->nPS, 1);
            break;
        case f_nDCs: rdar.readArray( "nDCs", &CSD->nDCs, 1);
            break;
        case f_nICb: rdar.readArray( "nICb", &CSD->nICb, 1);
            break;
        case f_nDCb: rdar.readArray( "nDCb", &CSD->nDCb, 1);
            break;
        case f_nPHb: rdar.readArray( "nPHb", &CSD->nPHb, 1);
            break;
        case f_nPSb: rdar.readArray( "nPSb", &CSD->nPSb, 1);
            break;
        case f_nTp: rdar.readArray( "nTp", &CSD->nTp, 1);
            break;
        case f_nPp: rdar.readArray( "nPp", &CSD->nPp, 1);
            break;
        case f_iGrd: rdar.readArray( "iGrd", &CSD->iGrd, 1);
            break;
        case f_fAalp: rdar.readArray( "fAalp", &CSD->nAalp, 1);
            break;
        case f_mLook: rdar.readArray( "mLook", &CSD->mLook, 1);
            break;
        }
        nfild = rdar.findNext();
    }

    // testing read
    std::string ret = rdar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from DataCH structure";
        Error( "Error", ret);
    }

    datach_realloc(CSD);
    //databr_free_internal(CNode);
    //databr_realloc();

    //dynamic data
    io_formats::TReadArrays<TIO>   rddar( 30, DataCH_dynamic_fields, in_format);

    if( CSD->iGrd  )
        rddar.setAlws( f_DD /*28 "DD"*/);

    // default set up
    for( ii=0; ii< CSD->nDCs*gridTP(CSD); ii++ )
    {
        if( CSD->DD ) CSD->DD[ii] = 0.;
    }
    for( ii=0; ii< CSD->nDC*gridTP(CSD); ii++ )
    {
        if( CSD->Cp0) CSD->Cp0[ii] = 0.;
        if( CSD->H0 ) CSD->H0[ii] = 0.;
        if( CSD->S0 ) CSD->S0[ii] = 0.;
        if( CSD->A0 ) CSD->A0[ii] = 0.;
        if( CSD->U0 ) CSD->U0[ii] = 0.;
        if( CSD->V0 ) CSD->V0[ii] = 0.;
    }
    CSD->Ttol = 0.1;
    CSD->Ptol = 10000;
    if( CSD->nIC == CSD->nICb )
    {
        rddar.setNoAlws( f_xic /*"xic"*/ );
        for( ii=0; ii< CSD->nICb; ii++ )
            CSD->xic[ii] = ii;
    }
    if(CSD->nDC == CSD->nDCb )
    {
        rddar.setNoAlws( f_xdc /*"xdc"*/);
        for( ii=0; ii< CSD->nDCb; ii++ )
            CSD->xdc[ii] = ii;
    }
    if(CSD->nPH == CSD->nPHb )
    {
        rddar.setNoAlws( f_xph /*"xph"*/);
        for( ii=0; ii< CSD->nPHb; ii++ )
            CSD->xph[ii] = ii;
    }

    nfild = rddar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f_xic: rddar.readArray( "xic", CSD->xic, CSD->nICb);
            break;
        case f_xdc: rddar.readArray( "xdc", CSD->xdc, CSD->nDCb);
            break;
        case f_xph: rddar.readArray( "xph", CSD->xph, CSD->nPHb);
            break;
        case f_ICNL: rddar.readArray( "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
            break;
        case f_ccIC: rddar.readArray( "ccIC", CSD->ccIC, CSD->nIC, 1 );
            break;
        case f_ICmm: rddar.readArray( "ICmm", CSD->ICmm, CSD->nIC);
            break;
        case f_DCNL: rddar.readArray( "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
            break;
        case f_ccDC: rddar.readArray( "ccDC", CSD->ccDC, CSD->nDC, 1 );
            break;
        case f_DCmm: rddar.readArray( "DCmm", CSD->DCmm, CSD->nDC);
            break;
        case f_PHNL: rddar.readArray( "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
            break;
        case f_ccPH: rddar.readArray( "ccPH", CSD->ccPH, CSD->nPH, 1 );
            break;
        case f_nDCinPH: rddar.readArray( "nDCinPH", CSD->nDCinPH, CSD->nPH);
            break;
        case f_A: rddar.readArray( "A", CSD->A, CSD->nDC*CSD->nIC );
            break;
        case f_Ttol: rddar.readArray( "Ttol", &CSD->Ttol, 1);
            break;
        case f_TKval: rddar.readArray( "TKval", CSD->TKval, CSD->nTp );
            break;
        case f_Ptol: rddar.readArray( "Ptol", &CSD->Ptol, 1);
            break;
        case f_Pval: rddar.readArray( "Pval", CSD->Pval, CSD->nPp );
            break;
        case f_denW: if( !CSD->denW )
                Error( "Error", "Array denW is not allocated in DCH!");
            rddar.readArray( "denW", CSD->denW, 5*gridTP(CSD) );
            break;
        case f_denWg: if( !CSD->denWg )
                Error( "Error", "Array denWg is not allocated in DCH!");
            rddar.readArray( "denWg", CSD->denWg, 5*gridTP(CSD) );
            break;
        case f_epsW: if( !CSD->epsW )
                Error( "Error", "Array epsW is not allocated in DCH!");
            rddar.readArray( "epsW", CSD->epsW,  5*gridTP(CSD) );
            break;
        case f_epsWg: if( !CSD->epsWg )
                Error( "Error", "Array epsWg is not allocated in DCH!");
            rddar.readArray( "epsWg", CSD->epsWg,  5*gridTP(CSD) );
            break;
        case f_V0: rddar.readArray( "V0", CSD->V0,  CSD->nDC*gridTP(CSD) );
            break;
        case f_G0: rddar.readArray( "G0", CSD->G0, CSD->nDC*gridTP(CSD) );
            break;
        case f_H0: if( !CSD->H0 )
                Error( "Error", "Array HO is not allocated in DCH!");
            rddar.readArray( "H0", CSD->H0,  CSD->nDC*gridTP(CSD));
            break;
        case f_S0: if( !CSD->S0 )
                Error( "Error", "Array S0 is not allocated in DCH!");
            rddar.readArray( "S0", CSD->S0,CSD->nDC*gridTP(CSD));
            break;
        case f_Cp0: if( !CSD->Cp0 )
                Error( "Error", "Array CpO is not allocated in DCH!");
            rddar.readArray( "Cp0", CSD->Cp0,CSD->nDC*gridTP(CSD) );
            break;
        case f_A0: if( !CSD->A0 )
                Error( "Error", "Array AO is not allocated in DCH!");
            rddar.readArray( "A0", CSD->A0, CSD->nDC*gridTP(CSD) );
            break;
        case f_U0: if( !CSD->U0 )
                Error( "Error", "Array UO is not allocated in DCH!");
            rddar.readArray( "U0", CSD->U0, CSD->nDC*gridTP(CSD) );
            break;
        case f_DD: if( !CSD->DD )
                Error( "Error", "Array DD is not allocated in DCH!");
            rddar.readArray( "DD", CSD->DD, CSD->nDCs*gridTP(CSD));
            break;
        case f_Psat: rddar.readArray( "Psat", CSD->Psat, CSD->nTp );
            break;
        }
        nfild = rddar.findNext();
    }

    // Set up flags
    if( CSD->ccPH[0] != PH_AQUEL )
    {
        rddar.setNoAlws( f_denW );
        rddar.setNoAlws( f_denWg );
        rddar.setNoAlws( f_epsW );
        rddar.setNoAlws( f_epsWg );
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
void datach_to_file(const DATACH* CSD, GemDataStream& ff)
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
    ff.writeArray( CSD->Psat,  CSD->nTp );
    ff.writeArray( CSD->Pval,  CSD->nPp );

    ff.writeArray( CSD->ccIC, CSD->nIC );
    ff.writeArray( CSD->ccDC, CSD->nDC );
    ff.writeArray( CSD->ccPH, CSD->nPH );

    if( CSD->ccPH[0] == PH_AQUEL )
    { ff.writeArray( CSD->denW,  5*gridTP(CSD) );
        ff.writeArray( CSD->denWg,  5*gridTP(CSD) );
        ff.writeArray( CSD->epsW, 5*gridTP(CSD) );
        ff.writeArray( CSD->epsWg, 5*gridTP(CSD) );
    }
    ff.writeArray( CSD->G0,  CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->V0,  CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->H0,  CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->S0, CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->Cp0, CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->A0, CSD->nDC*gridTP(CSD) );
    ff.writeArray( CSD->U0, CSD->nDC*gridTP(CSD) );
    if(  CSD->iGrd  )
        ff.writeArray( CSD->DD, CSD->nDCs*gridTP(CSD) );

    ff.writeArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
    ff.writeArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
    ff.writeArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );
}

// Reading DataCH structure from binary file
void datach_from_file(DATACH* CSD, GemDataStream& ff)
{
    // const data
    ff.readArray( &CSD->nIC, 14 );
    ff.readArray( &CSD->Ttol, 4 );

    datach_realloc(CSD);

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
    ff.readArray( CSD->Psat,  CSD->nTp );
    ff.readArray( CSD->Pval,  CSD->nPp );

    ff.readArray( CSD->ccIC, CSD->nIC );
    ff.readArray( CSD->ccDC, CSD->nDC );
    ff.readArray( CSD->ccPH, CSD->nPH );

    if( CSD->ccPH[0] == PH_AQUEL )
    {
        ff.readArray( CSD->denW,  5*gridTP(CSD) );
        ff.readArray( CSD->denWg,  5*gridTP(CSD) );
        ff.readArray( CSD->epsW, 5*gridTP(CSD) );
        ff.readArray( CSD->epsWg, 5*gridTP(CSD) );
    }
    ff.readArray( CSD->G0,  CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->V0,  CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->H0,  CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->S0, CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->Cp0, CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->A0, CSD->nDC*gridTP(CSD) );
    ff.readArray( CSD->U0, CSD->nDC*gridTP(CSD) );
    if(  CSD->iGrd  )
        ff.readArray( CSD->DD, CSD->nDCs*gridTP(CSD) );

    ff.readArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
    ff.readArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
    ff.readArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );
}

// writing DataBR to binary file
void databr_to_file(const DATACH* CSD, const DATABR* CNode, GemDataStream& ff)
{
    // const data
    ff.writeArray( &CNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
    if( CNode->NodeStatusFMT != No_nodearray )
        ff.writeArray( &CNode->TK, 32 );
    else
        ff.writeArray( &CNode->TK, 15 );
#else
    ff.writeArray( &CNode->TK, 15 );
#endif

    //dynamic data
    ff.writeArray( CNode->bIC, CSD->nICb );
    ff.writeArray( CNode->rMB, CSD->nICb );
    ff.writeArray( CNode->uIC, CSD->nICb );
    ff.writeArray( CNode->bSP, CSD->nICb );

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
    ff.writeArray( CNode->amru, CSD->nPSb );
    ff.writeArray( CNode->amrl, CSD->nPSb );
    ff.writeArray(  CNode->omPH, CSD->nPHb );
    //   datach_to_text_file();
    //   databr_to_text_file();
}

// Reading work dataBR structure from binary file
void databr_from_file(const DATACH* CSD, DATABR* CNode, GemDataStream& ff)
{
    databr_free_internal(CNode);
    databr_realloc(CSD, CNode);

    // const data
    ff.readArray( &CNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
    if( CNode->NodeStatusFMT != No_nodearray )
        ff.readArray( &CNode->TK, 32 );
    else
        ff.readArray( &CNode->TK, 15 );
#else
    ErrorIf(CNode->NodeStatusFMT != No_nodearray, ff.GetPath(),
            "Error reading work dataBR structure from binary file (No_nodearray)");
    ff.readArray( &CNode->TK, 15 );
#endif
    //dynamic data
    ff.readArray( CNode->bIC, CSD->nICb );
    ff.readArray( CNode->rMB, CSD->nICb );
    ff.readArray( CNode->uIC, CSD->nICb );
    ff.readArray( CNode->bSP, CSD->nICb );

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
    ff.readArray( CNode->amru, CSD->nPSb );
    ff.readArray( CNode->amrl, CSD->nPSb );
    ff.readArray(  CNode->omPH, CSD->nPHb );
}


#ifdef USE_NLOHMANNJSON
template void  databr_to_text_file<io_formats::NlohmannJsonWrite>(const DATACH* pCSD, const DATABR* pCNode, io_formats::NlohmannJsonWrite& out_format, bool, bool );
template void  databr_from_text_file<io_formats::NlohmannJsonRead>(const DATACH* pCSD, DATABR* pCNode, io_formats::NlohmannJsonRead& out_format );
template void  datach_to_text_file<io_formats::NlohmannJsonWrite>(const DATACH* pCSD, io_formats::NlohmannJsonWrite& out_format, bool, bool, bool);
template void  datach_from_text_file<io_formats::NlohmannJsonRead>(DATACH* pCSD, io_formats::NlohmannJsonRead& out_format, bool );
#else
template void  databr_to_text_file<io_formats::SimdJsonWrite>(const DATACH* pCSD, const DATABR* pCNode, io_formats::SimdJsonWrite& out_format, bool, bool );
template void  databr_from_text_file<io_formats::SimdJsonRead>(const DATACH* pCSD, DATABR* pCNode, io_formats::SimdJsonRead& out_format );
template void  datach_to_text_file<io_formats::SimdJsonWrite>(const DATACH* pCSD, io_formats::SimdJsonWrite& out_format, bool, bool, bool);
template void  datach_from_text_file<io_formats::SimdJsonRead>(DATACH* pCSD, io_formats::SimdJsonRead& out_format, bool );
#endif

template void  databr_to_text_file<io_formats::KeyValueWrite>(const DATACH* pCSD, const DATABR* pCNode, io_formats::KeyValueWrite& out_format, bool, bool);
template void  databr_from_text_file<io_formats::KeyValueRead>(const DATACH* pCSD, DATABR* pCNode, io_formats::KeyValueRead& out_format );
template void  datach_to_text_file<io_formats::KeyValueWrite>(const DATACH* pCSD, io_formats::KeyValueWrite& out_format, bool, bool, bool);
template void  datach_from_text_file<io_formats::KeyValueRead>(DATACH* pCSD, io_formats::KeyValueRead& out_format, bool);

}  // dbr_dch_api

