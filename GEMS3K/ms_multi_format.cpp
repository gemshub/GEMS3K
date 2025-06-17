//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi_format.cpp
/// Implementation of writing/reading IPM, DCH and DBR text I/O files
//
// Copyright (c) 2006-2022 S.Dmytriyeva,D.Kulik
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

#include "v_detail.h"
#include "io_template.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "io_keyvalue.h"
#include "ms_multi.h"
extern const std::string _GEMIPM_version_stamp;


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

//===================================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here. The fourth field contains
// the text of the comment for this data object, optionally written into the
// text-format output IPM file.
//
std::vector<io_formats::outField> MULTI_static_fields =  {  //8
    { "pa_PE" , 0 , 0, 0, "# PE: Flag for using electroneutrality condition in GEM IPM calculations (1 or 0)" },
    { "PV" ,    0 , 0, 0, "# PV: Flag for the volume balance constraint (on Vol IC) for indifferent equilibria at P_Sat (0 or 1)" },
    { "PSOL" ,  0 , 0, 0, "# PSOL: Total number of DCs in liquid hydrocarbon phases (0; reserved)" },
    { "PAalp" , 0 , 0, 0, "# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases " },
    { "PSigm" , 0 , 0, 0, "# PSigm: Flag for using (+) or ignoring (-) specific surface free energies  " },
    { "Lads" ,  0 , 0, 0, "# Lads: Total number of Dependent Components in sorption phases included into this system" },
    { "FIa" ,   0 , 0, 0, "# FIa: Number of sorption phases included in this system (0 if no sorption phases )" },
    { "FIat" ,  0 , 0, 0, "# FIat: Maximum number of surface types per adsorption phase (if FIa > 0, set FIat = 6)" }
};

std::vector<io_formats::outField> MULTI_dynamic_fields =  { //80
    // write/read dynamic (array) data to/from the text-format IPM file
    {  "sMod",  1 , 0, 0, "# sMod: Codes for TSolMod built-in  models of mixing in multicomponent phases [nPS*8]" },
    {  "LsMod", 1 , 0, 0, "\n# LsMod: Dimensions of TSolMod <IPxPH> and <PMc> data arrays [nPS*3]. In each row (for phase):"
       "\n# [0] number of interaction parameters (rows in <IPx>); [1] max. parameter order (columns in <IPx>);"
       "\n# [2] number of coefficients per interaction parameter in <PMc> array" },
    {  "LsMdc", 1 , 0, 0, "\n# LsMdc: Dimensions of TSolMod <DMc> and <MoiSN> arrays [nPS*3]: In each row (for phase):"
       "\n# [0] number of parameters per component; [1] 0; [2] 0. For multi-site (sublattice) models: "
       "\n#   [1] number of sublattices nS; [2] total number of moieties nM acting in sublattice sites" },
    {  "B",     1 , 0, 0, "# B: Full total bulk composition (vector b), moles [nIC] (will be partially re-written from DBR files)" },
    {  "DCCW",  0 , 0, 0, "# DCCW: internal DC class codes [nDC], will be reset automatically from DCH file content" },
    // DCCW is placeholder - something else can be used here, if needed
    {  "Pparc", 0 , 0, 0, "# Pparc: Partial pressures or fugacities of pure Dependent Components [nDC] (reserved)" },
    {  "fDQF",  0 , 0, 0, "\n# fDQF: DQF parameters of end members or pure gas fugacities, (J/mol/(RT) [nDC]" },
    {  "lnGmf", 0 , 0, 0, "\n# lnGmf: Natural logarithms of DC activity coefficients used at Simplex LP approximation only [nDC]" },
    {  "RLC",   0 , 0, 0, "# RLC: Code of metastability constraints for DCs {L U B (default)} [nDC]" },
    {  "RSC",   0 , 0, 0, "\n# RSC: Units of metastability/kinetic constraints for DCs {M} moles [nDC]" },
    {  "DLL",   0 , 0, 0, "\n# DLL: Lower metastability constraints on DC amounts <xDC>, moles [nDC] (default: 0)" },
    {  "DUL",   0 , 0, 0, "\n# DUL: Upper metastability constraints on DC amounts <xDC>, moles [nDC] (default: 1e6)" },
    {  "Aalp",  0 , 0, 0, "# Aalp: Specific surface areas of phases, m2/g [nPH]" },
    {  "Sigw",  0 , 0, 0, "\n# Sigw: Specific surface free energy for phase-water interface, J/m2 [nPH] (reserved)" },
    {  "Sigg",  0 , 0, 0, "\n# Sigg: Specific surface free energy for phase-gas interface, J/m2 (not yet used) [nPH]" },
    {  "YOF",   0 , 0, 0, "\n# YOF: Surface free energy parameter for phases in J/g (to accomodate for variable phase composition)  [nPH]" },
    {  "Nfsp",  0 , 0, 0, "\n# Nfsp: Fractions of the sorbent specific surface area allocated to surface types [nPS*6]" },
    {  "MASDT", 0 , 0, 0, "\n# MASDT: Total maximum site  density per surface type, mkmol/g [nPS*6]" },
    {  "C1",    0 , 0, 0, "\n# C1: Inner capacitance density parameter C1, F/m2 [nPS*6]" },
    {  "C2",    0 , 0, 0, "\n# C2: Outer capacitance density parameter C2, F/m2 [nPS*6]" },
    {  "C3",    0 , 0, 0, "\n# C3: Third capacitance density parameter C3, F/m2 [nPS*6]" },
    {  "pCh",   0 , 0, 0, "\n# pCh: Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [nPS*6]" },
    {  "SATX",  0 , 0, 0, "\n# SATX: Setup of surface sites and species (will be applied separately within each sorption phase) [Lads*4]"
       "\n# [0] surface type; [1] sorbent emd member; [2] surface site in surf. type; [3] surface EDL plane" },
    {  "MASDJ", 0, 0, 0, "\n# MASDJ: Parameters of surface species in surface complexation models [Lads*6]"
       "\n# [0] max site density mmol/g; [1] species charge allocated to 0 plane;"
       "\n# [2] species charge allocated to beta -or third plane; [3] Frumkin isotherm interaction parameter;"
       "\n# [4] denticity or coordination number CN; [5] reserved isoterm parameter" },
    {  "SCM",   0, 0, 0, "\n# SCM: Classifier of built-in electrostatic models applied to surface types in sorption phases [nPS*6]" },
    {  "SACT",  0, 0, 0, "\n# SACT: Classifier of applied SACT equations (isotherm corrections) [Lads]" },
    {  "DCads", 0, 0, 0, "\n# DCads: Classifier of DCs involved in sorption phases [Lads]" },
    // static: GEM IPM v3 numerical tolerances
    { "pa_DB" , 0 , 0, 0,  "# DB: Minimum amount of IC in the bulk composition, moles (except charge Zz) { 1e-17 }"},
    { "pa_DHB", 0 , 0, 0,  "\n# DHB: Maximum allowed relative mass balance residual for ICs { 1e-13 } " },
    { "pa_EPS", 0 , 0, 0,  "\n# EPS: Tolerance of the SolveSimplex() balance residual for ICs { 1e-10 } " },
    { "pa_DK",  0 , 0, 0,  "\n# DK: Tolerance for the Dikin's criterion of IPM convergence { 1e-6 } " },
    { "pa_DF" , 0 , 0, 0,  "\n# DF: Tolerance DF of the stability criterion for a lost phase to be inserted to mass balance { 0.01 } " },
    { "pa_DP",  0 , 0, 0,  "\n# DP: Maximal number of iterations in MassBalanceRefinement MBR() procedure { 130 }"  },
    { "pa_IIM", 0 , 0, 0,  "\n# IIM: Maximum allowed number of iterations in one main GEM IPM descent run { 7000 }" },
    { "pa_PD" , 0 , 0, 0,  "\n# PD: Mode of calculation of DC activity coefficients { 2 } " },
    { "pa_PRD" , 0 , 0, 0, "\n# PRD: Disable (0) or activate (-4 or less- max.dec.exp.for DC amount correction) SpeciationCleanup() { -5 }" },
    { "pa_AG" ,  0 , 0, 0, "\n# AG: Smoothing parameter 1 for non-ideal primal chemical potential increments (-1 to +1) { 1.0 }" },
    { "pa_DGC" , 0 , 0, 0, "\n# DGC: Smoothing parameter 2- exponent in smoothing function (-1 to +1) { 0 or 0.001 for adsorption }" },
    { "pa_PSM" , 0 , 0, 0, "\n# PSM: Level of diagnostic messages { 0- disabled (no ipmlog file); 1- default; 2-including warnings }" },
    { "pa_GAR" , 0 , 0, 0, "# GAR: Activity coefficient for major (M) species in solution phases at Simplex LP AIA { 1 }"  },
    { "pa_GAH" , 0 , 0, 0, "# GAH: Activity coefficient for minor (J) species in solution phases at Simplex LP AIA { 1000 }" },
    { "pa_DS",   0 , 0, 0, "\n# DS: Cutoff minimum amount of stable phase in GEM IPM primal solution, moles { 1e-20 }" },
    { "pa_XwMin" , 0 , 0, 0, "# XwMin: Cutoff mole amount of water-solvent for aqueous phase elimination { 1e-13 }" },
    { "pa_ScMin" , 0 , 0, 0, "# ScMin: Cutoff mole amount of solid sorbent for sorption phase elimination { 1e-13 }" },
    { "pa_DcMin" , 0 , 0, 0, "# DcMin: Cutoff mole amount for elimination of DC (species) in multi-component phase { 1e-33 }" },
    { "pa_PhMin" , 0 , 0, 0, "# PhMin: Cutoff mole amount for elimination of solution phases other than aqueous { 1e-20 }" },
    { "pa_ICmin" , 0 , 0, 0, "# ICmin: Cutoff effective molal ionic strength for calculation of aqueous activity coefficients { 1e-5 }" },
    { "pa_PC" ,   0 , 0, 0,  "\n# PC: Mode of Phase Selection: 1 old (Select-2), 2 new (PSSC), default { 2 }" },
    { "pa_DFM" ,  0 , 0, 0,  "# DFM: Tolerance for stability criterion for a phase to be eliminated from mass balance { 0.01 } " },
    { "pa_DFYw" , 0 , 0, 0,  "# DFYw: Insertion mole amount for water-solvent at Simplex()->MBR() bridge { 1e-5 }"},
    { "pa_DFYaq", 0 , 0, 0,  "# DFYaq: Insertion mole amount for aqueous species at Simplex()->MBR() bridge { 1e-5 }"  },
    { "pa_DFYid", 0 , 0, 0,  "\n# DFYid: Insertion mole amount for DCs of ideal solution phases at Simplex()->MBR() bridge { 1e-5 }" },
    { "pa_DFYr" , 0 , 0, 0,  "# DFYr: Insertion mole amount for major DCs in solution phases at Simplex()->MBR()bridge { 1e-5 }" },
    { "pa_DFYh" , 0 , 0, 0,  "# DFYh: Insertion mole amount for junior DCs in solution phases Simplex()->MBR() bridge{ 1e-5 }" },
    { "pa_DFYc" , 0 , 0, 0,  "# DFYc: Insertion mole amount for single-component phase at Simplex()->MBR() bridge { 1e-5 }" },
    { "pa_DFYs",  0 , 0, 0,  "# DFYs: Insertion mole amount for single-component phase in PSSC() algorithm { 1e-6 }" },
    { "pa_DW",    0 , 0, 0,  "# DW: Activate (1) or disable (0) error condition on maximum number of MBR() iterations DP { 1 }" },
    { "pa_DT",    0 , 0, 0,  "# DT: use DHB as relative maximum mass balance cutoff for all ICs (0, default); or for major ICs:"
      "\n# decimal exponent (<-6) applied to DHB cutoff; (1) use DHB also as an absolute cutoff { 1 }" },
    { "pa_GAS",   0 , 0, 0,  "\n# GAS: Threshold for primal-dual chemical potential difference used in SpeciationCleanup() { 0.001 }" },
    { "pa_DG",    0 , 0, 0,  "# Total number of moles used in internal re-scaling of the system (disabled if < 1e-4) { 1000 }" },
    { "pa_DNS" ,  0 , 0, 0,  "# DNS: Standard surface number density, nm-2 for calculating activity of surface species { 12.05 }" },
    { "pa_IEPS" , 0 , 0, 0,  "# IEPS: Tolerance for calculation of surface activity coefficient terms for surface species { 0.001 }" },
    { "pKin" ,    0 , 0, 0,  "\n# pKin: Flag for using metastability constraints on DC amounts in primal GEM solution { 1 } " },
    { "pa_DKIN" , 0 , 0, 0,  "# DKIN: Tolerance for non-trivial metastability constraints on DC amounts, moles { 1e-10 } " },
    { "mui" ,     0 , 0, 0,  "\n# mui: IC indices in parent RMULTS IC list (not used in standalone GEMS3K)" },
    { "muk" ,     0 , 0, 0,  "\n# muk: Phase indices in parent RMULTS Phase list (not used in standalone GEMS3K)" },
    { "muj" ,     0 , 0, 0,  "\n# muj: DC indices in parent RMULTS DC list (not used in standalone GEMS3K)" },
    { "pa_PLLG" , 0 , 0, 0,  "# pa_PLLG: Tolerance for checking divergence in IPM dual solution, 1 to 32001 { 30000 }, 0 disables" },
    { "tMin" ,    0 , 0, 0,  "# tMin: Type of thermodynamic potential to minimize (reserved)" },
    { "dcMod",    0 , 0, 0,  "# dcMod: Codes for PT corrections of DC thermodynamic data [nDC] (reserved)" },
    //TKinMet
    { "kMod",    0 , 0, 0,  "# kMod: Codes for built-in kinetic models [Fi*6]" },
    { "LsKin",    0 , 0, 0,  "# LsKin: number of parallel reactions; of species in activity products; of parameter coeffs in parallel reaction;\n"
      "# of parameters per species; parameter coefficients in As correction; of (separately considered) crystal faces or surface patches ( 1 to 4 ) [Fi][6]" },
    { "LsUpt",    0 , 0, 0,  "# LsUpt: number of uptake kinetics model parameters (coefficients) numpC[k]; (reserved)" },
    { "xICuC",    0 , 0, 0,  "# xICuC: Collected array of IC species indexes used in partition (fractionation) coefficients  ->L1[k]" },
    { "PfFact",    0 , 0, 0,  "# PfFact: form factors for phases (taken from TKinMet or set from TNode) [FI] (reserved)" },
    // TSorpMod stuff
    { "LsESmo",    0 , 0, 0,  "# LsESmo: number of EIL model layers; EIL params per layer; CD coefs per DC; reserved  [Fis][4]" },
    { "LsISmo",    0 , 0, 0,  "# LsISmo: number of surface sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of DC [Fis][4]" },
    // TSorpMod & TKinMet stuff
    { "SorMc",    0 , 0, 0,  "# SorMc: Phase-related kinetics and sorption model parameters: [Fis][16]" },
    // TSolMod stuff
    { "LsMdc2",    0 , 0, 0,  "# LsMdc2: [3*FIs] - number of DQF coeffs; reciprocal coeffs per end member" },
    { "LsPhl",    0 , 0, 0,  "# LsPhl: Number of phase links; number of link parameters; [Fi][2]" }
};


//===================================================================

/// Writing structure MULTI (GEM IPM work structure)
template<typename TIO>
void TMultiBase::to_text_file_gemipm( TIO& out_format, bool addMui,
                                      bool with_comments, bool brief_mode )
{
    const BASE_PARAM *pa_p = base_param();
    bool _comment = with_comments;
    char PAalp;
    char PSigm;
    get_PAalp_PSigm( PAalp, PSigm);

    out_format.put_head( GEMS3KGenerator::gen_ipm_name( out_format.set_name() ), "ipm");
    io_formats::TPrintArrays<TIO>  prar1( 8, MULTI_static_fields, out_format );
    io_formats::TPrintArrays<TIO>  prar( 80, MULTI_dynamic_fields, out_format );

    // set up array flags for permanent fields
    if( !( pm.FIs > 0 && pm.Ls > 0 ) )
    {
        prar.setNoAlws( (long int)(f_sMod ) );
        prar.setNoAlws( f_LsMod );
        prar.setNoAlws( f_LsMdc );
    }
    if( PSigm == S_OFF )
    {
        prar.setNoAlws( f_Sigw);
        prar.setNoAlws( f_Sigg);
    }
    if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
    { /* ADSORPTION AND ION EXCHANGE */
        prar.setNoAlws( f_Nfsp);
        prar.setNoAlws( f_MASDT);
        prar.setNoAlws( f_C1 );
        prar.setNoAlws( f_C2 );
        prar.setNoAlws( f_C3 );
        prar.setNoAlws( f_pCh );
        prar.setNoAlws( f_SATX );
        prar.setNoAlws( f_MASDJ );
        prar.setNoAlws( f_SCM );
        prar.setNoAlws( f_SACT );
        prar.setNoAlws( f_DCads );
    }

    if( _comment )
    {
        prar.writeComment( _comment, std::string("# ") + _GEMIPM_version_stamp );
        // << "# File: " << path << endl;
        prar.writeComment( _comment, "# Comments can be marked with # $ ; as the first character in the line");
        prar.writeComment( _comment, "# IPM text input file for the internal GEM IPM-3 kernel data");
        prar.writeComment( _comment, "# (should be read after the DCH file and before DBR files)\n");
        prar.writeComment( _comment, "# ID key of the initial chemical system definition");
    }

    prar.addField( "ID_key", char_array_to_string(pm.stkey, EQ_RKLEN));

    if( _comment )
        prar.writeComment( _comment, "\n## (1) Flags that affect memory allocation");

    if(!brief_mode || pa_p->PE != pa_p_.PE )
        prar1.writeField(f_pa_PE, pa_p->PE, _comment, false  );

    //   ff << "# Do not know if this stuff is really necessary" << endl;
    //   ff << "# 'GWAT'         55.50837344" << endl;
    //   ff << left << setw(12) << "<GWAT> " <<  right << setw(8) << pm.GWAT << endl;

    prar1.writeField(f_PV, pm.PV, _comment, brief_mode  );
    prar1.writeField(f_PSOL, pm.PSOL, _comment, brief_mode  );

    prar1.writeField(f_PAalp, PAalp, _comment, brief_mode  );
    prar1.writeField(f_PSigm, PSigm, _comment, brief_mode  );

    if( !brief_mode || pm.FIat > 0 || pm.Lads > 0 )
    {
        if( _comment )
            prar.writeComment( _comment, "## (2) Dimensionalities that affect memory allocation");
        prar1.writeField(f_Lads, pm.Lads, _comment, false  );
        prar1.writeField(f_FIa, pm.FIa, _comment, false  );
        prar1.writeField(f_FIat,  pm.FIat, _comment, false  );

        //   ff << left << setw(12) << "<sitNc> " <<  right << setw(8) << pm.sitNcat << endl;
        //   ff << left << setw(12) << "<sitNa> " <<  right << setw(8) << pm.sitNan << endl;
    } // brief_mode

    prar.writeComment( true, "\n<END_DIM>\n");

    // static data not affected by dimensionalities
    if( _comment )
    {
        prar.writeComment( _comment, "## (3) Numerical controls and tolerances of GEM IPM-3 kernel");
        prar.writeComment( _comment,"#      - Need to be changed only in special cases (see gems3k_ipm.html)");
    }

    if( !brief_mode || !essentiallyEqual(pa_p->DB, pa_p_.DB ) )
        prar.writeField(f_pa_DB, pa_p->DB, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->DHB, pa_p_.DHB ) )
        prar.writeField(f_pa_DHB, pa_p->DHB, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->EPS, pa_p_.EPS ) )
        prar.writeField(f_pa_EPS, pa_p->EPS, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->DK, pa_p_.DK ) )
        prar.writeField(f_pa_DK, pa_p->DK, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->DS, pa_p_.DS ) )
        prar.writeField(f_pa_DS,  pa_p->DS, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->DF, pa_p_.DF ) )
        prar.writeField(f_pa_DF, pa_p->DF, _comment, false  );
    if( !brief_mode || !essentiallyEqual(pa_p->DFM, pa_p_.DFM ) )
        prar.writeField(f_pa_DFM,  pa_p->DFM, _comment, false  );
    if(!brief_mode || pa_p->DP != pa_p_.DP )
        prar.writeField(f_pa_DP,  pa_p->DP, _comment, false  );
    if(!brief_mode || pa_p->IIM != pa_p_.IIM )
        prar.writeField(f_pa_IIM,  pa_p->IIM, _comment, false  );
    if(!brief_mode || pa_p->PD != pa_p_.PD )
        prar.writeField(f_pa_PD,  pa_p->PD, _comment, false  );
    if(!brief_mode || pa_p->PRD != pa_p_.PRD )
        prar.writeField(f_pa_PRD,  pa_p->PRD, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->AG, pa_p_.AG ) )
        prar.writeField(f_pa_AG,  pa_p->AG, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DGC, pa_p_.DGC ) )
        prar.writeField(f_pa_DGC,  pa_p->DGC, _comment, false  );
    if(!brief_mode || pa_p->PSM != pa_p_.PSM )
        prar.writeField(f_pa_PSM,  pa_p->PSM, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->GAR, pa_p_.GAR ) )
        prar.writeField(f_pa_GAR,  pa_p->GAR, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->GAH, pa_p_.GAH ) )
        prar.writeField(f_pa_GAH,  pa_p->GAH, _comment, false  );

    if(!brief_mode)
        if( _comment )
        {
            prar.writeComment( _comment, "\n# X*Min: Cutoff amounts for elimination of unstable species ans phases from mass balance");
        }

    if(!brief_mode || !essentiallyEqual(pa_p->XwMin, pa_p_.XwMin ) )
        prar.writeField(f_pa_XwMin,  pa_p->XwMin, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->ScMin, pa_p_.ScMin ) )
        prar.writeField(f_pa_ScMin,  pa_p->ScMin, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DcMin, pa_p_.DcMin ) )
        prar.writeField(f_pa_DcMin,  pa_p->DcMin, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->PhMin, pa_p_.PhMin ) )
        prar.writeField(f_pa_PhMin,  pa_p->PhMin, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->ICmin, pa_p_.ICmin ))
        prar.writeField(f_pa_ICmin,  pa_p->ICmin, _comment, false  );
    if(!brief_mode || pa_p->PC != pa_p_.PC )
        prar.writeField(f_pa_PC,  pa_p->PC, _comment, false  );

    if( _comment )
        prar.writeComment( _comment, "# DFY: Insertion mole amounts used after the LPP AIA and in PhaseSelection() algorithm\n");

    if(!brief_mode || !essentiallyEqual(pa_p->DFYw, pa_p_.DFYw ) )
        prar.writeField(f_pa_DFYw,  pa_p->DFYw, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYaq, pa_p_.DFYaq) )
        prar.writeField(f_pa_DFYaq,  pa_p->DFYaq, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYid, pa_p_.DFYid ) )
        prar.writeField(f_pa_DFYid,  pa_p->DFYid, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYr, pa_p_.DFYr ) )
        prar.writeField(f_pa_DFYr,  pa_p->DFYr, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYh, pa_p_.DFYh ) )
        prar.writeField(f_pa_DFYh,  pa_p->DFYh, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYc, pa_p_.DFYc ) )
        prar.writeField(f_pa_DFYc,  pa_p->DFYc, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DFYs, pa_p_.DFYs ) )
        prar.writeField(f_pa_DFYs,  pa_p->DFYs, _comment, false  );

    if( _comment )
        prar.writeComment( _comment, "# Tolerances and controls of the high-precision IPM-3 algorithm ");

    if(!brief_mode || pa_p->DW != pa_p_.DW )
        prar.writeField(f_pa_DW,  pa_p->DW, _comment, false  );
    if(!brief_mode || pa_p->DT != pa_p_.DT )
        prar.writeField(f_pa_DT,  pa_p->DT, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->GAS, pa_p_.GAS ) )
        prar.writeField(f_pa_GAS,  pa_p->GAS, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DG, pa_p_.DG ) )
        prar.writeField(f_pa_DG,  pa_p->DG, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->DNS, pa_p_.DNS ) )
        prar.writeField(f_pa_DNS, pa_p->DNS, _comment, false  );
    if(!brief_mode || !essentiallyEqual(pa_p->IEPS, pa_p_.IEPS) )
        prar.writeField(f_pa_IEPS, pa_p->IEPS, _comment, false  );
    prar.writeField(f_pKin, pm.PLIM, _comment, brief_mode  );
    if(!brief_mode || !essentiallyEqual(pa_p->DKIN, pa_p_.DKIN ) )
        prar.writeField(f_pa_DKIN, pa_p->DKIN, _comment, false  );
    if(!brief_mode || pa_p->PLLG != pa_p_.PLLG )
        prar.writeField(f_pa_PLLG, pa_p->PLLG, _comment, false  );
    if(!brief_mode || pm.tMin != G_TP_ )
        prar.writeField(f_tMin, pm.tMin, _comment, false  );

    //dynamic arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        if( _comment )
            prar.writeComment( _comment, "\n## (4) Initial data for multicomponent phases (see DCH file for dimension nPHs)");
        prar.writeArrayF(  f_sMod, pm.sMod[0], pm.FIs, 8L, _comment, brief_mode );

        long int LsModSum;
        long int LsIPxSum;
        long int LsMdcSum;
        long int LsMsnSum;
        long int LsSitSum;
        getLsModsum( LsModSum, LsIPxSum );
        getLsMdcsum( LsMdcSum, LsMsnSum, LsSitSum );

        prar.writeArray(  f_LsMod, pm.LsMod, pm.FIs*3, 3L, _comment, brief_mode);

        if(LsIPxSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# IPxPH: Index lists (in TSolMod convention) for interaction parameters of non-ideal solutions");
            prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
        }
        if(LsModSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# PMc: Tables (in TSolMod convention) of interaction parameter coefficients  for non-ideal solutions");
            prar.writeArray(  "PMc", pm.PMc,  LsModSum);
        }
        prar.writeArray(  f_LsMdc, pm.LsMdc, pm.FIs*3, 3L, _comment, brief_mode);
        if(LsMdcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# DMc: Tables (in TSolMod convention) of  parameter coefficients for dependent components");
            prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);
        }
        if(LsMsnSum )
        {
            if( _comment )
                prar.writeComment( _comment,  "\n# MoiSN:  end member moiety / site multiplicity number tables (in TSolMod convention) ");
            prar.writeArray(  "MoiSN", pm.MoiSN,  LsMsnSum);
        }
        long int DQFcSum, rcpcSum;
        getLsMdc2sum( DQFcSum, rcpcSum );
        prar.writeArray(  f_LsMdc2, pm.LsMdc2, pm.FIs*3, 3L, _comment, brief_mode);
        if(DQFcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# DQFc:  Collected array of DQF parameters for DCs in phases ");
            prar.writeArray(  "DQFc", pm.DQFc,  DQFcSum);
        }
        //   if(rcpcSum )
        //   {   if( _comment )
        //          prar.writeComment( _comment, "\n# rcpc:  Collected array of reciprocal parameters for DCs in phases ");
        //     prar.writeArray(  "rcpc", pm.rcpc,  rcpcSum);
        //   }
        long int PhLinSum, lPhcSum;
        getLsPhlsum( PhLinSum,lPhcSum );
        prar.writeArray(  f_LsPhl, pm.LsPhl, pm.FI*2, 2L, _comment, brief_mode);
        if(PhLinSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# PhLin:  indexes of linked phases and link type codes ");
            prar.writeArray(  "PhLin", &pm.PhLin[0][0], PhLinSum*2);
        }
        if(lPhcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# lPhc:  Collected array of phase link parameters ");
            prar.writeArray(  "lPhc", pm.lPhc,  lPhcSum);
        }
        prar.writeArray(  f_SorMc, pm.SorMc, pm.FIs*16, 16L, _comment, brief_mode);

        // TSorpMod stuff
        long int IsoCtSum, IsoScSum;
        long int IsoPcSum, xSMdSum;
        getLsISmosum( IsoCtSum,IsoScSum,IsoPcSum, xSMdSum );
        prar.writeArray(  f_LsISmo, pm.LsISmo, pm.FIs*4, 4L, _comment, brief_mode);
        if(xSMdSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# xSMd:  denticity of surface species per surface site (site allocation) ");
            prar.writeArray(  "xSMd", pm.xSMd, xSMdSum);
        }
        if(IsoPcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# IsoPc:  Collected isotherm coefficients per DC ");
            prar.writeArray(  "IsoPc", pm.IsoPc,  IsoPcSum);
        }
        if(IsoScSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# IsoSc:  Collected isotherm coeffs per site");
            prar.writeArray(  "IsoSc", pm.IsoSc, IsoScSum);
        }
        if(IsoCtSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# IsoCt:  Collected isotherm and SATC codes for surface site types");
            prar.writeArray(  "IsoCt", pm.IsoCt,  IsoCtSum, 1L);
        }
        long int EImcSum, mCDcSum;
        getLsESmosum( EImcSum, mCDcSum );
        prar.writeArray(  f_LsESmo, pm.LsESmo, pm.FIs*4, 4L, _comment, brief_mode);
        if(EImcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# EImc:  Collected EIL model coefficients");
            prar.writeArray(  "EImc", pm.EImc, EImcSum);
        }
        if(mCDcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# mCDc:  Collected CD EIL model coefficients per DC ");
            prar.writeArray(  "mCDc", pm.mCDc,  mCDcSum);
        }
        // TKinMet stuff
        prar.writeArrayF(  f_kMod, pm.kMod[0], pm.FI, 6L, _comment, brief_mode);
        long int xSKrCSum, ocPRkC_feSArC_Sum;
        long int rpConCSum, apConCSum, AscpCSum;
        getLsKinsum( xSKrCSum, ocPRkC_feSArC_Sum, rpConCSum, apConCSum, AscpCSum );
        prar.writeArray(  f_LsKin, pm.LsKin, pm.FI*6, 6L, _comment, brief_mode);
        if(xSKrCSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# xSKrC:  Collected array of aq/gas/sorption species indexes used in activity products");
            prar.writeArray(  "xSKrC", pm.xSKrC, xSKrCSum);
        }
        if(ocPRkC_feSArC_Sum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# ocPRkC:  Collected array of operation codes for kinetic parallel reaction terms");
            prar.writeArray(  "ocPRkC", &pm.ocPRkC[0][0],  ocPRkC_feSArC_Sum*2);
        }
        if(ocPRkC_feSArC_Sum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# feSArC:  Collected array of fractions of surface area related to parallel reactions");
            prar.writeArray(  "feSArC", pm.feSArC, ocPRkC_feSArC_Sum);
        }
        if(rpConCSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# rpConC:  Collected array of kinetic rate constants");
            prar.writeArray(  "rpConC", pm.rpConC,  rpConCSum);
        }
        if(apConCSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# apConC:  Collected array of parameters per species involved in activity product terms");
            prar.writeArray(  "apConC", pm.apConC, apConCSum);
        }
        if(AscpCSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# AscpC:  parameter coefficients of equation for correction of specific surface area");
            prar.writeArray(  "AscpC", pm.AscpC,  AscpCSum);
        }
        prar.writeArray(  f_PfFact, pm.PfFact, pm.FI, 1L, _comment, brief_mode);

        long int UMpcSum, xICuCSum;
        getLsUptsum( UMpcSum, xICuCSum );
        prar.writeArray(  f_LsUpt, pm.LsUpt, pm.FIs*2, 2L, _comment, brief_mode);
        if( UMpcSum )
        {
            if( _comment )
                prar.writeComment( _comment, "\n# UMpcC:  Collected array of uptake model coefficients");
            prar.writeArray(  "UMpcC", pm.UMpcC, UMpcSum);
        }
        if( pm.xICuC )
        {
            prar.writeArray(  f_xICuC, pm.xICuC, xICuCSum, _comment, brief_mode);
        }

    } // sMod

    if( _comment )
        prar.writeComment( _comment, "\n## (5) Data arrays which are provided neither in DCH nor in DBR files");
    prar.writeArray(  f_B, pm.B,  pm.N, -1L, _comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n# Initial data for DCs - see DATACH file for dimensions nDC, nDCs");
    prar.writeArray(  f_Pparc, pm.Pparc,  pm.L, -1L, _comment, brief_mode);
    //  ff << "\n\n# This is not necessary - can be calculated from G0 ???????????";
    // prar.writeArray(  "G0", pm.G0,  pm.L);
    prar.writeArray(  f_fDQF, pm.fDQF,  pm.L, -1L, _comment, brief_mode);
    prar.writeArray(  f_lnGmf, pm.lnGmf,  pm.L, -1L, _comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n# (6) Metastability constraints on DC amounts from above (DUL) and below (DLL)");
    prar.writeArrayF(  f_RLC, pm.RLC, pm.L, 1L, _comment, brief_mode );
    prar.writeArrayF(  f_RSC, pm.RSC, pm.L, 1L, _comment, brief_mode );
    prar.writeArray(  f_DLL, pm.DLL, pm.L, -1L, _comment, brief_mode);
    prar.writeArray(  f_DUL, pm.DUL,  pm.L, -1L, _comment, brief_mode);

    if( _comment )
        prar.writeComment( _comment, "\n# (7) Initial data for Phases\n");
    prar.writeArray(  f_Aalp, pm.Aalp,  pm.FI, -1L, _comment, brief_mode);
    if( PSigm != S_OFF )
    {
        prar.writeArray(  f_Sigw, pm.Sigw,  pm.FI, -1L, _comment, brief_mode);
        prar.writeArray(  f_Sigg, pm.Sigg,  pm.FI, -1L, _comment, brief_mode);
    }
    prar.writeArray(  f_YOF, pm.YOF,  pm.FI, -1L, _comment, brief_mode);

    if( pm.FIat > 0 &&  pm.FIs > 0 )
    { // ADSORPTION AND ION EXCHANGE
        if( _comment )
            prar.writeComment( _comment, "\n# (8) Initial data for sorption phases");

        prar.writeArray(  f_Nfsp, &pm.Nfsp[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_MASDT, &pm.MASDT[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_C1, &pm.XcapA[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_C2, &pm.XcapB[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_C3, &pm.XcapF[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_pCh, &pm.Xetaf[0][0], pm.FIs*pm.FIat, pm.FIat, _comment, brief_mode);
        prar.writeArray(  f_SATX, &pm.SATX[0][0], pm.Lads*4, 4L, _comment, brief_mode);
        prar.writeArray(  f_MASDJ, &pm.MASDJ[0][0], pm.Lads*DFCN, (long int)DFCN, _comment, brief_mode);
        prar.writeArrayF(  f_SCM, pm.SCM[0], pm.FIs, pm.FIat, _comment, brief_mode );
        prar.writeArrayF(  f_SACT, pm.SATT, pm.Lads, 1L, _comment, brief_mode );
        prar.writeArrayF(  f_DCads, pm.DCC3, pm.Lads, 1L, _comment, brief_mode );
    }

    //if(!brief_mode || prar.getAlws("dcMod" ))
    prar.writeArrayF(  f_dcMod, pm.dcMod[0], pm.L, 6L, _comment, brief_mode );

    /*
   outArray( ff, "Vol", pm.Vol,  pm.L);
   outArray( ff, "G0", pm.G0,  pm.L);
   outArray( ff, "PUL", pm.PUL,  pm.L);
   outArray( ff, "PLL", pm.PLL,  pm.L);
   outArray( ff, "lnGam", pm.lnGam,  pm.L);
   outArray( ff, "F0", pm.F0,  pm.L);
*/

    if( addMui && !brief_mode )
    {
        prar.writeArray(  f_mui, pm.mui,  pm.N, -1L, _comment, brief_mode);
        prar.writeArray(  f_muk, pm.muk,  pm.FI, -1L, _comment, brief_mode);
        prar.writeArray(  f_muj, pm.muj,  pm.L, -1L, _comment, brief_mode);
    }

    out_format.dump(  _comment );

}

/// Reading structure MULTI (GEM IPM work structure)
template<typename TIO>
void TMultiBase::from_text_file_gemipm( TIO& in_format,  DATACH  *dCH )
{
    BASE_PARAM *pa_p = base_param();
    long int ii, nfild;
    size_t len;

    //static values
    char PAalp;
    char PSigm;

    set_def();
    //mem_set( &pm.N, 0, 39*sizeof(long int));
    //mem_set( &pm.TC, 0, 55*sizeof(double));
    // get sizes from DATACH
    pm.TC = pm.TCc = 25.;
    pm.T = pm.Tc =298.15;
    pm.P = pm.Pc = 1.;
    pm.N = pm.NR = dCH->nIC;
    pm.L = dCH->nDC;
    pm.FI = dCH->nPH;
    pm.FIs = dCH->nPS;
    //
    pm.Ls = 0; //dCH->nDCs;
    for( ii=0; ii<dCH->nPS; ii++)
    {
        pm.Ls += dCH->nDCinPH[ii];
        if( dCH->ccPH[ii] == 'a' )
            pm.LO = pm.Ls-1;
        if( dCH->ccPH[ii] == 'g' || dCH->ccPH[ii] == 'p' || dCH->ccPH[ii] == 'f')
            pm.PG = dCH->nDCinPH[ii];
    }

    // setup default constants
    base_param()->PE = 1;
    pm.E = 1;
    pm.PV = 0;
    pm.PSOL = 0;
    PAalp = '+';
    PSigm = '+';
    pm.Lads = 0;
    pm.FIa = 0;
    pm.FIat = 0; //6
    pm.PLIM  = 1;

    // static arrays
    io_formats::TReadArrays<TIO> rdar( 8, MULTI_static_fields, in_format);
    if( !in_format.skip_line() ) // Skip line without <ID_key> in old format
    {
      rdar.readNext( "ID_key");
      rdar.readArray( "ID_key", pm.stkey,  1, EQ_RKLEN);
    }

    nfild = rdar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f_pa_PE: rdar.readArray("pa_PE" , &pa_p->PE, 1);
            pm.E = pa_p->PE;
            break;
        case f_PV: rdar.readArray("PV" , &pm.PV, 1);
            break;
        case f_PSOL: rdar.readArray("PSOL" , &pm.PSOL, 1);
            break;
        case f_PAalp: rdar.readArray("PAalp" , &PAalp, 1, 1);
            break;
        case f_PSigm: rdar.readArray("PSigm" , &PSigm, 1, 1);
            break;
        case f_Lads: rdar.readArray("Lads" , &pm.Lads, 1);
            break;
        case f_FIa: rdar.readArray("FIa" , &pm.FIa, 1);
            break;
        case f_FIat: rdar.readArray("FIat" , &pm.FIat, 1);
            break;
        }
        nfild = rdar.findNext();
    }

    // testing read
    std::string ret = rdar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from the MULTI structure";
        Error( "Error", ret);
    }

    PAalp_ = PAalp;
    PSigm_ = PSigm;
    multi_realloc( PAalp, PSigm );

    // get dynamic data from DATACH file
    for( ii=0; ii<dCH->nPH; ii++)
        pm.L1[ii] = dCH->nDCinPH[ii];

    for( ii=0; ii<dCH->nIC*dCH->nDC; ii++)
        pm.A[ii] = dCH->A[ii];

    if( pm.EZ )
    { long int iZ=-1;
        for(  ii=0; ii<dCH->nDC; ii++ )
            if( dCH->ccIC[ii] == IC_CHARGE )
                break;
        if( ii< dCH->nDC )
        { iZ = ii;
            for( ii=0; ii<dCH->nDC; ii++)
                pm.EZ[ii] = pm.A[pm.N*ii+iZ];
        }
    }

    for( ii=0; ii< dCH->nIC; ii++ )
    { pm.Awt[ii]  = dCH->ICmm[ii]*1e3;
        fillValue(pm.SB[ii], ' ', MaxICN );
        len = strlen(dCH->ICNL[ii]);
        //len = min(  len,MaxICN);
        copyValues( pm.SB[ii], dCH->ICNL[ii], std::min<size_t>(len,MAXICNAME));
        pm.SB[ii][MaxICN] = dCH->ccIC[ii];
        pm.ICC[ii] =  dCH->ccIC[ii];
    }

    if( std::fabs(dCH->DCmm[0]) < 1e-32 )  // Restore DCmm if skipped from the DCH file
        for( long int jj=0; jj< dCH->nDC; jj++ )  // Added by DK on 03.03.2007
        {
            dCH->DCmm[jj] = 0.0;
            for( ii=0; ii< dCH->nIC; ii++ )
                dCH->DCmm[jj] += dCH->ICmm[ii]*dCH->A[jj*dCH->nIC+ii];
        }

    for( ii=0; ii< dCH->nDC; ii++ )
    {
        pm.MM[ii] = dCH->DCmm[ii]*1e3;
        pm.DCC[ii] = dCH->ccDC[ii];
        len =strlen(dCH->DCNL[ii]);
        //len = min(  len,MaxDCN);
        copyValues( pm.SM[ii], dCH->DCNL[ii], std::min<size_t>(len,MAXDCNAME) );
    }

    for( ii=0; ii< dCH->nPH; ii++ )
    {
        len =strlen(dCH->PHNL[ii]);
        //len = min(  len,MaxPHN);
        fillValue( pm.SF[ii], ' ', MAXPHNAME+MAXSYMB );
        copyValues( pm.SF[ii]+MAXSYMB, dCH->PHNL[ii], std::min<size_t>(len,MAXPHNAME) );
        pm.SF[ii][0] = dCH->ccPH[ii];
        pm.PHC[ii] = dCH->ccPH[ii];
    }

    // !!!!  copyValues( pm.DCCW, dCH->ccDCW, dCH->nDC);
    // set up DCCW
    ConvertDCC();

    //dynamic data
    io_formats::TReadArrays<TIO>   rddar( 80, MULTI_dynamic_fields, in_format);

    // set up array flags for permanent fields

    if( !( pm.FIs > 0 && pm.Ls > 0 ) )
    {
        rddar.setNoAlws( (long int)(f_sMod ));
        rddar.setNoAlws( f_LsMod );
        rddar.setNoAlws( f_LsMdc );
    }
    if( PSigm == S_OFF )
    {
        rddar.setNoAlws( f_Sigw );
        rddar.setNoAlws( f_Sigg );
    }
    if( !( pm.FIat > 0 &&  pm.FIs > 0 ) )
    { /* ADSORPTION AND ION EXCHANGE */
        rddar.setNoAlws( f_Nfsp );
        rddar.setNoAlws( f_MASDT );
        rddar.setNoAlws( f_C1 );
        rddar.setNoAlws( f_C2 );
        rddar.setNoAlws( f_C3 );
        rddar.setNoAlws( f_pCh );
        rddar.setNoAlws( f_SATX );
        rddar.setNoAlws( f_MASDJ );
        rddar.setNoAlws( f_SCM );
        rddar.setNoAlws( f_SACT );
        rddar.setNoAlws( f_DCads );
    }

    // Read dynamic arrays
    nfild = rddar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        { case f_sMod: if( !pm.sMod )
                Error( "Error", "Array sMod is not used in this problem");
            rddar.readArray( "sMod" , pm.sMod[0], pm.FIs, 8 );
            break;
        case f_LsMod:{ if( !pm.LsMod )
                Error( "Error", "Array LsMod is not used in this problem");
            rddar.readArray( "LsMod" , pm.LsMod, pm.FIs*3) ;
            long int LsModSum;
            long int LsIPxSum;
            getLsModsum( LsModSum, LsIPxSum );
            if(LsIPxSum )
            {
                rddar.readNext( "IPxPH");
                alloc_IPx(LsIPxSum);
                rddar.readArray( "IPxPH", pm.IPx,  LsIPxSum);
            }
            if(LsModSum )
            {
                rddar.readNext( "PMc");
                alloc_PMc(LsModSum);
                rddar.readArray( "PMc", pm.PMc,  LsModSum);
            }
            break;
        }
        case f_LsMdc: { if( !pm.LsMdc )
                Error( "Error", "Array LsMdc not used in this problem");
            rddar.readArray( "LsMdc" , pm.LsMdc, pm.FIs*3 );
            long int LsMdcSum;
            long int LsMsnSum;
            long int LsSitSum;
            getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
            if(LsMdcSum )
            {
                rddar.readNext( "DMc");
                alloc_DMc(LsMdcSum);
                rddar.readArray( "DMc", pm.DMc,  LsMdcSum);
            }
            if(LsMsnSum )
            {
                rddar.readNext( "MoiSN");
                alloc_MoiSN(LsMsnSum);
                alloc_SitFr(LsSitSum);
                fillValue( pm.SitFr, 0., LsSitSum );
                rddar.readArray( "MoiSN", pm.MoiSN,  LsMsnSum);
            }
            break;
        }
        case f_LsMdc2:
        {
            if( !pm.LsMdc2 )
                Error( "Error", "Array LsMdc2 not used in this problem");
            rddar.readArray(  "LsMdc2", pm.LsMdc2, pm.FIs*3);
            long int DQFcSum, rcpcSum;
            getLsMdc2sum( DQFcSum, rcpcSum );
            if(DQFcSum )
            {
                rddar.readNext( "DQFc");
                alloc_DQFc(DQFcSum);
                rddar.readArray(  "DQFc", pm.DQFc,  DQFcSum);
            }
            break;
        }
        case f_LsPhl:
        { if( !pm.LsPhl )
                Error( "Error", "Array LsPhl not used in this problem");
            rddar.readArray(  "LsPhl",  pm.LsPhl, pm.FI*2);
            long int PhLinSum, lPhcSum;
            getLsPhlsum( PhLinSum,lPhcSum );

            if(PhLinSum )
            {
                rddar.readNext( "PhLin");
                alloc_PhLin(PhLinSum);
                rddar.readArray(  "PhLin", &pm.PhLin[0][0], PhLinSum*2);
            }
            if(lPhcSum )
            {
                rddar.readNext( "lPhc");
                alloc_lPhc(lPhcSum);
                rddar.readArray(  "lPhc", pm.lPhc,  lPhcSum);
            }
            break;
        }
        case f_SorMc:
            rddar.readArray(  "SorMc", pm.SorMc, pm.FIs*16);
            break;
            // TSorpMod stuff
        case f_LsISmo:
        { if( !pm.LsISmo )
                Error( "Error", "Array LsISmo not used in this problem");
            rddar.readArray(  "LsISmo",  pm.LsISmo, pm.FIs*4);

            long int IsoCtSum, IsoScSum;
            long int IsoPcSum, xSMdSum;
            getLsISmosum( IsoCtSum,IsoScSum,IsoPcSum, xSMdSum );

            if(xSMdSum )
            {
                rddar.readNext( "xSMd");
                alloc_xSMd(xSMdSum);
                rddar.readArray(  "xSMd", pm.xSMd, xSMdSum);
            }
            if(IsoPcSum )
            {
                rddar.readNext( "IsoPc");
                alloc_IsoPc(IsoPcSum);
                rddar.readArray(  "IsoPc", pm.IsoPc,  IsoPcSum);
            }
            if(IsoScSum )
            {
                rddar.readNext( "IsoSc");
                alloc_IsoSc(IsoScSum);
                rddar.readArray(  "IsoSc", pm.IsoSc, IsoScSum);
            }
            if(IsoCtSum )
            {
                rddar.readNext( "IsoCt");
                alloc_IsoCt(IsoCtSum);
                rddar.readArray(  "IsoCt", pm.IsoCt,  IsoCtSum, 1L);
            }
            break;
        }
        case f_LsESmo:
        {
            if( !pm.LsESmo )
                Error( "Error", "Array LsESmo not used in this problem");
            rddar.readArray(  "LsESmo",  pm.LsESmo, pm.FIs*4);
            long int EImcSum, mCDcSum;
            getLsESmosum( EImcSum, mCDcSum );

            if(EImcSum )
            {
                rddar.readNext( "EImc");
                alloc_EImc(EImcSum);
                rddar.readArray(  "EImc", pm.EImc, EImcSum);
            }
            if(mCDcSum )
            {
                rddar.readNext( "mCDc");
                alloc_mCDc( mCDcSum );
                rddar.readArray(  "mCDc", pm.mCDc,  mCDcSum);
            }
            break;
        }
            // TKinMet stuff
        case f_kMod:
            rddar.readArray(  "kMod", pm.kMod[0], pm.FI, 6L);
            break;
        case f_LsKin:
        {
            if( !pm.LsKin )
                Error( "Error", "Array LsKin not used in this problem");
            rddar.readArray(  "LsKin",  pm.LsKin, pm.FI*6);

            long int xSKrCSum, ocPRkC_feSArC_Sum;
            long int rpConCSum, apConCSum, AscpCSum;
            getLsKinsum( xSKrCSum, ocPRkC_feSArC_Sum, rpConCSum, apConCSum, AscpCSum );
            if(xSKrCSum )
            {
                rddar.readNext( "xSKrC");
                alloc_xSKrC(xSKrCSum);
                rddar.readArray(  "xSKrC", pm.xSKrC, xSKrCSum);
            }
            if(ocPRkC_feSArC_Sum )
            {
                rddar.readNext( "ocPRkC");
                alloc_ocPRkC(ocPRkC_feSArC_Sum);
                alloc_feSArC(ocPRkC_feSArC_Sum);
                rddar.readArray(  "ocPRkC", &pm.ocPRkC[0][0],  ocPRkC_feSArC_Sum*2);
                rddar.readNext( "feSArC");
                rddar.readArray(  "feSArC", pm.feSArC, ocPRkC_feSArC_Sum);
            }
            if(rpConCSum )
            {
                rddar.readNext( "rpConC");
                alloc_rpConC(rpConCSum);
                rddar.readArray(  "rpConC", pm.rpConC,  rpConCSum);
            }
            if(apConCSum )
            {
                rddar.readNext( "apConC");
                alloc_apConC(apConCSum);
                rddar.readArray(  "apConC", pm.apConC, apConCSum);
            }
            if(AscpCSum )
            {
                rddar.readNext( "AscpC");
                alloc_AscpC(AscpCSum);
                rddar.readArray(  "AscpC", pm.AscpC,  AscpCSum);
            }
            break;
        }
        case f_LsUpt:
        {
            if( !pm.LsUpt )
                Error( "Error", "Array LsUpt not used in this problem");
            rddar.readArray(  "LsUpt",  pm.LsUpt, pm.FIs*2);

            long int UMpcSum, xICuCSum;
            getLsUptsum( UMpcSum, xICuCSum );
            if(UMpcSum )
            {
                rddar.readNext( "UMpcC");
                alloc_UMpcC(UMpcSum);
                rddar.readArray(  "UMpcC", pm.UMpcC, UMpcSum);
            }
            break;
        }
        case f_PfFact:  rddar.readArray(  "PfFact", pm.PfFact, pm.FI );
            break;
        case f_xICuC:
        {
            long int UMpcSum, xICuCSum;
            getLsUptsum( UMpcSum, xICuCSum );
            alloc_xICuC(xICuCSum);
            rddar.readArray(  "xICuC", pm.xICuC, xICuCSum );
            break;
        }
        case f_B: rddar.readArray( "B", pm.B,  pm.N);
            break;
        case f_DCCW: rddar.readArray( "DCCW", pm.DCCW,  pm.L, 1);
            break;
        case f_Pparc: rddar.readArray( "Pparc", pm.Pparc,  pm.L);
            break;
        case f_fDQF: rddar.readArray( "fDQF", pm.fDQF,  pm.L);
            break;
        case f_lnGmf: rddar.readArray( "lnGmf", pm.lnGmf,  pm.L);
            break;
        case f_RLC: rddar.readArray( "RLC", pm.RLC, pm.L, 1 );
            break;
        case f_RSC: rddar.readArray( "RSC", pm.RSC, pm.L, 1 );
            break;
        case f_DLL: rddar.readArray( "DLL", pm.DLL,  pm.L);
            break;
        case f_DUL: rddar.readArray( "DUL", pm.DUL,  pm.L);
            break;
        case f_Aalp: rddar.readArray( "Aalp", pm.Aalp,  pm.FI);
            break;
        case f_Sigw: if( !pm.Sigw )
                Error( "Error", "Array Sigw not used in this problem");
            rddar.readArray( "Sigw", pm.Sigw,  pm.FI);
            break;
        case f_Sigg: if( !pm.Sigg )
                Error( "Error", "Array Sigg not used in this problem");
            rddar.readArray( "Sigg", pm.Sigg,  pm.FI);
            break;
        case f_YOF: rddar.readArray( "YOF", pm.YOF,  pm.FI);
            break;
        case f_Nfsp: if( !pm.Nfsp )
                Error( "Error", "Array Nfsp not used in this problem");
            rddar.readArray( "Nfsp", &pm.Nfsp[0][0], pm.FIs*pm.FIat);
            break;
        case f_MASDT: if( !pm.MASDT )
                Error( "Error", "Array MASDT not used in this problem");
            rddar.readArray( "MASDT", &pm.MASDT[0][0], pm.FIs*pm.FIat);
            break;
        case f_C1: if( !pm.XcapA )
                Error( "Error", "Array XcapA not used in this problem");
            rddar.readArray( "C1", &pm.XcapA[0][0], pm.FIs*pm.FIat);
            break;
        case f_C2: if( !pm.XcapB )
                Error( "Error", "Array XcapB not used in this problem");
            rddar.readArray( "C2", &pm.XcapB[0][0], pm.FIs*pm.FIat);
            break;
        case f_C3: if( !pm.XcapF )
                Error( "Error", "Array XcapF not used in this problem");
            rddar.readArray( "C3", &pm.XcapF[0][0], pm.FIs*pm.FIat);
            break;
        case f_pCh: if( !pm.Xetaf )
                Error( "Error", "Array Xetaf not used in this problem");
            rddar.readArray( "pCh", &pm.Xetaf[0][0], pm.FIs*pm.FIat);
            break;
        case f_SATX: if( !pm.SATX )
                Error( "Error", "Array SATX not used in this problem");
            rddar.readArray( "SATX", &pm.SATX[0][0], pm.Lads*4);
            break;
        case f_MASDJ: if( !pm.MASDJ )
                Error( "Error", "Array MASDJ not used in this problem");
            rddar.readArray( "MASDJ", &pm.MASDJ[0][0], pm.Lads*DFCN);
            break;
        case f_SCM: if( !pm.SCM )
                Error( "Error", "Array SCM not used in this problem");
            rddar.readArray( "SCM", pm.SCM[0], pm.FIs, pm.FIat );
            break;
        case f_SACT: if( !pm.SATT )
                Error( "Error", "Array SATT not used in this problem");
            rddar.readArray( "SACT", pm.SATT, pm.Lads, 1 );
            break;
        case f_DCads: if( !pm.DCC3 )
                Error( "Error", "Array DCC3 not used in this problem");
            rddar.readArray( "DCads", pm.DCC3, pm.Lads, 1 );
            break;
        case f_pa_DB: rddar.readArray( "pa_DB" , &pa_p->DB, 1);
            break;
        case f_pa_DHB: rddar.readArray("pa_DHB", &pa_p->DHB, 1);
            break;
        case f_pa_EPS: rddar.readArray("pa_EPS" , &pa_p->EPS, 1);
            break;
        case f_pa_DK: rddar.readArray("pa_DK" , &pa_p->DK, 1);
            break;
        case f_pa_DF: rddar.readArray("pa_DF" , &pa_p->DF, 1);
            break;
        case f_pa_DP: rddar.readArray("pa_DP", &pa_p->DP, 1);
            break;
        case f_pa_IIM: rddar.readArray("pa_IIM", &pa_p->IIM, 1);
            break;
        case f_pa_PD: rddar.readArray("pa_PD" , &pa_p->PD, 1);
            break;
        case f_pa_PRD: rddar.readArray("pa_PRD" , &pa_p->PRD, 1);
            break;
        case f_pa_AG: rddar.readArray("pa_AG" , &pa_p->AG, 1);
            break;
        case f_pa_DGC: rddar.readArray("pa_DGC" , &pa_p->DGC, 1);
            break;
        case f_pa_PSM: rddar.readArray("pa_PSM" , &pa_p->PSM, 1);
            break;
        case f_pa_GAR: rddar.readArray("pa_GAR" , &pa_p->GAR, 1);
            break;
        case f_pa_GAH: rddar.readArray("pa_GAH" , &pa_p->GAH, 1);
            break;
        case f_pa_DS: rddar.readArray("pa_DS", &pa_p->DS, 1);
            break;
        case f_pa_XwMin: rddar.readArray("pa_XwMin" , &pa_p->XwMin, 1);
            break;
        case f_pa_ScMin: rddar.readArray("pa_ScMin" , &pa_p->ScMin, 1);
            break;
        case f_pa_DcMin: rddar.readArray("pa_DcMin" , &pa_p->DcMin, 1);
            break;
        case f_pa_PhMin: rddar.readArray("pa_PhMin" , &pa_p->PhMin, 1);
            break;
        case f_pa_ICmin: rddar.readArray("pa_ICmin" , &pa_p->ICmin, 1);
            break;
        case f_pa_PC: rddar.readArray("pa_PC" , &pa_p->PC, 1);
            break;
        case f_pa_DFM: rddar.readArray("pa_DFM" , &pa_p->DFM, 1);
            break;
        case f_pa_DFYw: rddar.readArray("pa_DFYw" , &pa_p->DFYw, 1);
            break;
        case f_pa_DFYaq: rddar.readArray("pa_DFYaq" , &pa_p->DFYaq, 1);
            break;
        case f_pa_DFYid: rddar.readArray("pa_DFYid" , &pa_p->DFYid, 1);
            break;
        case f_pa_DFYr: rddar.readArray("pa_DFYr" , &pa_p->DFYr, 1);
            break;
        case f_pa_DFYh: rddar.readArray("pa_DFYh" , &pa_p->DFYh, 1);
            break;
        case f_pa_DFYc: rddar.readArray("pa_DFYc" , &pa_p->DFYc, 1);
            break;
        case f_pa_DFYs: rddar.readArray("pa_DFYs", &pa_p->DFYs, 1);
            break;
        case f_pa_DW: rddar.readArray("pa_DW", &pa_p->DW , 1);
            break;
        case f_pa_DT: rddar.readArray("pa_DT", &pa_p->DT , 1);
            break;
        case f_pa_GAS: rddar.readArray("pa_GAS", &pa_p->GAS, 1);
            break;
        case f_pa_DG: rddar.readArray("pa_DG" , &pa_p->DG, 1);
            break;
        case f_pa_DNS: rddar.readArray("pa_DNS" , &pa_p->DNS, 1);
            break;
        case f_pa_IEPS: rddar.readArray("pa_IEPS" , &pa_p->IEPS, 1);
            break;
        case f_pKin: rddar.readArray("pKin" , &pm.PLIM, 1);
            break;
        case f_pa_DKIN: rddar.readArray("pa_DKIN" , &pa_p->DKIN, 1);
            break;
        case f_mui: rddar.readArray("mui" , pm.mui, pm.N);
            break;
        case f_muk: rddar.readArray("muk" , pm.muk, pm.FI);
            break;
        case f_muj: rddar.readArray("muj" , pm.muj, pm.L);
            break;
        case f_pa_PLLG: rddar.readArray("pa_PLLG" , &pa_p->PLLG, 1);
            break;
        case f_tMin: rddar.readArray("tMin" , &pm.tMin, 1);
            break;
        case f_dcMod:   rddar.readArray( "dcMod" , pm.dcMod[0], pm.L, 6 );
            break;
        }
        nfild = rddar.findNext();
    }
    // testing read
    ret = rddar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from the MULTY structure";
        Error( "Error", ret);
    }
}


/// Writes Multi to a json/key-value string
/// \param brief_mode - Do not write data items that contain only default values
/// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string TMultiBase::gemipm_to_string( bool addMui, const std::string& test_set_name, bool with_comments, bool brief_mode )
{
    std::stringstream ss;
    write_ipm_format_stream( ss, GEMS3KGenerator::default_type_f, addMui, with_comments, brief_mode, test_set_name );
    return ss.str();
}

/// Reads Multi structure from a json/key-value string
bool TMultiBase::gemipm_from_string( const std::string& data,  DATACH  *dCH, const std::string& test_set_name )
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_ipm_format_stream( ss, GEMS3KGenerator::default_type_f, dCH, test_set_name );
    return true;
}

void  TMultiBase::read_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes  type_f, DATACH  *dCH, const std::string& test_set_name )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonRead in_format( stream, test_set_name, "ipm" );
        from_text_file_gemipm( in_format, dCH );
    }
#else
    {
        io_formats::SimdJsonRead in_format( stream, test_set_name, "ipm");
        from_text_file_gemipm( in_format, dCH );
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueRead in_format( stream );
        from_text_file_gemipm( in_format, dCH );
    }
        break;
    }
}

void  TMultiBase::write_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f,
                                           bool addMui, bool with_comments, bool brief_mode, const std::string& test_set_name )
{
    switch( type_f )
    {
    case GEMS3KGenerator::f_binary:
        break;
    case GEMS3KGenerator::f_json:
    case GEMS3KGenerator::f_thermofun:
#ifdef USE_NLOHMANNJSON
    {
        io_formats::NlohmannJsonWrite out_format( stream, test_set_name );
        to_text_file_gemipm( out_format, addMui, with_comments, brief_mode );
    }
#else
    {
        io_formats::SimdJsonWrite out_format( stream, test_set_name, with_comments );
        to_text_file_gemipm( out_format, addMui, with_comments, brief_mode );
    }
#endif
        break;
    case GEMS3KGenerator::f_key_value:
    case GEMS3KGenerator::f_kv_thermofun:
    {
        io_formats::KeyValueWrite out_format( stream );
        to_text_file_gemipm( out_format, addMui, with_comments, brief_mode );
    }
        break;
    }
}


#ifdef USE_NLOHMANNJSON
template void TMultiBase::from_text_file_gemipm<io_formats::NlohmannJsonRead>( io_formats::NlohmannJsonRead& in_format,  DATACH  *dCH );
template void TMultiBase::to_text_file_gemipm<io_formats::NlohmannJsonWrite>( io_formats::NlohmannJsonWrite& out_format, bool addMui, bool with_comments, bool brief_mode );
#else
template void TMultiBase::from_text_file_gemipm<io_formats::SimdJsonRead>( io_formats::SimdJsonRead& in_format,  DATACH  *dCH );
template void TMultiBase::to_text_file_gemipm<io_formats::SimdJsonWrite>( io_formats::SimdJsonWrite& out_format, bool addMui, bool with_comments, bool brief_mode );
#endif
template void TMultiBase::from_text_file_gemipm<io_formats::KeyValueRead>( io_formats::KeyValueRead& in_format,  DATACH  *dCH );
template void TMultiBase::to_text_file_gemipm<io_formats::KeyValueWrite>( io_formats::KeyValueWrite& out_format, bool addMui, bool with_comments, bool brief_mode );

//=============================================================================
// ms_multi_format.cpp

