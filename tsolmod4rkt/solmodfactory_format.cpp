//-------------------------------------------------------------------
/// \file solmodfactory_format.cpp
///
/// Implementation of writing/reading IPM text I/O files
//
// Copyright (c) 2023 S.Dmytriyeva,D.Kulik
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

#include "solmodfactory.h"
#include "v_detail.h"
#include "io_template.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "io_keyvalue.h"
#include "jsonconfig.h"


const char *_GEMIPM_version_stamp = " GEMS3K v.3.9.6 c.9a8c970";  // interim version, need merge for release v.4.0.0

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
   // TSolMod extern const double bar_to_Pa, m3_to_cm3,  kg_to_g;
   { "LsMdc2",    0 , 0, 0,  "# LsMdc2: [3*FIs] - number of DQF coeffs; reciprocal coeffs per end member" },
   { "LsPhl",    0 , 0, 0,  "# LsPhl: Number of phase links; number of link parameters; [Fi][2]" }
 };


/// Reading structure MULTI (GEM IPM work structure)
template<typename TIO>
void SolModFactory::from_text_file_gemipm( TIO& in_format,  DATACH  *dCH )
{
    long int ii, nfild;
    size_t len;

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
    pm.E = 1;
    pm.PV = 0;
    pm.PSOL = 0;
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
        case f_pa_PE: rdar.readArray("pa_PE" , &pm.E, 1);
            break;
        case f_PV: rdar.readArray("PV" , &pm.PV, 1);
            break;
        case f_PSOL: rdar.readArray("PSOL" , &pm.PSOL, 1);
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
    {
        ret += " - fields must be read from the MULTI structure";
        Error( "Error", ret);
    }

    multi_realloc();

    // get dynamic data from DATACH file
    for( ii=0; ii<dCH->nPH; ii++)
        pm.L1[ii] = dCH->nDCinPH[ii];

    for( ii=0; ii<dCH->nIC*dCH->nDC; ii++)
        pm.A[ii] = dCH->A[ii];

    if( pm.EZ )
    {
        long int iZ=-1;
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
    {
        pm.Awt[ii]  = dCH->ICmm[ii]*1e3;
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
    if( true /*PSigm == S_OFF*/ )
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
        {
        case f_sMod:
            if( !pm.sMod )
                Error( "Error", "Array sMod is not used in this problem");
            rddar.readArray( "sMod" , pm.sMod[0], pm.FIs, 8 );
            break;
        case f_LsMod:  {
            if( !pm.LsMod )
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
        case f_LsMdc: {
            if( !pm.LsMdc )
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
        case f_LsPhl: {
            if( !pm.LsPhl )
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
        case f_YOF: rddar.readArray( "YOF", pm.YOF,  pm.FI);
            break;
        case f_pKin: rddar.readArray("pKin" , &pm.PLIM, 1);
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

/// Reads Multi structure from a json/key-value string
bool SolModFactory::gemipm_from_string( const std::string& data,  DATACH  *dCH, const std::string& test_set_name )
{
    if( data.empty() )
        return false;

    std::stringstream ss;
    ss.str(data);
    read_ipm_format_stream( ss, GEMS3KGenerator::default_type_f, dCH, test_set_name );
    return true;
}

void  SolModFactory::read_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes  type_f, DATACH  *dCH, const std::string& test_set_name )
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

void SolModFactory::getLsModsum( long int& LsModSum, long int& LsIPxSum )
{
    LsModSum = 0;
    LsIPxSum = 0;
    for(long int i=0; i<pm.FIs; i++)
    {
        LsModSum += (pm.LsMod[i*3]*pm.LsMod[i*3+2]);
        LsIPxSum += (pm.LsMod[i*3]*pm.LsMod[i*3+1]);
    }
}

void SolModFactory::getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum )
{
    LsMdcSum = 0;
    LsMsnSum = 0;
    LsSitSum = 0;

    for(long int i=0; i<pm.FIs; i++)
    {
        LsMdcSum += (pm.LsMdc[i*3]*pm.L1[i]);
        LsMsnSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]*pm.L1[i]);
        LsSitSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]);
    }
}

// dimensions from LsPhl array
void SolModFactory::getLsPhlsum( long int& PhLinSum,long int& lPhcSum )
{
    PhLinSum = 0;
    lPhcSum = 0;

    for(long int i=0; i<pm.FI; i++)
    {
        PhLinSum += (pm.LsPhl[i*2]);
        lPhcSum += (/*pm.LsPhl[i*2]**/pm.LsPhl[i*2+1]);

    }
}

// dimensions from LsMdc2 array
void SolModFactory::getLsMdc2sum( long int& DQFcSum,long int& rcpcSum )
{
    DQFcSum = 0;
    rcpcSum = 0;

    for(long int i=0; i<pm.FIs; i++)
    {
        DQFcSum += (pm.LsMdc2[i*3]*pm.L1[i]);
        //       rcpcSum += (pm.LsMdc2[i*3+1]*pm.L1[i]);
    }
}

//---------------------------------------------------------//
/// Writing structure MULTI ( free format file  )
void SolModFactory::to_text_file( const std::string& path, bool append )
{
    std::ios::openmode mod = std::ios::out;
    if( append )
        mod = std::ios::out|std::ios::app;
    std::fstream ff(GemsSettings::with_directory(path), mod );
    ErrorIf( !ff.good() , path, "Fileopen error");

    io_formats::KeyValueWrite out_format( ff );
    out_format.put_head( "", "ipm");
    io_formats::TPrintArrays<io_formats::KeyValueWrite>  prar( 0, {}, out_format );

    if( append )
        prar.writeComment( true,"\nNext record" );
    prar.writeComment( true, char_array_to_string(pm.stkey, EQ_RKLEN)+"\n" );
    //  TProfil::pm->pa.p.write(ff);

    //prar.writeArray( "Short_PARAM",  &base_param()->PC, 10L );
    //prar.writeArray( "Double_PARAM",  &base_param()->DG, 28L );
    prar.writeArray( "Short_Const",  &pm.N, 39L );
    prar.writeArray(  "Double_Const",  &pm.TC, 53, 20 );
    // prar.writeArray(  "Add_Double_Const",  &pm.XwMinM, 12, 20 );
    prar.writeArray(  "EpsW", pm.epsW, 5);
    prar.writeArray(  "EpsWg", pm.epsWg, 5);
    prar.writeArray(  "DenW", pm.denW, 5);
    prar.writeArray(  "DenWg", pm.denWg, 5);
    prar.writeComment( true, std::string("Error Code ")+ pm.errorCode);
    prar.writeComment( true, std::string("Error Message") + pm.errorBuf);

    //dynamic values

    // Part 1
    prar.writeArray(  "L1", pm.L1,  pm.FI);
    prar.writeArray(  "DUL", pm.DUL,  pm.L);
    prar.writeArray(  "DLL", pm.DLL,  pm.L);
    prar.writeArray(  "Vol", pm.Vol,  pm.L);
    if(pm.V0)
        prar.writeArray("V0",pm.V0, pm.L);

    prar.writeArray(  "Pparc", pm.Pparc,  pm.L);
    prar.writeArray(  "MM", pm.MM,  pm.L);
    prar.writeArray(  "Awt", pm.Awt, pm.N);
    prar.writeArray(  "A", pm.A,  pm.N*pm.L);
    prar.writeArray(  "G", pm.G,  pm.L);
    prar.writeArray(  "G0", pm.G0,  pm.L);
    prar.writeArray(  "lnGam", pm.lnGam,  pm.L);
    prar.writeArray(  "lnGmo", pm.lnGmo,  pm.L);
    prar.writeArray(  "B", pm.B,  pm.N);
    
    prar.writeArray(  "XF", pm.XF,  pm.FI);
    prar.writeArray(  "X", pm.X,  pm.L);
    prar.writeArray(  "YOF", pm.YOF,  pm.FI);
    prar.writeArray(  "lnGmM", pm.lnGmM,  pm.L);
    prar.writeArray(  "fDQF", pm.fDQF,  pm.L);
    prar.writeArray(  "FVOL", pm.FVOL,  pm.FI);
    prar.writeArray(  "FWGT", pm.FWGT,  pm.FI);

    if( pm.L > 0 )
    {
        prar.writeArray(  "Wx", pm.Wx,  pm.L);
        prar.writeArray(  "Gamma", pm.Gamma,  pm.L);
        prar.writeArray(  "lnGmf", pm.lnGmf,  pm.L);
    }

    // Part 2  not always required arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        prar.writeArray(  "XFA", pm.XFA,  pm.FIs);
    }
    if( pm.LO > 1 )
    {
        prar.writeArray(  "Y_m", pm.Y_m,  pm.L);
    }
    if( pm.E )
    {
        prar.writeArray(  "EZ", pm.EZ,  pm.L);
    }

    // Part 3  new Phase definition
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        prar.writeArray(  "sMod", &pm.sMod[0][0], pm.FIs,8L);
        prar.writeArray(  "LsMod", pm.LsMod, pm.FIs*3);
        long int LsModSum;
        long int LsIPxSum;
        getLsModsum( LsModSum, LsIPxSum );
        prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
        prar.writeArray(  "PMc", pm.PMc,  LsModSum);
        long int LsMdcSum;
        long int LsMsnSum;
        long int LsSitSum;
        prar.writeArray(  "LsMdc", pm.LsMdc, pm.FIs*3);
        getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
        prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);
        prar.writeArray(  "MoiSN", pm.MoiSN,  LsMsnSum);
        prar.writeArray(  "SitFr", pm.SitFr,  LsSitSum);
        long int DQFcSum, rcpcSum;
        getLsMdc2sum( DQFcSum, rcpcSum );
        prar.writeArray(  "LsMdc2", pm.LsMdc2, pm.FIs*3);
        prar.writeArray(  "DQFc", pm.DQFc,  DQFcSum);
        //      prar.writeArray(  "rcpc", pm.rcpc,  rcpcSum);
        long int PhLinSum, lPhcSum;
        getLsPhlsum( PhLinSum,lPhcSum );
        prar.writeArray(  "LsPhl", pm.LsPhl, pm.FI*2);
        prar.writeArray(  "PhLin", &pm.PhLin[0][0], PhLinSum*2);
        prar.writeArray(  "lPhc", pm.lPhc,  lPhcSum);

        prar.writeArray(  "lnDQFt", pm.lnDQFt, pm.Ls);
        prar.writeArray(  "lnRcpt", pm.lnRcpt, pm.Ls);
        prar.writeArray(  "lnExet", pm.lnExet, pm.Ls);
        prar.writeArray(  "lnCnft", pm.lnCnft, pm.Ls);
    }

    if(pm.H0)
        prar.writeArray("H0",pm.H0, pm.L);
    if(pm.A0)
        prar.writeArray("A0",pm.A0, pm.L);
    if(pm.U0)
        prar.writeArray("U0",pm.U0, pm.L);
    if(pm.S0)
        prar.writeArray("S0",pm.S0, pm.L);
    if(pm.Cp0)
        prar.writeArray("Cp0",pm.Cp0, pm.L);

    prar.writeArray(  "VPh", &pm.VPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "GPh", &pm.GPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "HPh", &pm.HPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "SPh", &pm.SPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "CPh", &pm.CPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "APh", &pm.APh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "UPh", &pm.UPh[0][0], pm.FIs*MIXPHPROPS);
}

void SolModFactory::unpackDataBr( bool uPrimalSol )
{
    long int ii;
    pmp->kTau = CNode->Tm;
    pmp->kdT = CNode->dt;

    pmp->TCc = CNode->TK-C_to_K;
    pmp->Tc = CNode->TK;
    pmp->Pc  = CNode->P/bar_to_Pa;
    pmp->VXc = CNode->Vs/1.e-6; // from cm3 to m3
    // Obligatory arrays - always unpacked!
    for( ii=0; ii<CSD->nDCb; ii++ ) {
        pmp->DUL[ CSD->xdc[ii] ] = CNode->dul[ii];
        pmp->DLL[ CSD->xdc[ii] ] = CNode->dll[ii];
        if( pmp->DUL[ CSD->xdc[ii] ] < pmp->DLL[ CSD->xdc[ii] ] )   {
            Error("unpackDataBr", std::string("Upper kinetic restriction less than the lower one for DC&RC")
                  +char_array_to_string( pmp->SM[CSD->xdc[ii]], MAXDCNAME));
        }
    }
    for( ii=0; ii<CSD->nICb; ii++ )  {
        pmp->B[ CSD->xic[ii] ] = CNode->bIC[ii];
        if( ii < CSD->nICb-1 && pmp->B[ CSD->xic[ii] ] <  1e-17 /*multi_ptr()->base_param()->DB*/ ) {
            Error("unpackDataBr", std::string("Bulk mole amount of IC ")+
                  char_array_to_string(pmp->SB[CSD->xic[ii]], 6)+" is "+
                    std::to_string(pmp->B[ CSD->xic[ii] ])+" - out of range" );
        }
    }
    //    for( ii=0; ii<CSD->nPHb; ii++ ) {
    //        if( CSD->nAalp > 0 )
    //            pmp->Aalp[ CSD->xph[ii] ] = CNode->aPH[ii]/kg_to_g;
    //        pmp->Falp[ CSD->xph[ii] ] = CNode->omPH[ii];
    //    }

    for( ii=0; ii<CSD->nPHb; ii++ ) {
        pmp->XF[ CSD->xph[ii] ] = CNode->xPH[ii];
    }

    if( !uPrimalSol ) {
        //  Using primal solution retained in the MULTI structure instead -
        ; // the primal solution data from the DATABR structure are not unpacked

        //   pmp->IT = 0;
    }
    else {   // Unpacking primal solution provided in the node DATABR structure
        pmp->IT = CNode->IterDone; // ?  pmp->ITF+pmp->IT;
        //pmp->IT = 0;
        pmp->MBX = CNode->Ms;
        pmp->IC = CNode->IC;
        pmp->Eh = CNode->Eh;
        for( ii=0; ii<CSD->nDCb; ii++ ) {
            pmp->X[ CSD->xdc[ii] ] = CNode->xDC[ii];
        }

        //        for( ii=0; ii<CSD->nPSb; ii++ )  {
        //            pmp->PUL[ CSD->xph[ii] ] = CNode->amru[ii];
        //            pmp->PLL[ CSD->xph[ii] ] = CNode->amrl[ii];
        //        }

        for( ii=0; ii<CSD->nPHb; ii++ )  {
            pmp->XF[ CSD->xph[ii] ] = Ph_Moles(ii);
            pmp->FVOL[ CSD->xph[ii] ] = Ph_Volume(ii)*m3_to_cm3;
            pmp->FWGT[ CSD->xph[ii] ] = Ph_Mass(ii)*kg_to_g;
        }

        //        for( long int k=0; k<CSD->nPSb; k++ )
        //            for(long int i=0; i<CSD->nICb; i++ ) {
        //                long int dbr_ndx= (k*CSD->nICb)+i,
        //                        mul_ndx = ( CSD->xph[k]*CSD->nIC )+ CSD->xic[i];
        //                pmp->BF[ mul_ndx ] = CNode->bPS[dbr_ndx];
        //            }

        for( ii=0; ii<CSD->nPSb; ii++ ) {
            pmp->XFA[ CSD->xph[ii] ] = CNode->xPA[ii];
        }

        for( ii=0; ii<CSD->nDCb; ii++ )  {
            pmp->Gamma[ CSD->xdc[ii] ] = CNode->gam[ii];
        }

        //        long int jb, je = 0;
        //        for( long int k=0; k<pmp->FIs; k++ ) { // loop on solution phases
        //            jb = je;
        //            je += pmp->L1[ k ];
        //            // Load activity coeffs for phases-solutions
        //            for( ii=jb; ii<je; ii++ )  {
        //                pmp->lnGam[ii] =  multi_ptr()->PhaseSpecificGamma( ii, jb, je, k, 1L );
        //            }
        //        }
    }
    //  End
}

//Retrieves the current phase volume in m3 ( xph is DBR phase index) in the reactive sub-system.
// Works both for multicomponent and for single-component phases. Returns 0.0 if the phase mole amount is zero.
double  SolModFactory::Ph_Volume( const long int xBR ) const
{
    double vol;
    if( xBR < CSD->nPSb )
    {    // Phase-solution
        vol = CNode->vPS[xBR];
        // Perhaps not yet accounting for the volume of mixing!
    }
    else
    {
        long int xdc = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
        vol = pm.V0[xdc];//  /1e5
        vol *= CNode->xDC[DC_xCH_to_xDB(xdc)];
    }
    return vol;
}

//Retrieves the current phase amount in moles (xph is DBR phase index) in the reactive sub-system.
double  SolModFactory::Ph_Moles( const long int xBR ) const
{
    double mol;
    mol = CNode->xPH[xBR];

    return mol;
}

// Retrieves the phase mass in kg (xph is DBR phase index).
// Works for multicomponent and for single-component phases. Returns 0.0 if phase amount is zero.
double  SolModFactory::Ph_Mass( const long int xBR ) const
{
    double mass;
    if( xBR < CSD->nPSb )
        mass = CNode->mPS[xBR];
    else
    {
        long int xDC = Phx_to_DCx( Ph_xDB_to_xCH( xBR ));
        mass = CNode->xDC[ DC_xCH_to_xDB(xDC) ] * CSD->DCmm[xDC];
    }
    return mass;
}

// Converts the DC DCH index into the DC DBR index
// or returns -1 if this DC is not used in the data bridge
long int SolModFactory::DC_xCH_to_xDB( const long int xCH ) const
{
    for(long int ii = 0; ii<CSD->nDCb; ii++ )
        if( CSD->xdc[ii] == xCH )
            return ii;
    return -1;
}

// Returns the DCH index of the first DC belonging to the phase with DCH index Phx
long int  SolModFactory::Phx_to_DCx( const long int Phx ) const
{
    long int k, DCx = 0;
    for( k=0; k<CSD->nPHb; k++ )
    {
        if( k == Phx )
            break;
        DCx += CSD->nDCinPH[ k];
    }
    return DCx;
}


#ifdef USE_NLOHMANNJSON
template void SolModFactory::from_text_file_gemipm<io_formats::NlohmannJsonRead>( io_formats::NlohmannJsonRead& in_format,  DATACH  *dCH );
#else
template void SolModFactory::from_text_file_gemipm<io_formats::SimdJsonRead>( io_formats::SimdJsonRead& in_format,  DATACH  *dCH );
#endif
template void SolModFactory::from_text_file_gemipm<io_formats::KeyValueRead>( io_formats::KeyValueRead& in_format,  DATACH  *dCH );

//--------------------- end of tsolmod_multi_format.cpp ---------------------------


