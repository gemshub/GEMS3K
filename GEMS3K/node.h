//-------------------------------------------------------------------
// $Id$
/// \file node.h
/// Declaration of TNode class that implements a simple C/C++ interface
/// between GEMS3K and another code.
//
/// \class TNode node.h
/// Implements a simple C/C++ interface between GEM IPM and FMT codes.
/// Works with DATACH and work DATABR structures without using
/// the TNodearray class.
//
// Copyright (c) 2006-2021 S.Dmytriyeva, D.Kulik, G.Kosakowski, G.D.Miron, A.Leal
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
/// \mainpage GEMS3K Solver of GeoChemical Equilibria and its TNode class interface.
///
/// GEMS3K (formerly GEMIPM2K) is a C/C++ code implementing the efficient numerical kernel
/// IPM-3 of the GEM-Selektor v.3 package for geochemical  thermodynamic modeling of complex
/// heterogeneous multicomponent-multiphase  systems. GEMS3K results from substantial
/// improvements of convex programming Gibbs energy minimization algorithms achieved since
/// 2000, when development and support of GEMS was taken over by LES in Paul Scherrer Institut
/// (since 2008 jointly with IGP ETHZ) by the GEMS Development Team, currently consisting of
/// D.Kulik (lead), G.D.Miron, A.Leal, S.Dmytrieva, S.Nichenko, G. Kosakowski, A.Yapparova.
///
/// Standalone variant of the GEMS3K code can be coupled to reactive mass transport simulation
/// codes, also those running on high-performance computers. Input files (in text format) for
/// GEMS3K can be exported with a few mouse-clicks from the GEM-Selektor v.3  code, or prepared
/// manually using a simple ASCII text editor. Data exchange with the mass transport part of
/// the coupled code can be implemented in computer memory using TNode class functions.
///
/// The standalone GEMS3K code is licensed as the open-source software in order to promote its
/// broad application in hydrothermal-/ waste geochemistry and related research communities.
/// Other potential areas of GEMS3K application include coupled parameter-fitting codes
/// and phase diagram tools.
///
/// Copyright (C) 2012 GEMS Development Team
/// Available on web at http://gems.web.psi.ch/GEMS3K

#ifndef NODE_H
#define NODE_H

// #include "allan_ipm.h"
#include "datach_api.h"
#include "activities.h"
#include "v_detail.h"
#include "v_service.h"

#include "ms_multi.h"

#ifdef USE_THERMOFUN
#include "ThermoFun/ThermoFun.h"
#endif

class TActivity;
class TKinetics;

extern const double bar_to_Pa,
m3_to_cm3,
kg_to_g;

const long int
MaxICnameLength =      6,      // IC name length
MaxDCnameLength =      16,     // DC name length
MaxPHnameLength =      16;     // PH name length

/// \class TNode (GEMS3K kernel)
/// Implements a simple C/C++ interface between GEM IPM and FMT codes.
/// Works with DATACH and work DATABR structures without using
/// the TNodearray class.
class TNode
{
    friend class TNodeArray;

    /// Place for the *dbr. I/O file name
    std::string dbr_file_name;
    /// Last error message logged to the ipmlog file
    std::string ipmlog_error;

    /// The same GEMS3K input set (currently GEM-Selektor asks for the "set" name as .lst file name)
    std::string current_input_set_name = "";
    /// The same GEMS3K output set (currently GEM-Selektor asks for the "set" name as .lst file name)
    std::string current_output_set_name = "";
    /// false - Ensures the re-reading of the system properties into GEM IMP data structure
    bool load_thermodynamic_data = true; // false; // internal value

protected:

    MULTI* pmm;  ///< \protected Pointer to GEM IPM work data structure (ms_multi.h)
    DATACH* CSD;  ///< Pointer to chemical system data structure CSD (DATACH)
    DATABR* CNode;  ///< Pointer to a work node data bridge structure (node)

#ifdef USE_THERMOFUN
    std::unique_ptr<ThermoFun::ThermoEngine> thermo_engine;
#endif
    std::string thermo_json_string;

    // These pointers are only used in standalone GEMS3K programs
    TMultiBase* multi_base = nullptr;     // GEM IPM3 implementation class
    std::shared_ptr<TMultiBase> internal_multi;
    //    TAllan *ipm;       // Allan's GEM IPM implementation class
    // more speciation algorithms classes, when provided
    std::shared_ptr<TActivity> atp;   // Activity term class
    std::shared_ptr<TKinetics> kip;  // MW reaction kinetics class
    ///< used for exchanging input data and results between FMT and GEM IPM
    ACTIVITY* AiP; ///< Pointer to DC activities in phases and related properties


    // These four values are set by the last GEM_run() call
    double CalcTime;  ///< \protected GEMIPM2 calculation time, s
    long int
    PrecLoops,    ///< \protected Number of performed IPM-2 precision refinement loops
    NumIterFIA,   ///< \protected Total Number of performed FIA entry iterations
    NumIterIPM;   ///< \protected Total Number of performed IPM main iterations

    TMultiBase* multi_ptr() const
    {
        return multi_base;
    }

public:

    /// Default logger for TNode class
    static std::shared_ptr<spdlog::logger> node_logger;
    /// Default logger for ipmLogFile
    static std::shared_ptr<spdlog::logger> ipmlog_file;

    /// Constructor of the class instance in memory for standalone GEMS3K or coupled program
    TNode();
    TNode( const TNode& otherNode );

    /// Destructor
    virtual ~TNode();

    // Typical sequence for using TNode class ----------------------------------
    /// (1)
    /// Initialization of GEM IPM2 data structures in coupled RMT-GEM programs
    ///  that use GEMS3K module. Also reads in the IPM, DCH and DBR text input files
    ///  in key-value, json or binary format. Parameters:
    ///  ipmfiles_lst_name - name of a text file that contains:
    ///    " -j | -t |-b <DCH_DAT file name> <IPM_DAT file name> [<>] <dataBR file name>
    ///    or " -f <DCH_DAT file name> <IPM_DAT file name> <ThermoFun JSON format file> <dataBR file name>
    ///  dbfiles_lst_name - name of a text file that contains:
    ///    <dataBR  file name1>, ... , <dataBR file nameN> "
    ///    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
    ///    exist in the same directory where the ipmfiles_lst_name file is located.
    ///    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
    ///    and must contain valid initial chemical systems (of the same structure
    ///    as described in the DCH_DAT file) to set up the initial state of the FMT
    ///    node array.
    ///  If -t flag or no flag is specified then all data files must be in key-value text
    ///    (ASCII) format (and file names must have .dat extension);
    ///  If -j and -f flag is specified then all data files must be in JSON format (and file names
    ///    must have .json extension);
    ///  If -f flag is specified then the use of ThermoFun along with GEMS3K in place of the interpolation
    ///  of lookup arrays for standard thermodynamic data for substances;
    ///  if -b flag is specified then all data files are assumed to be binary (little-endian)
    ///    files.
    ///  @returns:
    ///    0: OK;
    ///    1: GEM IPM read file or input file format error;
    ///   -1: System error (e.g. memory allocation).
    long int  GEM_init( const char *ipmfiles_lst_name );

    /// (1a)
    /// Initialization of GEM IPM3 data structures in coupled programs that use GEMS3K module.
    /// Also reads the input data from the IPM, DCH and one DBR JSON input strings
    /// (e.g. exported from GEM-Selektor or retrieved from JSON database).
    /// Parameters:
    ///  @param dch_json -  DATACH - the Data for CHemistry data structure as a json/key-value string
    ///  @param ipm_json -  Parameters and settings for GEMS3K IPM-3 algorithm as a json/key-value string
    ///  @param dbr_json -  DATABR - the data bridge structure as a json/key-value string
    ///  @param fun_json -  ThermoFun data structure as a json string
    ///  @returns:
    ///    0: OK;
    ///    1: GEM IPM read file or input file format error;
    ///   -1: System error (e.g. memory allocation).
    long int  GEM_init( std::string dch_json, std::string ipm_json,
                        std::string dbr_json, std::string fun_json = "" );

    /// String i/o functions
    /// Writes CSD (DATACH structure) to a json/key-value string
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    std::string datach_to_string( bool with_comments = true, bool brief_mode = false ) const
    {
        return dbr_dch_api::datach_to_string(current_output_set_name, CSD, with_comments, brief_mode);
    }
    /// Reads CSD (DATACH structure) from a json/key-value string
    bool datach_from_string( std::string data )
    {
        return dbr_dch_api::datach_from_string(current_input_set_name, CSD, data);
    }
    /// Writes work node (DATABR structure) to a json/key-value string
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    std::string databr_to_string( bool with_comments = true, bool brief_mode = false ) const
    {
#ifndef NODEARRAYLEVEL
        CNode->NodeStatusFMT = No_nodearray;
#endif
        return dbr_dch_api::databr_to_string(current_output_set_name, CSD, CNode, with_comments, brief_mode);
    }
    /// Reads work node (DATABR structure) from a json/key-value string
    bool databr_from_string( std::string data )
    {
        return dbr_dch_api::databr_from_string(current_input_set_name, CSD, CNode, data);
    }
    /// Writes Multi to a json/key-value string
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    std::string gemipm_to_string( bool addMui, bool with_comments = true, bool brief_mode = false )
    {
        return  multi_ptr()->gemipm_to_string( addMui, current_output_set_name, with_comments, brief_mode );
    }

    /// Return current system name
    std::string system_id() const;

    /// Return current set name
    const std::string& input_set_name() const
    {
        return current_input_set_name;
    }

    /// Return current set name
    const std::string& output_set_name() const
    {
        return current_output_set_name;
    }

    /// Return code of error in IPM
    std::string code_error_IPM() const
    {
        return char_array_to_string(pmm->errorCode, 100);
    }
    /// Returndescription of error in IPM
    std::string description_error_IPM() const
    {
        return char_array_to_string(pmm->errorBuf, 1024);
    }

    //  Calls for direct coupling of a FMT code with GEMS3K

    /// (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
    ///  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
    ///   chemical data for all nodes used in the mass-transport model.
    void GEM_restore_MT(
            long int  &p_NodeHandle,   ///< Node identification handle
            long int  &p_NodeStatusCH, ///< Node status code;  see typedef NODECODECH
            //                                           GEM input output  FMT control
            double &p_TK,      ///< Temperature T, Kelvin                       +       -      -
            double &p_P,      ///< Pressure P,  Pa                              +       -      -
            double &p_Vs,     ///< Volume V of reactive subsystem,  m3         (+)      -      +
            double &p_Ms,     ///< Mass of reactive subsystem, kg               -       -      +
            double *p_bIC,    ///< Bulk mole amounts of IC  [nICb]              +       -      -
            double *p_dul,    ///< Upper restrictions to amounts of DC [nDCb]   +       -      -
            double *p_dll,    ///< Lower restrictions to amounts of DC [nDCb]   +       -      -
            double *p_aPH     ///< Specific surface areas of phases,m2/kg[nPHb] +       -      -
            );

    /// (6s) Passes (copies) the GEMS3K input data from the work instance of DATABR structure with TKinMet phases.
    ///  In addition, copies the AMRs for solution phases, and phase stability indices (to keep for kinetics).
    ///  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
    ///   chemical data for all nodes used in the mass-transport model.
    void GEM_restore_MT(
            long int  &p_NodeHandle,   ///< Node identification handle
            long int  &p_NodeStatusCH, ///< Node status code;  see typedef NODECODECH
            //                                           GEM input output  FMT control
            double &p_TK,      ///< Temperature T, Kelvin                       +       -      -
            double &p_P,      ///< Pressure P,  Pa                              +       -      -
            double &p_Vs,     ///< Volume V of reactive subsystem,  m3         (+)      -      +
            double &p_Ms,     ///< Mass of reactive subsystem, kg               -       -      +
            double *p_bIC,    ///< Bulk mole amounts of IC  [nICb]              +       -      -
            double *p_dul,    ///< Upper restrictions to amounts of DC [nDCb]   +       -      -
            double *p_dll,    ///< Lower restrictions to amounts of DC [nDCb]   +       -      -
            double *p_aPH,    ///< Specific surface areas of phases,m2/kg[nPHb] +       -      -
            double *p_omPH,    ///< Stability indices of phases,log10 scale [nPHb] (+)     +      -
            double *p_amru,      ///< Upper AMR to masses of sol. phases [nPSb]   +       -      -
            double *p_amrl       ///< Lower AMR to masses of sol. phases [nPSb]   +       -      -

            );

    /// (8) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
    ///     This call is usually preceeding the GEM_run() call
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double p_Vs,     ///< Volume V of reactive subsystem, m3               -       -      +
            double p_Ms,     ///< Mass of reactive subsystem, kg                   -       -      +
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll,   ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            double *p_aPH    ///< Specific surface areas of phases, m2/kg [nPHb]   +       -      -
            );

    /// (8a) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
    /// This overloaded variant uses the xDC speciation vector for setting the
    /// new bulk chemical composition to be used in the next GEM_run() calculation.
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double p_Vs,     ///< Volume V of reactive subsystem, m3               -       -      +
            double p_Ms,     ///< Mass of reactive subsystem, kg                   -       -      +
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll,   ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            double *p_asPH,  ///< Specific surface areas of phases, m2/kg [nPHb]   +       -      -
            double *p_xDC    ///< Mole amounts of DCs [nDCb] - will be convoluted
            ///< and added to the bIC GEM input vector (if full speciation
            ///< and not just increments then p_bIC vector must be zeroed off -
            ///< it will be calculated from p_xDC and stoichiometry matrix A
            );

    /// (8b) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
    /// In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
    /// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
    /// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
    /// structure contains a full list of Dependent Components used in GEM IPM2 calculations.
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double p_Vs,     ///< Volume V of reactive subsystem, m3               -       -      +
            double p_Ms,     ///< Mass of reactive subsystem, kg                   -       -      +
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll,   ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            double *p_asPH,  ///< Specific surface areas of phases, m2/kg [nPHb]   +       -      -
            double *p_xDC,   ///< Mole amounts of DCs [nDCb] - old primal soln.     +      -      -
            double *p_gam    ///< DC activity coefficients [nDCb] - old primal s.   +      -      -
            );

    /// (8c) Loads the minimum of GEMS3K input data for a given mass-transport node into the work instance
    ///   of DATABR structure. This call is usually preceeding the GEM_run() call.
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll    ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            );

    /// (8d) Loads the minimum of GEMS3K input data for a given mass-transport node into the work instance of
    ///   the DATABR structure. In addition, loads specific surface areas of phases,
    ///   stability indexes of phases (kept from previous time point), and AMRs for phases-solutions.
    ///   This call is usually preceeding the GEM_run() call.
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll,   ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            double *p_asPH,  ///< Specific surface areas of phases, m2/kg [nPHb]   +       -      -
            double *p_omPH,   ///< Stability indices of phases,log10 scale [nPHb]  (+)      +      -
            double *p_amru,   ///< Upper AMR to masses of sol. phases [nPSb]        +       -      -
            double *p_amrl    ///< Lower AMR to masses of sol. phases [nPSb]        +       -      -
            );

    /// (8e) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure,
    /// with stability indexes of phases (kept from previous time point), and AMRs for phases-solutions.
    /// In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
    /// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
    /// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
    /// structure contains a full list of Dependent Components used in GEM IPM2 calculations.
    void GEM_from_MT(
            long int  p_NodeHandle,   ///< Node identification handle
            long int  p_NodeStatusCH, ///< Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
            //                                              GEM input output  FMT control
            double p_TK,     ///< Temperature T, Kelvin                            +       -      -
            double p_P,      ///< Pressure P, Pa                                   +       -      -
            double p_Vs,     ///< Volume V of reactive subsystem, m3               -       -      +
            double p_Ms,     ///< Mass of reactive subsystem, kg                   -       -      +
            double *p_bIC,   ///< Bulk mole amounts of IC [nICb]                   +       -      -
            double *p_dul,   ///< Upper restrictions to amounts of DC [nDCb]       +       -      -
            double *p_dll,   ///< Lower restrictions to amounts of DC [nDCb]       +       -      -
            double *p_asPH,  ///< Specific surface areas of phases, m2/kg [nPHb]   +       -      -
            double *p_omPH,   ///< Stability indices of phases,log10 scale [nPHb]  (+)      +      -
            double *p_amru,   ///< Upper AMR to masses of sol. phases [nPSb]       +       -      -
            double *p_amrl,    ///< Lower AMR to masses of sol. phases [nPSb]      +       -      -
            double *p_xDC,   ///< Mole amounts of DCs [nDCb] - old primal soln.    +       -      -
            double *p_gam    ///< DC activity coefficients [nDCb] - old primal s.  +       -      -
            );

    /// (9) Optional, for passing the current mass transport time and time step into the work
    /// DATABR structure (for using it in TKinMet, or tracing/debugging, or in writing DBR files for nodes)
    /// This call should be used instead of obsolete GEM_set_MT() (provided below for compatibility with older codes)
    //
    void GEM_from_MT_time(
            //   long int  NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
            //   long int  NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
            double p_Tm,      ///< Actual total simulation time, s               +       -      -
            double p_dt       ///< Actual time step duration, s                  +       -      -
            );

    void GEM_set_MT(  // misleading name of the method - use instead GEM_from_MT_time(), see above
                      //   long int  NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
                      //   long int  NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
                      double p_Tm,      ///< Actual total simulation time, s               +       -      -
                      double p_dt       ///< Actual time step duration, s                  +       -      -
                      );

    // Sets time step (for kinetics) if step >= 0 or gets it if step < 0
    // direct access to MULTI - needs reimplementation for TKinetics in TNode
    long int GEM_step_MT( const long int step )
    {
        if( step >= 0 )
            pmm->ITau = step;
        return pmm->ITau;
    }


    /// (5) Reads another DBR file (with input system composition, T,P etc.) \ . The DBR file must be compatible with
    /// the currently loaded IPM and DCH files (see description  of GEM_init() function call).
    /// \param fname     Null-terminated (C) string containing a full path to the input DBR disk file.
    /// \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2)
    /// \return  0  if successful; 1 if input file(s) has not found been or is corrupt; -1 if internal memory allocation error occurred.
    long int GEM_read_dbr( const char* fname, GEMS3KGenerator::IOModes type_f );

    /// (5j) Reads another DBR object (with input system composition, T,P etc.) from JSON string \ .
    /// It must be compatible with the currently loaded IPM and DCH files
    ///  (see descriptions of GEM_init() methods).
    /// \param dbr_json  String containing a JSON document for the input DBR object.
    /// \param check_dch_compatibility  If true, forces checking the compatibility of DBR object with active DCH/IPM
    /// \return  0  if successful; 1 if input JSON string is empty or corrupt;-1 if internal memory allocation error occurred;
    ///          2 if checking the not compatibility of DBR object with active DCH/IPM.
    long int GEM_read_dbr( std::string dbr_json, const bool check_dch_compatibility = true );

    /// (2) Main call for GEM IPM calculations using the input bulk composition, temperature, pressure
    ///   and metastability constraints provided in the work instance of DATABR structure.
    ///   Actual calculation will be performed only when dBR->NodeStatusCH == NEED_GEM_SIA (5) or dBR->NodeStatusCH = NEED_GEM_AIA (1).
    ///   By other values of NodeStatusCH, no calculation will be performed and the status will remain unchanged.
    ///  In "smart initial approximation" (SIA) mode, the program can automatically switch into the "automatic initial
    ///  approximation" (AIA) mode and return  OK_GEM_AIA instead of OK_GEM_SIA.
    ///  \param uPrimalSol  flag to define the mode of GEM smart initial approximation
    ///                     (only if dBR->NodeStatusCH = NEED_GEM_SIA has been set before GEM_run() call).
    ///                     false  (0) -  use speciation and activity coefficients from previous GEM_run() calculation
    ///                     true  (1)  -  use speciation provided in the DATABR memory structure (e.g. after reading the DBR file)
    ///  \return NodeStatusCH  (the same as set in dBR->NodeStatusCH). Possible values (see "databr.h" file for the full list)
    long int  GEM_run( bool uPrimalSol );   // calls GEM for a work node

    /// Returns GEMIPM2 calculation time in seconds elapsed during the last call of GEM_run() - can be used for monitoring
    ///                      the performance of calculations.
    /// \return double number, may contain 0.0 if the calculation time is less than the internal time resolution of C/C++ function
    double GEM_CalcTime(long int& NumK2, long int& NumIterFIA, long int& NumIterIPM);

    /// To obtain the number of GEM IPM2 iterations performed during the last call of GEM_run() e.g. for monitoring the
    /// performance of GEMS3K in AIA or SIA modes, or for problem diagnostics.
    /// Parameters:  long int variables per reference (must be allocated before calling GEM_Iterations(), previous values will be lost. See Return values.
    /// \return  Function Total number of EFD + IPM iterations from the last call to GEM_run()
    /// \param PrecLoops   Number of performed IPM-2 precision refinement loops
    /// \param NumIterFIA  Total number of performed MBR() iterations to obtain a feasible initial approximation for the IPM algorithm.
    /// \param NumIterIPM  Total number of performed IPM main descent algorithm iterations.
    long int GEM_Iterations( long int& PrecLoops, long int& NumIterFIA, long int& NumIterIPM );

    /// (3) Writes the contents of the work instance of the DATABR structure into a disk file with path name  fname.
    ///   \param fname         null-terminated (C) string containing a full path to the DBR disk file to be written.
    ///                 NULL  - the disk file name path stored in the  dbr_file_name  field of the TNode class instance
    ///                 will be used, extended with ".out".  Usually the dbr_file_name field contains the path to the last input DBR file.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    ///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
    ///                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
    ///                 If   false (0), comments will not be written;
    ///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
    ///                                   with_comments=off as "condensed JSON"
    ///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
    void  GEM_write_dbr( const char* fname, GEMS3KGenerator::IOModes type_f,
                         bool with_comments = true, bool brief_mode = false);

    /// (3j) Writes the contents of the node work instance (the DATABR structure) into provided string in JSON format.
    ///   \param dbr_string  reference to a string where to write JSON string from the current DBR object.
    long int  GEM_write_dbr( std::string& dbr_json );

    /// (4) Produces a formatted text file with detailed contents (scalars and arrays) of the GEM IPM work structure.
    /// This call is useful when GEM_run() returns with a NodeStatusCH value indicating a GEM calculation error
    /// (see  above).  Another use is for a detailed comparison of a test system calculation after the version upgrade of GEMS3K.
    /// \param fname   null-terminated (C) string containing a full path to the disk file to be written.
    ///                NULL  - the disk file name path stored in the  dbr_file_name  field of the TNode class instance will be used,
    ///                extended with ".dump.out".  Usually the dbr_file_name field contains the path to the last input DBR file.
    void  GEM_print_ipm( const char* fname );

    /// The export to ThermoFun JSON format file should include all IComp, DComp and ReacDC records
    /// from the project database, not just the records needed for a particular system
    /// (where some elements, DComps or ReacDCs can be switched off) as done in preparation of DCH lookup arrays.
    ///  \param stream     stream to output json file
    virtual void  write_ThermoFun_format_stream(std::iostream&  stream, bool compact);

    /// (7)  Retrieves the GEMIPM2 chemical speciation calculation results from the work DATABR structure instance
    ///   into memory provided by the mass transport part. Dimensions and order of elements in the arrays must correspond
    ///   to those in currently existing DATACH memory structure.
    void GEM_to_MT(
            long int &p_NodeHandle,    ///< Node identification handle
            long int &p_NodeStatusCH,  ///< Node status code (changed after GEM calculation); see typedef NODECODECH
            long int &p_IterDone,      ///< Number of iterations performed in the last GEM IPM calculation
            //                                                  GEM input output  FMT control
            // Chemical scalar variables
            double &p_Vs,    ///< Total volume V of reactive subsystem at given P,T, m3    -      -      +     +
            double &p_Ms,    ///< Total mass of the reactive subsystem, kg                 -      -      +     +
            double &p_Gs,    ///< Total Gibbs energy of the reactive subsystem, J          -      -      +     +
            double &p_Hs,    ///< Total enthalpy of reactive subsystem, J (reserved)       -      -      +     +
            double &p_IC,    ///< Effective aqueous ionic strength, molal                  -      -      +     +
            double &p_pH,    ///< pH of aqueous solution                                   -      -      +     +
            double &p_pe,    ///< pe of aqueous solution                                   -      -      +     +
            double &p_Eh,    ///< Eh of aqueous solution, V                                -      -      +     +
            // Dynamic data - dimensions see in DATACH.H structure
            double  *p_rMB,  ///< Mole balance residuals for Independent Components [nICb] -      -       +     +
            double  *p_uIC,  ///< Dual solution: IC chemical potentials, mol/mol [nICb]    -      -       +     +
            double  *p_xDC,  ///< Primal solution: DC mole amounts  [nDCb]                 -      -       +     +
            double  *p_gam,  ///< External activity coefficients of DC [nDCb]              -      -       +     +
            double  *p_xPH,  ///< Total mole amounts of all phases [nPHb]                  -      -       +     +
            double  *p_vPS,  ///< Total volumes of multicomponent phases, m3   [nPSb]      -      -       +     +
            double  *p_mPS,  ///< Total mass of multicomponent phase (carrier),kg [nPSb]   -      -       +     +
            double  *p_bPS,  ///< Bulk compositions of multicomponent phases  [nPSb][nICb] -      -       +     +
            double  *p_xPA,  ///< Amount of carrier in a multicomponent asymmetric phase[nPSb]-    -      +     +
            double  *p_aPH,  ///< Calculated surface areas of phases (m2) [nPHb]           -      -       +     +
            double  *p_bSP   ///< Bulk composition of all solids, moles [nICb]             -      -       +     +
            );

    /// (7a) Optional, to check if the time step in the work DATABR structure was o.k. for TKinMet calculations,
    ///  compared with the time step p_dt given before the GEM calculation. Checks the criteria for the validity
    ///  of time step. If time step was acceptable by a TKinMet model used, returns the actual time step after
    ///  copying (changed) AMRs into p_dul and p_dll vectors, as well as (changed) specific surface areas of
    ///  some (kinetically controlled) phases. Otherwise, returns a (smaller) suggested time step, while the
    ///  p_dul, p_pll, and p_asPH vectors remain unchanged.
    ///  Returns 0 or a negative number (unchanged p_dul and p_dll), if TKinMet calculations failed.
    //
    double GEM_to_MT_time(
            double p_dt,       ///< Actual time step, s                                     -       -     (+)   (+)
            double *p_dul,    ///< Upper AMR restrictions to amounts of DC [nDCb]          -       -      +     -
            double *p_dll,    ///< Lower AMR restrictions to amounts of DC [nDCb]          -       -      +     -
            double *p_amru,    ///< Upper AMR to masses of solution phases [nPSb]           -       -      +     -
            double *p_amrl,    ///< Lower AMR to masses of solution phases [nPSb]           -       -      +     -
            double *p_asPH    ///< Specific surface areas of phases m2/kg  [nPHb]          -       -      +     -
            );


    // Access methods for direct or protected manipulation of CSD and DBR data

    DATACH* pCSD() const  /// Get the pointer to chemical system definition data structure
    {     return CSD;   }

    DATABR* pCNode() const  /// Get pointer to work node data structure
    /// usage on the level of TNodearray is not recommended !
    {        return CNode;     }

    TActivity* pActiv() const  /// Get pointer to TActivity class instance
    { return atp.get(); }

    // These methods get contents of fields in the work node structure
    double cTC() const     /// Get current node Temperature T, Celsius
    {  return CNode->TK-C_to_K;   }

    // These methods get contents of fields in the work node structure
    double cTK() const     /// Get current node Temperature T, Kelvin
    {  return CNode->TK;   }

    double cP() const     /// Get current node Pressure P, Pa
    {        return CNode->P;   }

    double cMs() const     /// Get current node mass in kg (reactive part)
    {        return CNode->Ms;   }

    double cVs() const     /// Get current node volume in m3 (reactive part)
    {        return CNode->Vs;   }

    /// Set current node identification handle to value of \param jj
    void setNodeHandle( long int jj )
    {      CNode->NodeHandle = jj;  }

    // Useful methods facilitating the communication between DataCH (or FMT)
    // and DataBR (or node) data structures for components and phases
    // (i.e. between the chemical system definition and the node)

    //DM 25.02.2014 (call from gemsfit)
    void Access_GEM_IMP_init();
    long get_sizeTSolMod();
    /// Return a pointer to a phase (TSolMod) with index xPH
    void *get_ptrTSolMod( int xPH ) const;

    /// Returns DCH index of IC given the IC Name string (null-terminated)
    /// or -1 if no such name was found in the DATACH IC name list
    long int IC_name_to_xCH( const char *Name ) const;

    /// Returns DCH index of DC given the DC Name string
    /// or -1 if no such name was found in the DATACH DC name list
    long int DC_name_to_xCH( const char *Name ) const;

    /// Returns DC Name string given the DCH index of DC, check MaxDCnameLength
    /// or -1 if no such name was found in the DATACH DC name list
    char* xCH_to_DC_name( int xCH ) const
    {return CSD->DCNL[xCH];}

    char* xCH_to_PH_name( int xCH ) const
    {return CSD->PHNL[xCH];}

    /// Returns IC Name string given the ICH index of IC, check MaxICnameLength
    /// or -1 if no such name was found in the DATACH IC name list
    char* xCH_to_IC_name( int xCH ) const
    {return CSD->ICNL[xCH];}

    /// Returns the class codes of phase given the ICH index of PH, check check MaxPHnameLength
    /// or -1 if no such name was found in the DATACH ccPH list
    char xCH_to_ccPH( int xCH ) const
    {return CSD->ccPH[xCH];}

    /// Returns the number of interaction parameters. Max parameter order (cols in IPx),
    /// and number of coefficients per parameter in PMc table [3*FIs]
    long int* Get_LsMod ( ) const
    {return pmm->LsMod;}

    /// Returns the Dimensions of TSolMod <DMc> and <MoiSN> arrays [nPS*3]: In each row (for phase):
    /// number of parameters per component; [1] 0; [2] 0. For multi-site (sublattice) models:
    /// [1] number of sublattices nS; [2] total number of moieties nM acting in sublattice sites
    long int* Get_LsMdc ( ) const
    {return pmm->LsMdc;}

    /// Returns DCH index of Phase given the Phase Name string
    /// or -1 if no such name was found in the DATACH Phase name list
    long int Ph_name_to_xCH( const char *Name ) const;

    /// Returns DBR index of IC given the IC Name string
    /// or -1 if no such name was found in the DATACH IC name list
    inline long int IC_name_to_xDB( const char *Name ) const
    { return IC_xCH_to_xDB( IC_name_to_xCH( Name ) ); }

    /// Returns DBR index of DC given the DC Name string
    /// or -1 if no such name was found in the DATACH DC name list
    inline long int DC_name_to_xDB( const char *Name ) const
    { return DC_xCH_to_xDB( DC_name_to_xCH( Name ) ); }

    /// Returns DBR index of Phase given the Phase Name string
    /// or -1 if no such name was found in the DATACH Phase name list
    inline long int Ph_name_to_xDB( const char *Name ) const
    { return Ph_xCH_to_xDB( Ph_name_to_xCH( Name ) ); }

    /// Converts the IC DCH index into the IC DBR index
    /// or returns -1 if this IC is not used in the data bridge
    long int IC_xCH_to_xDB( const long int xCH ) const;

    /// Converts the DC DCH index into the DC DBR index
    /// or returns -1 if this DC is not used in the data bridge
    long int DC_xCH_to_xDB( const long int xCH ) const;

    /// Converts the Phase DCH index into the Phase DBR index
    /// or returns -1 if this Phase is not used in the data bridge
    long int Ph_xCH_to_xDB( const long int xCH ) const;

    /// Converts the IC DBR index into the IC DCH index
    inline long int IC_xDB_to_xCH( const long int xBR ) const
    { return CSD->xic[xBR]; }

    /// Converts the DC DBR index into the DC DCH index
    inline long int DC_xDB_to_xCH( const long int xBR ) const
    { return CSD->xdc[xBR]; }

    /// Converts the Phase DBR index into the Phase DCH index
    inline long int Ph_xDB_to_xCH( const long int xBR ) const
    { return CSD->xph[xBR]; }

    /// Returns the DCH index of the first DC belonging to the phase with DCH index Phx
    long int Phx_to_DCx( const long int Phx ) const;

    /// Returns the DCH index of the first DC belonging to the phase with DCH index Phx,
    /// plus returns through the nDCinPh (reference) parameter the number of DCs included into this phase
    long int  PhtoDC_DCH( const long int Phx, long int& nDCinPh ) const;

    /// Returns the DCH index of the Phase to which the Dependent Component with index xCH belongs
    long int  DCtoPh_DCH( const long int xCH ) const;

    ///  Returns the DBR index of the first DC belonging to the phase with DBR index Phx,
    /// plus returns through the nDCinPh (reference) parameter the number of DCs included into DBR for this phase
    long int  PhtoDC_DBR( const long int Phx, long int& nDCinPh ) const;

    /// Returns the DBR index of the Phase to which the  Dependent Component with index xBR belongs
    long int  DCtoPh_DBR( const long int xBR ) const;

    // Data exchange methods between GEMIPM and work node DATABR structure
    // Are called inside of GEM_run()
    virtual void packDataBr();   ///<  Packs GEMIPM calculation results into work node structure
    void unpackDataBr( bool uPrimalSol ); ///<  unpacks work DATABR content into GEMIPM data structure

    // Access to interpolated thermodynamic data from DCH structure
    /// Checks if given temperature T (K) and pressure P (Pa) fit within the interpolation
    /// intervals of the DATACH lookup arrays (returns true) or not (returns false)
    bool  check_TP( double T, double P ) const;

    /// Tests TK as a grid point for the interpolation of thermodynamic data.
    /// \return index in the lookup grid array or -1  if it is not a grid point
    long int  check_grid_T( double TK ) const
    {
        return dbr_dch_api::check_grid_T(CSD, TK);
    }

    /// Tests P as a grid point for the interpolation of thermodynamic data.
    /// \return index in the lookup grid array or -1 if it is not a grid point
    long int  check_grid_P( double P ) const
    {
        return dbr_dch_api::check_grid_P(CSD, P);
    }

    /// Tests T (K) and P (Pa) as a grid point for the interpolation of thermodynamic data using DATACH
    /// lookup arrays. \return -1L if interpolation is needed, or 1D index of the lookup array element
    /// if TK and P fit within the respective tolerances.
    long int  check_grid_TP(double T, double P) const
    {
        return dbr_dch_api::check_grid_TP(CSD, T, P);
    }

    /// Returns number of temperature and  pressure grid points for one dependent component
    inline long int gridTP() const
    {
        return dbr_dch_api::gridTP(CSD);
    }

    /// Returns 1 if a Psat value corresponding to the temperature of interest was found in the GEMS3K input file
    double get_Ppa_sat( double Tk );

    /// Returns index of Tk point - Psat point pair
    long int get_grid_index_Ppa_sat( double Tk );

    double Get_Psat (double Tk);

    /// Sets new molar Gibbs energy G0(P,TK) value for Dependent Component for the lookup grid point (P,TK)
    /// in the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \param new_G0 in J/mol;
    /// \return 0
    double Set_DC_G0(const long int xCH, const double P, const double TK, const double new_G0 );

    /// Retrieves (interpolated) molar Gibbs energy G0(P,TK) value for Dependent Component
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \param norm defines in wnich units the value is returned: false - in J/mol; true (default) - in mol/mol
    /// \return G0(P,TK) or 7777777., if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_G0(const long int xCH, const double P, const double TK,  bool norm=true) const;

    /// Retrieves (interpolated, if necessary) molar volume V0(P,TK) value for Dependent Component (in J/Pa)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return V0(P,TK) (in J/Pa) or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_V0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) molar enthalpy H0(P,TK) value for Dependent Component (in J/mol)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return H0(P,TK) (in J/mol) or 7777777., if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_H0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) absolute molar enropy S0(P,TK) value for Dependent Component (in J/K/mol)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return S0(P,TK) (in J/K/mol) or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_S0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) constant-pressure heat capacity Cp0(P,TK) value for Dependent Component (in J/K/mol)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return Cp0(P,TK) (in J/K/mol) or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_Cp0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) Helmholtz energy  of Dependent Component (in J/mol)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return Helmholtz energy (in J/mol) or 7777777., if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_A0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) Internal energy of  Dependent Component (in J/mol)
    /// from the DATACH structure.
    /// \param xCH is the DC DCH index
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return Internal energy (in J/mol) or 7777777., if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DC_U0(const long int xCH, const double P, const double TK) const;

    /// Retrieves (interpolated) density and its derivatives of liquid water at (P,TK) from the DATACH structure or 0.0,
    /// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
    /// \param P refers to the pressure in Pascal
    /// \param TK refers to the temperature in Kelvin
    /// \param DensAW contains the density of water (at P and Tk) and its temperature and pressure derivatives
    void DensArrayH2Ow( const double P, const double TK, std::vector<double>& DensAW );

    /// Retrieves (interpolated) dielectric constant and its derivatives of liquid water at (P,TK) from the DATACH structure or 0.0,
    /// if TK (temperature, Kelvin) or P (pressure, Pa) parameters go beyond the valid lookup array intervals or tolerances.
    /// \param P refers to the pressure in Pascal
    /// \param TK refers to the temperature in Kelvin
    /// \param EpsAW contains the permittivity of water (at P and Tk) and its temperature and pressure derivatives
    void EpsArrayH2Ow( const double P, const double TK, std::vector<double>& EpsAW );


    /// Retrieves (interpolated) dielectric constant of liquid water at (P,TK) from the DATACH structure.
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return dielectric constant or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double EpsH2Ow(const double P, const double TK);

    /// Retrieves (interpolated) density of liquid water (in kg/m3) at (P,TK) from the DATACH structure.
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return density or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DenH2Ow(const double P, const double TK);

    /// Retrieves (interpolated) viscosity of liquid water (in kg/m3) at (P,TK) from the DATACH structure.
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return viscosity or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double VisH2Ow(const double P, const double TK);

    /// Retrieves (interpolated) dielectric constant of H2O vapor at (P,TK) from the DATACH structure.
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return dielectric constant of H2O vapor or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double EpsH2Og(const double P, const double TK);

    /// Retrieves (interpolated) density of H2O vapor (in kg/m3) at (P,TK) from the DATACH structure.
    /// \param P pressure, Pa
    /// \param TK temperature, Kelvin
    /// \return density of H2O vapor or 0.0, if TK or P  go beyond the valid lookup array intervals or tolerances.
    double DenH2Og(const double P, const double TK);

    /// Retrieves the current phase volume in m3 in the reactive sub-system.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase volume in m3 or 0.0, if the phase mole amount is very close to zero.
    double  Ph_Volume( const long int xBR ) const;

    /// Retrieves the current phase total Gibbs energy in J in the reactive sub-system.
    /// Chemical potential of phase (J/mol) is Gibbs energy divided by moles of the phase.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase total gibbs energy in J or 0.0, if the phase mole amount is very close to zero.
    double  Ph_GibbsEnergy( const long int xBR ) const;

    /// Retrieves the current phase enthalpy in J in the reactive sub-system.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase enthalpy in J or 0.0, if the phase mole amount is very close to zero.
    double  Ph_Enthalpy( const long int xBR ) const;

    /// Retrieves the current phase entropy in J/K in the reactive sub-system.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase entropy in J/K or 0.0, if the phase mole amount is very close to zero.
    double  Ph_Entropy( const long int xBR ) const;

    /// Retrieves the current phase heat capacity Cp in J/K in the reactive sub-system.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase heat capacity Cp in J/K or 0.0, if the phase mole amount is very close to zero.
    double  Ph_HeatCapacityCp( const long int xBR ) const;

    /// Retrieves the current phase amount (in equilibrium) in moles in the reactive sub-system.
    /// Works both for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the current phase volume in moles or 0.0, if the phase mole amount is very close to zero.
    double  Ph_Moles( const long int xBR ) const; // Added for consistency on Oct 1, 2020
    double  Ph_Mole( const long int xBR ) const; // obsolete, not recommended to use in new coupled codes

    /// Retrieves the phase mass in kg.
    /// Works for multicomponent and for single-component phases.
    /// \param xph is DBR phase index
    /// \return the phase mass in kg or 0.0, if the phase mole amount is very close to zero.
    double  Ph_Mass( const long int xBR ) const;

    /// Retrieves the phase stability (saturation) index (lgOmega).
    /// Works for multicomponent and for single-component phases.
    /// \param xBR is DBR phase index
    /// \return the phase saturation index in log10 scale (also for phases in zero amount).
    double Ph_SatInd(const long int xBR ) const;

    /// Retrieval of the phase bulk composition into memory indicated by  ARout.
    /// This function works for multicomponent and for single-component phases
    /// \param xph is DBR phase index
    /// \param ARout is array of at least [dCH->nICb elements] or ARout = NULL
    /// \return pointer to ARout which may also be  allocated inside of Ph_BC()
    /// in the case if parameter ARout = NULL is specified;
    /// to avoid a memory leak, you will have to free this memory wherever appropriate.
    double* Ph_BC( const long int xph, double* ARout = nullptr ) const;

    /// Retrieves total dissolved aqueous molality of Independent Component with DBR index xic.
    /// \param xic is IC DBR index
    /// \return total dissolved aqueous molality or 0.0,
    /// if there is no water in the node or no aqueous phase in DATACH.
    double Get_mIC( const long xic ) const;

    /// Retrieves pH of the aqueous solution
    double Get_pH( );

    /// Retrieves pe of the aqueous solution
    double Get_pe( );

    /// Retrieves Eh of the aqueous solution
    double Get_Eh( );

    /// Retrieves IC (effective molal ionic strength of aqueous electrolyte) of the aqueous solution
    double Get_IC( );

    /// Sets the TK in the work DATABR structure.
    /// \param TK is the temperature value
    void Set_TK(const double TK)
    {  CNode->TK = TK;  }

    /// Sets the P in the work DATABR structure.
    /// \param P is the presure value
    void Set_P(const double P)
    {  CNode->P = P;  }

    /// Retrieves the pressure P (Pa) in the current (work) node
    inline double Get_P( ) const
    {  return CNode->P;  }

    /// Retrieves the temperature T_K (Kelvin) in the current (work) node
    inline double Get_TK( ) const
    {  return CNode->TK;  }

    /// Sets the amount of IC  in the bIC input vector of the work DATABR structure.
    /// \param xic is IC DBR index
    /// \param bIC is amount of IC
    void Set_bIC( const long int xic, const double bIC)
    {  CNode->bIC[xic] = bIC;  }

    /// Retrieves the current dual chemical potential of IC in normalized scale.
    /// \param xic is IC DBR index
    inline double Get_uIC(const long int xic) const
    {  return CNode->uIC[xic];  }

    /// Retrieves the current amount of Independent Component.
    /// \param xic is IC DBR index
    inline double Get_bIC(const long int xic) const
    {  return CNode->bIC[xic];  }

    /// Sets the metastability constraint from below to the amount of DC
    /// in the dll vector of the work DATABR structure.
    /// \param xdc is DC DBR index
    inline void Set_dll( const long int xdc, const double dll)
    {  CNode->dll[xdc] = dll;  }

    /// Sets the metastability constraint from above to the amount of DC
    /// in the dul vector of the work DATABR structure.
    /// \param xdc is DC DBR index
    inline void Set_dul( const long int xdc, const double dul)
    {  CNode->dul[xdc] = dul;  }

    /// Sets the amount of DC in the xDC vector of the work DATABR structure.
    /// \param xdc is DC DBR index
    void Set_nDC( const long int xdc, const double nDC)
    {  CNode->xDC[xdc] = nDC;  }

    /// Retrieves the current mole amount of Dependent Component.
    /// \param xdc is DC DBR index
    inline double Get_nDC(const long int xdc) const
    {  return CNode->xDC[xdc];  }

    /// Retrieval of (dual-thermodynamic) chemical potential of the DC.
    /// \param xdc is DC DBR index
    /// \param norm defines the scale: if true (1) then in mol/mol, otherwise in J/mol
    double Get_muDC( const long int xDC, bool norm=true ) const;

    /// Retrieval of (dual-thermodynamic) activity of the DC.
    /// \param xdc is DC DBR index
    /// \param scale if true then activity is returned, if false then log10(activity)
    double Get_aDC( const long int xdc, bool scale=true ) const;

    /// Retrieves concentration of Dependent Component in its phase
    /// in the respective concentration scale. For aqueous species, molality is returned;
    /// for gas species, mole fraction not partial pressure; for surface complexes - molality;
    /// for species in other phases - mole fraction.
    /// \param xdc is DC DBR index
    /// \return 0.0, if DC has zero amount.
    double Get_cDC( const long int xdc ) const;

    /// Retrieves the activity coefficient of Dependent Component
    /// in its phase in the respective scale.
    /// \param xdc is DC DBR index
    /// \return 1.0, if DC has zero amount.
    inline double Get_gDC(const long int xdc) const
    {  return ( noZero( CNode->xDC[xdc] ) ? CNode->gam[xdc]: 1.0);  }

    /// Retrieves the molar mass of Dependent Component in kg/mol.
    /// \param xdc is DC DBR index
    inline double DCmm( const long int xdc ) const
    { return CSD->DCmm[ CSD->xdc[xdc]]; }

    /// Retrieves the molar mass of Independent Component in kg/mol.
    /// \param xic is IC DBR index
    inline double ICmm( const long int xic ) const
    { return CSD->ICmm[ CSD->xic[xic]]; }

    /// Retrieves the stoichiometry coefficient a[xdc][xic] of IC in the formula of DC.
    /// \param xdc is DC DBR index
    /// \param xic is IC DBR index
    inline double DCaJI( const long int xdc, const long int xic) const
    { return CSD->A[ CSD->xic[xic] + CSD->xdc[xdc] * CSD->nIC ]; }

    // These methods can only be used for the current work node (direct access to GEM IPM data)

    /// Sets the total amount of Independent Component.
    /// Also amount of ICs not included into DATABR list can be retrieved.
    /// Internal re-scaling to mass of the system is applied.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is IC DCH index
    inline void Set_IC_b( const double b_val, const long int xCH)
    { pmm->B[xCH] = b_val; }

    // used in GEMSFIT
    /// Sets the mLook Mode of lookup-interpolation: 0 interpolation (on nTp*nPp grid).
    /// \param mLook is 0 or 1
    void Set_mLook(const double mLook)
    {  CSD->mLook = static_cast<long>(mLook);  load_thermodynamic_data = false;}

    /// Sets the value of the interaction parameter.
    /// Internal re-scaling to mass of the system is applied.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xPMC is the index of the interaction parameter
    inline void Set_PMc( const double PMc_val, const long int xPMc)
    { pmm->PMc[xPMc] = PMc_val; load_thermodynamic_data = false; }

    /// Gets the value of the interaction parameter.
    inline void Get_PMc( double &PMc_val, const long int xPMc)
    {  PMc_val = pmm->PMc[xPMc]; load_thermodynamic_data = false; }

    /// Gets code of the aquesous solution model
    void Get_sMod( int ndx, std::string &sMod);

    /// Sets the value of the phase component parameter.
    /// Internal re-scaling to mass of the system is applied.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xDMC is the index of the interaction parameter
    inline void Set_DMc( const double DMc_val, const long int xDMc)
    { pmm->DMc[xDMc] = DMc_val; load_thermodynamic_data = false; }

    /// Gets the value of the interaction parameter.
    inline void Get_DMc( double &DMc_val, const long int xDMc)
    {  DMc_val = pmm->DMc[xDMc];  }



    /// Retrieves the current total amount of Independent Component.
    /// Also amount of ICs not included into DATABR list can be retrieved.
    /// Internal re-scaling to mass of the system is applied
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is IC DCH index
    inline double IC_b(const long int xCH) const
    { return pmm->B[xCH]; }

    /// Retrieves the current mole amount of DC directly from GEM IPM work structure.
    /// Also amount of DCs not included into DATABR list can be retrieved.
    /// Internal re-scaling to mass of the system is applied.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    inline double DC_n(const long int xCH) const
    {  return pmm->X[xCH]; }

    /// Retrieves the current (dual-thermodynamic) activity of DC
    /// directly from GEM IPM work structure. Also activity of a DC not included into DATABR list
    /// can be retrieved. If DC has zero amount, its dual-thermodynamic activity is returned anyway.
    /// For single condensed phase component, this value has a meaning of the saturation index,
    /// also in the presence of metastability constraint(s).
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    double DC_a(const long int xCH);

    // GEMSFIT access functions
    /// Functions for accessing parameters of mixing and properties of phase components used in TSolMod class

    /// Retrieves indices of origin in TSolMod composite arrays for a solution phase of interest index_phase.
    /// \param IN: index_phase is the DCH index of phase of interest.
    /// \param OUT: ipaIPx, ipaIPc, ipaDCc are origin indices of this phase in aIPx, aIPc and aDCc arrays, respectively.
    void Get_IPc_IPx_DCc_indices( long int &ipaIPx, long int &ipaIPc, long int &ipaDCc, const long int &index_phase );

    /// Retrieves dimensions of TSolMod array for a solution phase of interest index_phase.
    /// \param IN: index_phase is the DCH index of phase of interest.
    /// \param OUT: NPar, NPcoef, MaxOrd, NComp, NP_DC, are number of interaction parameters, number of coefficients per parameter,
    /// \param   maximum parameter order (i.e. row length in aIPx), number of components in the phase, and number of coefficients
    /// \param   per component, respectively.
    void Get_NPar_NPcoef_MaxOrd_NComp_NP_DC( long int &NPar, long int &NPcoef, long int &MaxOrd,
                                             long int &NComp, long int &NP_DC, const long int &index_phase );

    /// Gets values of the aIPc array (of interaction parameter coefficients) for the solution phase of interest index_phase.
    /// \param IN: ipaIPc is the origin index (of the first element) of the aIPc array; index_phase is the DCH index of phase of interest.
    /// \param OUT: returns vaIPc - vector with the contents of the aIPc sub-array.
    void Get_aIPc ( std::vector<double> &vaIPc, const long int &ipaIPc, const long int &index_phase );

    /// Gets values of the aIPx list array (of indexes of interacting moieties or components) for the solution phase of interest index_phase.
    /// \param IN: ipaIPx is the origin index (of the first element) of the aIPx array; index_phase is the DCH index of phase of interest.
    /// \param OUT: returns vaIPx - vector with the contents of the aIPx sub-array.
    void Get_aIPx ( std::vector<long int> &vaIPx,   const long int &ipaIPx, const long &index_phase );

    /// Gets values of the aDCc array (of components property coefficients) for the solution phase of interest index_phase.
    /// \param IN: ipaDCc is the origin index (of the first element) of the aDCc array; index_phase is the DCH index of phase of interest.
    /// \param OUT: returns vaDCc - vector with the contents of the aDCc sub-array.
    void Get_aDCc ( std::vector<double> &vaDCc, const long &ipaDCc, const long &index_phase );

    /// Sets values of the aIPc array (of interaction parameter coefficients) for the solution phase of interest index_phase.
    /// \param IN: vaIPc - vector with the contents of the aIPc sub-array to be set; ipaIPc is the origin index (of the first element)
    /// \param     of the aIPc array; index_phase is the DCH index of phase of interest.
    void Set_aIPc ( const std::vector<double> vaIPc, const long int &ipaIPc, const long &index_phase );

    /// Sets values of the aDCc array (of components property coefficients) for the solution phase of interest index_phase.
    /// \param IN: vaDCc - vector with the contents of the aDCc sub-array to be set. ipaDCc is the origin index (of the first element)
    /// \param of the aDCc array; index_phase is the DCH index of phase of interest.
    void Set_aDCc ( const std::vector<double> vaDCc, const long &ipaDCc, const long &index_phase );

    /// These methods set contents of fields in the work node structure
    /// Direct access to set temperature T_K in the current (work) node
    void Set_Tk   ( const double &T_k );

    /// Direct access to set pressure (P_b given in bar) in the current (work) node
    void Set_Pb   ( const double &P_b );
    // End GEMSFIT access functions

    /// Retrieves the current concentration of Dependent Component in its
    /// phase directly from the GEM IPM work structure. Also the activity of a DC not included into
    /// DATABR list can be retrieved. For aqueous species, molality is returned; for gas species,
    /// partial pressure; for surface complexes - density in mol/m2; for species in other phases -
    /// mole fraction. If DC has zero amount, the function returns 0.0.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    double DC_c(const long int xCH);

    /// Retrieves the current activity coefficient of DC in its plase
    /// directly from GEM IPM work structure. Also activity coefficient of a DC not included
    /// into DATABR list can be retrieved. If DC has zero amount, this function returns 1.0.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    inline double DC_g(const long int xCH) const
    {  return pmm->Gamma[xCH];  }

    /// Retrieves the natural logarithm of the internal activity coefficient of species at DCH index xCH
    /// \param xCH index of species DCH
    inline double DC_lng( const long int xCH ) const
    {  return pmm->lnGam[ xCH ]; }

    /// Retrieves the current (dual-thermodynamic) chemical potential of DC
    /// directly from GEM IPM work structure. Also for any DC not included into DATABR or having zero amount.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    /// \param norm defines in wnich units the chemical potential value is returned:
    ///             false - in J/mol; true (default) - in mol/mol
    double DC_mu(const long int xCH, bool norm=true);

    /// Retrieves the standard chemical potential of DC directly
    /// from GEM IPM work structure at current pressure and temperature.
    /// Also for any DC not included into DATABR or having zero amount.
    /// These methods can only be used for the current work node (direct access to GEM IPM data)
    /// \param xCH is DC DCH index
    /// \param norm defines in which units the chemical potential value is returned:
    ///      false - in J/mol; true (default) - in mol/mol
    double DC_mu0(const long int xCH, bool norm=true);

    /// Writes work node (DATABR structure) to a text VTK file
    virtual void databr_to_vtk( std::fstream& ff, const char*name, double time, long int  cycle,
                                long int nFilds, long int (*Flds)[2]);

    /// Get the last error message logged to the ipmlog file
    const std::string& ipmLogError() const {
        return ipmlog_error;
    }

    /// Clear error message
    void clearipmLogError() {
        ipmlog_error.clear();
    }

    /// Generate thermodynamic data from ThermoEngine
    bool load_all_thermodynamic_from_thermo( double TK, double PPa );

    /// Clear thermodynamic data from ThermoEngine
    void clear_ThermoEngine();

    /// Read ThermoEngine
    bool load_ThermoEngine(const std::string& thermo_file_or_string);


    /// Deletes fields of DATABR structure indicated by data_BR_
    /// and sets the pointer data_BR_ to NULL
    DATABR* databr_free(DATABR* data_BR_);

protected:

    void allocMemory();
    void freeMemory();

    void databr_realloc()
    {
        dbr_dch_api::databr_realloc(CSD, CNode);
    }

    ///  Reads the contents of the work instance of the DATABR structure from a disk file with path name dbr_file.
    ///   \param dbr_file  string containing a full path to the DBR disk file to be read.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    void  read_dbr_format_file( const std::string& dbr_file, GEMS3KGenerator::IOModes type_f );

    /// Writes the contents of the work instance of the DATABR structure into a disk file with path name  dbr_file.
    ///   \param dbr_file  string containing a full path to the DBR disk file to be read.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    ///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
    ///                           in the DBR file. If set to true (1), the comments will be written for all data entries (default).
    ///                           If   false (0), comments will not be written.
    ///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
    ///                                   with_comments=off as "condensed JSON"
    ///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
    void  write_dbr_format_file( const std::string& dbr_file, GEMS3KGenerator::IOModes type_f, bool with_comments, bool brief_mode );


    // Methods to perform output to vtk files

    /// Prints data of CNode data object element with handle nfild and index ndx into a VTK file
    /// referenced by ff
    void databr_element_to_vtk( std::fstream& ff, DATABR *CNode_, long int nfild, long int ndx );

    /// Prints name of CNode data object element with handle nfild and index ndx into a VTK file
    /// referenced by ff
    void databr_name_to_vtk( std::fstream& ff, long int nfild, long int ndx, long int ndx2=0 );

    /// Prints size of data with handle nfild with dimensions nel and nel2 into a VTK file
    /// referenced by ff
    void databr_size_to_vtk(  long int nfild, long int& nel, long int& nel2 );

    /// Prints header of VTK data with name, time, loop # and xyz coordinates set to 1
    /// into a VTK file referenced by ff
    void databr_head_to_vtk( std::fstream& ff, const char*name, double time, long cycle,
                             long nx = 1, long ny = 1, long nz = 1 );

    // Copying CSD and CNode data structures for TNodeArray parallelization
    /// Copy CSD (DATACH structure) data from other structure.
    void datach_copy( DATACH* otherCSD );
    /// Reads node (work DATABR structure) data from other DBR.
    void databr_copy( DATABR* otherCNode );

    /// Test to reload thermodynamic data from grid
    void CheckMtparam();
    virtual void init_into_gems3k();

public:

    // Added by AL and DK in 2014-2018 as an alternative (more generic for the chemical system) Activity API

    /// Initializes and fills out a new TActivity class instance
    /// from already existing CSD, work GEMS3K and CNode data structures
    void InitCopyActivities( DATACH* CSD, MULTI* pmm, DATABR* CNode );

    // Generic access methods to contents of the new TActivity class

    // Handling of standard thermodynamic data for Dependent Components included in the system

    /// Sets temperature T (in units of K)
    void setTemperature(const double T);

    // setS pressure p (in units of Pa)
    void setPressure(const double P);

    /// Computes standard Gibbs energies g0 for currently set T,P for all DComps in J/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateStandardGibbsEnergies();

    /// Computes standard molar volumes V0 for currently set T,P for all DComps in m3/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateStandardVolumes();

    /// Computes standard enthalpies h0 for currently set T,P for all DComps in J/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateStandardEnthalpies();

    /// Computes standard absolute entropies S0 for currently set T,P for all DComps in J/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateStandardEntropies();

    /// Computes standard heat capacities Cp0 for currently set T,P for all DComps in J/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateStandardHeatCapacities();

    /// Computes all above standard properties for currently set T,P for all DComps in J/mol
    /// (so far only P,T interpolation) from the lookup arrays in CSD
    void updateThermoData( bool norm = true );

    // Handling speciation, concentrations, activities, and other properties of (mixed) phases
    //   and their components

    /// Set the speciation vector or parimal solution (called X or Y in GEMS3K (in units of moles)
    ///   by copying from a provided vector n
    void setSpeciation(const double* n);

    /// (Re)Computes concentrations of DComps in all phases of the currently loaded system
    /// (usually on a current iteration of GEM)
    void updateConcentrations();

    /// initialize models of mixing in multicomponent phases before the GEM run
    void initActivityCoefficients();

    /// (Re)Computes activity coefficients for all DComps in all phases
    ///  (usually on a current iteration of GEM)
    void updateActivityCoefficients();

    /// (Re)Computes integral properties of all (multicpmponent?) phases
    ///   usually after GEM run is finished (TBD)
    void getIntegralPhaseProperties();

    /// (Re)Computes primal chemical potentials of all DComps in all phases
    ///  (usually on a current iteration of GEM)
    void updateChemicalPotentials();

    ///  (Re)Computes and returns total Gibbs energy of the system
    ///  (usually on a current iteration of GEM)
    double updateTotalGibbsEnergy();

    // compute primal activities
    void updateActivities();

    /// Updates all above chemical properties of DComps and phases
    /// after changed speciation and/or T,P or upon GEM iteration ?
    void updateChemicalData();

    // Interface to kinetics and metastability controls (using TKinetics class)  TBD

    // void setSpeciesUpperAMRs( const double* nu );

    // void setSpeciesLowerAMRs( const double* nl );

    // void setPhasesUpperAMRs( const double* nfu );

    // void setPhasesLowerAMRs( const double* nfl );

    // long int updateKineticsMetastability( long int LinkMode );


};

// Redo into a function with interpolation
// Diffusion coefficient of dependent component with node DBr index ICx
// #define nodeCH_DD( DCx )    ( TNode::na->pCSD()->DD[
//                              TNode::na->pCSD()->xDC[(DCx)]] )


#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// _node_h_

