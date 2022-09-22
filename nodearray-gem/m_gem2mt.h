//-------------------------------------------------------------------
// $Id: m_gem2mt.h 1373 2009-07-22 12:25:22Z gems $
// To be finalized in Version 3.10 (2021)
// Declaration of GEM2MT class, config and calculation functions
//
// Copyright (C) 2005,2011,2021 D.Kulik, S. Dmytriyeva
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Selektor GUI GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the GPL v.3 license
//
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _m_gem2mt_h_
#define _m_gem2mt_h_

#include "particlearray.h"

namespace  io_formats {
class TRWArrays;
}

// internal
enum grr_constants { // std::string len for graph
    MAXAXISNAME = 9,
    MAXGSNAME = 70,
    MAXGRNAME = 16
};

const unsigned int  MAXFORMULA =     81,
                    V_SD_RKLEN = 32,
                    V_SD_VALEN = 24;
const int  MAXIDNAME = 12;
const   int MAXFORMUNITDT=     40;

enum pe_valind { /* index control */
    START_, STOP_, STEP_
};

typedef enum {

    GS_INDEF   = '0',
    GS_GOING    = '1',
    GS_DONE    = '2',
    GS_ERR     = '3',

    AS_INDEF   = '0',
    AS_READY   = '1',
    AS_RUN     = '2',
    AS_DONE    = '3'


} GS_AS_CLASSES;

const long int MT_RKLEN = 80,
               SIZE_HYDP = 7;

typedef struct
{ // Description of GEM2MT data structure (prototype of 04 Dec. 2007)
  // Record key format: the same as Proces
  char
   PunE,          // Units of energy   { j;  J c C N reserved }
   PunV,          // Units of volume   { j;  c L a reserved }
   PunP,          //  Units of pressure  { b;  B p P A reserved }
   PunT,          // Units of temperature  { C; K F reserved }

// Allocation and setup flags (14)
   PvICi,    // Use IC quantities for initial system compositions? { + * - }
   PvAUi,    // Use formula units for initial sub-system compositions? { + * - }
   PvMSt,    // Use math script for start setup (+ -)?
   PvMSg,    // Use math script for graphic presentation (+ -)?
   PvEF,     // Use empirical data for graphics  (+ -)?
   PvPGD,    // Use phase groups definitions (+ -)?
   PvFDL,    // Use flux definition list (+ -)
   PvSFL,    // Use source fluxes and elemental stoichiometries for them? (+ -)
   PvGrid,   // Use array of grid point locations? (+ -)
   PvDDc,    //  Use diffusion coefficients for DC - DDc vector (+ -)
   PvDIc,    //  Use diffusion coefficients for IC - DIc vector (+ -)
   PvDCH,    //  Select ICs, DCs and phases to be exchanged via DATABR file (take all, if unchecked) (+ -)?
   PvnVTK,   //  Use selected fields to VTK format (+ -)?
   PvMSc,    // Use math script for control on time steps (+ -)?

     // Controls on operation (14)
   PsMode,  // Code of GEM2MT mode of operation { S F B A W D C T }
//  Status and control flags (+ -)
   gStat,  // GEM2MT generation status: 0 -indefinite; 1 on-going generation run;
        //  2 - done generation run; 3 - generation run error (+ 5: the same in stepwise mode)
   iStat,  // GEM2MT iteration status:  0 - indefinite; 1 ready for analysis;
        //  2 - analysis run; 3 - analysis done; ( + 5: the same using stepwise mode)
   PsSYd,  // Save generated SysEq records to data base (+ -)
   PsSIA,    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
   PsSdat, //  Save DataCH and inital DataBR files as text json files (j|f), as key value files (t|+|o) or binary (b|-)
   PsSdef, //  Do not write data items that contain only default values (+ -)
   PsScom, //  Write files with comments for all data entries ( text mode ) or as "pretty JSON"  ( json mode ) (+ -)
   PsTPai,  //  Create T,P values in Tval, Pval using iterator(+) or enter(-)
   PsTPpath,  //  Disable  Interpolation ( pre-defined TP paths ) (+ -)
   PsMO,     // Use non stop debug output from nodes (+ -)?
   PsVTK,    // Use non stop output from nodes to VTK format files (+ -)?
   PsMPh,   //  Type of transport in box models ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
   PsSmode, //  Using stepwise mode(+ -)?

   name[MAXFORMULA],  //  Name of GEM2MT task
   notes[MAXFORMULA];  //  Comments


  long int  // input L (32+6)
   nC,   // nQ - number of local equilibrium compartments (nodes, boxes, volumes)
   // xC, yC, zC  numbers of nodes along x, y, z coordinates
   nIV,  // number of initial variants of the chemical system, nIV <= nC
   nPG,  // number of mobile phase groups MGP (0 or >= 1)
   nFD,  // total number of MGP flux definitions, incl elemental source fluxes (0 or >=1 )
   nSFD,   // number of elemental source flux definitions (0 or >= 1 )
   nEl,   // number of electrolytes for setting up electrolyte diffusion coefficients in mDEl vector
   nPTypes,     // res Number of allocated particle types (< 20 ? )
   nProps,      // res Number of particle statistic properties (for monitoring) >= anPTypes
   Lbi,  // Lb - number of formula units to set compositions in initial variants
   Nsd,  // N of references to data sources
   Nqpt, // Number of elements in the script work array qpi for transport
   Nqpg, // Number of elements in the script work array qpc for graphics

   Nb,   // N - number of independent components (set automatically from RMults)
   FIb,  // N - number of phases (set automatically from RMults)
   Lb,   // N - number of dependent components in multycomponent phases (set automatically from RMults)
   bTau, // Time point for the simulation break (Tau[0] at start)   // not used
   ntM,  // Maximum allowed number of time iteration steps (default 1000)
   nYS,  // number of plots (columns in the yt array)
   nE,   // Total number of experiments points (in xEt, yEt) to plot over
   nYE,  // number of experimental parameters (columns in the yEt array)
   nPai,  // Number of P points in MTP interpolation array in DataCH ( 1 to 10 )
   nTai,  // Number of T points in MTP interpolation array in DataCH ( 1 to 20 )
   Lsf,       // number of DCs in phases-solutions in Multi (DATACH) for setting box-fluxes
  // These dimensionalities define sizes of dynamic data in DATABR structure!!!
  // Needed to reduce on storage demand for data bridge instances (nodes)!
  // Connection occurs through xIC, xPH and xDC lists!
    nICb,       // number of stoichiometry units (<= nIC) used in the data bridge
    nDCb,      	// number of DC (chemical species, <= nDC) used in the data bridge
    nPHb,     	// number of phases (<= nPH) used in the data bridge
    Nf,       // nICb number of ICs in  (DATABR) for setting box-fluxes
    FIf,      // nPHb number of phases in (DATABR) for setting box-fluxes
    nVTKfld, //  Number of selected fields to VTK format
nRes1,  //   reserved
nRes2,  //   reserved
nRes3,  //   reserved

// iterators for generating syseq record keys for initial system variants
   tmi[3],   // SYSTEM CSD definition #: start, end, step (initial)
   NVi[3]; // Restrictions variant #: start, end, step
  short  axisType[6];  // axis graph type, background(3), graph type, reserved

  double
  // Input for compositions of initial systems
   Msysb, // Masses (kg) and volumes (L) for initial systems: Ms (total mass, normalize)
   Vsysb, // Vs (total volume of the object, for volume concentrations)
   Mwatb, // M(H2O) (mass of water-solvent for molalities)
   Maqb,  // Maq (mass of aqueous solution for ppm etc.)
   Vaqb,  // Vaq (volume of aqueous solution for molarities)
   Pgb,   // Pg (pressure in gas, for partial pressures)
   Tmolb, // MOL total mole amount for basis sub-system composition calculations
   WmCb,  // mole fraction of the carrier DC (e.g. sorbent or solvent)
   Asur,  // Specific surface area of the sorbent (for adsorbed species)
// "ADpar" data object (11 doubles)
   tf,     // time step reduction factor (usually 1)
   vol_in, // initial total node volume (m3),  was column length (m)  -> to sizeLc[0]
   fVel,   // Advection/diffusion mass transport: initial fluid advection velocity (m/sec)
   eps_in, // initial node effective porosity (0 < eps < 1), usually 1
   Km_in,  // initial effective permeability, m2, usually 1
   al_in,  // initial specific longitudinal dispersivity (m), usually 1e-3
   Dif_in, // initial general diffusivity (m2/sec), usually 1e-9
   nto_in, // initial tortuosity factor, usually 1
   cdv,     // cutoff for IC amount differences in the node between time steps (mol, 1e-9)
   cez,     // cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12)
   ADrs1   // reserved
  ;

  double
     // Iterators for MTP interpolation array in DataCH
     Pai[4],  // Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
     Tai[4],  // Temperature T, C: start, end, increment for MTP array in DataCH , Ttol
     Tau[3],  // Physical time iterator
     sizeLc[3];   // spatial dimensions of the medium, also define topology of nodes
      // graphics
   float
      size[2][4]; // Graph axis scale for the region and the fragment


  // Inital systems
  long int
  // These lists of indices connect the DATABR arrays with this structure
     *xIC,   // ICNL indices in DATABR IC vectors [nICb]
     *xDC,   // DCNL indices in DATABR DC list [nDCb]
     *xPH,   // PHNL indices in DATABR phase vectors [nPHb]
     *NPmean,       // array of initial mean particle type numbers per node ( size: nPTypes )
     *nPmin,        // minimum average total number of particles of each type per one node nPTypes
     *nPmax; // maximum average total number of particles of each type per one node nPTypes
  long int  (*xVTKfld)[2];     //  list of selected fields to VTK format

  double (*PTVm)[5];// Pressure P, bar for initial systems (values within Pai range) [nV]
                    // Temperature T, C for initial systems (values within Pai range)
                    // Volume of the system (L), usually zeros changed to m3 on 20.11.2021 by DK
                    // Initial mass (kg) (to normalize composition of initial system in nodes)
                    // Actual mass of initial system (kg) - to set initial node masses
   double
     *Bn,    //  [nIV][Nb] Table of bulk compositions of initial systems
     *CIb,   // [nIV][Nb] Table of quantity/concentration of IC in initial systems
     *CAb,    // [nIV][Lbi] Table of quantity/concentration of formulae for initial systems

     *Tval,   // discrete values of T [nTai] in grid arrays in DataCH
     *Pval;   // discrete values of P [nPai]

  // Init of nodes
    long int (*DiCp)[2];     // array of indexes of initial system variants for distributing to nodes [nC]
    long int (*ParTD)[6];  // array of particle type definitions at t0 or after interruption nPTypes
    double  (*HydP)[SIZE_HYDP]; // [nC][6] hydraulic parameters for nodes in mass transport model
                                 //  value order to be described
    double (*StaP)[4]; // [nC][4] Pressure P, bar for nodes (values within Pai range)
                   // Temperature T, C for nodes (values within Pai range)
                   // Volume of the nodes in m3 (usually zeros)
                   // Initial mass of the node at time step 0 (in kg)
    double (*grid)[3];  // Array of grid point locations, size is nC

  // Fluxes
    long int (*FDLi)[2]; //[nFD][2] Indexes of nodes where this flux begins and ends
              // negative value means one-side flux (source or sink)
             // for source fluxes, -2 means "source flux stoichiometry with index 1
             // line in the BSF table", and so on
double  (*FDLf)[4]; // [nFD][4] Part of the flux defnition list (flux order, flux rate, MGP quantities)
    double
     *PGT,  // Quantities of phases in MGP [FIf][nPG]
     *BSF,   // [nSFD][Nf] table of bulk compositions of elemental fluxes
              //  More to be added here for seq reactors?
     *MB,    // [nC][Nf]  Table of current IC masses in the boxes (in kg)
     *dMB    // [nC][Nf]  Table of current derivatives dM/dTau for ICs in boxes (in kg/s)
     ;

 // graphics
    double
     *xEt,    // Abscissa for experimental points [nXE]
     *yEt,    // Ordinates for experimental points to plot [nXE, nYE]
     *qpi,   //  [Nqpi] Work array for initial systems math script
     *qpc,   //  [Nqpc] Work array for mass transport math script,
     *xt,    //  Abscissa for sampled data [nS]
     *yt;    //  Ordinates for sampled data [nS][nYS]


 //optional
  double
   *DDc,  //  [Lsf] diffusion coefficients for DC
   *DIc,  //  [Nf] diffusion coefficients for IC
   *DEl   //  [nE] diffusion coefficients for electrolyte salts
    ;
// reserved before, used since 17.Nov.2021 for GEM2MT B or F sim output
  double
   *BM,  // arr1 [nC] Current masses of the nodes (boxes) (in kg)
   *BdM, // [nC] Current time derivatives of mass of the nodes (boxes) (in kg/s)
   *BmgpM, // [nC][nPG] Current masses of MGPs in nodes/boxes (in kg), only outgoing
   *FmgpJ, // [nFD] Current MGP fluxes (in kg/s), only outgoing relative to boxes
   *arr5;  // [nSFD] Current source fluxes (defined in the BSF table), kg/s
  char
   *tExpr,  // Math script text for calculation of mass transport
   *gExpr,  // Math script text for calculation of data sampling and plotting
//
   *CIclb, // [Nb] Units of IC quantity/concentration for initial systems compositions
   *AUcln, // [Lbi] Units of setting UDF quantities for initial system compositions
   *UMGP;  // [FIf] units for setting phase quantities in MGP (see PGT )

  char (*sdref)[V_SD_RKLEN]; // "List of bibl. refs to data sources" [0:Nsd-1]
  char (*sdval)[V_SD_VALEN];  // "Parameters taken from the respective data sources"[0:Nsd-1]
  char (*nam_i)[MAXIDNAME]; // [nIV][12] id names of initial systems
  char (*for_i)[MAXFORMUNITDT]; // [Lbi][40] formulae for setting initial system compositions
  char (*for_e)[MAXFORMUNITDT]; // [nE][40] formulae for diffusing dissolved electrolytes
  char (*stld)[EQ_RKLEN]; // List of SysEq record keys for initial systems [nIV]
  char (*FDLid)[MAXSYMB]; // [nFD] ID of fluxes
  char (*FDLop)[MAXSYMB]; // [nFD] Operation codes (letters): flux order,  type codes
  char (*FDLmp)[MAXSYMB]; // [nFD] ID of MGP to move in this flux (if starts with letter)
                             // Otherwise BSF row index of elemental flux (if 0,1,2,...) 
  char (*MGPid)[MAXSYMB]; // [nPG] ID list of mobile phase groups

  char (*SBM)[MAXICNAME+MAXSYMB];  // Keys (names) of IC

  //  graphics
  char  xNames[MAXAXISNAME],        // Abscissa name
        yNames[MAXAXISNAME];      // Ordinate name
  char  (*lNam)[MAXGRNAME];        // List of ID of lines on Graph [nYS]
  char  (*lNamE)[MAXGRNAME];       // List of ID of lines of empirical data [nYE]


/* Work arrays */
 double
   *An,  // [Lbi][N] stoich matrix for formula units from the for_i list
   *Ae   // [nE][N] stoich matrix for diffusing electrolytes (not used for now)
   ;
 double
   *gfc,  // [nC][nPG][Nf] Array of element partition coefficients between given MGP quantity and its source box
   *yfb,  // [nC][nPG][Nf] Array of MGP element bulk compositions in boxes at current time point (moles)
   (*tt)[9];
 char sykey[EQ_RKLEN+10],   // Key of currently processed SysEq record
   *etext,              // internal
   *tprn;              // internal

 //work data indices/counters
 long int
   ctm,    // current CSD #
   cnv,    //  current restriction variant #
   qc,     // current index of the node/box/compartment ( 0 to nC-1 )
   kv,     // current index of the initial system variant (0 to nIV-1 )
   jqc,    // script c-style index (= qc-1) for transport
   jqs,    // script c-style index (= qc-1) for graphics
   jt,     // current index of sampled point (for sampling scripts)
   jdd,   // current index of diffusing DC
   jdi,   // current index of diffusing IC
   ide,   // current index of diffusing electrolyte
   ct,    // actual time iterator
   qf     // current index of flux (in the flux table)
   ;

 double
   cT, // current value of T
   cP, // current value of P
   cV, // current value of V
   cTau, // current physical time
   dTau, // current time step value
   oTau, // old time step value
   dx,   // node distance [L/nC]
   TimeGEM, // pure GEM runtime, in seconds
   ref3,
   ref4
  ;

   char timep[16], TCp[16], Pp[16], NVp[16], Bnamep[16];
}
GEM2MT;

// Pointers to DataCH and DataBR (current) fields in memory ?
// They are located in MULTI

// Current GEM2MT
class TGEM2MT
{
  GEM2MT mt[1];

  std::shared_ptr<TNodeArray> na = nullptr;       // pointer to nodearray class instance
  TParticleArray* pa_mt = nullptr;       // pointer to TParticleArray class instance

    std::string pathVTK;
    std::string nameVTK;
    std::string prefixVTK;

    void logProfilePhMol( FILE* logfile, int inode )
    {
        if( pa_mt )
            pa_mt->logProfilePhMol( logfile, inode );
    }

protected:

    void AllocNa();

    void keyTest( const char *key );
    void Expr_analyze( int obj_num );
    void CalcPoint( int nPoint );
    bool test_sizes();
    void SelectNodeStructures( bool select_all );
    void init_arrays( bool mode );
    void calc_eqstat( bool startSys );
    void outMulti();
    void mt_next();
    void mt_reset();
    void gen_task( bool startSys );
    void make_A( long int siz_, char (*for_)[MAXFORMUNITDT] );
    void Bn_Calc();
    void gen_TPval();
    void CalcStartScript();
    void CalcControlScript();

    void  copyNodeArrays();
    void  NewNodeArray();
    void  putHydP( DATABRPTR* C0 );
    void  LinkNode0(  long int nNode );
    void  LinkNode1(  long int nNode );
    void  LinkCSD(  long int nNode );
    void  allocNodeWork();
    void  freeNodeWork();

    void  CalcGraph();
    long int CheckPIAinNodes1D( char mode,
              long int start_node = 0, long int end_node = 1000 );

    bool  CalcIPM( char mode, long int start_node = 0,
         long int end_node = 1000, FILE* diffile = nullptr );

    void  MassTransAdvecStart();
    void  MassTransAdvecStep( bool ComponentMode = true );
    void  MassTransCraNicStart();
    void  MassTransCraNicStep( bool ComponentMode = true );
    void  MassTransParticleStart();
    void  MassTransParticleStep( bool ComponentMode = true );
    bool Trans1D( char mode  );  // return true if canceled

   // for box flux ODE integration sims with variable time step (TBD)
    long int MaxIter,  // max number of iterations
         nfcn,      // number of functional estimates
         nstep,     // number of steps
         naccept,   // number of permissible steps
         nrejct;    // number of unpermissible steps
    double *x = nullptr;
    double *dx = nullptr;
    double *tv = nullptr;

    // Flow-through box-flux transport simulations
    void  BoxFluxTransportStart();
    void  FlowThroughBoxFluxStep();

    double BoxMasses( long int q );
    void ComposMGPinBox( long int q );
//    void DisCoefMGPinBox( long int q );
    void  dMBZeroOff(  double *dm );
    double MassICinHfe( long int fe, long int i );
    double MassHfe(long int fe);
    double MassICinMGP( long int q, long int f, long int i);
    double MassMGP(long int q, long int f  );
    double dMBfluxDir( long int q, long int i, double *dm, double fRate, double sign = 1.);
    long int LookUpXMGP( const char* MGPid );
    void  dMBflux( long int kk, double *dm );
    void  BoxComposUpdate( long int q );
    void  BoxesBCupdate();  // was CalcNodeFlux()
    void  CalcMGPdata();
    // Calculate new box states for tcur = x
    bool BoxEqStatesUpdate( long int Ni,long int pr, double x, double step );
    bool CalcSeqReacModel( char mode ); // Calculation of S-mode sequential reactors model
    bool CalcBoxFluxModel( char mode ); // integrate boxes with fluxes of Mobile Groups of Phases

    // calculate 1-step for the system of ODEs
    void Solut( double *m, double *dm, double t );
    // internal point j calculation
    void MIDEX( long int j, double t, double h );
    // ODE integration from t_begin to t_end with time step step (step can be reduced inside)
    // returns current (possibly reduced) step value or negative value in case of error
   double INTEG( double eps, double step, double t_begin, double t_end );

    double PrintPoint( long int nPoint, FILE* diffile = nullptr, FILE* logfile = nullptr, FILE* ph_file = nullptr);
    
public:

    static TGEM2MT* pm;
    
    GEM2MT *mtp;

    std::shared_ptr<TNodeArray> nodeArray()
    { return na; }
    explicit TGEM2MT( size_t nrt );

    ~TGEM2MT();

    const char* GetName() const
    {
        return "GEM2MT";
    }


    void set_def(int q);
    void mem_kill(int q);
    void mem_new(int q);

    double Reduce_Conc( char UNITP, double Xe, double DCmw, double Vm,
        double R1, double Msys, double Mwat, double Vaq, double Maq, double Vsys );

    // write/read gem2mt structure
    int ReadTask( const char *gem2mt_in1, const char *vtk_dir );
    int ReadTaskString( const std::string json_string );
    int WriteTask( const char *unsp_in1 );

    int MassTransInit( const char *lst_f_name, const char *dbr_lst_f_name );
    int MassTransStringInit(const std::string& dch_json, const std::string& ipm_json,
                            const std::vector<std::string>& dbr_json);
    void RecCalc();

    // for separate
    void checkAlws(io_formats::TRWArrays&  prar1, io_formats::TRWArrays&  prar) const;
    template<typename TIO>
    void to_text_file( TIO& out_format, bool with_comments, bool brief_mode ) const;
    template<typename TIO>
    void from_text_file(TIO& ff);

    bool userCancel;
    bool stepWise;
    bool calcFinished;
    std::string Vmessage;

   class UserCancelException {};
   bool internalCalc();
   void savePoint();

   GEMS3KGenerator GEMS3k_generator();
};

enum gem2mt_inernal {
  RMT_MODE_T = 'T',   // T Preparation of initial DataCH and DataBR files for external coupled RMT modeling using GEMS3K

  RMT_MODE_B = 'B',    // B Box-flux generic RMT simulation (with MGP or BSF fluxes)
  RMT_MODE_S = 'S',    // S Sequential reactors RMT simulation ( for moving fluids or solids )
  RMT_MODE_F = 'F',    // F Flow-through reactors RMT simulation ( with MGP fluxes )

  RMT_MODE_A = 'A',    // A 1D advection (finite-difference) coupled RMT model
  RMT_MODE_C = 'C',    // A 1D advection (Crank-Nicolson) coupled RMT model
  RMT_MODE_W = 'W',    // W random-walk advection-diffusion coupled RMT model
  RMT_MODE_D = 'D',    // D 1D diffusion (finite-volume) coupled RMT model (TBD)

  // Type of transport in box models ( Added 25.11.2011 DK )
   MGP_TT_UNDEF = '0',  // Undefined
   MGP_TT_AQS   = '1',  // Aqueous phase only (classic flush or flow)
   MGP_TT_GASF  = '2',  // Gaseous fluid phase only
   MGP_TT_AQGF  = '3',  // Aqueous and gaseous phases (all-fluids flush)
   MGP_TT_SOLID = '4'   // Total solid only (dialysis bags)

};

typedef enum {  /// Field index into outField structure
     f_PvPGD=0,  f_PvFDL,  f_PvSFL,   f_PvGrid,  f_PvDDc,
     f_PvDIc,  f_PvnVTK, f_PsMode,  f_PsSIA,  f_PsSdat,
     f_PsSdef, f_PsScom,  f_PsMO,   f_PsVTK,   f_PsMPh,
     f_nC, f_nIV,  f_nMGP,  f_nFD, f_nSFD,
     f_nEl, f_nPTypes, f_nProps, f_Nsd, f_bTau,
     f_ntM, f_nVTKfld, f_nPai, f_nTai, f_Lsf,
     f_Nf, f_FIf, f_Tau, f_sizeLc,   f_InpSys,
    f_Vsysb,  f_Mwatb,  f_Maqb,  f_Vaqb,  f_Pgb,
    f_Tmolb,  f_WmCb,  f_Asur,  ff_tf,  ff_Vt,
    ff_vp,  ff_eps,  ff_Km,  ff_al,  ff_Dif,
    ff_nto, ff_cdv,  ff_cez, f_Name, f_Note,
    f_mtWrkS, f_mtWrkF
} GEM2MT_STATIC_FIELDS;

typedef enum {  /// Field index into outField structure
     f_SDref=0, f__SDval, f__DiCp, f__FDLi, f__xFlds,
     f__mDDc,  f__mDIc,  f__mDEl,  f__HydP,   f__BSF,
     f__MB,   f__dMB, f__FDLf,   f__PGT, f__nam_i,
     f__for_e,   f__FDLid,  f__FDLop,  f__FDLmp,  f__MGPid,
     f__UMGP,  f__mGrid,  f__NPmean,  f__nPmin,  f__nPmax,
    f__ParTD
} GEM2MT_DYNAMIC_FIELDS;

#endif //_m_gem2mt_h_
