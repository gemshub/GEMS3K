//-------------------------------------------------------------------
// $Id: ms_gem2mt.h 693 2006-03-30 17:03:55Z gems $
// To be finalized in Version 3.0 (2006)
// Declaration of GEM2MT class (separete mode), config and calculation functions
//
// Copyright (C) 2005 D.Kulik, S. Dmytriyeva
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _ms_gem2mt_h_
#define _ms_gem2mt_h_

#include "particlearray.h"

//============================================================================
// internal
enum grr_constants { // gstring len for graph
    MAXAXISNAME = 9,
    MAXGSNAME = 70,
    MAXGRNAME = 7
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
    AS_DONE    = '3',


} GS_AS_CLASSES;

#define SIZE_HYDP 7

typedef struct
{ // Description of GEM2MT data structure (prototype of 04 Dec. 2007)
  // Record key format: the same as Proces
  char
   PunE,          // Units of energy   { j;  J c C N reserved }
   PunV,          // Units of volume   { j;  c L a reserved }
   PunP,          //  Units of pressure  { b;  B p P A reserved }
   PunT,          // Units of temperature  { C; K F reserved }

// Allocation flags
   PvICi,    // Use IC quantities for initial system compositions? { + * - }
   PvAUi,    // Use formula units for initial sub-system compositions? { + * - }
   PvMO,     // Use non stop debug output for nodes (+ -)?
   PvMSg,    // Use math script for graphic presentation (+ -)?
   PvEF,      // Use empirical data for graphics  (+ -)?
   PvPGD,    // Use phase groups definitions (+ -)?
   PvFDL,    // Use flux definition list (+ -)
   PvSFL,    // Use source fluxes and elemental stoichiometries for them? (+ -)
   PvGrid,    // Use array of grid point locations? (+ -)
   PvRes,    // reserved (+ -)

     // Controls on operation
   PsMode,  // Code of GEM2MT mode of operation { S F A D T }
//  Status and control flags (+ -)
   gStat,  // GEM2MT generation status: 0 -indefinite; 1 on-going generation run;
        //  2 - done generation run; 3 - generation run error (+ 5: the same in stepwise mode)
   iStat,  // GEM2MT iteration status:  0 - indefinite; 1 ready for analysis;
        //  2 - analysis run; 3 - analysis done; ( + 5: the same using stepwise mode)
   PsSYd,  // Save generated SysEq records to data base (+ -)
   PsSdat, //  Save DataCH and inital DataBR files as text files (+) or binary (-)
   PsDDc,  //  Use diffusion coefficients for DC - DDc vector (+ -)
   PsDIc,  //  Use diffusion coefficients for IC - DIc vector (+ -)
   PsTPai,  //  Create T,P values in Tval, Pval using iterator(+) or enter(-)

   name[MAXFORMULA],  //  Name of GEM2MT task
   notes[MAXFORMULA]  //  Comments
    ;
  short
// input I
   nC,   // nQ - number of local equilibrium compartments (nodes)
   xC, yC, zC,  // numbers of nodes along x, y, z coordinates
   nIV,  // number of initial variants of the chemical system, nIV <= nC
   nPG,  // number of mobile phase groups (0 or >1)
   nFD,  // number of MPG flux definitions (0 or >1)
   nSFD,   // number of source flux definitions (0 or < nFD )
   nEl, // number of electrolytes for setting up electrolyte diffusion coefficients in mDEl vector
   nPTypes,     // res Number of allocated particle types (< 20 ? )
   nProps,      // res Number of particle statistic properties (for monitoring) >= anPTypes
   Lbi,  // Lb - number of formula units to set compositions in initial variants
   Nsd,  // N of references to data sources
   Nqpt, // Number of elements in the script work array qpi for transport
   Nqpg, // Number of elements in the script work array qpc for graphics
   Nb,   // N - number of independent components (set automatically from RMults)
   FIb,   // N - number of phases (set automatically from RMults)
   Lb,   // N - number of dependent components in multycomponent phases (set automatically from RMults)
   bTau, // Time point for the simulation break (Tau[0] at start)
   ntM, // Maximum allowed number of time iteration steps (default 1000)
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

// iterators for generating syseq record keys for initial system variants
   tmi[3],   // SYSTEM CSD definition #: start, end, step (initial)
   NVi[3],    // Restrictions variant #: start, end, step
   axisType[6],  // axis graph type, background(3), graph type, reserved
// These lists of indices connect the DATABR arrays with this structure
    *xIC,   // ICNL indices in DATABR IC vectors [nICb]
    *xDC,   // DCNL indices in DATABR DC list [nDCb]
    *xPH,   // PHNL indices in DATABR phase vectors [nPHb]
    *NPmean,       // array of initial mean particle type numbers per node ( size: nPTypes )
    *nPmin,        // minimum average total number of particles of each type per one node nPTypes
    *nPmax;        // maximum average total number of particles of each type per one node nPTypes
 short (*DiCp)[2];     // array of indexes of initial system variants for
                       // distributing to nodes [nC]
 short (*ParTD)[6]; // array of particle type definitions at t0 or after interruption nPTypes

 short (*FDLi)[2]; //[nFD][2] Indexes of nodes where this flux begins and ends
              // negative value means one-side flux (source or sink)
         // for source fluxes, -2 means "source flux stoichiometry with index 1
         // line in the BSF table", and so on
   ;
  float        // input
   *Pi,    // Pressure P, bar for initial systems (values within Pai range) [nIV]
   *Ti,    // Temperature T, C for initial systems (values within Pai range) [nIV]
   *Vi;    // Volume of the system (L) [nIV] usually zeros
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
   fVel,   // Advection/diffusion mass transport: initial fluid advection velocity (m/sec)
   cLen,   // column length (m)
   tf,     // time step reduction factor
   cdv,    // cutoff factor for differences (1e-9)
   cez,    // cutoff factor for minimal amounts of IC in node bulk compositions (1e-12)
   al_in,  // initial value of longitudinal dispersivity (m), usually 1e-3
   Dif_in, // initial general diffusivity (m2/sec), usually 1e-9
   ADrs3,
   ADrs4,
   ADrs5,
   ADrs6
  ;
  float
// Iterators for MTP interpolation array in DataCH
   Pai[4],  // Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
   Tai[4],  // Temperature T, C: start, end, increment for MTP array in DataCH , Ttol
   Tau[3],  // Physical time iterator
// graphics
   size[2][4], // Graph axis scale for the region and the fragment
   sizeLc[3],   // spatial dimensions of the medium defines topology of nodes
   *xEt,    // Abscissa for experimental points [nXE]
   *yEt,       // Ordinates for experimental points to plot [nXE, nYE]
   *DDc,  //  [Lsf] diffusion coefficients for DC
   *DIc,  //  [Nf] diffusion coefficients for IC
   *DEl   //  [nE] diffusion coefficients for electrolyte salts
    ;
float (*grid)[3];      // Array of grid point locations, size is nC

 double
   *Bn,    //  [nIV][Nb] Table of bulk compositions of initial systems
   *qpi,   //  [Nqpi] Work array for initial systems math script
   *qpc,    //  [Nqpc] Work array for mass transport math script,
   *xt,    //  Abscissa for sampled data [nS]
   *yt,     //  Ordinates for sampled data [nS][nYS]
   *BSF,    // [nSFD][Nf] table of bulk compositions of source fluxes
            //  More to be added here for seq reactors?
   *MB,  // [nC] [Nf] column of current masses of boxes (in kg)
   *dMB // [nC][Nf]  Table of current derivatives dM for elements in reservoirs
    ;
 double  (*HydP)[SIZE_HYDP]; // [nC][6] hydraulic parameters for nodes in mass transport model
   //  value order to be described
 float
   *Tval,   // discrete values of T [nTai] in grid arrays in DataCH
   *Pval,   // discrete values of P [nPai]

   *CIb, // [nIV][Nb] Table of quantity/concentration of IC in initial systems
   *CAb, // [nIV][Lbi] Table of quantity/concentration of formulae for initial systems
   *PGT  // Quantities of phases in MPG [FIf][nPG]
    ;
 float  (*FDLf)[4]; // [nFD][4] Part of the flux defnition list (flux order, flux rate, MPG quantities)
 char
   *tExpr,  // Math script text for calculation of mass transport
   *gExpr,  // Math script text for calculation of data sampling and plotting
//
   *CIclb, // [Nb] Units of IC quantity/concentration for initial systems compositions
   *AUcln, // [Lbi] Units of setting UDF quantities for initial system compositions
   *UMPG;  // [FIf] units for setting phase quantities in MPG (see PGT )

  char (*sdref)[V_SD_RKLEN]; // "List of bibl. refs to data sources" [0:Nsd-1]
  char (*sdval)[V_SD_VALEN];  // "Parameters taken from the respective data sources"[0:Nsd-1]
  char (*nam_i)[MAXIDNAME]; // [nIV][12] id names of initial systems
  char (*for_i)[MAXFORMUNITDT]; // [Lbi][40] formulae for setting initial system compositions
  char (*for_e)[MAXFORMUNITDT]; // [nE][40] formulae for diffusing dissolved electrolytes
  char (*stld)[EQ_RKLEN]; // List of SysEq record keys for initial systems [nIV]
  char (*FDLid)[MAXSYMB]; // [nFD] ID of fluxes
  char (*FDLop)[MAXSYMB]; // [nFD] Operation codes (letters) flux type codes
  char (*FDLmp)[MAXSYMB]; // [nFD] ID of MPG to move in this flux  dim changed!
  char (*MPGid)[MAXSYMB]; // [nPG] ID list of mobile phase groups
//
  char (*SBM)[MAXICNAME+MAXSYMB];  // Keys (names) of IC
//  graphics
  char  xNames[MAXAXISNAME],        // Abscissa name
    yNames[MAXAXISNAME];      // Ordinate name
  char  (*lNam)[MAXGRNAME];        // List of ID of lines on Graph [nYS]
  char  (*lNamE)[MAXGRNAME];       // List of ID of lines of empirical data [nYE]


/* Work arrays */
 float
   *An,  // [Lbi][N] stoich matrix for formula units from the for_i list
    *Ae  // [nE][N] stoich matrix for for diffusing electrolytes
   ;
 double
 *gfc,  // [nC][nPG][Nf] Array of element partition coefficients between MPG and its source box
 *yfb;  // [nC][nPG][Nf] Array of MPG bulk compositions at current time point 
   ;
   char sykey[EQ_RKLEN+10],    // Key of currently processed SysEq record
   *etext,              // internal
   *tprn;              // internal
//work data
 short
   ctm,    // current CSD #
   cnv,    //  current restriction variant #
   qc,     // current index of the compartment ( 1 to nC )
   kv,     // current index of the initial system variant (1 to nIV )
   jqc,    // script c-style index (= qc-1) for transport
   jqs,    // script c-style index (= qc-1) for graphics
   jt,     // current index of sampled point (for sampling scripts)
   jdd,   // current index of diffusing DC
   jdi,   // current index of diffusing IC
   ide,   // current index of diffusing electrolyte
   ct,  // actual time iterator
   rei5
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

    TNodeArray* na;       // pointer to nodearray class instance
    TParticleArray* pa;       // pointer to TParticleArray class instance

protected:

   void Alloc();
   void Free();


    void  copyNodeArrays();
    int CheckPIAinNodes1D( char mode, int start_node = 0, int end_node = 1000 );
    bool  CalcIPM( char mode, int start_node = 0,
         int end_node = 1000, FILE* diffile = NULL );
    void  MassTransAdvecStart();
    void  MassTransAdvecStep( bool ComponentMode = true );
    void  MassTransParticleStart();
    void  MassTransParticleStep( bool ComponentMode = true );


    
    // for box floor model
    int MaxIter,  // max number of iterations
         nfcn,      // number of functional estimates
         nstep,     // number of steps
         naccept,   // number of permissible steps
         nrejct;    // number of unpermissible steps
    double *x;
    double *dx;
    double *tv;

    double (*tt)[9];

   
//    bool CalcBoxModel( char mode ); // calculate Mobile Phase-Group Flows
    int LookUpXMGP( const char* MGPid );
    // calculate 1-step from system of equation 
    void Solut( double *m, double *dm, double t );
    // Calculate new reservuir states for tcur = x
    bool CalcNewStates( int Ni,int pr, double x, double step );
    // internal point j calculation
    void MIDEX( int j, double t, double h );
    void INTEG( double eps, double& step, double t_begin, double t_end );

    
public:

	// write/read unspace structure	   
	void to_text_file( fstream& ff, bool with_comments );
	void from_text_file(fstream& ff);
//    void result_to_text_file( fstream& ff, bool with_comments );
	int ReadTask( const char *unsp_in1 );
	int WriteTask( const char *unsp_in1 );
	
//   int MassTransSetUp( const char *gem2mt_in1 );
   int MassTransInit( const char *chbr_in1 );
   bool Trans1D( char mode  );  // return true if canceled
   bool CalcBoxModel( char mode ); // calculate Mobile Phase-Group Flows

    static TGEM2MT* pm;

    GEM2MT *mtp;

    TGEM2MT();
    ~TGEM2MT();
    void RecCalc( );

    const char* GetName() const
    {
        return "GEM2MT";
    }

};

enum gem2mt_inernal {
//              A_CHI = 1,  A_GAM =2,
  GMT_MODE_S = 'S',   // Preparation of initial DataCH and DataBR files for external coupled RMT modeling using GEMIPM2K
  GMT_MODE_F = 'F',   // F Flux-box RMT scoping model (MPG fluxes, sequential reactors etc.)
  GMT_MODE_A = 'A',   // A 1D advection (numerical) coupled RMT scoping model
  GMT_MODE_D = 'D',   // D 1D diffusion (numerical) coupled RMT scoping model
  GMT_MODE_T = 'T',    // T 1D advection-diffusion coupled RMT scoping model
  GMT_MODE_W = 'W',    // W random-walk advection-diffusion coupled RMT scoping model
  GMT_MODE_V = 'V'    // V const volume coupled RMT scoping model

};


#endif //_ms_gem2mt_h_
