//-------------------------------------------------------------------
// $Id: m_unspace.h 861 2007-02-19 11:29:37Z gems $
// Under construction, to be included in version 3.0 (2008)
// Declaration of standalone variant of TUnSpace class and its methods
// 
// Attention: This header differs from m_unspace.h 
//   The ms_unspace.h and m)unspace.h headers cannot be used together!
//
// Copyright (C) 1995,2004,2008 K.Chudnenko, S.Dmytrieva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _ms_unspace_h_
#define _ms_unspace_h_

#include "nodearray.h"


//============================================================================
// internal definitions
//
enum grr_constants { // gstring len for graph
    MAXAXISNAME = 9,
    MAXGSNAME = 70,
    MAXGRNAME = 7
};

const unsigned int  MAXFORMULA =     81,
                      MAXPHSYMB = 8,
                      PH_RKLEN = 48,
                      V_SD_RKLEN = 32,
                      V_SD_VALEN = 24;
const int  MAXIDNAME = 12;
const   int MAXFORMUNITDT=     40;


const int UNSP_RKLEN = 80,
          NAME_SIZE = 10,
          UNSP_SIZE1 = 10,
          UNSP_SIZE2 = 8,
          PARNAME_SIZE = 20;

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


typedef struct
{ // Description of UnSpace data structure
    char
    //  the same record key as in Proces
    name[MAXFORMULA],      // Full name of UnSpace task
    notes[MAXFORMULA],     // Comments
//////*?*/    xNames[MAXAXISNAME], // Abscissa name
//////*?*/    yNames[MAXAXISNAME], // Ordinate name

/////    PunE,          // Units of energy   { j;  J c C N reserved }
/////   PunV,          // Units of volume   { j;  c L a reserved }
/////   PunP,          //  Units of pressure  { b;  B p P A reserved }
/////    PunT,          // Units of temperature  { C; K F reserved }

//  Status and control flags (+-)
    Gstat,   // UnSpace generation status: 0 -indefinite; 1 on-going generation run;
             //  2 - done generation run; 3 - generation run error (+ 5: the same in stepwise mode)
    Astat,   // UnSpace analysis status:  0 - indefinite; 1 ready for analysis;
             //  2 - analysis run; 3 - analysis done; ( + 5: the same using stepwise mode)

// Flags that control construction and sampling of UnSpace
    PsUnInt,  //  ( % 0 A 1 ) scale of uncertainty intervals for G S V: % or 0 - in %; A or 1- absolute
    PsUnInB,  //  ( % A L) scale for uncertainty intervals for B (% or 0 - in %; A or 1 - abs; L or 2 - log10)
    PsUnFltI, // ( A L ) scale for IC molality interval filters (A absolute, L log scale)
    PsUnFltD, // ( A L C) scale for DC activity/fugacity interval filters (A absolute, L log scale C mol fraction)
    PsGen[7], // (+ 1 on; - 0 off ) int zond[6]; probe type  [0]-G; [1]-S; [2]-b; [3]-T; [4]-P; [5]-V [6] pGam
/////    PsSY,     // Save generated SysEq records to data base (+ -)
/////    PsEqn,    // Will math script be specified in this UnSpace definition (+ -)
/////    PsGraph,    // Will graphics be specified in this UnSpace definition (+ -)

// Plags that control input of criteria and analysis of sampled results
    Pa_f_pha,   //  (+ -)   flag of input for the filter on phase association
    Pa_f_mol,   //  (+ -)   flag of input for filters on total IC molality
    Pa_f_fug,   //  (+ -)   flag of input for the filters on fugacity/activity
    Pa_f_mfr,   //  (+ -)   flag of input for the filters on mole fraction
    Pa_f_pH,    //  (+ -)   flag of input for the filter on pH
    Pa_f_Eh,    //  (+ -)   flag of input for the filter on Eh
    Pa_Adapt,   //  code for adaptive loops for G(298): 0-no; >=1-yes (max number of loops)
    Pa_OF,      //  code of function to construct the payoff matrix (A or 0; B or 1; C or 2; D or 3; E or 4; F or 5)
    Pa_Crit,    //  (0 1 2 3 4) OSP criterion: 0- PA frequency statistics, 1- Laplace(quantile),
                   //  2-Homeniuk(quantile) 3-Laplace (single point), 4-Homeniuk (single point)
    Pa_Zcp,     // mode of Laplace function calculation: + mean of abs.values; - abs. mean of values

// Allocation flags
    PvPOM,     //  (+ -) allocation flag for payoff matrix (if Q <= 1001)    + - + -
    PvPOR,     //  (+ -) allocation flag for payoff matrix row (if Q > 1001) + + - -
/////    PvGs,      //  (+ -) allocation flag for G0-related vectors
/////    PvSs,      //  (+ -) allocation flag for S0-related vectors
/////    PvVs,      //  (+ -) allocation flag for V0-related vectors
/////    PvPgam,    //  (+ -) allocation flag for ParGamma-related vectors (reserved)
    PvSi,      //  (+ -) flag of mode of calculation of chemical potentials in payoff function
    Pa_f_IC;    //  (+ -)   flag of input for the filter on Independent Components


// Dimensionalities related to the UnSpace problem (Project/System indexation)
    short
    N,         //   Number of independent components   (from mup->N)
    L,         //   Number of dependent components     (from mup->L)
    Ls,        //   Total number of DC in multi-comp.phases (from mup->Ls)
    Fi,        //   Number of phases                   (from mup->Fi)
    nPG,       //   Total number of uncertain input parameters (for all unspace groups)
/////    Nsd,       //   N of data source references (default 0)
// input I
    Q,     //  0 < Q < 1001 input number of sample GEM calculations to be generated
    qQ,    //  number of points in quantile
    ob,    //  number of solutions ( as optimal) that went trough superposition of all filters
// input II
    nG,    //  total number of UnSpace groups (coordinates, nGB + nGN + nGR, calc. automatically)
    nGB,   //  number of UnSpace parameter groups for sampling by Belov's grid
    nGN,   //  number of UnSpace groups used in Monte Carlo with normal distribution
    nGR,   //  number of UnSpace groups used in Monte Carlo with uniform distribution
    NgT,   //  Unspace group index for uncertain temperature T; 0- T is not uncertain
    NgP,   //  Unspace group index for uncertain pressure P; 0- P not uncertain
    NgV,   //  Unspace group index for uncertain volume constraint V, 0 - V not uncertain (reserved)

    // output
    Lapl,  //  index of optimal sample point selected by the Laplace criterion (within 0 to Q-1)
    Hurw,  //  index of optimal sample point selected by the Hurtvitz criterion
    Wald,  //  index of optimal sample point selected by the Wald criterion
    Homen, //  index of optimal sample point selected by the Homenyuk criterion
    nl,    //  number of sample points with the same phase assemblage as that selected by the Laplace criterion
    nh,    //  number of sample points with the same phase assemblage as that selected by the Hurtvitz criterion
    nw,    //  number of sample points with the same phase assemblage as that selected by the Wald criterion
    nHom,  //  number of sample points with the same phase assemblage as that selected by the Homenyuk criterion
    nPhA,  //  < Q  Number of different phase assemblages in all sampled points (variants)
    t,     // current index t of payoff matrix row (calculated sampled point) (<Q)
    q,     // current index q of payoff matrix column (input data variant) (<Q)
    g,     // current index of unspace group (< nG)
    i,     // current index of independent component (<N)
    j,     // current index of dependent component  (<L)
    k,     // current index of phase ( 0 to Fi-1)  (<Fi)
    ka;    // current index of phase assemblage  (<nPhA)

//////*?*/    dimEF[2],    // Dimensions of array of empirical data
//////*?*/    dimXY[2],    // Dimensions of data sampler tables: col.1 - N of records;

//////*?*/    axisType[6];  // axis graph type, background(3), graph type, reserved

// input
float
/*?*/    T[2],      //  Adapted/Initial Deterministic value of T (Celsius)
/*?*/    IntT[2],   //  Adapted/Initial 0 or uncertainty half-interval for T (C)
    Tc,     //  Current value of T
/*?*/   P[2],       //  Adapted/Initial Deterministic value of P (bar)
/*?*/   IntP[2],   // Adapted/Initial 0 or uncertainty half-interval for P (bar)
    Pc,     //  Current value of P
/*?*/    V[2],      //  Adapted/Initial Deterministic value of V (in cm3? L?)
/*?*/    IntV[2],   //  Adapted/Initial 0 or uncertainty half-interval for V (in cm3? L?)
    Vc,     //  Current value of V
// input intervals for filtering sampled GEM solution variants
    pH_lo,   // pH lower filter limit (default 0)
    pH_up,   // pH upper filter limit (default 14)
    Eh_lo,   // Eh lower filter limit (default -1)
    Eh_up,   // Eh upper filter limit (default 1)
    IC_lo,   // Ionic strength lower limit (default 0)
    IC_up,   // Ionic strength lower limit (default 3)
// Quantile level (size)
    quan_lev; // quantile level (default 0.05)
//////*?*/  size[2][4]; // Graph axis scale for region and fragment

// calculated
double
    CrL,     // Laplace single criterion value (from analysing payoff matrix)
    CrH,     // Hurtvitz single criterion value (from analysing payoff matrix)
    CrW,     // Wald single criterion value (from analysing payoff matrix)
    CrHom;   // Homenuk  single criterion value (from analysing payoff matrix)

// Arrays:
// Input
short
   *PbD,   //  [nG] generation codes for unspace groups: 0- Belov; 1-uniform; 2-normal
   *NgLg,  //  [L]  0-deterministic DC; value >0 is index of a group of varying G to which this DC is assigned
   *NgGam, //  [Ls] 0-deterministic DC; value >0 is index of a group varying gamma params (reserved)
   *NgLs,  //  [L]  0-deterministic DC; value >0 is index of a group varying S
   *NgLv,  //  [L]  0-deterministic DC; value >0 is index of a group varying V
   *NgNb,  //  [N]  0-deterministic DC; value >0 is index of a group varying b to which this IC is assigned
   *f_PhA; //  [N]  input phase indexes of known equilibrium phase association for a primary filter (Pa_f_pha)
                    //  in ascending order
// input uncertainty half-intervals; 0 means deterministic DC or IC
/*?*/  float (*IntLg)[2]; //  [L][2]  adapted/initial  uncertainty half-intervals for G variation (flag PsUnInt defines units)
/*?*/  float (*IntGam)[2]; //  [Ls][2]  adapted/initial uncertainty half-intervals for gamma params (reserved)
/*?*/  float (*IntLs)[2];  //  [L][2]  adapted/initial uncertainty half-intervals for S (reserved)
/*?*/  float (*IntLv)[2];  //  [L][2]  adapted/initial uncertainty half-intervals for V variation (flag PsUnInt defines units)
/*?*/  float (*IntNb)[2];  //  [N][2]  adapted/initial uncertainty half-intervals for b variation (flag PsUnInB defines units)
// input filters for sampled solution variants
float
   *m_t_lo, //  [N]  total IC molality filter - lower limits (default 0/-20) (units defined in PsUnFltI flag)
   *m_t_up, //  [N]  total IC molality filter - upper limits (default 10/+1) (units defined in PsUnFltI flag)
   *fug_lo, //  [Ls] fugacity/activity filter - lower limits (default 0/-20) (units defined in PsUnFltD flag)
   *fug_up, //  [Ls] fugacity/activity filter - upper limits (default 10/+1) (units defined in PsUnFltD flag)
   *ncp;    //  [Q][nG] table of normalised coordinates of points in uncertainty space (>= 0 <= 1)

// Input II (will be done automatically)
/*?*/ float   (*Gs)[2];     //  [L][2]  adapted/initial  copy of deterministic values of G298 for DCs
/*?*/ float   (*Ss)[2];     //  [L][2]  adapted/initial  copy of deterministic values of S298 for DCs (reserved)
/*?*/ float   (*Vs)[2];     //  [L][2]  adapted/initial  copy of deterministic values of V298 for DCs
/*?*/ float   (*GAMs)[2];   //  [Ls][2]  adapted/initial copy of deterministic Par Gamma (reserved)
/*?*/ double   (*Bs)[2];    //  [N][2]  adapted/initial  copy of deterministic values of bulk composition vector b

// collected GEM sample calculation output for creation of the payoff matrix
double
        *vG,     //   [Q][L]  G0 values used in sample input data variants (indexes q j)
/*?*/   *vB,     //   [Q][N]  B values used in sample input data variants (indexes q j)
/*?*/   *vS,     //   [Q][L]  S values used in sample input data variants (indexes q j)
/*?*/   *vmV,     //   [Q][L]  mV values used in sample input data variants (indexes q j)
/*?*/   *vNidP,     //   [Q][Ls]  GAMs values used in sample input data variants (indexes q j)

   *vY,     //   [Q][L]  x (DC mole amounts) from sample GEM solution variants (indexes t j)
   *vYF,    //   [Q][FI] XF (phase mole amounts) from sample GEM solution variants (indexes t k)
   *vGam,   //   [Q][L]  lnGam values for DC from sample GEM solution variants (indexes t j)
   *vMol,   //   [Q][N]  m_t values of total IC molality from sample GEM solution variants (indexes t i)
               // for filtering/statistics
   *vU,     //   [Q][N]  u values of dual chemical potentials from sample GEM solution variants (indexes t i)
   *vFug;   //   [Q][Ls] lga values of log10 fugacity/activity from sample GEM solution variants (indexes t j)
//////*?*/   *x0,   // Vector of abscissa dimXY[1]
//////*?*/   *y0;  // Sampled data array [Q] [dimXY[1]]
float
//////*?*/   *xE, *yE,         // Input scale data xS , yS
   *vT,     //   [Q]   TC vector of temperatures (C) from sample input data variants (index q)
   *vP,     //   [Q]   P vector of pressures (bar) from sample input data variants (index q)
/*?*/   *vV,     //   [Q]   value of V (in cm3? L?) (index q)
   *OVB,    //   [nGB+1]  pseudo-random numbers for Belov algorithm
   *OVR,    //   [nGR+1]  pseudo-random numbers by Monte Carlo with uniform distribution
   *OVN;    //   [nGN+1]  pseudo-random numbers by Monte Carlo with normal distribution

   float (*vpH)[3];    //   [Q][3]   vector of pH, Eh, IC  values from sample GEM solution variants (index t)

/*?*/float  (*quanCv)[4]; // [qQ][4] indices of sample GEM variants taken into quantile Laplace,
                     // Hurtvitz, Wald, Homenyuk (columns )
short
/*?*/  (*quanCx)[4]; // [qQ][4] Values taken into quantile Laplace,
                     // Hurtvitz, Wald, Homenyuk (columns )
// Work arrays
short
   *PhAndx,  //  [Q][N] -> [nPhA][N]  Table of system indices of phases to identify Phase Assemblages
   *sv,      //  [Q]   sv[bp->QT]  Indices of Phase Assemblage ids in sampled variants (PhAndx)
   *PhNum;   //  [nPhA]  vector of number of sample GEM solution variants containing a phase association
float
   *PhAfreq, //  [nPhA]  frequency of sample GEM solution variants containing a phase association (PhNum*100/Q)
// payoff matrix
   *pmr,     //  pointer to the beginning  of current line in payoff matrix ( PvROM/PvPOR )
   *POM,     //  [Q][Q]  payoff matrix     optional (see PvPOM flag)
   *POR;     //  [Q]  payoff matrix row    optional (see PvPOR flag)

double
   *Zcp,     //  [Q]  mean over rows of pay-off matrix for Laplace criterion
   *Zmin,    //  [Q]  minima in rows of pay-off matrix for Hurwitz criterion
   *Zmax,    //  [Q]  maxima in rows row pay-off matrix  for Hurwitz criterion
   *ZmaxAbs, //  [Q]  abs(max)values in rows of pay-off matrix for Wald criterion
   *Hom,     //  [Q]  array of calculated values for Homenyuk criterion
   *Prob;    //  [Q]  array of probabilities for calculating Homenuk criterion
   // array for statistics output
 double (*UnIC)  [UNSP_SIZE1]; //  [N][UNSP_SIZE1] statistics over independent components
 double (*UgDC)  [UNSP_SIZE1]; //  [nPG][UNSP_SIZE1] statistics set 1 over dependent components
 double (*UaDC)  [UNSP_SIZE1]; //  [Ls][UNSP_SIZE1] statistics set 2 over dependent components
 double (*UnDCA) [UNSP_SIZE2]; // [nPG][UNSP_SIZE2] statistics set 1 over dependent components

// input GEMS
 char (*UnICn)[NAME_SIZE]; //  Names of columns for output array UnIC
 char (*UgDCn)[NAME_SIZE]; //  Names of columns for output array UgDC & UaDc
/*?*/ char (*UaDCn)[NAME_SIZE]; //  Names of columns for output array UgDC & UaDc
 char  (*UnDCAn)[NAME_SIZE]; //  Names of columns for output array UnDCA

 char  (*PhAID)[8];   // [nPhA] list of phase assemblage ID names (def. 0001 0002 ... )
 char  (*PhAlst)[80]; // [nPhA] list of phase assemblage titles (made of Phase names)
///// char
/////    *Expr,     // Text with math script equations (params for activity coeffs ??? reserved )
//////*?*/    *ExprGraph, // Text with math script equations (params for activity coeffs ??? reserved )
/////    stkey[EQ_RKLEN+10],  //  currently loaded (processed) SysEq key
/////    timep[16], TCp[16], Pp[16], NVp[16], Bnamep[16],
/////    *tprn;               //  internal

//////*?*/ char (*lNam)[MAXGRNAME];   // List of ID of lines on Graph
/*?*/ char (*ParNames)[PARNAME_SIZE];   // [nPG]  List of insertain input parameters names

    char (*SGp)[MAXPHNAME];    // List of UnSpace group names [kG]
/////    char (*stl)[EQ_RKLEN];     // List of generated SysEq records [Q]
/////    char (*sdref)[V_SD_RKLEN]; // List of SDref keys to data sources [Nsd]
/////    char (*sdval)[V_SD_VALEN]; // Comments to data sources [Nsd]

    // work data
float
  *A;   // work stoichiometry matrix in project/system indexation

}
UNSPACE;



// Current UnSpace instances
class TUnSpace 
//: public TCModule
{
    UNSPACE us[2];

// !!! internaal data extracted fron GEMS modelling project
//    RMULTS* mup;
   short mup_Laq;     // LO  - of DC in aqueous phase
   short mup_Pg;      // PG  - of DC in gaseous phase
   char (*mup_SF)[MAXPHNAME+MAXSYMB];// List of PHASE definition keys [0:Fi-1]             DB   
   short  *mup_Ll;          // L1 vector, shows a number of DC included to each phase [0:Fi-1] DB

//    SYSTEM *syu; // syu = TProfil::pm->syp;
//   char  *syu_Dcl;  // DC selection switches { + * - } [0:mu.L-1]
//   float *syu_Guns;  //  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
//   float *syu_Vuns;  //  mu.L work vector of uncertainty space increments to tp->Vm
   double  *syu_B;  // Vector b of bulk chemical composition of the system, moles [0:mu.N-1]
   
//    MTPARM *tpp;
//   double *tpp_G; // Partial molar(molal) Gibbs energy g(TP) (always), J/mole 
   float *tpp_S;    // Partial molar(molal) entropy s(TP), J/mole/K
//   double *tpp_Vm;   // Partial molar(molal) volume Vm(TP) (always), J/bar
    
//    SYSTEM *syu; // syu = TProfil::pm->syp;
      MULTI *pmu; // pmu = TProfil::pm->pmp;

      TNode* na;       // pointer to nodearray class instance

protected:

// For separate mode

// memory allocation	
	 void Alloc();
	 void Free();
     void setup_defaults();

	   
//    void keyTest( const char *key );

    // internal methods
    void unsp_eqkey();           // calculate EqStat
    void setPhaseAssemb();       // set phase groups lists & counters
    void work_dyn_kill();
    void phase_lists_new();
    void work_dyn_new();
    void nG_dyn_new();
    short  near_prime_number( short num );
    bool test_sizes();
    void set_def_data_to_arrays( bool mode );   // must be changed with DK
    void init_generation();
    int calc_nPG();
    double value_nPG( int line, int q );
    void adapt_nPG( int line, double new_val, double new_int );
    void build_nPG_list();
    void init_analyse();

   //calc part big modules ( prosledit` ispol`zovanie, reorganizovat`)
    // kirpichiki dlya RecCalc
    void analyseArrays();
    void AdapG();                    // output to array  UnDCA
    void buildTestedArrays();

    //test switches for calculation flags & types, not calculate all  !!!!!!

    // must be changed
     int filters( int k );  // add new filters : f_PhA ...
     void out_QT( int Ngr );  // output to arrays UnIC, UgDC, UaDC

   // must be changed for new structure phase list (nPhA, Phndx ...)
      short kol_in_sol( int j );
      int kolgrup();
      int kvant_parag( int type /*short *mas*/ );
      void kvant_index(float per,int N,double *mas, short type /**ind*/);
      int sel_parg( short *sv );

    // calculation part  ( inernal from Kostya, not to be changed )
    void Un_criteria();             // calculates payoff matrix & criteria
    void  NexT(int J );              // Next point in space uncertainty
    void UNIFORM0( int reg );
    void BELOV(int QT, int NG, float *OVB);
    int  WL(int i,int j,int QT);
    double ePO( int i, int j );
    double ePO1( int i,int j );
    double ePO2( int i, int j );
    double ePO3( int i,int j );
    double ePO4( int i, int j );
    double ePO5( int i,int j );
    double g_u( int j, double *U,float *A,
                      double lgGAM,double x,double Xw);
    double g_uw( int j, double *U, float *A,
       double lgGAM,double x,double xjw, double Xw);
    double norrand(double *x);
    double ravrand(double *x);
    float PNTBEL(int J,int QT,int OV);

    // graphics
//    void text_analyze( int nObj );
//    void calc_graph();




public:

	// write/read the UnSpace structure data	   
	void to_text_file( fstream& ff, bool with_comments );
	void from_text_file(fstream& ff);
    void result_to_text_file( fstream& ff, bool with_comments );
    int TaskSystemInit( const char *chbr_in1 );
	int ReadTask( const char *unsp_in1 );
	int WriteTask( const char *unsp_in1 );

	
    static TUnSpace* pu;

    UNSPACE *usp;

    TUnSpace();
    ~TUnSpace();

    const char* GetName() const
    {
        return "UnSpace";
    }



//    void RecInput( const char *key );
//    void MakeQuery();
//    int RecBuild( const char *key, int mode = VF_UNDEF );
    void CalcTask( const char *key );


};

typedef enum {
  // codes of distribution law for UnSpace groups  formerly PRB_CLASSES
    UNSP_DL_BELOV    = 'b',     /* belov */
    UNSP_DL_UNIFORM  = 'u',     /* uniform */
    UNSP_DL_NORMAL   = 'n',      /* normal */

    UNSP_UN_PRS    = '%',
    UNSP_UN_ABS    = 'A',
    UNSP_UN_LOG    = 'L',
    UNSP_UN_CMOL   = 'C',

/*  GS_INDEF   = '0',
    GS_GOING    = '1',
    GS_DONE    = '2',
    GS_ERR     = '3',
    AS_INDEF   = '0',
    AS_READY   = '1',
    AS_RUN     = '2',
    AS_DONE    = '3',
*/

    UNSP_OF_A      = 'A',
    UNSP_OF_B      = 'B',
    UNSP_OF_C      = 'C',
    UNSP_OF_D      = 'D',
    UNSP_OF_E      = 'E',
    UNSP_OF_F      = 'F',

    UNSP_CRIT_PA       =  '0',
    UNSP_CRIT_LAPL_QAN  = '1',
    UNSP_CRIT_HUR_QAN  =  '2',
    UNSP_CRIT_WALD_QAN  = '3',
    UNSP_CRIT_HOME_QAN  = '4'

  //

} UNSP_CLASSES;



#endif  // _ms_unspace_h_

