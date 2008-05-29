//-------------------------------------------------------------------
// $Id: node.h 1066 2008-05-16 14:16:59Z gems $
//
// TNode class - implements a simple C/C++ interface
// between GEM IPM and FMT codes
// Works with DATACH and work DATABR structures
// without using the Tnodearray class
//
// (c) 2006,2008 S.Dmytriyeva, D.Kulik
//
// This file is part of GEMIPM2K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization
// developed in the Laboratory for Waste Management, 
//   Paul Scherrer Institute

// This file may be distributed under the licence terms
// defined in GEMIPM2K.QAL
//
// See also http://gems.web.psi.ch/
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#ifndef _node_h_
#define _node_h_

#include "m_param.h"
#include "datach.h"
#include "databr.h"

#ifndef IPMGEMPLUGIN
class QWidget;
#endif

class TNode
{
  gstring dbr_file_name;  // place for the *dbr. I/O file name 
  
protected:
   MULTI* pmm;  // Pointer to GEM IPM work data structure (see ms_multi.h)

#ifdef IPMGEMPLUGIN
       // These pointers are only used in standalone GEMIPM2K programs
    TMulti* multi;
    TProfil* profil;
#endif

    DATACH* CSD;  // Pointer to chemical system data structure CSD (DATACH)
    DATABR* CNode;  // Pointer to a work node data bridge structure (node)
         // used for exchanging input data and results between FMT and GEM IPM

    // These four values are set by the last GEM_run() call
    double CalcTime;  // GEMIPM2 calculation time, s 
    int PrecLoops,    // Number of performed IPM-2 precision refinement loops
        NumIterFIA,   // Total Number of performed FIA entry iterations 
        NumIterIPM;   // Total Number of performed IPM main iterations
    
    // Checks if given Tc and P fit within the interpolation intervals
    bool  check_TP( double& Tc, double& P );

    // Tests Tc as a grid point for the interpolation of thermodynamic data
    // Returns index in the lookup grid array or -1  if it is not a grid point
    int  check_grid_T( double& Tc );

    // Tests P as a grid point for the interpolation of thermodynamic data
    // Return index in the lookup grid array or -1 if it is not a grid point
    int  check_grid_P( double& P );

    void allocMemory();
    void freeMemory();

   // Functions that maintain DATACH and DATABR memory allocation
    void datach_realloc();
    void datach_free();
    void databr_realloc();

    // deletes fields of DATABR structure indicated by data_BR_
    // and sets the pointer data_BR_ to NULL
    DATABR* databr_free( DATABR* data_BR_ =0);

    // Binary i/o functions
    // including file i/o using GemDataStream class (with account for endianness)
      // writes CSD (DATACH structure) to a binary DCH file
    void datach_to_file( GemDataStream& ff );
      // reads CSD (DATACH structure) from a binary DCH file
    void datach_from_file( GemDataStream& ff );
      // writes node (work DATABR structure) to a binary DBR file
    void databr_to_file( GemDataStream& ff );
      // reads node (work DATABR structure) from a binary DBR file
    void databr_from_file( GemDataStream& ff );

    // Text i/o functions
      // writes CSD (DATACH structure) to a text DCH file
    void datach_to_text_file( fstream& ff, bool with_comments = true );
      // reads CSD (DATACH structure) from a text DCH file
    void datach_from_text_file( fstream& ff);
      // writes work node (DATABR structure) to a text DBR file
    void databr_to_text_file(fstream& ff, bool with_comments = true );
      // reads work node (DATABR structure) from a text DBR file
    void databr_from_text_file(fstream& ff );

    // virtual functions for interaction with TNodeArray class (not used at TNode level)
    virtual void  setNodeArray( int , int*  ) { }
    virtual void  checkNodeArray( int, int*, const char* ) { }
    virtual int nNodes()  const // virtual call for interaction with TNodeArray class
    { return 1; }

#ifndef IPMGEMPLUGIN
    // Integration in GEMS-PSI GUI environment
    // Prepares and writes DCH and DBR files for reading into the coupled code
    void makeStartDataChBR(
         TCIntArray& selIC, TCIntArray& selDC, TCIntArray& selPH,
         short nTp_, short nPp_, float Ttol_, float Ptol_,
         float *Tai, float *Pai );

    // Creates lookup arrays for interpolation of thermodynamic data 
    void G0_V0_H0_Cp0_DD_arrays(); // to be written into DCH file

    // Virtual function for interaction with TNodeArray class
    virtual void  setNodeArray( gstring& , int , bool ) { }
#endif

public:

static TNode* na;   // static pointer to this TNode class instance

#ifndef IPMGEMPLUGIN
   TNode( MULTI *apm );   // constructor for integration in GEMS environment
#else

  TNode();      // constructor for standalone GEMIPM2K or coupled program
#endif

  virtual ~TNode();      // destructor

// Typical sequence for using TNode class ----------------------------------
// (1)
// For separate coupled FMT-GEM programs that use GEMIPM2K module
// Reads in the MULTI, DATACH and optionally DATABR files prepared
// in text format from GEMS and fills out nodes in node arrays
// according to distribution vector nodeTypes ( only for compatibility
// with TNodeArray class; If TnodeArray is not used then nodeTypes = 0
// must be set).
// Optional parameter getNodT1 defines whether the DATABR (node) files
// will be read (if set to true or 1). In this case, on the TNode level,
// only the contents of last file (in the ipmfiles_lst list) will be
// accessible because all DATABR files are read into a single work DATABR
// structure. On the level of TNodeArray, the initial node contents
// from DATABR files will be distributed among nodes in T1 node array
// according to distribution list nodeTypes.
//
    int  GEM_init( const char *ipmfiles_lst_name,
                   int *nodeTypes = 0, bool getNodT1 = false);

#ifdef IPMGEMPLUGIN
//  Calls for direct coupling of a FMT code with GEMIPM2K

// (2)
// Passes (copies) the GEM input data from an already loaded DATABR
//  work structure into parameters provided by the FMT part
//
   void GEM_restore_MT(
    short  &p_NodeHandle,   // Node identification handle
    short  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                      //                                    GEM input output  FMT control
    double &p_TC,      // Temperature T, C                       +       -      -
    double &p_P,      // Pressure P, bar                         +       -      -
    double &p_Vs,     // Volume V of reactive subsystem, cm3    (+)      -      +
    double &p_Ms,     // Mass of reactive subsystem, kg          -       -      +
    double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
    double *p_dul,    // upper kinetic restrictions [nDCb]       +       -      -
    double *p_dll,    // lower kinetic restrictions [nDCb]       +       -      -
    double *p_aPH     // Specific surface areas of phases (m2/g) +       -      -
   );

// (3)
// Loads GEM input data  provided in parameters by the FMT part
// into the DATABR work structure for the subsequent GEM calculation
//
   void GEM_from_MT(
    short  p_NodeHandle,   // Node identification handle
    short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                     //                                     GEM input output  FMT control
    double p_TC,     // Temperature T, C                        +       -      -
    double p_P,      // Pressure P, bar                         +       -      -
    double p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
    double p_Ms,     // Mass of reactive subsystem, kg          -       -      +
    double *p_bIC,    // bulk mole amounts of IC [nICb]         +       -      -
    double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
    double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
    double *p_aPH  // Specific surface areas of phases (m2/g)   +       -      -
   );

// Overload - uses also xDC vector for bulk composition input to GEM
// added by DK on 09.07.2007
void GEM_from_MT(
 short  p_NodeHandle,   // Node identification handle
 short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                  //                                     GEM input output  FMT control
 double p_TC,     // Temperature T, C                        +       -      -
 double p_P,      // Pressure P, bar                         +       -      -
 double p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
 double p_Ms,     // Mass of reactive subsystem, kg          -       -      +
 double *p_bIC,    // bulk mole amounts of IC [nICb]         +       -      -
 double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
 double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
 double *p_aPH,  // Specific surface areas of phases (m2/g)   +       -      -
 double *p_xDC  // Optional: mole amounts of DCs [nDCb] - will be convoluted
                 // and added to the bIC GEM input vector
);

// Overload - uses xDC and gam vectors as old primal solution for the node
// in GEM IPM2 input when NEED_GEM_SIA flag is set for calculation
// Important! This variant works only when DATACH contains a full list of DCs
// with passed through the DATABR structure.
// added by DK on 17.09.2007
void GEM_from_MT(
 short  p_NodeHandle,   // Node identification handle
 short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                  //                                     GEM input output  FMT control
 double p_TC,     // Temperature T, C                         +      -      -
 double p_P,      // Pressure P, bar                          +      -      -
 double p_Vs,     // Volume V of reactive subsystem, cm3      -      -      +
 double p_Ms,     // Mass of reactive subsystem, kg           -      -      +
 double *p_bIC,    // bulk mole amounts of IC [nICb]          +      -      -
 double *p_dul,   // upper kinetic restrictions [nDCb]        +      -      -
 double *p_dll,   // lower kinetic restrictions [nDCb]        +      -      -
 double *p_aPH,  // Specific surface areas of phases (m2/g)   +      -      -
 double *p_xDC,  // Amounts of DCs [nDCb] - old primal soln.  +      -      -
 double *p_gam   // DC activity coeffs [nDCb] - old primal s. +      -      -
);

#endif

// (3 alternative)
// Reads work node (DATABR structure) from a file path name fname
// Parameter binary_f defines if the file is in binary format (true or 1)
// or in text format (false or 0, default)
//
   int GEM_read_dbr( const char* fname, bool binary_f=false );

// (4)
// Main call for GEM IPM calculation, returns p_NodeStatusCH value
// see databr.h for p_NodeStatusCH flag values
// Before calling GEM_run(), make sure that the node data are
// loaded using GEM_from_MT() call; after calling GEM_run(),
// check the return code and retrieve chemical speciation etc.
// using the GEM_to_MT() call
//
   int  GEM_run( bool uPrimalSol );   // calls GEM for a work node

// Returns GEMIPM2 calculation time in sec after the last call to GEM_run()
   double GEM_CalcTime();

// Returns total number of FIA + IPM iterations after the last call to GEM_run()
// More detailed info is returned via parameters by reference:
//    PrecLoops:  Number of performed IPM-2 precision refinement loops
//    NumIterFIA: Total Number of performed FIA entry iterations
//    NumIterIPM: Total Number of performed IPM main iterations
   int GEM_Iterations( int& PrecLoops, int& NumIterFIA, int& NumIterIPM ); 
   
// (5) For interruption/debugging
// Writes work node (DATABR structure) into a file path name fname
// Parameter binary_f defines if the file is to be written in binary
// format (true or 1, good for interruption of coupled modeling task
// if called in loop for each node), or in text format
// (false or 0, default). Parameter with_comments, if true, tells that 
// the text file will be written with comments for all data entries. 
//
   void  GEM_write_dbr( const char* fname,  bool binary_f=false, bool with_comments = true);

// (5a) For detailed examination of GEM work data structure:
// writes GEMIPM internal MULTI data structure into text file
// path name fname in debugging format (different from MULTI input format).
// This file cannot be read back with GEM_init()!
//
   void  GEM_print_ipm( const char* fname );

#ifdef IPMGEMPLUGIN
// (6)
// Copies GEM calculation results into parameters provided by the
// FMT part (dimensions and order of elements in arrays must correspond
// to those in currently existing DATACH structure )
//
   void GEM_to_MT(
   short &p_NodeHandle,    // Node identification handle
   short &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
   short &p_IterDone,      // Number of iterations performed by GEM IPM
                         //                                     GEM input output  FMT control
    // Chemical scalar variables
    double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
    double &p_Ms,    // Mass of reactive subsystem, kg          -      -      +     +
    double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
    double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
    double &p_IC,    // Effective molal aq ionic strength       -      -      +     +
    double &p_pH,    // pH of aqueous solution                  -      -      +     +
    double &p_pe,    // pe of aqueous solution                  -      -      +     +
    double &p_Eh,    // Eh of aqueous solution, V               -      -      +     +
    // Dynamic data - dimensions see in DATACH.H structure
    double  *p_rMB,  // MB Residuals from GEM IPM [nICb]         -      -       +     +
    double  *p_uIC,   // IC chemical potentials (mol/mol)[nICb]  -      -       +     +
    double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]  -      -       +     +
    double  *p_gam,    // activity coeffs of DC [nDCb]           -      -       +     +
    double  *p_xPH,  // total mole amounts of phases [nPHb]      -      -       +     +
    double  *p_vPS,  // phase volume, cm3/mol        [nPSb]      -      -       +     +
    double  *p_mPS,  // phase (carrier) mass, g      [nPSb]      -      -       +     +
    double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]   -      -    +     +
    double  *p_xPA  // amount of carrier in phases  [nPSb] ??       -      -    +     +
  );

#endif

// Access methods for direct or protected manipulation of CSD and DBR data
//
    DATACH* pCSD() const  // get pointer to chemical system data structure
    {     return CSD;   }

    DATABR* pCNode() const  // get pointer to work node data structure
    {        return CNode;
    }  // usage on the level of Tnodearray is not recommended !

    // These methods get contents of fields in the work node structure
    double cTC() const     // get current Temperature T, C
    {  return CNode->TC;   }

    double cP() const     // get current Pressure P, bar
    {        return CNode->P;   }

    // Setting node identification handle
    void setNodeHandle( int jj )
    {      CNode->NodeHandle = (short)jj;  }

// Useful methods facilitating the communication between DataCH (or FMT)
// and DataBR (or node) data structures for components and phases
// (i.e. between the chemical system definition and the node)

   // Returns DCH index of IC given the IC Name string (null-terminated)
   // or -1 if no such name was found in the DATACH IC name list
   int IC_name_to_xCH( const char *Name );

   // Returns DCH index of DC given the DC Name string
   // or -1 if no such name was found in the DATACH DC name list
   int DC_name_to_xCH( const char *Name );

   // Returns DCH index of Phase given the Phase Name string
   // or -1 if no such name was found in the DATACH Phase name list
   int Ph_name_to_xCH( const char *Name );

   // Returns DBR index of IC given the IC Name string
   // or -1 if no such name was found in the DATACH IC name list
   inline int IC_name_to_xDB( const char *Name )
   { return IC_xCH_to_xDB( IC_name_to_xCH( Name ) ); }

   // Returns DBR index of DC given the DC Name string
   // or -1 if no such name was found in the DATACH DC name list
   inline int DC_name_to_xDB( const char *Name )
   { return DC_xCH_to_xDB( DC_name_to_xCH( Name ) ); }

   // Returns DBR index of Phase given the Phase Name string
   // or -1 if no such name was found in the DATACH Phase name list
   inline int Ph_name_to_xDB( const char *Name )
   { return Ph_xCH_to_xDB( Ph_name_to_xCH( Name ) ); }

   // Converts the IC DCH index into the IC DBR index
   // or returns -1 if this IC is not used in the data bridge
   int IC_xCH_to_xDB( const int xCH );

   // Converts the DC DCH index into the DC DBR index
   // or returns -1 if this DC is not used in the data bridge
   int DC_xCH_to_xDB( const int xCH );

   // Converts the Phase DCH index into the Phase DBR index
   // or returns -1 if this Phase is not used in the data bridge
   int Ph_xCH_to_xDB( const int xCH );

   // Converts the IC DBR index into the IC DCH index
   inline int IC_xDB_to_xCH( const int xBR )
   { return CSD->xIC[xBR]; }

   // Converts the DC DBR index into the DC DCH index
   inline int DC_xDB_to_xCH( const int xBR )
   { return CSD->xDC[xBR]; }

   // Converts the Phase DBR index into the Phase DCH index
   inline int Ph_xDB_to_xCH( const int xBR )
   { return CSD->xPH[xBR]; }

   // Converts the Phase DCH index into the DC DCH index (for pure phases)
    int Phx_to_DCx( const int Phx );

   // Converts the Phase DCH index into the DC DCH index (1-st)
   // returns into nDCinPh number of DC included into Phx phase
    int  PhtoDC_DCH( const int Phx, int& nDCinPh );

   // Converts the Phase DBR index into the DC DBR index (1-st selected )
   // returns into nDCinPh number of DC selected into Phx phase
    int  PhtoDC_DBR( const int Phx, int& nDCinPh );

    // Data exchange methods between GEMIPM and work node DATABR structure
    // Are called inside of GEM_run()
    void packDataBr();   //  packs GEMIPM calculation results into work node structure
    void unpackDataBr( bool uPrimalSol ); //  unpacks work DATABR content into GEMIPM data structure

    // Access to interpolated thermodynamic data from DCH structure
    // Test Tc and P as grid point for the interpolation of thermodynamic data
    // Return index in grid matrix or -1
     int  check_grid_TP(  double& Tc, double& P );
    // Access to interpolated G0 from DCH structure ( xCH is the DC DCH index)
     double  DC_G0_TP( const int xCH, double& Tc, double& P );
    // Access to interpolated V0 from DCH structure ( xCH is the DC DCH index)
     double  DC_V0_TP( const int xCH, double& Tc, double& P );

// To be provided - access to interpolated thermodynamic data from DCH structure
//  DC_H0_TP
//  DC_S0_TP
//  DC_Cp0_TP
//  DC_DD_TP

     // Retrieval of Phase Volume ( xBR is DBR phase index)
      double  Ph_Volume( const int xBR );
     // Retrieval of Phase mass ( xBR is DBR phase index)
      double  Ph_Mass( const int xBR );
     // Retrieval of Phase composition ( xBR is DBR phase index)
      void  Ph_BC( const int xBR, double *ARout=0 );


#ifndef IPMGEMPLUGIN
// These calls are used only inside the GEMS-PSI GEM2MT module

    // Makes start DATACH and DATABR data using GEMS internal data (MULTI and other)
    // interaction variant (the user must select ICs, DCs and phases to be included
    // in DATABR lists)
    void MakeNodeStructures( QWidget* par, bool select_all,
             float *Tai, float *Pai, short nTp_ = 1 ,
             short nPp_ = 1 , float Ttol_ = 1., float Ptol_ =1. );

    // Overloaded variant - takes lists of ICs, DCs and phases according to
    // already existing index vectors axIC, axDC, axPH (with anICb, anDCb,
    // anPHb, respectively)
    void MakeNodeStructures(  short anICb, short anDCb,  short anPHb,
             short* axIC, short* axDC,  short* axPH,
             float* Tai, float* Pai,  short nTp_,
             short nPp_, float Ttol_, float Ptol_  );

#endif
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Data direct access macroses (for FMT programs written in C++)
// Work on both sides of assignment - use with caution!
//
// Molar mass of Independent Component with node DATABRIDGE index ICx
#define nodeCH_ICmm( ICx )  (  TNode::na->pCSD()->ICmm[ \
                               TNode::na->pCSD()->xIC[(ICx)]] )

// Molar mass of Dependent Component with node DATABRIDGE index DCx
#define nodeCH_DCmm( DCx )  (  TNode::na->pCSD()->DCmm[ \
                               TNode::na->pCSD()->xDC[(DCx)]] )

// Redo into a function with interpolation
// Diffusion coefficient of dependent component with node DBr index ICx
// #define nodeCH_DD( DCx )    ( TNode::na->pCSD()->DD[
//                              TNode::na->pCSD()->xDC[(DCx)]] )

// stoichiometry coefficient A[j][i] of IC with node DATABRIDGE index ICx
// in the formula of DC with node index DCx
#define nodeCH_A( DCx, ICx )  ( (double)(TNode::na->pCSD()->A[ \
                                 (TNode::na->pCSD()->xIC[(ICx)])+ \
                                 (TNode::na->pCSD()->xDC[(DCx)]) * \
                                  TNode::na->pCSD()->nIC]) )

// more will be added soon!

#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// _node_h_

