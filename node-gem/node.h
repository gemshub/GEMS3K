//-------------------------------------------------------------------
// TNode class - implements a simple C/C++ interface
// between GEM IPM and FMT codes
// Works with DATACH and work DATABR structures
// without using the nodearray class
//
// Written by S.Dmytriyeva,  D.Kulik
// Copyright (C) 2006 S.Dmytriyeva, D.Kulik
//
// This file is part of GEMIPM2K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization

// This file may be distributed under the licence terms
// defined in GEMIPM2K.QAL
//
// See also http://les.web.psi.ch/Software/GEMS-PSI
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

protected:

   MULTI* pmm;  // Pointer to GEMIPM work data structure (see ms_multi.h)

#ifdef IPMGEMPLUGIN
               // This is used in the isolated GEMIPM2K module for coupled codes
    TMulti* multi;
    TProfil* profil;

#endif

    DATACH* CSD;     // Pointer to chemical system data structure CSD (DATACH)

    DATABR* CNode;   // Pointer to a work node data bridge structure (node)
      // used for exchanging input data and results between FMT and GEMIPM

    void check_TP();  // Checks if given T and P fit within interpolation intervals
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
    // including file i/o using GemDataStream class (different endianness)
      // writes CSD (DATACH structure) to a binary file
    void datach_to_file( GemDataStream& ff );
      // reads CSD (DATACH structure) from a binary file
    void datach_from_file( GemDataStream& ff );
      // writes node (work DATABR structure) to a binary file
    void databr_to_file( GemDataStream& ff );
      // reads node (work DATABR structure) from a binary file
    void databr_from_file( GemDataStream& ff );

    // Text i/o functions
      // writes CSD (DATACH structure) to a text file
    void datach_to_text_file( fstream& ff );
      // reads CSD (DATACH structure) from a text file
    void datach_from_text_file( fstream& ff);
      // writes work node (DATABR structure) to a text file
    void databr_to_text_file(fstream& ff );
      // reads work node (DATABR structure) from a text file
    void databr_from_text_file(fstream& ff );

    // virtual functions for interaction with nodearray class (not used at TNode level)
    virtual void  setNodeArray( int , int*  ) { }
    virtual void  checkNodeArray( int, int*, const char* ) { }
    virtual int nNodes()  const // virtual call for interaction with nodearray class
    { return 1; }

#ifndef IPMGEMPLUGIN
    // Integration in GEMS
    // - prepares DATACH and DATABR files for reading into the coupled code
    void makeStartDataChBR(
         TCIntArray& selIC, TCIntArray& selDC, TCIntArray& selPH,
         short nTp_, short nPp_, float Ttol_, float Ptol_,
         float *Tai, float *Pai );

    // creates arrays of thermodynamic data for interpolation
    void G0_V0_H0_Cp0_DD_arrays(); // which are written into DATACH file

    // virtual function for interaction with nodearray class
    virtual void  setNodeArray( gstring& , int , bool ) { }
#endif

public:

static TNode* na;   // static pointer to this class

#ifndef IPMGEMPLUGIN
   TNode( MULTI *apm );   // constructor for integration in GEMS
#else

  TNode();      // constructor for GEMIPM2K
#endif

    virtual ~TNode();      // destructor

// Typical method call sequence ----------------------------------------------
// (1)
    // For separate coupled FMT-GEM programs that use GEMIPM2K module
    // Reads in the MULTI, DATACH and DATABR files prepared from GEMS
    // and fills out nodes in node arrays according to distribution vector
    // nodeTypes ( only for TNodeArray ). If TnodeArray is not used then
    // nodeTypes = 0 must be set
    int  GEM_init( const char *ipmfiles_lst_name,
                   int *nodeTypes = 0, bool getNodT1 = false);

#ifdef IPMGEMPLUGIN
//  Calls for direct coupling of a FMT code with GEMIPM2K

// (2)
// Loads GEM input data  provided in parameters by the FMT part
// into the DATABR work structure for the subsequent GEM calculation
   void GEM_from_MT(
    short  p_NodeHandle,   // Node identification handle
    short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                     //                                     GEM input output  FMT control
    double p_T,      // Temperature T, K                        +       -      -
    double p_P,      // Pressure P, bar                         +       -      -
    double p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
    double p_Ms,     // Mass of reactive subsystem, kg          -       -      +
    double *p_bIC,    // bulk mole amounts of IC [nICb]         +       -      -
    double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
    double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
    double *p_aPH  // Specific surface areas of phases (m2/g)   +       -      -
   );

// (3)
// Passes (copies) the GEM input data from an already loaded DATABR
//  work structure into parameters provided by the FMT part
   void GEM_restore_MT(
    short  &p_NodeHandle,   // Node identification handle
    short  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                      //                                    GEM input output  FMT control
    double &p_T,      // Temperature T, K                        +       -      -
    double &p_P,      // Pressure P, bar                         +       -      -
    double &p_Vs,     // Volume V of reactive subsystem, cm3    (+)      -      +
    double &p_Ms,     // Mass of reactive subsystem, kg          -       -      +
    double *p_bIC,    // bulk mole amounts of IC [nICb]          +       -      -
    double *p_dul,    // upper kinetic restrictions [nDCb]       +       -      -
    double *p_dll,    // lower kinetic restrictions [nDCb]       +       -      -
    double *p_aPH     // Specific surface areas of phases (m2/g) +       -      -
   );

#endif

   // (2 alternative)
   // reads work node (DATABR structure) from a file
   int  GEM_read_dbr( bool binary_f, const char *fname );

   // (4)
   // Main call for GEM IPM calculation, returns p_NodeStatusCH value
   int  GEM_run();   // calls GEM for a work node

   // (5 debugging/interruption)
   // For examining GEM calculation results:
   // Prints current multi and/or work node structures to files with
   // names given in the parameter list (if any parameter is NULL
   // then writing the respective file is skipped)
   void  GEM_printf( const char* multi_file, const char* databr_text,
                         const char* databr_bin );

#ifdef IPMGEMPLUGIN
//  Calls for direct coupling of an FMT code with GEMIPM2K
   // (5)
   // Copies GEM calculation results into parameters provided by the
   // FMT part (dimensions and order of elements in arrays must correspond
   // to those in currently existing DATACH structure )
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
    double cT() const     // get current Temperature T, K
    {        return CNode->T;   }

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
   int IC_name_to_x( const char *Name );

   // Returns DCH index of DC given the DC Name string
   // or -1 if no such name was found in the DATACH DC name list
   int DC_name_to_x( const char *Name );

   // Returns DCH index of Phase given the Phase Name string
   // or -1 if no such name was found in the DATACH Phase name list
   int Ph_name_to_x( const char *Name );

   // Returns DBR index of IC given the IC Name string
   // or -1 if no such name was found in the DATACH IC name list
   inline int IC_name_to_xDB( const char *Name )
   { return IC_xCH_to_xDB( IC_name_to_x( Name ) ); }

   // Returns DBR index of DC given the DC Name string
   // or -1 if no such name was found in the DATACH DC name list
   inline int DC_name_to_xDB( const char *Name )
   { return DC_xCH_to_xDB( DC_name_to_x( Name ) ); }

   // Returns DBR index of Phase given the Phase Name string
   // or -1 if no such name was found in the DATACH Phase name list
   inline int Ph_name_to_xDB( const char *Name )
   { return Ph_xCH_to_xDB( Ph_name_to_x( Name ) ); }

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
   inline int IC_xBR_to_xCH( const int xBR )
   { return CSD->xIC[xBR]; }

   // Converts the DC DBR index into the DC DCH index
   inline int DC_xBR_to_xCH( const int xBR )
   { return CSD->xDC[xBR]; }

   // Converts the Phase DBR index into the Phase DCH index
   inline int Ph_xBR_to_xCH( const int xBR )
   { return CSD->xPH[xBR]; }

   // Converts the Phase DCH index into the DC DCH index (for pure phases)
    int Phx_to_DCx( const int Phx );

    // Data exchange methods between GEMIPM and work node DATABR structure
    // Are called inside of GEM_run()
    void packDataBr();      //  packs GEMIPM calculation results into work node structure
    void unpackDataBr();    //  unpacks work node structure into GEMIPM data structure

// To be provided - access to interpolated thermodynamic data from DCH structure
//  G0TP
//  V0TP
//  H0TP
//  S0TP
// Cp0TP
//  DDTP

#ifndef IPMGEMPLUGIN
// These calls are used only inside of the GEMS-PSI GEM2MT module

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
// Molar mass of Independent Component with node index ICx
#define nodeCH_ICmm( ICx )  (  TNode::na->pCSD()->ICmm[ \
                               TNode::na->pCSD()->xIC[(ICx)]] )

// Molar mass of Dependent Component with node index DCx
#define nodeCH_DCmm( DCx )  (  TNode::na->pCSD()->DCmm[ \
                               TNode::na->pCSD()->xDC[(DCx)]] )

// Redo into a function with interpolation
// Diffusion coefficient of dependent component with node index ICx
// #define nodeCH_DD( DCx )    ( TNode::na->pCSD()->DD[
//                              TNode::na->pCSD()->xDC[(DCx)]] )

// stoichiometry coefficient A[j][i] of IC with node index ICx
// in the formula of DC with node index DCx
#define nodeCH_A( DCx, ICx )  ( (double)(TNode::na->pCSD()->A[ \
                                 (TNode::na->pCSD()->xIC[(ICx)])+ \
                                 (TNode::na->pCSD()->xDC[(DCx)]) * \
                                  TNode::na->pCSD()->nIC]) )

// more will be added soon!

#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// _node_h_

