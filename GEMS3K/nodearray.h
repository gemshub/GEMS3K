//-------------------------------------------------------------------
// $Id$
/// \file nodearray.h
/// Contains declaration of TNodeArray class implementing an advanced
/// interface for development of coupled codes involving GEMS3K.
//
/// \class TNodeArray nodearray.h
/// Implements an advanced (level 2) C/C++ interface with GEMS3K for the
/// development of coupled reactive transport codes.
/// Works with DATACH and an array of DATABR structures; uses TNode class
//
// Copyright (C) 2006-2017 S.Dmytriyeva, D.Kulik
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
//

#ifndef _nodearray_h_
#define _nodearray_h_

#include <functional>
#include "node.h"

/// The function is executed when ProcessProgress
using  ProcessProgressFunction = std::function<bool( const gstring& message, long point )>;

// These structures are needed for implementation of Random Walk and
// similar particle-based transport algorithms
enum  PTCODE /// Codes of particle type
{
   DISSOLVED = 20,
   ADVECTIVE = 21,
   DIFFUSIVE = 22,
   COLLOID = 23
};

struct  LOCATION /// Location (coordinates) of a point in space
                 /// for implementation of particle transport algorithms
{  double x,
        y,
        z;

  LOCATION():
   x(0.), y (0.), z(0.) {}

  LOCATION( double ax, double ay=0., double az=0. ):
     x(ax), y (ay), z(az) {}

  LOCATION( const LOCATION& loc ):
        x(loc.x), y (loc.y), z(loc.z) {}

//   const LOCATION& operator= (const LOCATION& loc);
//   {
//      x = loc.x;
//      y = loc.y;
//      z = loc.z;
//    }
};

/// Struct to setup GEMRun properties
struct TestModeGEMParam
{
    char mode;        ///   mode can be NEED_GEM_PIA (smart algorithm) or NEED_GEM_AIA
    char useSIA;      /// Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    long int step;    /// extern actual time iterator
    double cdv;       /// cutoff for IC amount differences in the node between time steps (mol, 1e-9)
    double cez;       /// cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12)

    TestModeGEMParam( char amode, char  auseSIA,
                      long int astep,   double acdv,  double acez ):
       mode(amode), useSIA(auseSIA), step(astep), cdv(acdv), cez(acez) {}
};


class TParticleArray;

// Definition of TNodeArray class
class TNodeArray
{
protected:

    TNode calcNode;

    DATABR* (*NodT0);  ///< array of nodes for previous time point
    DATABR* (*NodT1);  ///< array of nodes for current time point

    long int anNodes;  ///< Number of allocated nodes
    long int sizeN;	   ///< Number of nodes along x direction
    long int sizeM;    ///< Number of nodes along y direction
    long int sizeK;    ///< Number of nodes along z direction

    LOCATION size;     ///< spatial dimensions of the medium ( x, 0, 0 - 1D; x,y,0 - 2D; x,0,z - 2D; x,y,z - 3D )
                       ///< defines topology of nodes (N of grid points per node): 1D- 2; 2D- 4; 3D- 8 )
                       ///< relative to coordinate origin (0,0,0) units
    LOCATION* grid;    ///< Array of grid point locations, size is anNodes

    char* tcNode;      ///< Node type codes (see databr.h), size anNodes
    bool* iaNode;      ///< GEM IA status for all nodes (true: NEED_GEM_AIA, false: NEED_GEM_SIA)
    
    void allocMemory();
    void freeMemory();

   /// Prototypes of functions to manage location of particles
   /// within nodes relative to the whole grid of the node walls
   LOCATION getGrid( long int iN, long int jN, long int kN ) const;

   /// Test if the location cxyz resides in the node ( ii,jj,kk )
   bool isLocationInNode( long int ii, long int jj, long int kk, LOCATION cxyz ) const;
   /// Test if the location cxyz resides in the node with absolute index iNode
   bool isLocationInNode( long int iNode, LOCATION cxyz ) const;

   /// Copies data from the work DATABR structure into the node ndx in
   /// the node arrays NodT0 and NodT1  (as specified in nodeTypes array)
   void  setNodeArray( long int ndx, long int* nodeTypes  );

   /// Test setup of the boundary condition for all nodes in the task
   void  checkNodeArray( long int i, long int* nodeTypes, const char*  datachbr_file );

   ///  Copies data for a node ndx from the array of nodes anyNodeArray that
   /// contains nNodes into the work node data bridge structure
   ///  (implementation thread-safe)
   void CopyWorkNodeFromArray( TNode& wrkNode, long int ndx, long int nNodes, DATABRPTR* anyNodeArray );

   ///  Moves work node data to the ndx element of the node array anyNodeArray
   /// that has nNodes. Previous contents of the ndx element will be lost,
   /// work node will be allocated new and will contain no data
   ///  (implementation thread-safe)
   void MoveWorkNodeToArray( TNode& wrkNode, long int ndx, long int nNodes, DATABRPTR* anyNodeArray );

   /// New Stuff--------------------------------------------------------------

   /// Here we compare this node for current time and for previous time to test need to recalculate equilibrium
   /// in this node because its vector b has changed
   /// Zeroing charge off in bulk composition
   bool NeedGEMS( TNode& wrkNode, const TestModeGEMParam& modeParam, DATABR* C0, DATABR* C1  );

   ///  Testing indicators for IA in the ii node for smart algorithm
   long int SmartMode( const TestModeGEMParam& modeParam, long int ii,  bool* piaN  );

   /// Build bad GEM result message
   gstring ErrorGEMsMessage( long int RetCode,  long int ii, long int step  );

   ///  Here we do a GEM calculation in box ii (implementation thread-safe)
   bool CalcIPM_Node(  const TestModeGEMParam& modeParam, TNode& wrkNode,
                       long int ii, DATABRPTR* C0, DATABRPTR* C1, bool* iaN, FILE* diffile  );

   // alloc new memory
    DATABR * allocNewDBR( TNode& wrkNode)
    {
        DATABR *node = new DATABR;
        wrkNode.databr_reset( node, 1 );
        wrkNode.databr_realloc(node);
        return node;
    }

    // end of new stuff -------------------------------------------------------

public:

  static TNodeArray* na;   ///< static pointer to this class

  const TNode& getCalcNode()
  { return calcNode;}

  DATACH* pCSD() const  /// Get the pointer to chemical system definition data structure
  {  return calcNode.pCSD();   }

  /// Returns number of temperature and  pressure grid points for one dependent component
  inline long int gridTP() const
  {  return calcNode.gridTP();   }


  /// Retrieves the stoichiometry coefficient a[xdc][xic] of IC in the formula of DC.
  /// \param xdc is DC DBR index
  /// \param xic is IC DBR index
  inline double DCaJI( const long int xdc, const long int xic) const
  { return calcNode.DCaJI( xdc, xic); }

  /// Retrieves the molar mass of Independent Component in kg/mol.
  /// \param xic is IC DBR index
  inline double ICmm( const long int xic ) const
  { return calcNode.ICmm(  xic); }

  /// Retrieves the molar mass of Dependent Component in kg/mol.
  /// \param xdc is DC DBR index
  inline double DCmm( const long int xdc ) const
  { return calcNode.DCmm( xdc); }

  /// Converts the Phase DBR index into the Phase DCH index
  inline long int Ph_xDB_to_xCH( const long int xBR ) const
  { return calcNode.Ph_xDB_to_xCH( xBR ); }

#ifndef IPMGEMPLUGIN
// These calls are used only inside of GEMS-PSI GEM2MT module

   /// Constructor for integration in GEM2MT module of GEMS-PSI
   TNodeArray( long int nNodes, MULTI *apm );

   /// Constructor that uses 3D node arrangement
   TNodeArray( long int asizeN, long int asizeM, long int asizeK, MULTI *apm );

  /// Prints MULTI, DATACH and DATABR files structure prepared from GEMS.
  /// Prints files for separate coupled FMT-GEM programs that use GEMS3K module
  /// or if putNodT1 == true  as a break point for the running FMT calculation
  /// \param nIV - Number of allocated nodes
  /// \param bin_mode - Write IPM, DCH and DBR files in binary mode ( false - txt mode)
  /// \param brief_mode - Do not write data items that contain only default values
  /// \param with_comments -Write files with comments for all data entries ( in text mode)
  /// \param addMui - Print internal indices in RMULTS to IPM file for reading into Gems back
  gstring PutGEM2MTFiles(  QWidget* par, long int nIV,
       bool bin_mode = false, bool brief_mode = false, bool with_comments = true,
       bool putNodT1=false, bool addMui=false );

  /// Generate MULTI, DATACH and DATABR files structure prepared from GEMS.
  /// Prints files for separate coupled FMT-GEM programs that use GEMS3K module
  /// or if putNodT1 == true  as a break point for the running FMT calculation
  /// \param filepath - IPM work structure file path&name
  /// \param message - callback message function
  /// \param nIV - Number of allocated nodes
  /// \param bin_mode - Write IPM, DCH and DBR files in binary mode ( false - txt mode)
  /// \param brief_mode - Do not write data items that contain only default values
  /// \param with_comments -Write files with comments for all data entries ( in text mode)
  /// \param addMui - Print internal indices in RMULTS to IPM file for reading into Gems back
  /// \return DATABR list file name
  gstring genGEMS3KInputFiles(  const gstring& filepath, ProcessProgressFunction message,
          long int nIV, bool bin_mode, bool brief_mode, bool with_comments,
          bool putNodT1, bool addMui );

   /// Reads DATABR files saved by GEMS as a break point of the FMT calculation.
   /// Copying data from work DATABR structure into the node array NodT0
   /// and read DATABR structure into the node array NodT1 from file dbr_file
   void  setNodeArray( gstring& dbr_file, long int ndx, bool binary_f );

   /// Overloaded variant - takes lists of ICs, DCs and phases according to
   /// already existing index vectors axIC, axDC, axPH (with anICb, anDCb,
   /// anPHb, respectively)
   void InitCalcNodeStructures(  long int anICb, long int anDCb,  long int anPHb,
               long int* axIC, long int* axDC,  long int* axPH, bool no_interpolat,
            double* Tai, double* Pai,  long int nTp_,
            long int nPp_, double Ttol_, double Ptol_  )
   {
       calcNode.MakeNodeStructures(  anICb, anDCb,  anPHb,
                axIC, axDC,  axPH, no_interpolat,
                Tai,  Pai,  nTp_,  nPp_, Ttol_,  Ptol_  ) ;
   }

   /// Makes start DATACH and DATABR data using GEMS internal data (MULTI and other)
   /// interaction variant. The user must select ICs, DCs and phases to be included
   /// in DATABR lists
   /// Lookup arays from iterators
   void MakeNodeStructuresOne( QWidget* par, bool select_all,
       double Tai[4], double Pai[4]  )
   {
     calcNode.MakeNodeStructures( par, select_all,  Tai, Pai  );
     // setup dataBR and NodeT0 data
     NodT0[0] = allocNewDBR( calcNode);
     NodT1[0] = allocNewDBR( calcNode);
     MoveWorkNodeToArray( calcNode, 0, 1, na->pNodT0() );
     MoveWorkNodeToArray( calcNode, 0, 1, na->pNodT1() );
   }

   /// Makes start DATACH and DATABR data using GEMS internal data (MULTI and other)
   /// interaction variant. The user must select ICs, DCs and phases to be included
   /// in DATABR lists
   void MakeNodeStructuresOne( QWidget* par, bool select_all,bool no_interpolat,
            double *Tai, double *Pai, long int nTp_ = 1 ,
            long int nPp_ = 1 , double Ttol_ = 1., double Ptol_ =1. )
   {
     calcNode.MakeNodeStructures( par, select_all, no_interpolat,
                Tai, Pai, nTp_, nPp_, Ttol_, Ptol_ );
     NodT0[0] = allocNewDBR( calcNode);
     NodT1[0] = allocNewDBR( calcNode);
     MoveWorkNodeToArray( calcNode, 0, 1, na->pNodT0() );
     MoveWorkNodeToArray( calcNode, 0, 1, na->pNodT1() );
   }

   const TNode& LinkToNode( long int ndx, long int nNodes, DATABRPTR* anyNodeArray )
   {
     CopyWorkNodeFromArray( calcNode, ndx, nNodes, anyNodeArray );
     return calcNode;
   }

   void SaveToNode( long int ndx, long int nNodes, DATABRPTR* anyNodeArray )
   {
       // Save databr
       calcNode.packDataBr();
       if( !NodT0[ndx] )
         NodT0[ndx] = allocNewDBR( calcNode);
       if( !NodT1[ndx] )
         NodT1[ndx] = allocNewDBR( calcNode);
       MoveWorkNodeToArray( calcNode, ndx, nNodes, anyNodeArray );
  }

   DATABR *reallocDBR( long int ndx, long int nNodes, DATABRPTR* anyNodeArray)
   {
       if( ndx < 0 || ndx>= nNodes )
         return nullptr;
       // free old memory
       if( anyNodeArray[ndx] )
       {
            anyNodeArray[ndx] = calcNode.databr_free( anyNodeArray[ndx] );
            delete[] anyNodeArray[ndx];
       }
      anyNodeArray[ndx] = allocNewDBR( calcNode );
      if( !NodT1[ndx] )
        NodT1[ndx] = allocNewDBR( calcNode);
      return anyNodeArray[ndx];
   }

#else
// Used in GEMIPM2 standalone module only
   TNodeArray( long int nNod );   ///< Constructors for 1D arrangement of nodes
   TNodeArray( long int asizeN, long int asizeM, long int asizeK );
   ///< Constructor that uses 3D node arrangement
#endif

   /// Makes one absolute node index from three spatial coordinate indexes
   inline long int iNode( long int indN, long int indM, long int indK ) const
     { return  (( indK * sizeM + indM  ) * sizeN + indN);  }

   /// Get i index along N (x axis) from the absolute index ndx
   inline long int indN( long int ndx ) const
    { return  (ndx % sizeN);  }

   /// Get j index along M (y axis) from the absolute index ndx
   inline long int indM( long int ndx ) const
    {
	   long int j = (ndx - ndx % sizeN);
         j /=  sizeN;
     return  (j % sizeM);
    }

    /// Get k index along K (z axis) from the absolute index ndx
   inline long int indK( long int ndx ) const
    {
	   long int k = ndx - ndx % sizeN;
          k /= sizeN;
          k = k - k % sizeM;
      return  k/sizeM;
    }

    ~TNodeArray();      /// Destructor

    long int nNodes() const  /// Get total number of nodes in the node array
    { return anNodes; }

    long int SizeN() const  /// Get number of nodes in N direction (along x coordinate)
    { return sizeN; }

    long int SizeM() const  /// Get number of nodes in M direction (along y coordinate)
    { return sizeM; }

    long int SizeK() const  /// Get number of nodes in K direction (along z coordinate)
    { return sizeK; }

    DATABRPTR* pNodT0() const /// Get pointer to array of nodes for the previous time point
    { return NodT0; }

    DATABRPTR* pNodT1() const /// Get pointer to array of nodes for the current time point
    { return NodT1; }

    bool* piaNode() const /// Get pointer to IA switches for nodes
    { return iaNode; }
    
    char* ptcNode() const /// Get pointer to boundary condition codes for nodes
    { return tcNode; }
    
    /// New Stuff--------------------------------------------------------------

    ///  Here we do a GEM calculation in box ii
    bool CalcIPM_One(  const TestModeGEMParam& modeParam, long int ii, FILE* diffile )
    {
        return CalcIPM_Node(  modeParam, calcNode, ii, pNodT0(), pNodT1(), piaNode(), diffile );
    }

    ///  Here we do a GEM calculation in boxes from  start_node to end_node
    bool CalcIPM_List( const TestModeGEMParam& modeParam, long int start_node, long int end_node, FILE* diffile );

    ///
    /// Initialization of GEM IPM3 data structures in coupled programs
    /// that use GEMS3K module. Also reads in the IPM, DCH and one or many DBR text input files.
    ///  \param ipmfiles_lst_name  pointer to a null-terminated C string with a path to a text file
    ///                      containing the list of names of  GEMS3K input files.
    ///                      Example: file "test.lst" with a content:    -t "dch.dat" "ipm.dat" "dbr-0.dat"
    ///                      (-t  tells that input files are in text format)
    ///  \param dbrfiles_lst_name  a pointer to a null-terminated C string with a path to a text file
    ///                      containing the list of comma-separated names of  DBR input files.
    ///                      Example: file "test-dbr.lst" with a content:    "dbr-0.dat" , "dbr-1.dat" , "dbr-2.dat"
    ///  \param nodeTypes    the initial node contents
    ///                      from DATABR files will be distributed among nodes in array according to the
    ///                      distribution index list nodeTypes
    ///  \param getNodT1     used only when reading multiple DBR files after the modeling
    ///                      task interruption  in GEM-Selektor
    ///  \return 0  if successful; 1 if input file(s) were not found or corrupt;
    ///                      -1 if internal memory allocation error occurred.
    long int  GEM_init( const char* ipmfiles_lst_name,
             const char* dbrfiles_lst_name, long int* nodeTypes, bool getNodT1);

    // end of new stuff -------------------------------------------------------

    /// Calls GEM IPM calculation for a node with absolute index ndx
    long int RunGEM( TNode& wrkNode,  long int  iNode, long int Mode, DATABRPTR* nodeArray );
    
    /// Calls GEM IPM calculation for a selected group of nodes of TNodeArray (that have nodeFlag = 1) 
    /// in a loop with an optional openmp parallelization
    void RunGEM( long int Mode, int nNodes, DATABRPTR* nodeArray, long int* nodeFlags, long int* retCodes );

    /// Calls GEM IPM for one node with three indexes (along x,y,z)
    long int  RunGEM( TNode& wrkNode, long int indN, long int indM, long int indK,
                      long int Mode, DATABRPTR* nodeArray  )
    { return RunGEM( wrkNode, iNode( indN, indM, indK ), Mode, nodeArray ); }
        // (both calls clean the work node DATABR structure)

    /// Initialization of TNodeArray data structures. Reads in the DBR text input files and
    /// copying data from work DATABR structure into the node array
    ///  \param dbrfiles_lst_name  pointer to a null-terminated C string with a path to a text file
    ///                      containing the list of names of  DBR input files.
    ///                      Example: file "test-dbr.lst" with a content:    "dbr-0.dat" , "dbr-1.dat" , "dbr-2.dat"
    ///  \param nodeTypes    the initial node contents from DATABR files will be distributed among nodes in array
    ///                      according to the distribution list nodeTypes
    ///  \param getNodT1     optional parameter used only when reading multiple DBR files after modeling
    ///                      task interruption  in GEM-Selektor
    void  InitNodeArray( const char *dbrfiles_lst_name, long int *nodeTypes, bool getNodT1, bool binary_f  );

   //---------------------------------------------------------
   // Methods for working with node arrays (access to data from DBR)
   /// Calculate phase (carrier) mass, kg  of single component phase
   double get_mPH( long int ia, long int nodex, long int PHx );
   /// Calculate phase volume, cm3  of single component phase
   double get_vPH( long int ia, long int nodex, long int PHx );
   /// Calculate bulk compositions  of single component phase
   double get_bPH( long int ia, long int nodex, long int PHx, long int IC );


   //---------------------------------------------------------
   // Methods for working with node arrays

    /// Copies a node from the node array arr_From to the same place in the
    /// node array arr_To. Previous contents of the ndx element in arr_To
    /// will be lost. Uses the work node structure which will be newly
    /// allocated and contain no data afterwards
    void CopyNodeFromTo( long int ndx, long int nNodes, DATABRPTR* arr_From,
         DATABRPTR* arr_To );

    //---------------------------------------------------------
    // Data collection for monitoring differences
    // formatted writing into text file that must be already open 
    //
    /// Prints difference increments in all nodes (cells) for step t (time point at)
    void logDiffsIC( FILE* diffile, long int t, double at, long int nx, long int every_t );

    /// Prints dissolved elemental molarities in all cells for time point t / at
    void logProfileAqIC( FILE* logfile, long int t, double at, long int nx, long int every_t );

    /// Prints total elemental amounts in all cells for time point t / at
    void logProfileTotIC( FILE* logfile, long int t, double at, long int nx, long int every_t );

    /// Prints amounts of phases in all cells for time point t / at
    void logProfilePhMol( FILE* logfile, TParticleArray* pa, long int t, double at, long int nx, long int every_t );
    
    /// Prints volumes of phases in all cells for time point t / at
    void logProfilePhVol( FILE* logfile, long int t, double at, long int nx, long int every_t );
    
    /// Prints dissolved species molarities in all cells for time point t / at
    void logProfileAqDC( FILE* logfile, long int t, double at, long int nx, long int every_t );

    //---------------------------------------------------------
    // Working with the node grid (mainly used in Random Walk algorithms)

    ///  Set grid coordinate array use predefined array aGrid
    /// or set up regular scale
     void SetGrid( double aSize[3], double (*aGrid)[3] = nullptr );

     /// Finds a node absolute index for the current
     /// point location. Uses grid coordinate array grid[].
     /// Performance-important functions to be used e.g. in particle tracking methods
     long int FindNodeFromLocation( LOCATION cxyz, long int old_node = -1 ) const;

     /// Get 3D sizes for node (  from cxyz[0] - to cxyz[1] )
     void GetNodeSizes( long int ndx, LOCATION cxyz[2] );

     /// Get 3D location for node (  from cxyz[0] - to cxyz[1] )
     LOCATION& GetNodeLocation( long int ndx )
     { return grid[ndx]; }

     /// Get 3D size of the whole region
     LOCATION& GetSize()
     { return size; }

     /// Get full mass particle type in the node ndx
     double GetNodeMass( long int ndx, char type, char tcode, unsigned char ips );

     /// Move a mass m_v from node ndx_from to node ind_to, for particle type
     void MoveParticleMass( long int ndx_from, long int ind_to, char type, char ComponentMode, 
    		 char tcode, unsigned char ips, double m_v );

     /// Writes work node (DATABR structure) to a text VTK file
     void databr_to_vtk( fstream& ff, const char*name, double time, long cycle,
                               long int nFields=0, long int (*Flds)[2]=nullptr);

};

//Node scalar data access macros (added in devPhase), work both as get() and put()

// Temperature T (Kelvin) in T0 node
#define node0_TK( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->TK)
// Temperature T (Kelvin) in T1 node
#define node1_TK( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->TK)

// Pressure P (Pa) in T0 node
#define node0_P( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->P)
// Pressure P (Pa) in T1 node
#define node1_P( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->P)

// Volume V of reactive subsystem (m3) in T0 node
#define node0_Vs( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Vs)
// Volume V of reactive subsystem (m3) in T1 node
#define node1_Vs( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Vs)

// Volume of inert subsystem (m3) in T0 node
#define node0_Vi( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Vi)
// Volume of inert subsystem (m3) in T1 node
#define node1_Vi( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Vi)

// Mass of reactive subsystem (kg) in T0 node
#define node0_Ms( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Ms)
// Mass of reactive subsystem (kg) in T1 node
#define node1_Ms( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Ms)

// Mass of inert subsystem (kg) in T0 node
#define node0_Mi( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Mi)
// Mass of inert subsystem (kg) in T1 node
#define node1_Mi( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Mi)

/*
Gs,     ///< Total Gibbs energy of the reactive subsystem (J/RT) (norm)  -      -      +     +
Hs, 	///< Total enthalpy of reactive subsystem (J) (reserved)         -      -      +     +
Hi,     ///< Total enthalpy of inert subsystem (J) (reserved)            +      -      -     +
*/

// Effective aqueous ionic strength (molal) in T0 node
#define node0_IC( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->IC)
// Effective aqueous ionic strength (molal) in T1 node
#define node1_IC( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->IC)

// pH of aqueous solution in the activity scale (-log10 molal) in T0 node
#define node0_pH( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->pH)
// pH of aqueous solution in the activity scale (-log10 molal) in T1 node
#define node1_pH( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->pH)

// pe of aqueous solution in the activity scale (-log10 molal) in T0 node
#define node0_pe( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->pe)
// pe of aqueous solution in the activity scale (-log10 molal) in T1 node
#define node1_pe( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->pe)

// Eh of aqueous solution (V) in T0 node
#define node0_Eh( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Eh)
// Eh of aqueous solution (V) in T1 node
#define node1_Eh( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Eh)

// Actual total simulation time (s) in T0 node
#define node0_Tm( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Tm)
// Actual total simulation time (s) in T1 node
#define node1_Tm( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Tm)

// Actual time step (s) in T0 node (needed for TKinMet, can change in GEM calculation)
#define node0_dt( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->dt)
// Actual time step (s) in T1 node (needed for TKinMet, can change in GEM calculation)
#define node1_dt( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->dt)

#ifdef NODEARRAYLEVEL
// \section  FMT variables (units or dimensionsless) - to be used for storing them
//  at the nodearray level, normally not used in the single-node FMT-GEM coupling

// General diffusivity of disolved matter (m2/s) in T0 node
#define node0_Dif( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Dif)
// General diffusivity of disolved matter (m2/s) in T1 node
#define node1_Dif( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Dif)

// Total volume of the node (m3) (Vs + Vi) in T0 node
#define node0_Vt( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Vt)
// Total volume of the node (m3) (Vs + Vi) in T1 node
#define node1_Vt( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Vt)

// Advection velocity (in pores)  (m/s) in T0 node
#define node0_vp( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->vp)
// Advection velocity (in pores)  (m/s) in T1 node
#define node1_vp( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->vp)

// Effective (actual) porosity normalized to 1 in T0 node
#define node0_eps( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->eps)
// Effective (actual) porosity normalized to 1 in T1 node
#define node1_eps( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->eps)

// Actual permeability (m2) in T0 node
#define node0_Km( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Km)
// Actual permeability (m2) in T1 node
#define node1_Km( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Km)

// Actual Darcy`s constant (m2/s) in T0 node
#define node0_Kf( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Kf)
// Actual Darcy`s constant (m2/s) in T1 node
#define node1_Kf( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Kf)

// Specific storage coefficient, dimensionless in T0 node
#define node0_S( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->S)
// Specific storage coefficient, dimensionless in T1 node
#define node1_S( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->S)

// Transmissivity, m2/s in T0 node
#define node0_Tr( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->Tr)
// Transmissivity, m2/s in T1 node
#define node1_Tr( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->Tr)

// Actual hydraulic head (hydraulic potential) (m) in T0 node
#define node0_h( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->h)
// Actual hydraulic head (hydraulic potential) (m) in T1 node
#define node1_h( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->h)

// Actual carrier density for density-driven flow (kg/m3) in T0 node
#define node0_rho( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->rho)
// Actual carrier density for density-driven flow (kg/m3) in T1 node
#define node1_rho( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->rho)

// Specific longitudinal dispersivity of porous media (m) in T0 node
#define node0_al( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->al)
// Specific longitudinal dispersivity of porous media (m) in T1 node
#define node1_al( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->al)

// Specific transversal dispersivity of porous media (m) in T0 node
#define node0_at( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->at)
// Specific transversal dispersivity of porous media (m) in T1 node
#define node1_at( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->at)

// Specific vertical dispersivity of porous media (m) in T0 node
#define node0_av( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->av)
// Specific vertical dispersivity of porous media (m) in T1 node
#define node1_av( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->av)

// Hydraulic longitudinal dispersivity (m2/s) in T0 node
#define node0_hDl( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->hDl)
// Hydraulic longitudinal dispersivity (m2/s) in T1 node
#define node1_hDl( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->hDl)

// Hydraulic transversal dispersivity (m2/s) in T0 node
#define node0_hDt( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->hDt)
// Hydraulic transversal dispersivity (m2/s) in T1 node
#define node1_hDt( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->hDt)

// Hydraulic vertical dispersivity (m2/s) in T0 node
#define node0_hDv( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->hDv)
// Hydraulic vertical dispersivity (m2/s) in T1 node
#define node1_hDv( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->hDv)

// Tortuosity factor (dimensionless) in T0 node
#define node0_nto( nodex ) (TNodeArray::na->pNodT0()[(nodex)]->nto)
// Tortuosity factor (dimensionless) in T1 node
#define node1_nto( nodex ) (TNodeArray::na->pNodT1()[(nodex)]->nto)

#endif

//IC node data access macros

#define node0_bIC( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->bIC[(ICx)])
#define node1_bIC( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->bIC[(ICx)])
#define node0_rMB( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->rMB[(ICx)])
#define node1_rMB( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->rMB[(ICx)])
#define node0_uIC( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->uIC[(ICx)])
#define node1_uIC( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->uIC[(ICx)])
// equilibrium bulk composition of solid part of the system, moles from T0 node
#define node0_bSP( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->bSP[(ICx)])
// equilibrium bulk composition of solid part of the system, moles from T1 node
#define node1_bSP( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->bSP[(ICx)])

//DC node data access macros

  // amount of DC with index DCx from T0 node with index nodex
#define node0_xDC( nodex, DCx ) (TNodeArray::na->pNodT0()[(nodex)]->xDC[(DCx)])
  // amount of DC with index DCx from T1 node with index nodex
#define node1_xDC( nodex, DCx ) (TNodeArray::na->pNodT1()[(nodex)]->xDC[(DCx)])

  // activity coefficient of DC with index DCx from T0 node with index nodex
#define node0_gam( nodex, DCx ) (TNodeArray::na->pNodT0()[(nodex)]->gam[(DCx)])
  // activity coefficient of DC with index DCx from T1 node with index nodex
#define node1_gam( nodex, DCx ) (TNodeArray::na->pNodT1()[(nodex)]->gam[(DCx)])

  // upper constraint on amount of DC with index DCx from T0 node with index nodex
#define node0_dul( nodex, DCx ) (TNodeArray::na->pNodT0()[(nodex)]->dul[(DCx)])
  // upper constraint on amount of DC with index DCx from T1 node with index nodex
#define node1_dul( nodex, DCx ) (TNodeArray::na->pNodT1()[(nodex)]->dul[(DCx)])

  // lower constraint on amount of DC with index DCx from T0 node with index nodex
#define node0_dll( nodex, DCx ) (TNodeArray::na->pNodT0()[(nodex)]->dll[(DCx)])
  // lower constraint on amount of DC with index DCx from T1 node with index nodex
#define node1_dll( nodex, DCx ) (TNodeArray::na->pNodT1()[(nodex)]->dll[(DCx)])

//Phase node data access macros
  // amount of phase with index PHx from T0 node with index nodex
#define node0_xPH( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->xPH[(PHx)])
  // amount of phase with index PHx from T1 node with index nodex
#define node1_xPH( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->xPH[(PHx)])

// amount of phase with index PHx from T0 node with index nodex
#define node0_omPH( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->omPH[(PHx)])
// amount of phase with index PHx from T1 node with index nodex
#define node1_omPH( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->omPH[(PHx)])

  // volume of multicomponent phase with index PHx from T0 node with index nodex
#define node0_vPS( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->vPS[(PHx)])
  // volume of multicomponent phase with index PHx from T1 node with index nodex
#define node1_vPS( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->vPS[(PHx)])

  // volume of single-component phase with index PHx from T0 node with index nodex
#define node0_vPH( nodex, PHx ) (TNodeArray::na->get_vPH( 0, (nodex), (PHx)))
  // volume of single-component phase with index PHx from T1 node with index nodex
#define node1_vPH( nodex, PHx ) (TNodeArray::na->get_vPH( 1, (nodex), (PHx)))

  // mass of multicomponent phase with index PHx from T0 node with index nodex
#define node0_mPS( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->mPS[(PHx)])
  // mass of multicomponent phase with index PHx from T1 node with index nodex
#define node1_mPS( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->mPS[(PHx)])

  // mass of single-component phase with index PHx from T0 node with index nodex
#define node0_mPH( nodex, PHx )  (TNodeArray::na->get_mPH( 0, (nodex), (PHx)))
  // mass of single-component phase with index PHx from T1 node with index nodex
#define node1_mPH( nodex, PHx )  (TNodeArray::na->get_mPH( 1, (nodex), (PHx)))

// Upper AMRs on amounts of multi-component phases (mol) with index PHx from T0 node with index nodex
#define node0_amru( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->amru[(PHx)])
// Upper AMRs on amounts of multi-component phases (mol) with index PHx from T1 node with index nodex
#define node1_amru( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->amru[(PHx)])

// Lower AMRs on amounts of multi-component phases (mol) with index PHx from T0 node with index nodex
#define node0_amrl( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->amrl[(PHx)])
// Lower AMRs on amounts of multi-component phases (mol) with index PHx from T1 node with index nodex
#define node1_amrl( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->amrl[(PHx)])

  // amount of solvent/sorbent in phase with index PHx from T0 node with index nodex
#define node0_xPA( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->xPA[(PHx)])
  // amount of solvent/sorbent in phase with index PHx from T1 node with index nodex
#define node1_xPA( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->xPA[(PHx)])   

// Phase compositions node data access macros
// amount of independent component ICx in multi-component phase PHx in T0 node nodex
#define node0_bPS( nodex, PHx, ICx ) ( TNodeArray::na->pNodT0()[(nodex)]->bPS[ \
                                       (PHx)*TNodeArray::na->pCSD()->nICb+(ICx)])
// amount of independent component ICx in multi-component phase PHx in T1 node nodex
#define node1_bPS( nodex, PHx, ICx ) ( TNodeArray::na->pNodT1()[(nodex)]->bPS[ \
                                       (PHx)*TNodeArray::na->pCSD()->nICb+(ICx)])

// amount of independent component ICx in single-component phase PHx in T0 node nodex
#define node0_bPH( nodex, PHx, ICx )  (TNodeArray::na->get_bPH( 0, (nodex), (PHx), (ICx)))
// amount of independent component ICx in single-component phase PHx in T1 node nodex
#define node1_bPH( nodex, PHx, ICx )  (TNodeArray::na->get_bPH( 1, (nodex), (PHx), (ICx)))

#endif   // _nodearray_h_

// end nodearray.h
