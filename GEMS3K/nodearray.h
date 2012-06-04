//-------------------------------------------------------------------
// $Id$
//
/// \class TNodeArray nodearray.h
/// Implements an advanced (level 2) C/C++ interface with
/// GEM IPM for development of coupled and GEM-FMT codes
/// Works with DATACH and work DATABR structures
/// Uses the TNode class
//
// Written by S.Dmytriyeva,  D.Kulik
/// \author Copyright (C) 2006,2008 S.Dmytriyeva, D.Kulik
//
// This file is part of GEMS3K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization

// This file may be distributed under the licence terms
// defined in GEMS3K.QAL
//
// See also http://gems.web.psi.ch/
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#ifndef _nodearray_h_
#define _nodearray_h_

#include "node.h"

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

  LOCATION( float ax, float ay=0., float az=0. ):
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

// Definition of TNodeArray class
class TNodeArray : public TNode
{
protected:

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

public:

  static TNodeArray* na;   ///< static pointer to this class

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

   /// Reads DATABR files saved by GEMS as a break point of the FMT calculation.
   /// Copying data from work DATABR structure into the node array NodT0
   /// and read DATABR structure into the node array NodT1 from file dbr_file
   void  setNodeArray( gstring& dbr_file, long int ndx, bool binary_f );

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
    
    /// Calls GEM IPM calculation for a node with absolute index ndx
    long int RunGEM( long int ndx, long int Mode );

    /// Calls GEM IPM for one node with three indexes (along x,y,z)
    long int  RunGEM( long int indN, long int indM, long int indK, long int Mode )
    { return RunGEM( iNode( indN, indM, indK ), Mode); }
        // (both calls clean the work node DATABR structure)

    /// Copies data from the work DATABR structure into the node ndx in
    /// the node arrays NodT0 and NodT1  (as specified in nodeTypes array)
    void  setNodeArray( long int ndx, long int* nodeTypes  );

   /// Test setup of the boundary condition for all nodes in the task
    void  checkNodeArray( long int i, long int* nodeTypes, const char*  datachbr_file );

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

    ///  Copies data for a node ndx from the array of nodes anyNodeArray that
    /// contains nNodes into the work node data bridge structure
    void CopyWorkNodeFromArray( long int ndx, long int nNodes, DATABRPTR* anyNodeArray );

    ///  Moves work node data to the ndx element of the node array anyNodeArray
    /// that has nNodes. Previous contents of the ndx element will be lost,
    /// work node will be allocated new and will contain no data
    void MoveWorkNodeToArray( long int ndx, long int nNodes, DATABRPTR* anyNodeArray );

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
    void logProfilePhMol( FILE* logfile, long int t, double at, long int nx, long int every_t );
    
    /// Prints volumes of phases in all cells for time point t / at
    void logProfilePhVol( FILE* logfile, long int t, double at, long int nx, long int every_t );
    
    /// Prints dissolved species molarities in all cells for time point t / at
    void logProfileAqDC( FILE* logfile, long int t, double at, long int nx, long int every_t );

    //---------------------------------------------------------
    // Working with the node grid (mainly used in Random Walk algorithms)

    ///  Set grid coordinate array use predefined array aGrid
    /// or set up regular scale
     void SetGrid( double aSize[3], double (*aGrid)[3] = 0 );

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
                               long int nFields=0, long int (*Flds)[2]=0);

};

//IC node data access macroses
#define node0_bIC( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->bIC[(ICx)])
#define node1_bIC( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->bIC[(ICx)])
#define node0_rMB( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->rMB[(ICx)])
#define node1_rMB( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->rMB[(ICx)])
#define node0_uIC( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->uIC[(ICx)])
#define node1_uIC( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->uIC[(ICx)])

//DC node data access macroses

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

//Phase node data access macroses
  // amount of phase with index PHx from T0 node with index nodex
#define node0_xPH( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->xPH[(PHx)])
  // amount of phase with index PHx from T1 node with index nodex
#define node1_xPH( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->xPH[(PHx)])

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


  // amount of solvent/sorbent in phase with index PHx from T0 node with index nodex
#define node0_xPA( nodex, PHx ) (TNodeArray::na->pNodT0()[(nodex)]->xPA[(PHx)])
  // amount of solvent/sorbent in phase with index PHx from T1 node with index nodex
#define node1_xPA( nodex, PHx ) (TNodeArray::na->pNodT1()[(nodex)]->xPA[(PHx)])

// equilibrium bulk composition of solid part of the system, moles from T0 node
#define node0_bSP( nodex, ICx ) (TNodeArray::na->pNodT0()[(nodex)]->bSP[(ICx)])
// equilibrium bulk composition of solid part of the system, moles from T1 node
#define node1_bSP( nodex, ICx ) (TNodeArray::na->pNodT1()[(nodex)]->bSP[(ICx)])


// Phase compositions node data access macroses
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
