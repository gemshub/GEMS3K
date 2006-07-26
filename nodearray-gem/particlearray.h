//-------------------------------------------------------------------
// $Id: particlearray.h 684 2005-11-23 11:19:27Z gems $
//
// C/C++ interface for particle tracking methods for FMT node array
// Working whith DATACH and DATABR structures
//
// (C) 2006 S.Dmytriyeva, D.Kulik, W.Pfingsten
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _particlearray_h_
#define _particlearray_h_

#include "nodearray.h"


enum  PMCODE // Codes of particle movement
{
  IMMOBILE_C_VOLUME = 0,
  MOBILE_C_VOLUME = 1,
  MOBILE_VIRTUAL = 2,
  IMMOBILE_C_MASS = 10,
  MOBILE_C_MASS = 11

};

// This structure describes a single transport particle
typedef struct
{
   char ptype;    // particle type index ( 1 to 255 )
   char mmode;    // particle movement mode (see enum PMCODE)
   char tcode;    // particle transport mechanism code (see enum PTCODE)
   unsigned char ips;  // DataBr index of phase or species to which this particle is connected

   LOCATION xyz;       // current particle location

   double m_v;    // mass or volume of the particle (depending on ptype and mmode)

   int node;      // (absolute) index of the parent node

   unsigned int res; // reserved

}  // 32 bytes
PARTICLE;


// This class keeps the whole array of particles and is linked to the TNodeArray class
class TParticleArray
{
    PARTICLE* ParT0;  // array of particles for current time point
    PARTICLE* ParT1;  // array of particles for previous time point

    int anParts;      // Number of allocated particles (in each array for T0 and T1) <= nNodes*nPmax
    int anPTypes;     // Number of allocated particle types (< 20 ? )
    int anProps;      // Number of particle statistic properties (for monitoring) >= anPTypes

    TNodeArray* nodes; // Pointer to TNodeArray class

// Data structures for particle statistics setup

    // data for particle array setup
    PARTICLE *ParTD;  // array of particle type definitions at t0 or after interruption anPTypes
    short *NPmean;       // array of initial mean particle type numbers per node ( size: anPTypes )

    short* nPmin;        // minimum average total number of particles of each type per one node anPTypes
    short* nPmax;        // maximum average total number of particles of each type per one node anPTypes

    int nNodes;       // Total number of allocated nodes (for convenience, to be copied from nodearray class)
    short* NPnum;       // array of particle numbers in nodes ( size: nNodes * (anPTypes))

   // internal setup variables
    int  ndxCsource; // Cauchy source index

    // Current step data for particle statistics
    int cpTypes;      // current number of particle types >= 1 and <= anPTypes;

    int* cParts;       // current total number of particles per type (>= nPmin*nNodes and <= anParts;)

    int cpx;          // current particle index            ??
    int cptx;         // int current particle type index   ??

    int* (*NPlist);   // list of particle indexes coming into each node (at T1) size: nNodes * nPmax
    float* (*NPstat); // array of particle statistic properties in nodes (nNodes * anProps)
                      // for monitoring transport

//  Time variables
    double tau_max;   // maximum physical time (sec), origin 0

    double dt_min;    // minimum time step length
    double dt_max;    // maximum time step length

    double dt;        // current time step (dt_min <= dt <= dt_max)
    double dt_prev;   // previous time step

    double t0;        // current time point
    double t1;        // next time point, t1 = t0 + dt

    int tstep;        // current time step (for t0 time point)

// More data fields for interpolation, smoothing, etc. to add




// Methods

// internal methods

   void freeMemory();
   LOCATION setPointInNode( LOCATION nodeSize[2] );
   double InterpolationVp_hDl_1D( int px,double& vp, double& al, double& Dif );

  // Important for masstransport step
  // Calculation of new particle locations
  int DisplaceParticle( int px, double t0, double t1 );
  // Walk (transport step) for particle px between nodes
  int MoveParticleBetweenNodes( int px, double t0, double t1 );
  // call to the whole Random Walk method time step over all particles and nodes
 // returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
 // GEM was called before this function
  int RandomWalkIteration( int Mode, double t0, double t1 );
  // call to the whole FiniteCell Walk method time step over all particles and nodes
  // returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
  int FCellWalkIteration( int Mode, double t0, double t1 );

// not use functions
  void particles_to_text_file( fstream& ff );    // writes particle array(s) to a text file
  void particles_from_text_file( fstream& ff);   // reads particle array(s) from a text file
  void PGlists_to_text_file(fstream& ff );     // writes work node (DATABR structure) to text file
  void PGlists_from_text_file(fstream& ff );   // reads work node (DATABR structure) from text file


public:

   TParticleArray( short nPTypes, short nProps,
           short *aNPmean, short (*aParTD)[6],
           short *anPmin, short *anPmax, TNodeArray* aNodes );
// destructor
  ~TParticleArray();

   int nParts() const { return anParts; }  // returns total number of particles (in each time row)

   PARTICLE* pPART0( int px ) const  // get pointer to PARTICLE data structure (T0) with index Indx
      { return &ParT0[ px ]; }

   PARTICLE* pPART1( int px ) const  // get pointer to PARTICLE data structure (T1) with index Indx
      { return &ParT1[ px ]; }

   LOCATION* GetPartLocationT0( int px ) const  // returns pointer to the particle px location at time T0
    {  return &(ParT0[px].xyz); }

   LOCATION* GetPartLocationT1( int px ) const  // returns pointer to the particle px location at time T1
    {   return &(ParT1[px].xyz); }

   void SetPartLocationT0( int px, LOCATION& cxyz )  // sets particle px location at time T0
   {   ParT0[px].xyz = cxyz;  }

   void SetPartLocationT1( int px, LOCATION& cxyz )  // sets particle px location at time T1
  {    ParT1[px].xyz = cxyz;  }

  void  setUpCounters();
  void CopyfromT1toT0();  // Copy resalts of ParT1 step to ParT0
  void ParticleArrayInit();  // Particle array initialization
   // stub call for coupled mass transport calculation
  int GEMCOTAC( int Mode, double t0, double t1 );

//#ifdef IPMGEMPLUGIN

  static TParticleArray* pa;

//#endif

  int nPTypes() const
   { return anPTypes; }    // Number of allocated particle types (< 20 ? )
  short getNPnum( int iNode, int iType ) const // particle numbers in node
   { return NPnum[iType + iNode*anPTypes];   }
};

#endif   // _particlearray_h_

