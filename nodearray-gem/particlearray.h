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
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _particlearray_h_
#define _particlearray_h_

#include "nodearray.h"

// Random numbers ==========================================================
double randuni(double& x); // uniform
double randnorm(double& x); // normal
double ran2(long int & idum);  // uniform between 1 and 0
double ran3(long int & idum);  // uniform between 1 and 0

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

   long int node;      // (absolute) index of the parent node

   unsigned long int res; // reserved

}  // 32 bytes
PARTICLE;


// This class keeps the whole array of particles and is linked to the TNodeArray class
class TParticleArray
{
    PARTICLE* ParT0;  // array of particles for current time point
    PARTICLE* ParT1;  // array of particles for previous time point

    long int anParts;      // Number of allocated particles (in each array for T0 and T1) <= nNodes*nPmax
    long int anPTypes;     // Number of allocated particle types (< 20 ? )
    long int anProps;      // Number of particle statistic properties (for monitoring) >= anPTypes

    TNodeArray* nodes; // Pointer to TNodeArray class

// Data structures for particle statistics setup

    // data for particle array setup
    PARTICLE *ParTD;  // array of particle type definitions at t0 or after interruption anPTypes
    long int* NPmean;       // array of initial mean particle type numbers per node ( size: anPTypes )

    long int* nPmin;        // minimum average total number of particles of each type per one node anPTypes
    long int* nPmax;        // maximum average total number of particles of each type per one node anPTypes

    long int nNodes;       // Total number of allocated nodes (for convenience, to be copied from nodearray class)
    long int* NPnum;       // array of particle numbers in nodes ( size: nNodes * (anPTypes))

   // internal setup variables
    long int  ndxCsource; // Cauchy source index

    // Current step data for particle statistics
    long int cpTypes;      // current number of particle types >= 1 and <= anPTypes;

    long int* cParts;       // current total number of particles per type (>= nPmin*nNodes and <= anParts;)

    long int cpx;          // current particle index            ??
    long int cptx;         // int current particle type index   ??

    long int* (*NPlist);   // list of particle indexes coming into each node (at T1) size: nNodes * nPmax
    double* (*NPstat); // array of particle statistic properties in nodes (nNodes * anProps)
                      // for monitoring transport

//  Time variables
    double tau_max;   // maximum physical time (sec), origin 0

    double dt_min;    // minimum time step length
    double dt_max;    // maximum time step length

    double dt;        // current time step (dt_min <= dt <= dt_max)
    double dt_prev;   // previous time step

    double t0;        // current time point
    double t1;        // next time point, t1 = t0 + dt

    long int tstep;        // current time step (for t0 time point)

// More data fields for interpolation, smoothing, etc. to add


// Methods

// internal methods

   void freeMemory();
   LOCATION setPointInNode( LOCATION nodeSize[2] );
   double InterpolationVp_hDl_1D( long int px,double& vp, double& al, double& Dif, double& Dpm );

  // Important for masstransport step
  // Calculation of new particle locations
  long int DisplaceParticle( long int px, double t0, double t1 );
  // Walk (transport step) for particle px between nodes
  long int MoveParticleBetweenNodes( long int px, bool CompMode, double t0, double t1 );
  // call to the whole Random Walk method time step over all particles and nodes
 // returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
 // GEM was called before this function
  long int RandomWalkIteration( long int Mode, bool CompMode, double t0, double t1 );
  // call to the whole FiniteCell Walk method time step over all particles and nodes
  // returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
  long int FCellWalkIteration( long int Mode, bool CompMode, double t0, double t1 );

// not use functions
  void particles_to_text_file( std::fstream& ff );    // writes particle array(s) to a text file
  void particles_from_text_file( std::fstream& ff);   // reads particle array(s) from a text file
  void PGlists_to_text_file(std::fstream& ff );     // writes work node (DATABR structure) to text file
  void PGlists_from_text_file(std::fstream& ff );   // reads work node (DATABR structure) from text file


public:

   TParticleArray( long int nPTypes, long int nProps,
           long int *aNPmean, long int (*aParTD)[6],
           long int *anPmin, long int *anPmax, TNodeArray* aNodes );
// destructor
  ~TParticleArray();

   long int nParts() const { return anParts; }  // returns total number of particles (in each time row)

   PARTICLE* pPART0( long int px ) const  // get pointer to PARTICLE data structure (T0) with index Indx
      { return &ParT0[ px ]; }

   PARTICLE* pPART1( long int px ) const  // get pointer to PARTICLE data structure (T1) with index Indx
      { return &ParT1[ px ]; }

   LOCATION* GetPartLocationT0( long int px ) const  // returns pointer to the particle px location at time T0
    {  return &(ParT0[px].xyz); }

   LOCATION* GetPartLocationT1( long int px ) const  // returns pointer to the particle px location at time T1
    {   return &(ParT1[px].xyz); }

   void SetPartLocationT0( long int px, LOCATION& cxyz )  // sets particle px location at time T0
   {   ParT0[px].xyz = cxyz;  }

   void SetPartLocationT1( long int px, LOCATION& cxyz )  // sets particle px location at time T1
  {    ParT1[px].xyz = cxyz;  }

  void  setUpCounters();
  void CopyfromT1toT0();  // Copy results of ParT1 step to ParT0
  void ParticleArrayInit();  // Particle array initialization
   // stub call for coupled mass transport calculation
  long int GEMPARTRACK( long int Mode, bool ComponentMode, double t0, double t1 );

  void logProfilePhMol( FILE* logfile, int inode );

  long int nPTypes() const
   { return anPTypes; }    // Number of allocated particle types (< 20 ? )
  long getNPnum( long int iNode, long int iType ) const // particle numbers in node
   { return NPnum[iType + iNode*anPTypes];   }
};

#endif   // _particlearray_h_

