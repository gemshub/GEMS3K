//-------------------------------------------------------------------
// $Id: particlearray.h 684 2005-11-23 11:19:27Z gems $
//
// C/C++ interface for particle tracking methods for FMT node array
// Working whith DATACH and DATABR structures
//
// (C) 2006, 2019  S.Dmytriyeva, D.Kulik, W.Pfingsten
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
//
// This file may be distributed under the terms of GEMS3 Development
// Quality Assurance Licence (GEMS3.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include <cmath>
#include "particlearray.h"

static long int idum_static = -10000l;
static double Rand = -1;

TParticleArray::TParticleArray( long int nPTypes, long int nProps,
           long int *aNPmean,
           long int (*aParTD)[6],
           long int *anPmin, long int *anPmax,
           TNodeArray* aNodes ):
  ParT0(0), ParT1(0), anParts(0),
  anPTypes(nPTypes), anProps(nProps), nodes(aNodes)
{
  long int ii;
  nNodes = nodes->nNodes();

  NPmean = new long int[anPTypes];
  ParTD = new PARTICLE[anPTypes];
  nPmin = new long int[anPTypes];
  nPmax = new long int[anPTypes];

  for( ii=0; ii < anPTypes; ii++ )
  {
    NPmean[ii] = aNPmean[ii];
    nPmin[ii] = anPmin[ii];
    nPmax[ii] = anPmax[ii];

    ParTD[ii].ptype = (char)aParTD[ii][0];
    ParTD[ii].mmode = (char)aParTD[ii][1];
    ParTD[ii].tcode = (char)aParTD[ii][2];
    ParTD[ii].ips = (char)aParTD[ii][3];
  }

  anParts = 0;
  for( ii=0; ii < anPTypes; ii++ )
   anParts += NPmean[ii];
  anParts *= nNodes;

  ParT0 = new PARTICLE[anParts];
  ParT1 = new PARTICLE[anParts];
  NPnum = new long int[nNodes*anPTypes];
  ParticleArrayInit();

  cParts = 0;
  NPlist = 0;
  NPstat = 0;
}

void TParticleArray::freeMemory()
{
  delete[] NPmean;
  delete[] ParTD;
  delete[] nPmin;
  delete[] nPmax;
  delete[] ParT0;
  delete[] ParT1;
  delete[] NPnum;
}


TParticleArray::~TParticleArray()
{
   freeMemory();
}

// Particle array initialization
void TParticleArray::ParticleArrayInit()
{
    long int iNode, iType, k, cpx1;
  LOCATION nodeSize[2];

  cpx1 = 0;
  ndxCsource = 0;
  Rand = -1.;
  for( iNode=0; iNode < nNodes; iNode++ )
  {
    nodes->GetNodeSizes( iNode, nodeSize );
    for( iType=0; iType < anPTypes; iType++ )
    {
      NPnum[iNode*anPTypes+iType] = NPmean[iType];
//      double dd = (nodeSize[1].x-nodeSize[0].x)/NPnum[iNode*anPTypes+iType];
     for( k=0; k < NPmean[iType]; k++ )
      {
        ParT0[cpx1].ptype = iType;
        ParT0[cpx1].mmode = ParTD[iType].mmode;
        ParT0[cpx1].tcode = ParTD[iType].tcode;
        ParT0[cpx1].ips = ParTD[iType].ips;
        ParT0[cpx1].m_v = 0.;
        ParT0[cpx1].node = iNode;
//        ParT0[cpx1].xyz.x = nodeSize[0].x+k*dd;
        ParT0[cpx1].xyz = setPointInNode(nodeSize);
        ParT1[cpx1] = ParT0[cpx1];
        cpx1++;
      }
    }
  }
}

// Particle array initialization
void TParticleArray::setUpCounters()
{
	long int iNode;
  ndxCsource = 0;
  for( iNode=0; iNode < nNodes; iNode++ )
  {
    // set up Cauchy source index
    if( nodes->pNodT1()[iNode]->NodeTypeHY == NBC3source )
    { ndxCsource = iNode;
      break;
    }
  }
}

void TParticleArray::CopyfromT1toT0()  // Copy results of ParT1 step to ParT0
{
/*  fstream f_out("nods_particl.out", ios::out|ios::app  );
  for(long int iNode=0; iNode < nNodes; iNode++ )
  {
    for(long int iType=0; iType < anPTypes; iType++ )
       f_out << NPnum[iNode*anPTypes+iType] << " ";
    f_out << " || ";
  }
 f_out << endl;

*/
  for(long int k=0; k < anParts; k++ )
  {
    ParT0[k] = ParT1[k];
  }
}

// Returns a uniform random deviate between nodeSize[0] and nodeSize[1]
// for x, y, z coordinate
LOCATION TParticleArray::setPointInNode( LOCATION nodeSize[2] )
{
  LOCATION loc;
  if( nodes->SizeN() > 1 )
  {
    loc.x = randuni( Rand ); //  randnorm(Rand); // value from 0. to 1.
//    loc.x = ran3( idum ); //  ran2( idum ); // value from 0. to 1.
    loc.x = (loc.x*(nodeSize[1].x-nodeSize[0].x)+nodeSize[0].x);
  }
  else loc.x = 0;
  if( nodes->SizeM() > 1 )
  {
    loc.y = randnorm(Rand ); //  randuni( Rand); // value from 0. to 1.
//    loc.y = ran2( idum ); //  ran3( idum ); // value from 0. to 1.
    loc.y = (loc.y*(nodeSize[1].y-nodeSize[0].y)+nodeSize[0].y);
  }
  else loc.y = 0;
  if( nodes->SizeK() > 1 )
  {
    loc.z = randnorm(Rand ); //  randuni( Rand); // value from 0. to 1.
//    loc.z = ran2( idum ); //  ran3( idum ); // value from 0. to 1.
    loc.z = (loc.z*(nodeSize[1].z-nodeSize[0].z)+nodeSize[0].z);
  }
  else loc.z = 0;
  return loc;
}

// Implementation of interpolation for particle advection velocities and
// dispersivities between nodes ( in 1D case)
// px index of particle
double TParticleArray::InterpolationVp_hDl_1D( long int px,
   double& vp, double& al, double& Dif, double& Dpm )
{
  if( nodes->SizeM() > 1 ||  nodes->SizeK() > 1 )
     Error( "InterpolationVp_hDl_1D", "Error mode of interpolation." );

  DATABR *dbr1, *dbr2;    // nodes for interpolation
  long int nodInd1, nodInd2;  // number of nodes for interpolation
  double hDl, x1m, x2m, eps, nto; // middle-point coordinates in the nodes
  LOCATION nodeSize[2];

// set up location
  nodInd1 = ParT1[px].node;
  nodes->GetNodeSizes( nodInd1, nodeSize );
  x1m = nodeSize[0].x + (nodeSize[1].x-nodeSize[0].x)/2 ;
  if( ParT1[px].xyz.x <  x1m )
   nodInd2 = nodInd1-1;
  else
   nodInd2 = nodInd1+1;

// get new vp and hDl
  dbr1 = nodes->pNodT1()[nodInd1];  // nodes at current time point
  if( nodInd2 < 0 || nodInd2 >= nNodes )
  {
#ifdef NODEARRAYLEVEL
    vp = dbr1->vp;
    hDl = dbr1->hDl;
    al = dbr1->al;
    Dif = dbr1->Dif;
    eps = dbr1->eps;
    nto = dbr1->nto;
    Dpm = dbr1->Dif*dbr1->eps/dbr1->nto;
#endif
  }
  else
  {
    nodes->GetNodeSizes( nodInd2, nodeSize );
    x2m = nodeSize[0].x + (nodeSize[1].x-nodeSize[0].x)/2 ;
    dbr2 = nodes->pNodT1()[nodInd2];
    // vx = v1 - (v2-v1)/(x2-x1)*(x-x1);
    double d = (ParT1[px].xyz.x-x1m)/(x2m-x1m);
#ifdef NODEARRAYLEVEL
    vp = dbr1->vp;
    vp -= (dbr2->vp - dbr1->vp )*d;
    al = dbr1->al;
    al -= (dbr2->al - dbr1->al )*d;
    Dif = dbr1->Dif;
    Dif -= (dbr2->Dif - dbr1->Dif )*d;
    eps = dbr1->eps;
    eps -= (dbr2->eps - dbr1->eps )*d;
    nto = dbr1->nto;
    nto -= (dbr2->nto - dbr1->nto )*d;
    Dpm =  dbr1->Dif*dbr1->eps/dbr1->nto;    // fix?
    Dpm -= (dbr2->Dif*dbr2->eps/dbr2->nto - dbr1->Dif*dbr1->eps/dbr1->nto)*d;
//    Dpm = eps*Dif/nto;  // added account for tortuosity and porosity DK 9.05.19
//    hDl = al*vp+Dif;
    hDl = al*vp+Dpm;    // added DK 9.05.19
#endif
  }
  return hDl;
}

// Important for the mass-transport step
// Calculation of new particle locations
// Advective step, Brownian step, Dispersive step, Diffusion step
// px: index of particle
long int TParticleArray::DisplaceParticle( long int px, double /*t0*/, double /*t1*/ )
{
//  DATACH* ch = nodes->pCSD();       // DataCH structure
	long int nodInd = ParT1[px].node;
  double ds = 0.;
  double vp, hDl, al, Dif, Dpm;

  ErrorIf( nodInd < 0 , "DisplaceParticle", "Error index" );
//  DATABR* dbr = nodes->pNodT1()[nodInd];  // nodes at current time point

  switch( ParT1[px].mmode )
  {
   case IMMOBILE_C_VOLUME:
   case MOBILE_C_VOLUME:
   case MOBILE_VIRTUAL:
   case IMMOBILE_C_MASS:
   default:
                         break;

   case MOBILE_C_MASS:
        hDl = InterpolationVp_hDl_1D( px, vp, al, Dif, Dpm );
// vp = dbr->vp;     // testing without interpolation
// hDl = dbr->hDl;   // testing without interpolation
         if( hDl > 0)
//         ds = 2.*(randuni( idum )-0.5)*sqrt( 6.*hDl*vp*dt);
         ds = 2.*(ran3( idum_static )-0.5)* std::sqrt( 6.*hDl*dt);
         ParT1[px].xyz.x += vp*dt + ds;
                        break;
   }
  return 0;
}

// Walk (transport step) for particle px between nodes
// returns -1 if particle stays in the same node
// or an index of another node the particle enters
long int TParticleArray::MoveParticleBetweenNodes( long int px, bool CompMode, double /*t0*/, double /*t1*/ )
{
  long int old_node = ParT1[px].node;
  long int new_node = nodes->FindNodeFromLocation( ParT1[px].xyz, old_node );
  char type_ = ParT1[px].ptype;
  long int nodeType = nodes->pNodT1()[old_node]->NodeTypeHY;

  if( new_node == -1 )     // location behind region
  {
   //if( old_node != 0 && old_node != nodes->nNodes()-1 )
   // new_node = nodes->FindNodeFromLocation( ParT1[px].xyz  ); // all list
   // may be reflection (otrazhenie)
   // change xyz
  // may be discuss: sources, links, reflections . . .
  }

  if( old_node == new_node ) // particle left in old node
     return -1;

  // move mass of particle from old_node to new_node (if new_node==-1 mass may be lost)
  nodes->MoveParticleMass( old_node, new_node, type_, CompMode, 
            ParT1[px].tcode, ParT1[px].ips, ParT1[px].m_v );

  // check minimum/maximum particle number in node
  switch( nodeType )
  {
    case normal:    // normal reactive node
              break;
    case NBC3source:    // = 3, Cauchy source ( constant flux )
        if( new_node == -1 )
        {
// Only for 1D calculation !!! check for 2D and 3D
            double new_x = ParT1[px].xyz.x;
            if( new_x >= nodes->GetSize().x )
            {    new_x -= nodes->GetSize().x;
                 new_node = 0;
            }
            else // new_x < 0
            {    new_x  += nodes->GetSize().x;
                 new_node = nodes->nNodes()-1;
            }
            ParT1[px].xyz.x = new_x;
//            new_node = nodes->FindNodeFromLocation( ParT1[px].xyz, -1 );
        }
        break;
     case NBC3sink: // =-3, Cauchy sink (constant flux)
        if( new_node == -1 )
        {
// Only for 1D calculation !!! check for 2D and 3D
          double new_x = ParT1[px].xyz.x;
          if( new_x >= nodes->GetSize().x )
          {    new_x -= nodes->GetSize().x;
               new_node = 0;
          }
          else // new_x < 0
          {    new_x  += nodes->GetSize().x;
               new_node = nodes->nNodes()-1;
          }
          ParT1[px].xyz.x = new_x;
//            new_node = nodes->FindNodeFromLocation( ParT1[px].xyz, -1 );
        }
        break;
/*    if( new_node == -1 )
    {   LOCATION nodeSize[2];
        new_node = ndxCsource;
        nodes->GetNodeSizes( new_node, nodeSize );
// Only for 1D calculation !!! check for 2D and 3D
       double new_x = ParT1[px].xyz.x;
       if( new_x >= nodes->GetSize().x )
          new_x -= nodes->GetSize().x+nodeSize[0].x;
       else // new_x < 0
          new_x  += nodeSize[1].x;
       ParT1[px].xyz.x = new_x;
    }
*/
     case NBC2source: //    = 2, Neumann source ( constant gradient )
        break;
     case NBC2sink:   //    = -2, Neumann sink (constant gradient)
        break;
     case NBC1source: //    = 1, Dirichlet source ( constant concentration )
        break;
     case NBC1sink:   //    = -1, Dirichlet sink (constant concentration)
        break;
    default:     // Possibly error - wrong node type!
        break;
  }

  if( new_node == -1 )
  {
    char buff[300];
    sprintf( buff, " pxOld=%ld npxNew=%ld",  old_node, new_node  );
    Error("W003RWM",buff);
  }
  else
  {
    NPnum[new_node*anPTypes+type_]++;
    NPnum[old_node*anPTypes+type_]--;
    ParT1[px].node = new_node;
    // may be statistic
  }

  return 0;
}

// call to the whole Random Walk method time step over all particles and nodes
// returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
// GEM was called before this function
long int TParticleArray::RandomWalkIteration( long int /*Mode*/, bool CompMode, double t01, double t11 )
{

  long int iNode, iType, /*iRet = 0,*/ cpx1;
  double *mass, m_;

// set up new masses to particles after GEM calculations
  mass = new double[nNodes*anPTypes];
  for( iNode=0; iNode < nNodes; iNode++ )
    for( iType=0; iType < anPTypes; iType++ )
    {
       m_ = nodes->GetNodeMass( iNode, (char)iType,
            ParTD[iType].tcode, ParTD[iType].ips );
       mass[iNode*anPTypes+iType] = m_/ NPnum[iNode*anPTypes+iType];
  }
  for( cpx1=0; cpx1 < anParts; cpx1++ )
    ParT1[cpx1].m_v = mass[ParT1[cpx1].node*anPTypes+ParT1[cpx1].ptype];
  delete[] mass;

 // new  particle position xyz
  for( cpx1=0; cpx1 < anParts; cpx1++ )
     /*iRet =*/ DisplaceParticle( cpx1, t01, t11 );

 // Walk (transport step) for particles between nodes
  for( cpx1=0; cpx1 < anParts; cpx1++ )
    /* iRet =*/ MoveParticleBetweenNodes( cpx1, CompMode, t01, t11 );

  bool reset = false;
  for( iNode=0; iNode < nNodes; iNode++ )
    for( iType=0; iType < anPTypes; iType++ )
    {
      if( NPnum[iNode*anPTypes+iType] < nPmin[iType] )
      {
        char buff[300];

        sprintf( buff, " Node=%ld npNum=%ld npMin=%ld ",
        iNode, NPnum[iNode*anPTypes+iType], nPmin[iType] );
        Error("W005RWM",buff);
      // or  alternative
      // reset = true;   ParticleArrayInit - reset all particles
     }
    }
  if( reset )  // reset all particles as start point
    ParticleArrayInit();

  return 0;
}

// call to the whole FiniteCell Walk method time step over all particles and nodes
// returns 0 if time step is accepted; not 0 if rejected (another dt is needed)
long int TParticleArray::FCellWalkIteration( long int /*Mode*/, bool /* CompMode*/, double /*t0*/, double /*t1*/ )
{
  return 0;
}

 // stub call for coupled mass transport calculation
long int TParticleArray::GEMPARTRACK( long int Mode, bool ComponentMode, double t0_, double t1_ )
{
        long int iRet=0;

  t0 = t0_;
  t1 = t1_;
  dt = (t1-t0);

  switch( Mode )
  {
   case 'W':  iRet = RandomWalkIteration( Mode, ComponentMode, t0, t1 );
              break;
   case 'V':  iRet = FCellWalkIteration( Mode, ComponentMode, t0, t1 );
              break;
  }
  return iRet;
}

void TParticleArray::logProfilePhMol( FILE* logfile, int inode )
{
  long int npa;
  for( long int jp=0; jp < nPTypes(); jp++ )
       {
            npa = getNPnum( inode, jp);   // number of particles in the node
            fprintf( logfile, "%-8ld ", npa );
       }
  fprintf( logfile, "\n" );
}

//=========================================================================

// Random numbers (re-written from Numerical Recipes in C)
// with uniform distribution
//
double randuni (double& x)
{ double m35=34359738368., m36=68719476736., m37=137438953472.;
  double a=0.,b=1.;
  if( x < 0 ) // Initialize process
  {
	long int j;
    double R;
    j=rand();
    R = ceil(24359738368.*j/RAND_MAX + 10000000000.);
    if( approximatelyZero(fmod(R,2)) )
      R=R+1.;
    x = R;
  }
  x=x*5.;
  if(x>=m37) x=x-m37;
  if(x>=m36) x=x-m36;
  if(x>=m35) x=x-m35;
 return(x/m35*(b-a)+a);
}

// normal point
double randnorm(double& x)
{ double R1=0.;
  long int j;
  for(j=0;j<101;j++)
    R1+=randuni(x);
      R1=(R1-101./2.)/pow(101./12.,0.5);
           R1=1./6.*(R1-(-3.0));
           if(R1<0.) R1=0.;
           if(R1>1.) R1=1.;
           return(R1);
// return(1./9.*(R1-(-4.5)));
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS  1.2e-7
#define RNMX  (1.0-EPS)

// Long period (> 2'1018) random number generator of L'Ecuyer with
// Bays-Durham shuffle and added safeguards.
// Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
// the endpoint values). Call with idum a negative integer to initialize;
// thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less than 1.
// Modified from Numerical Recipes in C
//
double ran2(long int& idum)
{
   long int j;
   long int k;
   static long int idum2=123456789;
   static long int iy=0;
   static long int iv[NTAB];
   double temp;
   if (idum <= 0)
   { // Initialize.
      if ( -idum < 1)
         idum=1; // Be sure to prevent idum = 0.
      else
         idum = -idum;
      idum2= idum;
     for (j=NTAB+7;j>=0;j--) // Load the shuffle table (after 8 warm-ups).
     {
       k = idum/IQ1;
       idum = IA1*(idum-k*IQ1)-k*IR1;
       if (idum < 0)
         idum += IM1;
       if (j < NTAB)
         iv[j] = idum;
     }
   iy=iv[0];
   }
   k = idum/IQ1;                       // Start here when not initializing.
   idum = IA1 * (idum-k*IQ1) - k*IR1;  // Compute idum=(IA1*idum) % IM1 without
                                       // overflows by Schrage's  method.
   if ( idum < 0 )
      idum += IM1;
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;  // Compute idum2=(IA2*idum) % IM2 likewise.
   if (idum2 < 0)
      idum2 += IM2;
   j = iy/NDIV;                      // Will be in the range 0..NTAB-1.
   iy = iv[j]-idum2;                 // Here idum is shuffled, idum and idum2 are
   iv[j] = idum;                     // combined to generate output.
   if (iy < 1)
      iy += IMM1;
   if ( ( temp=AM*iy) > RNMX)
       return RNMX;         //  Because users don't expect endpoint values.
   return temp;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// According to Knuth, any large MBIG, and any smaller (but still large) MSEED
// can be substituted for the above values.
// Returns a uniform random deviate between 0.0 and 1.0.
// Set idum to any negative value to initialize or reinitialize the sequence.
// Modified From Numerical Recipes in C
//
double ran3(long int& idum)
{
  static long int inext,inextp;
  static long int ma[56];     // The value 56 (range ma[1..55]) is special and
  static long int iff=0;       // should not be modified; see  Knuth.
  long int mj,mk;
  long int i,ii,k;
  if ( idum < 0 || iff == 0) // Initialization.
  {
     iff=1;
     mj = labs( MSEED - labs(idum));  // Initialize ma[55] using the seed idum
     mj %= MBIG;                      // and the large number MSEED.
     ma[55] = mj;
     mk = 1;
     for (i=1;i<=54;i++)   // Now initialize the rest of the table,
     {  ii=(21*i) % 55;    // in a slightly random order,
        ma[ii]=mk;         // with numbers that are not especially random.
        mk=mj-mk;
        if (mk < MZ)
          mk += MBIG;
        mj=ma[ii];
      }
      for (k=1;k<=4;k++)    // We randomize them by 'warming up the generator'
        for(i=1;i<=55;i++)
        {
         ma[i] -= ma[1+(i+30) % 55];
         if (ma[i] < MZ)
            ma[i] += MBIG;
         }
       inext=0;     // Prepare indices for our first generated number.
       inextp=31;  //  The constant 31 is special; see Knuth.
       idum=1;
   }
   //  Here is where we start, except on initialization.
    if (++inext == 56)
       inext=1;           // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56)
       inextp=1;
    mj = ma[inext]-ma[inextp];  // Generate a new random number subtractively.
    if (mj < MZ)
         mj += MBIG;    //  Be sure that it is in range.
    ma[inext]=mj;       // Store it,
    return mj*FAC;      // and output the derived uniform deviate.
}


//-----------------------End of particlearray.cpp--------------------------

