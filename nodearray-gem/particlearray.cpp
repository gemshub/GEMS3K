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

#include "particlearray.h"
#include "num_methods.h"

static long idum = -10000l;
static double Rand = -1;

TParticleArray* TParticleArray::pa;


TParticleArray::TParticleArray( short nPTypes, short nProps,
           short *aNPmean,
           short (*aParTD)[6],
           short *anPmin, short *anPmax,
           TNodeArray* aNodes ):
  ParT0(0), ParT1(0), anParts(0),
  anPTypes(nPTypes), anProps(nProps), nodes(aNodes)
{
  int ii;
  nNodes = nodes->nNodes();

  NPmean = new short[anPTypes];
  ParTD = new PARTICLE[anPTypes];
  nPmin = new short[anPTypes];
  nPmax = new short[anPTypes];

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
  NPnum = new short[nNodes*anPTypes];
  ParticleArrayInit();

  cParts = 0;
  NPlist = 0;
  NPstat = 0;
  pa = this;
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
  int iNode, iType, k, cpx;
  LOCATION nodeSize[2];

  cpx = 0;
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
        ParT0[cpx].ptype = (char)iType;
        ParT0[cpx].mmode = ParTD[iType].mmode;
        ParT0[cpx].tcode = ParTD[iType].tcode;
        ParT0[cpx].ips = ParTD[iType].ips;
        ParT0[cpx].m_v = 0.;
        ParT0[cpx].node = iNode;
//        ParT0[cpx].xyz.x = nodeSize[0].x+k*dd;
        ParT0[cpx].xyz = setPointInNode(nodeSize);
        ParT1[cpx] = ParT0[cpx];
        cpx++;
      }
    }
  }
}

// Particle array initialization
void TParticleArray::setUpCounters()
{
  int iNode;
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

void TParticleArray::CopyfromT1toT0()  // Copy resalts of ParT1 step to ParT0
{
/*  fstream f_out("nods_particl.out", ios::out|ios::app  );
  for(int iNode=0; iNode < nNodes; iNode++ )
  {
    for(int iType=0; iType < anPTypes; iType++ )
       f_out << NPnum[iNode*anPTypes+iType] << " ";
    f_out << " || ";
  }
 f_out << endl;

*/
  for(int k=0; k < anParts; k++ )
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
    loc.x = (float)randuni( Rand ); //  randnorm(Rand); // value from 0. to 1.
//    loc.x = ran3( idum ); //  ran2( idum ); // value from 0. to 1.
    loc.x = (loc.x*(nodeSize[1].x-nodeSize[0].x)+nodeSize[0].x);
  }
  else loc.x = 0;
  if( nodes->SizeM() > 1 )
  {
    loc.y = (float)randnorm(Rand ); //  randuni( Rand); // value from 0. to 1.
//    loc.y = ran2( idum ); //  ran3( idum ); // value from 0. to 1.
    loc.y = (loc.y*(nodeSize[1].y-nodeSize[0].y)+nodeSize[0].y);
  }
  else loc.y = 0;
  if( nodes->SizeK() > 1 )
  {
    loc.z = (float)randnorm(Rand ); //  randuni( Rand); // value from 0. to 1.
//    loc.z = ran2( idum ); //  ran3( idum ); // value from 0. to 1.
    loc.z = (loc.z*(nodeSize[1].z-nodeSize[0].z)+nodeSize[0].z);
  }
  else loc.z = 0;
  return loc;
}

// Implementation of interpolation for particle advection velocities and
// dispersivities between nodes ( in 1D case)
// px index of particle
double TParticleArray::InterpolationVp_hDl_1D( int px,
   double& vp, double& al, double& Dif )
{
  if( nodes->SizeM() > 1 ||  nodes->SizeK() > 1 )
     Error( "InterpolationVp_hDl_1D", "Error mode of interpolation." );

  DATABR *dbr1, *dbr2;    // nodes for interpolation
  int nodInd1, nodInd2;  // number of nodes for interpolation
  double hDl, x1m, x2m;       // middle-point coordinates in the nodes
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
    vp = dbr1->vp;
    hDl = dbr1->hDl;
    al = dbr1->al;
    Dif = dbr1->Dif;
  }
  else
  {
    nodes->GetNodeSizes( nodInd2, nodeSize );
    x2m = nodeSize[0].x + (nodeSize[1].x-nodeSize[0].x)/2 ;
    dbr2 = nodes->pNodT1()[nodInd2];
    // vx = v1 - (v2-v1)/(x2-x1)*(x-x1);
    double d = (ParT1[px].xyz.x-x1m)/(x2m-x1m);
    vp = dbr1->vp;
    vp -= (dbr2->vp - dbr1->vp )*d;
    al = dbr1->al;
    al -= (dbr2->al - dbr1->al )*d;
    Dif = dbr1->Dif;
    Dif -= (dbr2->Dif - dbr1->Dif )*d;
    hDl = al*vp+Dif;
  }
  return hDl;
}

// Important for the mass-transport step
// Calculation of new particle locations
// Advective step, Brownian step, Dispersive step, Diffusion step
// px: index of particle
int TParticleArray::DisplaceParticle( int px, double /*t0*/, double /*t1*/ )
{
//  DATACH* ch = nodes->pCSD();       // DataCH structure
  int nodInd = ParT1[px].node;
  double ds = 0.;
  double vp, hDl, al, Dif;

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
        hDl = InterpolationVp_hDl_1D( px, vp, al, Dif );
// vp = dbr->vp;     // testing without interpolation
// hDl = dbr->hDl;   // testing without interpolation
         if( hDl > 0)
//         ds = 2.*(randuni( idum )-0.5)*sqrt( 6.*hDl*vp*dt);
         ds = 2.*(ran3( idum )-0.5)*sqrt( 6.*hDl*dt);
         ParT1[px].xyz.x += vp*dt + ds;
                        break;
   }

  return 0;
}

// Walk (transport step) for particle px between nodes
// returns -1 if particle stays in the same node
// or an index of another node the particle enters
int TParticleArray::MoveParticleBetweenNodes( int px, double /*t0*/, double /*t1*/ )
{
  int old_node = ParT1[px].node;
  int new_node = nodes->FindNodeFromLocation( ParT1[px].xyz, old_node );
  char type_ = ParT1[px].ptype;
  short nodeType = nodes->pNodT1()[old_node]->NodeTypeHY;

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

  // move mass of old_node to new_node (if new_node==-1 mass be lost)
  nodes->MoveParticleMass( old_node, new_node, type_,
            ParT1[px].tcode, ParT1[px].ips, ParT1[px].m_v );

  // check minimum/maximum particle number in node
  switch( nodeType )
  {
    case normal:
              break;
    case NBC3source:
    case NBC3sink:
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
             break;
  }

  if( new_node == -1 )
  {
    vstr buff(300);
    sprintf( buff, " pxOld=%d npxNew=%d",  old_node, new_node  );
    Error("W003RWM",buff.p);
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
int TParticleArray::RandomWalkIteration( int /*Mode*/, double t0, double t1 )
{

  int iNode, iType, iRet=0, cpx;
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
  for( cpx=0; cpx < anParts; cpx++ )
    ParT1[cpx].m_v = mass[ParT1[cpx].node*anPTypes+ParT1[cpx].ptype];
  delete[] mass;

// new  particle position xyz
  for( cpx=0; cpx < anParts; cpx++ )
     iRet = DisplaceParticle( cpx, t0, t1 );


 // Walk (transport step) for particles between nodes
  for( cpx=0; cpx < anParts; cpx++ )
     iRet = MoveParticleBetweenNodes( cpx, t0, t1 );

  bool reset = false;
  for( iNode=0; iNode < nNodes; iNode++ )
    for( iType=0; iType < anPTypes; iType++ )
    {
      if( NPnum[iNode*anPTypes+iType] < nPmin[iType] )
      {
        vstr buff(300);

        sprintf( buff, " Node=%d npNum=%d npMin=%d ",
        iNode, NPnum[iNode*anPTypes+iType], nPmin[iType] );
        Error("W005RWM",buff.p);
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
int TParticleArray::FCellWalkIteration( int /*Mode*/, double /*t0*/, double /*t1*/ )
{
  return 0;
}

 // stub call for coupled mass transport calculation
int TParticleArray::GEMCOTAC( int Mode, double t0_, double t1_ )
{
  int iRet;

  t0 = t0_;
  t1 = t1_;
  dt = (t1-t0);
  switch( Mode )
  {
   case 'W':  iRet = RandomWalkIteration( Mode, t0, t1 );
              break;
   case 'V':  iRet = FCellWalkIteration( Mode, t0, t1 );
              break;
  }
  return iRet;
}


//-----------------------End of particlearray.cpp--------------------------

