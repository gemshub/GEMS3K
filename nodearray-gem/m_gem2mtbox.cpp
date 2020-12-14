//-------------------------------------------------------------------
// $Id: m_gem2mtbox.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of TInteg/TGEM2MT classes, calculation functions
//
// Rewritten from C to C++ by S.Dmytriyeva  
// Copyright (C) 1995,2008,2011  S. Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Selektor GUI GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#include <cmath>
#include <cstdio>
#include <iomanip>


#ifndef IPMGEMPLUGIN

#include "m_gem2mt.h"
#include "m_compos.h"
#include "visor.h"
#include "stepwise.h"

#else

#include <ctime>
#include "m_gem2mt.h"
#include "nodearray.h"

#endif

#define dMB( q, i) ( dm[ (q)*mtp->Nf + (i)] ) 

#define MB( q, i)  ( m[ (q)*mtp->Nf + (i)] )

#define g(q,f,i)   ( mtp->gfc[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )
#define y(q,f,i)   ( mtp->yfb[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )

#define ord(kk)    ( ROUND(mtp->FDLf[kk][0]) )
#define v(kk)      ( (mtp->FDLf[kk][1]) )

#define H(f, i)    ( mtp->BSF[ (f) * mtp->Nf + ( i )] )


inline double TGEM2MT::MassICinMGP( long int q, long int f, long int i )
{
    double MassIC=0.;
    if(f >= 0)
       MassIC = y(q,f,i) * na->ICmm( i );
    return MassIC;
}

// Calculation of mass of f-th MGP in q-th box (in kg)
double TGEM2MT::MassMGP( long int q, long int f  )
{
    double Mass = 0.0;
    if ( f < 0 )
       Mass = 0.0;
    else {
       for(long int i=0; i<mtp->Nf; i++ )
          Mass += y(q,f,i) * na->ICmm( i );
    }  // Mass /= 1000.;
//    fstream f_log("ipmlog.txt", ios::out|ios::app );
//    f_log << " MGP= " << Mass << " q= " << q << " f= " << f << endl;
    return Mass;
}

//   Increment to mass differential for i-th IC in q-th box using rate
//
inline void TGEM2MT::dMBfluxDir( long int q, long int i, double *dm, double fRate, double sign )
{
    switch(mtp->DiCp[q][1])
    {                  // boundary condition node
       case NBC3source:  // 3: Cauchy source ( constant flux )
                 if( sign < 0 )
                      break;
                 [[fallthrough]];
       case NBC1source:  //1: Dirichlet source ( constant concentration )
       case NBC1sink:    // -1: Dirichlet sink
       case NBC2source:  // 2: Neumann source ( constant gradient )
       case NBC2sink:    // -2: Neumann sink
       case NBC3sink:    // -3: Cauchy sink
       case INIT_FUNK:   // 4: functional conditions (e.g. input time-depended functions)
       case normal: // normal node  // Checking phase assemblage
       default:    dMB(q,i) +=  fRate*sign;
                   break;
    }
}

// calculate kk-th flux between boxes with indexes in FDLi[kk]
//  !!!!!!!!!!!!!!!!!!!!!!!!! To be checked and revised!
//
void TGEM2MT::dMBflux( long int kk, double *m, double *dm, double /*t*/ )
{
  long int  q, p, f, fe=-1; int i;
  char FLXid[MAXSYMB+1], MGPid[MAXSYMB+1];
  long int FLXorder/*, FLXtype*/;
  double nbIC, fRate, fRateS, MGPo=0., MGPi=0., vkk=1., mqfi, gqfi, Hfei, MBqi;//, MBpi;
  bool sinkOut=false;

  fe = -1;
  q = mtp->FDLi[kk][0];    // index of outgoing box
  p = mtp->FDLi[kk][1];    // index of receiving box

  strncpy( MGPid, mtp->FDLmp[kk], MAXSYMB );
  MGPid[MAXSYMB] = 0;                           // MGP identifier
  f = LookUpXMGP( MGPid );                      // Getting MGP index

  if( f<0 ) // Elemental production flux?
  {
        sscanf( MGPid, "%ld", &fe );  // Reading BSF row index
        if( fe < 0 || fe >= mtp->nSFD )
           Error( "BOXFLUX:", "Wrong MGP identifier or BSF row index!" );
  }

  strncpy( FLXid, mtp->FDLid[kk], MAXSYMB );
  FLXid[MAXSYMB] = 0;      // Flux identifier
  FLXorder = ord(kk);      // flux order
  vkk = v(kk);             // flux rate constant
  MGPo = MassMGP( q, f  ); // calculating MGP masses: in outgoing box
  if( p >= 0 )
      MGPi = MassMGP( p, f  ); // in incoming box
  if( q >= 0 && f >= 0 )
  {            // Normal MGP flux from box q to box p
   // NB: Negative v(kk) means "production" in q box and "consumption" in p box
        if( p < 0 )
            sinkOut = true; // This is a sinkout flux from box q to nowhere (if v > 0)
                            //  or production of MGP in box q (if v < 0)
        switch( FLXorder )
        {
          case 0:  // Zero-order (absolute) flux;  v(kk) in 1/s units
                 for(i=0; i<mtp->Nf; i++ )
                 {
                    nbIC = node1_bIC(q, i); gqfi = g(q,f,i);
                    mqfi = gqfi * nbIC * na->ICmm( i );
                    fRate = mqfi * vkk;
                    if( mtp->PsMode != RMT_MODE_S || ( sinkOut && mtp->PsMode == RMT_MODE_S ) )
                    {
                        if( mtp->PsMode != RMT_MODE_S )
                            dMBfluxDir( q, i, dm, fRate, (-1.)); //    dMB(q,i) -=  fRate;
                        else { // remove the whole MGP stuff in S mode
                            fRateS = y(q,f,i) * na->ICmm( i ); // * vkk;
                            dMBfluxDir( q, i, dm, fRateS, (-1.));
                        }
                    }
                    if( !sinkOut /* && p != q */ )
                        dMBfluxDir( p, i, dm, fRate); //dMB(p,i) +=  fRate;
                 }
                 break;
          case 1:  // First-order to source box mass flux     v in 1/s/kg units
                 for(i=0; i<mtp->Nf; i++ )
                 {
                     gqfi = g(q,f,i); MBqi = MB(q,i); nbIC = node1_bIC(q, i);
                     fRate = gqfi * vkk * MBqi; //* * nbIC * na->ICmm( i ) ;
                     if( mtp->PsMode != RMT_MODE_S || ( sinkOut && mtp->PsMode == RMT_MODE_S ) )
                     {
                         if( mtp->PsMode != RMT_MODE_S )
                             dMBfluxDir( q, i, dm, fRate, (-1.)); //    dMB(q,i) -=  fRate;
                         else { // remove the whole MGP stuff in S mode
                             fRateS = y(q,f,i) * na->ICmm( i ); // * vkk;
                             dMBfluxDir( q, i, dm, fRateS, (-1.));
                         }
                     }
                     if( !sinkOut /* && p != q */ )
                         dMBfluxDir( p, i, dm, fRate); //dMB(p,i) +=  fRate;
                 }
                 break;
          case 2:  // Second-order to source and sink MGP mass flux v in 1/s/kg2 units
                         // (proportional to product of masses in both boxes)
                 if( sinkOut )
                     break;
                 for(i=0; i<mtp->Nf; i++ )
                 {
                    gqfi = g(q,f,i); MBqi = MB(q,i); //MBpi = MB(p,i);
                    fRate = gqfi * vkk * MGPo * MGPi; // * MBqi * MBpi;
                    dMBfluxDir( q, i, dm, fRate, (-1.)); // dMB(q,i) -=  fRate;
                    dMBfluxDir( p, i, dm, fRate); // dMB(p,i) +=  fRate;
                 }
                 break;
           case 3:  // was -1 First-order flux to MGP mass in receiving box v in 1/s/kg units
                 if( sinkOut )
                    break;
                 for(i=0; i<mtp->Nf; i++ )
                 {
                    gqfi = g(q,f,i);// MBpi = MB(p,i);
                    fRate = gqfi * vkk * MGPi; // * MBpi;
                    dMBfluxDir( q,i, dm, fRate, (-1.)); //dMB(q,i) -=  fRate;
                    dMBfluxDir( p,i, dm, fRate); //dMB(p,i) +=  fRate;
                 }
                 break;
          default:  // Other orders or flux types - to be done!
                 break;
        }
        ;
    }
   else {
      if( q >= 0 && f < 0 )
      {    // This is an elemental flux ( fe row in BSF table )
             // NB: Negative v(f) means "production" in q box and "consumption" in p box
         if( p < 0 )
             sinkOut = true; // This is an elemental sinkout flux from box q to nowhere (if v > 0)
                                  //  or elemental production in box q (if v < 0)
         switch( FLXorder )
         {
            case 0:  // Zero-order flux
                  for(i=0; i<mtp->Nf; i++ )
                  {
                     Hfei = H(fe,i);
                     fRate = Hfei * vkk;
                     dMBfluxDir( q, i, dm, fRate, (-1.)); //    dMB(q,i) -=  fRate;
                     if( !sinkOut && p != q )
                         dMBfluxDir( p,i, dm, fRate); //dMB(p,i) +=  fRate;
                   }
                   break;
             case 1:  // First-order to outgoing box flux
                   for(i=0; i<mtp->Nf; i++ )
                   {
                      Hfei = H(fe,i); MBqi = MB(q,i);
                      fRate = Hfei * vkk * MGPo; // * MBqi; // * MassICinMGP(q,f,i); // * MGPo;
                      dMBfluxDir( q, i, dm, fRate, (-1.)); //dMB(q,i) -=  fRate;
                      if( !sinkOut && p != q)
                          dMBfluxDir( p,i, dm, fRate); //dMB(p,i) +=  fRate;
                   }
                   break;
             case 2:  // Second-order to source and sink flux
                           // (proportional to product of MGP masses in both boxes)
                   if( sinkOut )
                         break;
                   for(i=0; i<mtp->Nf; i++ )
                   {
                      Hfei = H(fe,i); MBqi = MB(q,i); //MBpi = MB(p,i);
                      fRate = Hfei * vkk * MGPo * MGPi; // * MBqi * MBpi;
                      dMBfluxDir( q,i, dm, fRate, (-1.)); //dMB(q,i) -=  fRate;
                      if( p != q )
                          dMBfluxDir( p,i, dm, fRate); //dMB(p,i) +=  fRate;
                   }
                   break;
             case 3:  // was -1 First-order to receiver flux only
                   if( sinkOut )
                       break;
                   for(i=0; i<mtp->Nf; i++ )
                   {
                       Hfei = H(fe,i); //MBpi = MB(p,i);
                       fRate = Hfei * vkk * MGPi; // * MBpi;
                       dMBfluxDir( q,i, dm, fRate, (-1.)); //dMB(q,i) -=  fRate;
                       if( p != q )
                           dMBfluxDir( p,i, dm, fRate); //dMB(p,i) +=  fRate;
                   }
                   break;
             default:  // Other orders or flux types - nothing to be done yet!
                   break;
           }
        }
    }
}

// Cleaning up the table of IC mass differentials in boxes/nodes
void TGEM2MT::dMBZeroOff(  double *dm )
{
  long int q, i;
  // Zeroing IC mass differentials off
  for( q=0; q <mtp->nC; q++ )
        for(i =0; i< mtp->Nf; i++ )
      dMB(q,i) = 0.;
}

// calculate transport step at time cTau using the table of fluxes
void TGEM2MT::Solut( double *m, double *dm, double tau )
{ 
  long int kk;

  // Zeroing IC mass differentials off
  dMBZeroOff(  dm );

  for(kk=0; kk<mtp->nFD; kk++ )  // Looking through the list of fluxes 
  {
     dMBflux( kk, m, dm, tau );
  }   // kk
}

#undef dMB
#undef MB

#define Mb( q, i)  ( mtp->MB[(q)*mtp->Nf + (i)])
#define dMb( q, i)  (mtp->dMB[(q)*mtp->Nf + (i)])

// change bulk composition in the box q using changes dM*dTau over time step
void TGEM2MT::BoxComposUpdate( long int q )
{
   long int i;
// double n1bICo, n1bICn=0., dMbqi, ICmm, dTau;  // for debugging

   for(i =0; i< mtp->Nf; i++ )
   {
// n1bICo = node1_bIC(q, i); dMbqi = dMb( q, i); ICmm = na->ICmm( i ); dTau =  mtp->dTau;

       if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_AQUEL && i == mtp->Nf-1 )
       {
         node1_bIC(q, i) += dMb( q, i) * mtp->dTau;   // molar mass of charge is 0
// n1bICn = node1_bIC(q, i);
         node1_bIC(q, i) = 0.0;  // Provisorial - zeroing-off charge for the node
       }
       else {                                         // dTau is the time step (in s)
         node1_bIC(q, i) += dMb( q, i) / na->ICmm( i ) * mtp->dTau;
// n1bICn = node1_bIC(q, i);
         if( node1_bIC(q, i) < 1e-12 )
            node1_bIC(q, i) = 1e-12;    // preventing too small or negative amounts of ICs
       }
   }
}

// Update all box compositions at cTau using time step dTau
void TGEM2MT::BoxesBCupdate()
{
  long int q;
  for( q=0; q <mtp->nC; q++ )
  {
     BoxComposUpdate( q );
  }  // q
}

// Calculation of q-th box reactive IC masses in kg
// returns total mass of the box in kg
//
double TGEM2MT::BoxMasses( long int q )
{
   long int i;
   double BoxMass = 0.; //, Mbqi = 0.;
   for(i=0; i<mtp->Nf; i++ )
   {
       Mb( q, i) = node1_bIC( q, i ) * na->ICmm( i );
// Mbqi = Mb(q,i);    // debugging
       BoxMass += Mb(q,i);
   }
   return BoxMass;
}

// Calculates table of MGP bulk compositions for q-th box (now added phase conc. units)
// variant for simplified box-flux and sequential reactors transport models
// Differs from 'standard' module in interpretation of the mole amount of MPG:
//     1 'M' means the whole source amount (of aq, gas or solid phase) for S models
//     and just mole amount of phase for F and B models
//
void TGEM2MT::ComposMGPinBox( long int q  )
{
    double PhFract=0., Xincr=1., Xe, PHmw, Vm=1., R1=1., Msys=1.,
            Mwat=1., Vaq=1., Maq=1., Vsys=1., Mph=0., Mgas=0.;
    char UNITP = 'M';
    long int f, i, k, x_aq=-1, x_gf=-1, naqgf=0;

    for(f=0; f<mtp->nPG; f++ )
          for(i=0; i<mtp->Nf; i++ )
            y(q,f,i) = 0.0;
    Msys = TNodeArray::na->pNodT1()[q]->Ms; // in kg
    Vsys = TNodeArray::na->pNodT1()[q]->Vs *1000.; // in dm3

    // Index of aq phase in DBR
    if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_AQUEL )
        x_aq = 0;

    if( x_aq == 0 )
    {   if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_GASMIX ||
            na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_FLUID )
           x_gf = 1; }
    else if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_GASMIX ||
             na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_FLUID )
        x_gf = 0;

    for(f=0; f<mtp->nPG; f++ )
      for( k=0; k<mtp->FIf; k++)
      {
          Xe = mtp->PGT[k*mtp->nPG+f];  // given quantity of phase to add to MGP
          if( fabs(Xe)<1e-19 )
              continue;
          UNITP = mtp->UMGP[k];
          if( k < na->pCSD()->nPSb )
          {
             PHmw = node1_mPS( q, k )        // mass of solution phase kg
                  *1000./ node1_xPH( q, k ); // amount of phase, PHmw in g/mol
             Vm = node1_vPS( q, k ) *1e6  // molar volume of solution phase cm3/mol
                  / node1_xPH( q, k );
          }
          else {   // pure phase
              PHmw = node1_mPH( q, k ) * 1000. / node1_xPH( q, k );
              Vm = node1_vPH( q, k ) * 1e6 / node1_xPH( q, k ); // cm3/mol
          }
          Mph = node1_mPS( q, k );  // mass of phase in kg;
          R1 = node1_xPH( q, k );  // amount of phase in moles
          if( k == 0 && x_aq == 0 )
          {
             Mwat = node1_xPA( q, k )*18.0153/1000.;  // mass of water in kg
             Vaq = node1_vPS( q, k ) * 1000; // in dm3 assuming density 1 g/cm3
             Maq = node1_mPS( q, k );  // in kg
             R1 = Maq*1e3/PHmw; // approximate mole amount of aqueous phase
          }

          // Converting quantity Xe into amount (moles)
#ifndef IPMGEMPLUGIN
          Xincr = TCompos::pm->Reduce_Conc( UNITP, Xe, PHmw, Vm, R1, Msys,
                  Mwat, Vaq, Maq, Vsys );
#else
          Xincr = Reduce_Conc( UNITP, Xe, PHmw, Vm, R1, Msys,
                  Mwat, Vaq, Maq, Vsys );
#endif

          PhFract = 0.; // Initial fraction of real phase amount
          naqgf = 0;
          if( (x_aq + x_gf) == 1 )
              naqgf = 1;
          switch( mtp->PsMPh )
          {
           case MGP_TT_AQS:   // '1'   Assuming mole amount as fraction of the aq phase mass
                if( k == 0 )  //       and mass as 'true' mass in the 'S' mode
                {
                   if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL )
                   {
                       if( mtp->PsMode == RMT_MODE_S )
                           PhFract = Xincr; //  /R1;
                       else PhFract = Xincr/node1_xPH( q, k );
                   }
                   else if( UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM )
                       PhFract = Xincr * PHmw / 1e3 / Mph;
                }
              break;
           case MGP_TT_SOLID:  // '4'  TBD  assuming
              if( k > naqgf )  //   mole amount as fraction of the phase mass and mass as 'true' mass
              {
                 if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL )
                 {
                     if( mtp->PsMode == RMT_MODE_S )
                         PhFract = Xincr;  // /R1;
                     else PhFract = Xincr/node1_xPH( q, k );
                 }
                 else if( UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM )
                     PhFract = Xincr * PHmw /1e3 / ( Msys - Maq ); // Need Msolid and bICsolid vector in DATABR
              }
              break;
           case MGP_TT_AQGF:   // '3' Assuming mole amount as fraction of the aq or gas phase mass
              if( k == 0 || k == 1 )  //   and mass as 'true' mass
              {
                 if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL )
                 {
                     if( mtp->PsMode == RMT_MODE_S )
                         PhFract = Xincr; // /R1;
                     else PhFract = Xincr/node1_xPH( q, k );
                 }
                 else if( UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM )
                     PhFract = Xincr * PHmw / 1e3 / Mph;
              }
              break;
           case MGP_TT_GASF:   // '2'  TBD   assuming only gas/fluid and no aq phase is present?
              if( k == naqgf )  //    mole amount as fraction of the gas/fluid mass and mass as 'true' mass
              {
                 if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL )
                 {
                     if( mtp->PsMode == RMT_MODE_S )
                         PhFract = Xincr; // /R1;
                     else PhFract = Xincr/node1_xPH( q, k );
                 }
                 else if( UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM )
                     PhFract = Xincr * PHmw / 1e3 / Mgas;
              }
              break;
           default: break;
          }

          if( fabs(PhFract) > 0.0 )   // calculate IC compositions of MGP and MGP IC partition coeffs
          {
 //double yqfio, yqfin, n1bPHqki, gqfio, gqfin; // debugging
            if( k < na->pCSD()->nPSb )
                 for(i=0; i<mtp->Nf; i++ )
                 {
// yqfio = y(q,f,i); n1bPHqki = node1_bPS( q, k, i );
                    y(q,f,i) += node1_bPS( q, k, i );
// yqfin = y(q,f,i); gqfio = g(q,f,i);
                    g(q,f,i) = ( node1_bIC( q, i ) > 0. ? y(q,f,i)*PhFract/node1_bIC( q, i ): 0.);
// gqfin = g(q,f,i);
                 }
            else
                 for(i=0; i<mtp->Nf; i++ )
                 {
 // yqfio = y(q,f,i); n1bPHqki = node1_bPH( q, k, i );
                    y(q,f,i) += node1_bPH( q, k, i );
 // yqfin = y(q,f,i); gqfio = g(q,f,i);
                   g(q,f,i) = ( node1_bIC( q, i ) > 0. ? y(q,f,i)*PhFract/node1_bIC( q, i ): 0.);
 // gqfin = g(q,f,i);
                 }
            if( x_aq == 0 )
            {
                y(q,f,mtp->Nf-1) = 0.0;    // provisional - zeroing-off charge
                g(q,f,mtp->Nf-1) = 0.0;
            }
          }  // k
      }  // f
}

// Calculation of box masses, MGP compositions and MGP distribution coefficients in all boxes
//
void TGEM2MT::CalcMGPdata()
{
  long int q;
    // Calculation of current box reactive IC masses in kg
       for( q=0; q <mtp->nC; q++ )
           BoxMasses( q );
    // Calculation of MGP bulk compoisitions in boxes (in moles of ICs)
       for( q=0; q <mtp->nC; q++ )
           ComposMGPinBox( q );
 }

// Calculate new equilibrium states in the boxes for tcur = t
//  Ni 
//  pr
//  tcur - current time
//  step - current time step
//
bool
TGEM2MT::BoxEqStatesUpdate(  long int Ni, long int /*pr*/, double tcur, double step )
{
  bool iRet = true;
  FILE* diffile = nullptr;

  mtp->dTau = step;
  mtp->cTau = tcur;

#ifndef IPMGEMPLUGIN
  if( Ni >= 0)
  {
      char buf[300];
      //std::string Vmessage;
      sprintf(buf, "   step %ld; time %lg; dtime %lg  ", mtp->ct, mtp->cTau, mtp->dTau );
      Vmessage = "Simulating Reactive Transport in a Box-Flux chain: ";
      Vmessage += buf;
      Vmessage += ". Please, wait (may take time)...";

    if( mtp->PsSmode != S_OFF  )
    {
      STEP_POINT2();
    }
    else
      iRet = pVisor->Message( window(), GetName(),Vmessage.c_str(),
                           nstep, mtp->ntM );


   if( iRet )
         Error("GEM2MT Box-Flux model", "Cancelled by the user");
  }
#endif
  
 if( Ni >= 0)
 { // Update bulk compositions in boxes at current time point
   BoxesBCupdate();
 } 
 // Calculate new box equilibrium states at current time tcur
 if( mtp->PsSIA != S_ON )
     CalcIPM( NEED_GEM_AIA, 0, mtp->nC, diffile );
 else
     CalcIPM( NEED_GEM_SIA, 0, mtp->nC, diffile );
 
 if( Ni >= 0 )
 { // Here one has to compare old and new equilibrium phase assemblage
   // and pH/pe in all nodes and decide if the time step was Ok or it
   // should be decreased. If so then the nodes from C0 should be
   // copied to C1 (to be implemented)

   // Output of the results if step accepted
   if( mtp->PsMO != S_OFF  )
      PrintPoint( 1 );
 }

#ifndef IPMGEMPLUGIN
   // time step accepted - Copying nodes from C1 to C0 row
      pVisor->Update();
      CalcGraph();
#endif
   
  // copy node array for T0 into node array for T1
  mtp->oTau = mtp->cTau;
  copyNodeArrays();

  CalcMGPdata();

  return iRet;
}

// Initialization of box-flux transport calculations
//
void TGEM2MT::BoxFluxTransportStart()
{
    mtp->dTau = mtp->Tau[STEP_];;
    mtp->oTau =  mtp->Tau[START_];
    // mtp->cTau = mtp->Tau[START_];
    mtp->ct = 0;
    mtp->qf = 0;

#ifndef IPMGEMPLUGIN
    mtp->gfc = (double *)aObj[ o_mtgfc]->Alloc(  mtp->nC*mtp->nPG, mtp->Nf, D_);
    mtp->yfb = (double *)aObj[ o_mtyfb]->Alloc(  mtp->nC*mtp->nPG, mtp->Nf, D_);
    mtp->tt  = (double (*)[9])aObj[ o_mttt]->Alloc(  mtp->nC*mtp->Nf, 9, D_);
#else
    if( mtp->gfc )
            delete[] mtp->gfc;
    mtp->gfc = new double[mtp->nC * mtp->nPG * mtp->Nf];

    if( mtp->yfb )
            delete[] mtp->yfb;
    mtp->yfb = new double[mtp->nC * mtp->nPG * mtp->Nf];

    if( mtp->tt )
        delete[] mtp->tt;
    mtp->tt = new double[mtp->nC * mtp->Nf][9];
 #endif
    long int q, f, i;
       for( q=0; q <mtp->nC; q++ )
             for(i=0; i<mtp->Nf; i++ )
                 Mb( q, i) = 0.;

       for( q=0; q <mtp->nC; q++ )
           for(f=0; f<mtp->nPG; f++ )
               for(i=0; i<mtp->Nf; i++ )
               {
                   y(q,f,i) = 0.0;
                   g(q,f,i) = 0.0;
               }
}

void TGEM2MT::FlowThroughBoxFluxStep()
{
   mtp->ct += 1;
   mtp->oTau = mtp->cTau;
   mtp->cTau += mtp->dTau;
   // Calculate IC mass differentials
   Solut(  mtp->MB, mtp->dMB, mtp->cTau );
   // Update bulk compositions in boxes at cTau
   BoxesBCupdate();
}

// Returns MGP index from MGP identifier MGPid
//   or -1 if the identifier was not found in the MGP id list
long int TGEM2MT::LookUpXMGP( const char* MGPid )
{
        long int found = -1, smgpx = -1;
        // Check if the first character is 0 1 2 3 4 5 6 7 8 9
        switch(MGPid[0])
        {
            case '0': smgpx = 0; break;
            case '1': smgpx = 1; break;
            case '2': smgpx = 2; break;
            case '3': smgpx = 3; break;
            case '4': smgpx = 4; break;
            case '5': smgpx = 5; break;
            case '6': smgpx = 6; break;
            case '7': smgpx = 7; break;
            case '8': smgpx = 8; break;
            case '9': smgpx = 9; break;
            default: break;
        }
        // If so, this is index of elemental flux with composition from BSF table
        if(smgpx >= 0)
            return found;
        // Looking for a normal MGP index
        for( long int f=0; f < mtp->nPG; f++ )
        {
            if( strncmp( mtp->MGPid[f], MGPid, MAXSYMB ) )
                continue;
            found = f;
            break;
        }
        return found;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Calculate sequential-reactors model (complete moving of fluid)
// mode: SIA (0) or AIA (-1)
// Goes through boxes one-by-one with immediate re-equilibration
//
bool TGEM2MT::CalcSeqReacModel( char mode )
{

  try {
    //std::string Vmessage;
    long int p, i, kk, x_aq=-1, x_gf=-1;//, naqgf=1, lastp=0;
    bool iRet = false;
    clock_t outp_time = (clock_t)0;

    BoxFluxTransportStart();

#ifndef IPMGEMPLUGIN
    bool UseGraphMonitoring = false;
    char buf[300];


    if( mtp->PsSmode == S_OFF )
      if(  mtp->PvMSg != S_OFF && vfQuestion(window(),
             GetName(), "Use graphic monitoring?") )
        {
            RecordPlot( nullptr );
            UseGraphMonitoring = true;
        }
#endif

// na->CopyWorkNodeFromArray( 0, mtp->nC,  na->pNodT1() );
// na->GEM_write_dbr( "node0000.dat",  0, false );

    // In this mode, no start calculation of equilibria in all nodes!
    // Index of aq phase in DBR
    if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_AQUEL )
        x_aq = 0;
    if( x_aq == 0 )
    {   if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_GASMIX ||
            na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_FLUID )
           x_gf = 1; }
    else if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_GASMIX ||
             na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_FLUID )
        x_gf = 0;
    // here a distinction is needed between transport of aqueous, gas/fluid, and solids
    // Box 0 (source of aq fluid) is special case, not equilibrated
    switch( mtp->PsMPh )
    {
     case MGP_TT_AQS: // '1'
                     for(i=0; i<mtp->Nf; i++ )
                        node1_bPS( 0, x_aq, i ) = node1_bIC( 0, i );
                     TNodeArray::na->pNodT1()[0]->Ms = BoxMasses( 0 );  // in kg
                     TNodeArray::na->pNodT1()[0]->Vs = TNodeArray::na->pNodT1()[0]->Ms / 1e3; // in m3
                     node1_vPS( 0, 0 ) = TNodeArray::na->pNodT1()[0]->Vs;
                     node1_mPS( 0, 0 ) = TNodeArray::na->pNodT1()[0]->Ms;
                     node1_xPH( 0, 0 ) = node1_mPS( 0, 0 )*1000./18.0153;  // assuming it consists of H2O
                     node1_xPA( 0, 0 ) = node1_mPS( 0, 0 )*1000./18.0153;
        break;
     case MGP_TT_SOLID:  // '4'  TBD
                    /*naqgf = 0;
                    if( (x_aq + x_gf) == 1 )
                        naqgf = 1;*/

        break;
     case MGP_TT_AQGF:   // '3'  TBD
                     for(i=0; i<mtp->Nf; i++ )
                        node1_bPS( 0, x_aq, i ) = node1_bIC( 0, i )/2.;
                     for(i=0; i<mtp->Nf; i++ )
                        node1_bPS( 0, x_gf, i ) = node1_bIC( 0, i )/2.;
        break;
     case MGP_TT_GASF:   // '2'  TBD
                    for(i=0; i<mtp->Nf; i++ )
                        node1_bPS( 0, x_gf, i ) = node1_bIC( 0, i );
        break;
     default: break;
    }

    if( !na->CalcIPM_One( TestModeGEMParam(mode, mtp->PsSIA, mtp->ct, mtp->cdv, mtp->cez ), 0, 0 ) )
      iRet = false;  // Analysis of errors after GEM calculation?
    // Calculation of current box 0 reactive IC masses in kg
    BoxMasses( 0 );
    // Calculation of MGP bulk compositions in box 0 (in moles of ICs)
    ComposMGPinBox( 0 );

//sprintf(buf, "node_0000_wave_0000_noeq_dbr.dat" );
//na->CopyWorkNodeFromArray( 0, mtp->nC,  na->pNodT1() );
//na->GEM_write_dbr( buf, 0, false );

  //  This loop contains the overall transport time step (wave)
  do {

#ifndef IPMGEMPLUGIN

        if( mtp->ct > 0)
          CalcStartScript();

      sprintf(buf, "   step %ld; time %lg; dtime %lg  ", mtp->ct, mtp->cTau, mtp->dTau );
      Vmessage = "Simulating Transport through Sequential Reactors chain: ";
      Vmessage += buf;
      Vmessage += ". Please, wait (may take time)...";

      if( mtp->PsSmode != S_OFF  )
      {
          STEP_POINT2();
      }
      else
      {
          iRet = pVisor->Message( window(), GetName(),Vmessage.c_str(),
                                  mtp->ct, mtp->ntM, UseGraphMonitoring );
      }

      if( iRet )
             return iRet;// Error("GEM2MT SeqReac model", "Cancel by the user");
#endif

      // Set up new boxes states at cTau
     for( p = 1/*0*/; p < mtp->nC; p++ )
     {
         //if( p == mtp->nC-1 )
         //    lastp = p;  // for debugging
         for(i =0; i< mtp->Nf; i++ )
           dMb( p, i) = 0.;
         for(kk=0; kk < mtp->nFD-1; kk++ )  // Looking through the list of fluxes
         {
             if( mtp->FDLi[kk][1] == p ) // this flux comes into this box
             {
               if( mtp->ct > 0 ) // don't remove anything at 0-th wave!
               {
                  // Removing all this MGP from the p box first
                  long int oldp = mtp->FDLi[kk+1][1];
                  mtp->FDLi[kk+1][1] = -1;
                  dMBflux( kk+1, mtp->MB, mtp->dMB, mtp->cTau );
                  mtp->FDLi[kk+1][1] = oldp;
                  // mtp->FDLf[kk+1][0] = oldo;
               }
               // Adding the incoming MGP to p-th box
               dMBflux( kk, mtp->MB, mtp->dMB, mtp->cTau );
             }
         }   // kk

         // change bulk composition in the box p
         BoxComposUpdate( p );

//sprintf(buf, "node_%4.4d_wave_%4.4d_in_dbr.dat", p, mtp->ct );
//na->CopyWorkNodeFromArray( p, mtp->nC,  na->pNodT1() );
//na->GEM_write_dbr( buf,  0, false);

         // calculate equilibrium state in q-th box
         node1_Tm( p ) = mtp->cTau;
         node1_dt( p ) = mtp->dTau;
         if( !na->CalcIPM_One( TestModeGEMParam(mode, mtp->PsSIA, mtp->ct, mtp->cdv, mtp->cez ), p, 0 ) )
           iRet = false;  // Analysis of errors after GEM calculation?
         mtp->qc = p;
            //  iRet = false;
//         if( iRet == true )  // Error in GEM calculation
//             break;
//sprintf(buf, "node_%4.4d_wave_%4.4d_out_dbr.dat", p, mtp->ct );
//na->CopyWorkNodeFromArray( p, mtp->nC,  na->pNodT1() );
//na->GEM_write_dbr( buf,  0, false);
         // Calculation of current box reactive IC masses in kg
         BoxMasses( p );
         // Calculation of MGP bulk compositions in boxes (in moles of ICs)
         ComposMGPinBox( p );
       } // p

#ifndef IPMGEMPLUGIN
    // time step accepted - Copying nodes from C1 to C0 row
    pVisor->Update();
    CalcGraph();
#endif
       // copy node array T1 into node array T0
       copyNodeArrays();

        mtp->ct += 1;
        mtp->oTau = mtp->cTau;
        mtp->cTau += mtp->dTau;

        if( mtp->PsVTK != S_OFF )
            outp_time += PrintPoint( 0 );

   } while ( mtp->cTau < mtp->Tau[STOP_] && mtp->ct < mtp->ntM );

#ifndef IPMGEMPLUGIN
    pVisor->CloseMessage();
#endif

  }
  catch( TError& xcpt )
  {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), xcpt.title, xcpt.mess);
#else
       std::cerr << xcpt.title.c_str() << "  " <<  xcpt.mess.c_str() << std::endl;
#endif
       return 1;
  }
/*  catch(...)
  {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), "CalcSeqReacModel", "Unknown exception");
#else
       cerr << "CalcSeqReacModel Unknown exception" << endl;
#endif
      return -1;
  }*/
  return 0;
}

#undef dMb
#undef Mb

//Calculate generic box megasystem of 'B' type with MGP fluxes
//   mode - SIA or AIA
//
bool TGEM2MT::CalcBoxFluxModel( char /*mode*/ )
{
  try{
	 
    BoxFluxTransportStart();


#ifndef IPMGEMPLUGIN
    bool iRet = false;
    bool UseGraphMonitoring = false;

     if( mtp->PsSmode == S_OFF )
          if(  mtp->PvMSg != S_OFF && vfQuestion(window(),
             GetName(), "Use graphic monitoring?") )
        {
            RecordPlot( nullptr );
            UseGraphMonitoring = true;
        }
#endif

  // calculate inital states
  // Calculation of chemical equilibria in all nodes at the beginning
  // with the AIA initial approximation
  BoxEqStatesUpdate( -1,  0, mtp->cTau, mtp->dTau );

  //16/11/2011 we call this function in end of BoxEqStatesUpdate() ???
  //CalcMGPdata(); // Initial calculation of masses and MGP compositions
  // time iteration part
   nfcn = nstep = naccept = nrejct = 0;
   
#ifndef IPMGEMPLUGIN
       iRet = pVisor->Message( window(), GetName(),
           "Simulating Reactive Transport in a Box-Flux setup. "
           "Please, wait (may take time)...", nstep, mtp->ntM, UseGraphMonitoring );
       if( iRet )
        Error("GEM2MT generic box-flux model", "Cancelled by the user");
#endif
  INTEG( 1e-3, /*mtp->cdv,*/ mtp->dTau, mtp->Tau[START_], mtp->Tau[STOP_] );

#ifndef IPMGEMPLUGIN
    pVisor->CloseMessage();
#endif
    
  }
  catch( TError& xcpt )
   {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), xcpt.title, xcpt.mess);
#else
       std::cerr << xcpt.title.c_str() << "  " <<  xcpt.mess.c_str() << std::endl;
#endif
       return 1;
   }
/*  catch(...)
  {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), "CalcBoxFluxModel", "Unknown exception");
#else
       cerr << "CalcBoxFluxModel Unknown exception" << endl;
#endif
      return -1;
  }*/
  return 0;
}

//--------------------------------------------------------------------
// Integration process

//const long int NMAX = 800;
const long int KM = 8;
const double UROUND = 1.73e-18;
const double FAC1 = 2.e-2;
const double FAC2 = 4.0;
const double FAC3 = .9;
const double FAC4 = .8;
const double SAFE1 = .65;
const double SAFE2 = .94;
const double MAXSTEP = 1.7;
//const int MAXINTEGEXPR = 200;

//double *x;
//double *dx;
//double *tv;

//static double tt[9][ MAXINTEGEXPR ]; change 14/12/2007 [ MAXINTEGEXPR ][9]
static double hh[9];
static double w[9];
static double err;
static double epsd4;
static long int    nj[9]={2,4,6,8,10,12,14,16,18};
static double a1[9]={3.0,7.0,13.0,21.0,31.0,43.0,57.0,73.0,91.0};


// internal point j calculation
void TGEM2MT::MIDEX( long int j, double t, double h )
{
    double *z1=0, *z2=0, *dz=0;
    double hj, scal, fac, facmin, expo, ys, v1, v2;
    long int i,m,mm,l,n;
    try
    {
        n = mtp->nC*mtp->Nf;
        z1 = new double[n];
        z2 = new double[n];
        dz = new double[n];
        memset( z1, 0, sizeof(double)*n );
        memset( z2, 0, sizeof(double)*n );
        memset( dz, 0, sizeof(double)*n );

        hj = h / (double)nj[ j ];
        // 
        for( i=0; i < n; i++ )
        {
            z1[ i ] = x[ i ];
            z2[ i ] = x[ i ] + hj * dx[ i ];
        }
        // 
        m = nj[ j ] - 1;
        for( mm=0; mm < m; mm++ )
        {
            Solut(  z2, dz, t+hj*(double)(mm+1) );
            for( i=0; i < n; i++ )
            {
                ys = z1[ i ];
                z1[ i ] = z2[ i ];
                z2[ i ] = ys + hj * dz[ i ] * 2.;
            }
        }
        // 
        Solut(  z2, dz, t+h );
        for( i=0; i < n; i++ )
            mtp->tt[ i ][ j ] = ( z1[ i ] + z2[ i ] + hj * dz[ i ] ) / 2.;
        nfcn += nj[ j ];
        //
        if(j == 0)
            goto VEL;
        for( l=j; l >= 1; l-- )
        {
            fac = pow( (double)nj[ j ] / (double)nj[ l-1 ], 2.) - 1.;
            for( i=0; i < n; i++ )
                mtp->tt[ i ][ l-1 ] = mtp->tt[ i ][ l ] +
                                 ( mtp->tt[ i ][ l ] - mtp->tt[ i ][ l-1 ] ) / fac;
        }
        err = 0e0;
        for( i=0; i < n; i++ )
        { // 
            v2 = fabs( mtp->tt[ i ][ 0 ] );
            v1 = fabs( x[ i ] );
            v1 = std::max( v1, v2 );
            v2 = std::max( 1.e-6, UROUND / epsd4 );
            scal = std::max( v1, v2 );
            err += pow( (mtp->tt[ i ][ 0 ] - mtp->tt[ i ][ 1 ] ) / scal, 2. );
        }
        err = pow( err / (double)n, .5 );
        //
        expo = 1.e0 / (double)( 2 * ( j + 1 )-1 );
        facmin = pow( FAC1, expo );
        v1 = std::max( facmin, pow( err/epsd4, expo ) / SAFE2 );
        fac = std::min( FAC2/facmin, v1 );
        fac= 1.e0 / fac;
        hh[ j ] = std::min( h * fac, MAXSTEP );
        w[ j ] = a1[ j ] / hh[ j ];

VEL:
        delete[] z1;
        delete[] z2;
        delete[] dz;
    }
    catch( ... )
    {
        if( z1 ) delete[] z1;
        if( z2 ) delete[] z2;
        if( dz ) delete[] dz;
        Error( "INTEG ", "Error in MIDEX!");
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Subroutine INTEG - numerical integration of the system of n ordinary 
//    differential equations of the form 
//		 dxi/dt = f( x1,x2, ... ,xn,t )  i=1,...,n
//               xi( t0 )=xi0
//  over the time interval [t_begin, t_end] with minimal time step length 
//     step and precision eps 
//  Uses implicit method of rational extrapolation (??)
// 
void
TGEM2MT::INTEG( double eps, double& step, double t_begin, double t_end )
{
    double  t, h1, h, v1;
    long int     j,i,reject,last,kc=0,kopt,k;

    x = mtp->MB;
    dx = mtp->dMB;
    tv = &mtp->cTau;

    h = step;
    epsd4 = eps * SAFE1;
    v1 = -log10(eps)*.6+1.5;
    v1 = std::min( 8.0, v1 );
    k = std::max( 3.0, v1 ) - 1;
    t = t_begin;
    h1 = t_end-t;
    v1 = std::min( MAXSTEP, h1/2. );
    h = std::min( h, v1 );
    //BoxEqStatesUpdate( 0, k, t, h ); // 14/12/2007 ????? may be done before in calc
    err = w[ 0 ] = 0.0;
    reject = last = 0;   // false

    //
    while( fabs( h1 ) >= UROUND )
    {
        v1 = std::min( h1, MAXSTEP);
        h = std::min( h, v1 );
        if( h >= ( h1 - UROUND ) )  last = 1;      // true
        Solut(  x, dx, t );
        nfcn++;
        if (( nstep == 0 )||( last ))     // 1
        {
            nstep++;
            for( j=0; j <= k; j++ )
            {
                kc=j;
                MIDEX( j, t, h );
                if( ( j > 0 ) && ( err <= eps ) ) goto l60;
            }
            goto l55;
        }
        //
l30:
        nstep++;
        ErrorIf( nstep >= mtp->ntM /*MaxIter*/, "INTEG", "No convergence - too many iterations" ); // 14/12/2007 !!!!
        kc = k-1;
        for( j=0; j <= kc; j++ )
            MIDEX( j, t, h);
        // 
        if( !( k == 1 || reject ) )
        {
            if( err <= eps ) goto l60;
            if( ( err/eps ) >
                    ( pow( (double)(nj[k+1]*nj[k])/4., 2. ) ) )  goto l100;
        }
        MIDEX( k, t, h );
        kc = k;
        if( err <= eps ) goto l60;
        // 
l55:
        if( err/eps > pow( (double)(nj[k])/2.0, 2.0 ) ) goto l100;
        kc = k + 1;
        MIDEX( kc, t, h );
        if( err > eps ) goto l100;
        // 
l60:
        t += h;
        step = h;
//        mtp->cTau = t;
        for( i=0; i < mtp->nC*mtp->Nf; i++ )
            x[i] = mtp->tt[i][0];
        Solut(  x, dx, t );
        naccept++;
        mtp->ct++;
        BoxEqStatesUpdate( naccept, kc, t, h );
        if( mtp->PsVTK != S_OFF )
           PrintPoint( 0 );
        //
        if( kc == 1 )
        {
            kopt = 2;
            if( reject ) kopt = 1;
        }
        else if( kc <= k )
        {
            kopt = kc;
            if( w[kc-1] < w[kc]*FAC3 ) kopt = kc - 1;
            if( w[kc] < w[kc-1]*FAC3 )
                kopt = std::min( (kc+1) , (KM-1) );
        }
        else
        {
            kopt = kc-1;
            if( kc > 2 && w[kc-2] < w[kc-1]*FAC3 )
                kopt = kc - 2;
            if( w[kc] < w[kopt]*FAC3 )
                kopt = std::min( kc, (KM-1) );
        }
        // 
        if( reject )
        {
            k = std::min( kopt, kc );
            h = std::min( h, hh[ k ] );
            reject = 0;           // false
        }
        else
        {  // 
            if( kopt <= kc )
                h = hh[ kopt ];
            else if( kc < k && ( w[kc] < ( w[kc-1] * FAC4 ) ))
                h = hh[kc] * a1[kopt+1] / a1[kc];
            else h = hh[kc] * a1[kopt] / a1[kc];
            k = kopt;
        }
        h1=t_end-t;
    }     // end while
    *tv = t;
    step = h;
    return;
l100:
    k = std::min( k, kc );
    if( k > 1 && ( w[k-1] < w[k] * FAC3 ) )
        k--;
    nrejct++;
    h = hh[k];
    reject = 1;
    goto l30;
}

// --------------------- end of m_gem2mtbox.cpp ---------------------------


