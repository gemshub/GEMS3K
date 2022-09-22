//-------------------------------------------------------------------
// $Id: m_gem2mtbox.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of TInteg/TGEM2MT classes, calculation functions
//
// Rewritten from C to C++ by S.Dmytriyeva  
// Revised in 2019, 2021 by D.Kulik
// Copyright (C) 1995,2008,2011,2019,2021  S. Dmytriyeva, D.Kulik
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

#include "m_gem2mt.h"
#include "GEMS3K/nodearray.h"
#include "v_service.h"

#define dMB( q, i) ( dm[ (q)*mtp->Nf + (i)] )
#define MB( q, i)  ( m[ (q)*mtp->Nf + (i)] )

#define g(q,f,i)   ( mtp->gfc[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )
#define y(q,f,i)   ( mtp->yfb[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )

#define bmgpm(q,f)   ( mtp->BmgpM[ (q)*(mtp->nPG)+ (f)] )

#define ord(kk)    ( ROUND(mtp->FDLf[kk][0]) )
// #define v(kk)      ( (mtp->FDLf[kk][1]) )
#define v(kk)      ( (mtp->FDLf[kk][1])/mtp->tf )  // time scaling - need to test

#define H(f, i)    ( mtp->BSF[ (f) * mtp->Nf + ( i )] )


// Returns mass (kg) of IC i in H (BSF) fe row
//
double TGEM2MT::MassICinHfe( long int fe, long int i )
{
    double MassIC=0.;
    if(fe >= 0)
       MassIC = H(fe,i) * na->ICmm( i );  // na->ICmm(i) is already in kg/mol!
    return MassIC;
}

// Returns mass (kg) of the production/consumption fe row
//
double TGEM2MT::MassHfe(long int fe)
{
    double Mass = 0.0;
    if ( fe < 0 )
       Mass = 0.0;
    else {
       for(long int i=0; i<mtp->Nf; i++ )
          Mass += H(fe,i) * na->ICmm( i );  // in kg
    }
// fstream f_log("ipmlog.txt", ios::out|ios::app );
// f_log << " MGP= " << Mass << " q= " << q << " f= " << f << " ";
    return Mass;
}

// Returns mass (kg) of IC i in MGP f in box q
//
double TGEM2MT::MassICinMGP( long int q, long int f, long int i )
{
    double  MassIC=0.,
            yqfi = y(q,f,i),
            icmm = na->ICmm( i );   // in kg/mol
    if(f >= 0)
       MassIC = yqfi * icmm;
    return MassIC;
}

// Returns mass of f-th MGP in q-th box (in kg)
double TGEM2MT::MassMGP( long int q, long int f  )
{
    double Mass = 0.0;
    if ( f >= 0 )
    {
       for(long int i=0; i<mtp->Nf; i++ )
          Mass += MassICinMGP( q, f, i );  // to test
          // Mass += y(q,f,i) * na->ICmm( i );  // in kg
    }
// fstream f_log("ipmlog.txt", ios::out|ios::app );
// f_log << " MGP= " << Mass << " q= " << q << " f= " << f << " ";
    return Mass;
}

// Cleaning up the table of IC mass differentials in boxes/nodes
//
void TGEM2MT::dMBZeroOff(  double *dm )
{
  long int q, i;
  // Zeroing IC mass differentials off
  for( q=0; q <mtp->nC; q++ )
        for(i =0; i< mtp->Nf; i++ )
            dMB(q,i) = 0.;
}

//   Increment to mass differential for i-th IC in q-th box using the IC flux rate fRate
//      fRate: i-th element flux rate in kg/s
//      dm[]:  table of derivatives dm_i/dt in all boxes (in kg/s), output
//      sign: -1 outgoing flux     +1 incoming flux
// Returns increment to mass differential
//
double TGEM2MT::dMBfluxDir( long int q, long int i, double *dm, double fRate, double sign )
{
    switch(mtp->DiCp[q][1])
    {                  // boundary condition node
       case NBC3source:  // 3: Cauchy source ( constant flux )
                 if( sign < 0 )
                      break;
                 [[fallthrough]];
       case NBC1source:  //  1: Dirichlet source ( constant concentration )
       case NBC1sink:    // -1: Dirichlet sink
       case NBC2source:  //  2: Neumann source ( constant gradient )
       case NBC2sink:    // -2: Neumann sink
       case NBC3sink:    // -3: Cauchy sink
       case INIT_FUNK:   //  4: functional conditions (e.g. input time-depended functions)
       case normal: // normal node  // Checking phase assemblage
       default:    dMB(q,i) += fRate*sign;
                   break;
    }
//    double dmbqi = dMB(q,i);    // for debugging
    return dMB(q,i);
}

// Calculate kk-th flux between boxes with indexes in FDLi[kk] in kg/time
//   m  ( m(q,i) ) is the table of masses of element i in box q
//   dm (dm(q,i) ) is the table of time derivatives of masses of element i in box q
//   t is the current time in s (for debugging)
//   Checked and revised by DK on 30.08.2021
//
void TGEM2MT::dMBflux( long int kk, double *dm )
{
  long int  q, p, f, fe=-1; int i;
  char FLXid[MAXSYMB+1], MGPid[MAXSYMB+1];
  long int FLXorder/*, FLXtype*/;
  double fiRate, fiStepM = 0.0, MGPfq=0., MGPfp=0., vkk=1., dTau = mtp->dTau/mtp->tf;
  // mqfi, nbIC, MBo, Hfei, MBqi, MBpi;
  bool sinkOut=false;

  fe = -1;
  q = mtp->FDLi[kk][0];             // index of outgoing box
  p = mtp->FDLi[kk][1];             // index of receiving box
  double MBq = node1_Ms( q );       // Mass of box q
  // double MBp = node1_Ms( p );           // Mass of box p

  strncpy( MGPid, mtp->FDLmp[kk], MAXSYMB );
  MGPid[MAXSYMB] = 0;                           // MGP identifier
  f = LookUpXMGP( MGPid );                      // Getting MGP index

  if( f < 0 ) // is this an elemental production/sink flux BSF?
  {
        sscanf( MGPid, "%ld", &fe );  // Reading the BSF row index
        if( fe < 0 || fe >= mtp->nSFD )
           Error( "BOXFLUX:", "Wrong MGP identifier or BSF row index!" );
  }

  strncpy( FLXid, mtp->FDLid[kk], MAXSYMB );
  FLXid[MAXSYMB] = 0;               // Flux identifier
  FLXorder = ord(kk);               // flux order
  vkk = v(kk);                      // flux rate constant
  bmgpm(q,f) = MGPfq = MassMGP( q, f  ); // calculating MGP masses: in outgoing box (kg)
  if( p >= 0 )
      MGPfp = MassMGP( p, f  );     // MGP mass in incoming box (kg)
  mtp->FmgpJ[kk] = 0.0;  // cleaning before summation of total MGP mass flux over elements

  if( q >= 0 && f >= 0 )
  {                                 // Normal MGP flux from box q to box p
      // NB: Negative v(kk) means "production" in q box and "consumption" in p box
      if( p < 0 )
      {
          sinkOut = true;       // This is a sinkout flux from box q to nowhere (if v > 0)
                                //    or production of MGP in box q (if v < 0)
          MGPfp = MGPfq;        // conditional assumption
      }
      if(MGPfq >= 1e-19)
      {  // This MGP has a significant mass
        for(i=0; i<mtp->Nf; i++ )
        {
         double icmm = na->ICmm( i );       // mol mass of element i, kg/mol
         double yqfi = y(q,f,i);            // moles of element i in MGP f in box q
         double mqfi = yqfi * icmm;         // mass of element i in mgp f in box q  kg
         double fqfi = 0.0;
         if(MGPfq > 1e-19)
            fqfi = mqfi/MGPfq;          // mass fraction of element i relative to MGP f mass in box q
         fiRate = 0.0;
// For debugging:
//         double icamqi = node1_bIC(q, i);   // amount of element i in box q
//         double icampi = node1_bIC(p, i);   // amount of element i in box p
//         double mqi = icmm * icamqi;        // mass of element i in box q   kg
//         double mpi = icmm * icampi;        // mass of element i in box q   kg
//         double ypfi = y(p,f,i);            // moles of element i in MGP f in box p
//         double gqfi = g(q,f,i);            // partition coeff of element i in MGP in box q
//         double gpfi = g(p,f,i);            // partition coeff of element i in MGP in box q
//         double mpfi = ypfi * icmm;         // mass of element i in mgp f in box p  kg
//         double fpfi = mpfi/MGPfp;          // mass fraction of element i relative to MGP f mass in box q

         switch( FLXorder )
         {
            case 0:  // Zero-order (absolute) flux of element i;  v(kk) is MGP flux rate in kg/s units
               fiRate = fqfi * vkk;                 // Flux of IC i in MGP f from box q in kg/s
               break;

            case 1:  // First-order to source box MGP mass flux;  v(kk) in 1/s units
               fiRate = fqfi * vkk * MGPfq;         // Flux of IC i in MGP f from box q in kg/s
               break;

            case 2:  // Second-order to source and sink MGP mass flux:  vkk in 1/s units
                         // (proportional to the SUM of MGP masses in both boxes)
               fiRate = fqfi * vkk * (MGPfq + MGPfp); // Flux of IC i in MGP f from box q in kg/s
               break;

            case 3:  // was -1 First-order flux to MGP mass in receiving box in 1/s units
               fiRate = fqfi * vkk * MGPfp;         // Flux of IC i in MGP f from box q in kg/s
               break;

            default:  // Other orders or flux types - TBA here!
               break;
        }

        // Convention at positive rate constant: dm/dt < 0 is outgoing, dm/dt > 0 is incoming flux
        fiStepM += dTau * dMBfluxDir( q, i, dm, fiRate, -1.0);  // outgoing flux of element i in box q over time step dt
        mtp->FmgpJ[kk] += fiRate*(-1.0); // Only outgoing fluxes are counted to show the user
        if( !sinkOut )
           dMBfluxDir( p, i, dm, fiRate, 1.0);     // incoming flux of element i into box p
        ;
       } // i
      }
    }
   else {
      if( q >= 0 && f < 0 && fe >= 0)
      {         // This is an elemental flux ( fe row in BSF table )
                // NB: Negative v(f) means "production" in q box and "consumption" in p box
          for(i=0; i<mtp->Nf; i++ )
          {                                  //  or elemental production flux in box q (if v < 0)
            fiRate = 0.0;
            double ghfei = MassICinHfe( fe, i )/ MassHfe( fe ); // mass partition coefficient of IC i in BS(fe) row

            switch( FLXorder )
            {
                case 0:     // Zero-order flux      vkk in kg/s
                   fiRate = ghfei * vkk;        // Flux of IC i in BSF fe from/to box q in kg/s
                   break;

                case 1:     // was -1 First-order flux to mass in box q in 1/s units
                   fiRate = ghfei * vkk * MBq;  // Flux of IC i in BSF fe from/to box q in kg/s
                   break;

                default:  // Other orders or flux types - nothing to be done yet!
                   break;
           }
           mtp->FmgpJ[kk] += dMBfluxDir( q, i, dm, fiRate, 1.0 );
          } // i
        }
    }
}

// calculate transport step at time cTau using the table of fluxes
void TGEM2MT::Solut( double *m, double *dm, double tau )
{ 
  // Zeroing IC mass differentials off
  dMBZeroOff(  dm );

  for(long int kk=0; kk<mtp->nFD; kk++ )  // Looking through the list of fluxes
  {
     dMBflux( kk, dm );
  }
}

#undef dMB
#undef MB

#define Mb( q, i)  ( mtp->MB[(q)*mtp->Nf + (i)])
#define dMb( q, i)  (mtp->dMB[(q)*mtp->Nf + (i)])

// change bulk composition in the box q using changes dM*dTau over time step
//
void TGEM2MT::BoxComposUpdate( long int q )
{
   long int i;
   double  dMbqi, ICmm, dTau;
   // double n1bICo, n1bICn=0.;       // for debugging
   dTau =  mtp->dTau/mtp->tf;         // time differential in s (time step), rescaled
   mtp->BdM[q] = 0.0;
   for(i =0; i< mtp->Nf; i++ )
   {
//       n1bICo = node1_bIC(q, i);
       dMbqi = dMb( q, i);          // in kg/s
       mtp->BdM[q]  += dMbqi;       // Increment to total mass differential in the box
       ICmm = na->ICmm( i );        // in kg

       if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_AQUEL && i == mtp->Nf-1 )
       {
         node1_bIC(q, i) = 0.0;  // Provisorial - zeroing-off charge for the node
       }
       else {
           // normal element or ligand - dTau is the time step (in s)
         node1_bIC(q, i) += dMbqi/ICmm * dTau;
//         n1bICn = node1_bIC(q, i);
         if( node1_bIC(q, i) < 1e-12 )
            node1_bIC(q, i) = 1e-12;        // preventing too small or negative amounts of ICs
         Mb(q,i) = node1_bIC(q, i) * ICmm;  // Refreshing mass of IC i in box q
       }
   } // i
}

// Update all box compositions at time cTau using time step dTau
//
void TGEM2MT::BoxesBCupdate()
{
  // Debug trace output
  TNode::ipmlog_file->trace("Time t= {}   dt= {}", mtp->cTau/mtp->tf, mtp->dTau/mtp->tf);

  for( long int q=0; q <mtp->nC; q++ )
  {
    BoxComposUpdate( q );

    // Debug trace output
    double dMsq = 0.0;
    TNode::ipmlog_file->trace("box qc=  {}   mass= {}", q, node1_Ms( q ));
    for(long int i =0; i< mtp->Nf; i++ )
    {
       dMsq += dMb(q,i);
    }  // i
    mtp->BdM[q] = dMsq; // Alternatively calculated mass differential of the box
    TNode::ipmlog_file->trace("box qc=  {}   mass= {}", q, dMsq);
  }  // q
}

// Calculation of q-th box reactive IC masses in kg
// returns a total mass of the box in kg
//
double TGEM2MT::BoxMasses( long int q )
{
   long int i;
   double BoxMass = 0.;
//   double Mbqi = 0.;          // debugging
   for(i=0; i<mtp->Nf; i++ )
   {
       Mb( q, i ) = node1_bIC( q, i ) * na->ICmm( i );   // in kg
//       Mbqi = Mb(q,i);        // debugging
       BoxMass += Mb(q,i);
   }
   return BoxMass;
}

// Calculates a table of MGP bulk compositions for q-th box (now added phase conc. units)
//   This variant is for simplified box-flux and sequential reactors transport models
// Differs from B flux module in interpretation of the composition of MGP:
//     1 'M' or 1 'G' means the whole amount (of aq, gas or solid phase) for S models
//     and just 1 mole or 1 kg of phase for F and B models
//     1 'n' means the whole current amount of any phase for S, F and B models
//
void TGEM2MT::ComposMGPinBox( long int q  )
{
    double PhFract=0., Xincr=1., Xe, PHmw, Vm=1., R1=1., Msys=1.,
            Mwat=1., Vaq=1., Maq=1., Vsys=1.;
    // double Mph=0., Mgas=0.;  // debugging
    char UNITP = 'n';                   // mole fraction (of phase into MGP)
    long int f, i, k, x_aq=-1, x_gf=-1;
    // long int naqgf=0;          // debugging

    // Cleanind y(q,f,i) tables for this q-th box
    for(f=0; f<mtp->nPG; f++ )
    {
        for(i=0; i<mtp->Nf; i++ )
        {
           y(q,f,i) = 0.0;
           g(q,f,i) = 0.0;
        }
    }
    // Index of aq phase in DBR
    if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_AQUEL )
        x_aq = 0;
    else
        x_aq = -1;

    if( x_aq == 0 )
    {
        if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_GASMIX ||
            na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 1 )] == PH_FLUID )
           x_gf = 1;
    }
    else
    {
        if( na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_GASMIX ||
             na->pCSD()->ccPH[na->Ph_xDB_to_xCH( 0 )] == PH_FLUID )
        x_gf = 0;  // Index of gas phase
    }

    for(f=0; f<mtp->nPG; f++ )
    {                               // Loop over MGP
      for( k=0; k<mtp->FIf; k++)
      {                             // Loop over phases
          Xe = mtp->PGT[k*mtp->nPG+f];  // given quantity of phase to add to MGP
          if( fabs(Xe)<1e-19 )
              continue;
          R1 = node1_xPH( q, k );  // amount of phase in moles:: mole fraction is related to amount
          if( fabs(R1)<1e-19 )
              continue;
          UNITP = mtp->UMGP[k];    // one list for all MGPs (so far)
          if( k < na->pCSD()->nPSb )
          {                         // Phase-solution
             PHmw = node1_mPS( q, k ) * 1e3 / R1;   // mass of solution phase, g
             Vm = node1_vPS( q, k ) *1e6 / R1;      // molar volume of solution phase cm3/mol
//             Mph = node1_mPS( q, k );          // mass of phase in kg;
          }
          else {                    // pure phase
              PHmw = node1_mPH( q, k ) * 1e3 / R1;  // g/mol
              Vm = node1_vPH( q, k ) * 1e6 / R1;    // cm3/mol
//              Mph = node1_mPH( q, k );           // mass of phase in kg;
          }

          if( k == 0 && x_aq == 0 )
          {                         // This is aqueous phase
             Mwat = node1_xPA( q, k )*18.0153/1000.;    // mass of water in kg
             Vaq = node1_vPS( q, k ) * 1e3;             // in dm3 (L)
             Maq = node1_mPS( q, k );                   // in kg
             Msys = Maq;   // For aqueous phase, its total mass, volume are taken instead of that of the system
             Vsys = Vaq;   //  (provisionally)
          }
          else {            // phases other that aqueous - normalizing by mass of the whole system
            Msys = node1_Ms( q );       // in kg
            Vsys = node1_Vs( q )*1e3;   // in L
          }

          // Converting quantity Xe into amount (moles)
          Xincr = Reduce_Conc( UNITP, Xe, PHmw, Vm, R1, Msys,
                  Mwat, Vaq, Maq, Vsys );

          PhFract = 0.; // Initial fraction of real phase amount

          if( x_aq >= 0 && k == x_aq )  //    Aqueous phase
          {  //    This is aqueous phase
             if( mtp->PsMode == RMT_MODE_S )
             {
                 PhFract = 1.0;  // For sequential reactors S mode, assuming the whole aq phase is taken into MGP
             }
             else if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL ||
                      UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM ||
                      UNITP == CON_MOLFR || UNITP == CON_MOLPROC || UNITP == CON_pMOLFR ||
                      UNITP == CON_WTFR || UNITP == CON_WTPROC || UNITP == CON_PPM ||
                      UNITP == CON_VOLFR || UNITP == CON_VOLPROC || UNITP == CON_pVOLFR )
             {
                PhFract = Xincr / R1;
             }
             else {
                 PhFract = 0.0; // Other units are not allowed here (for the phase in MGP)
                 // - issue an error message!
             }
          }

          if( x_gf >= 0 && k == x_gf )   // Gas/fluid phase
          {
             if( mtp->PsMode == RMT_MODE_S )
             {
                 PhFract = 1.0;  // For sequential reactors S mode, assuming the whole gas fluid phase is taken into MGP
             }
             else if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL ||
                      UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM ||
                      UNITP == CON_MOLFR || UNITP == CON_MOLPROC || UNITP == CON_pMOLFR ||
                      UNITP == CON_WTFR || UNITP == CON_WTPROC || UNITP == CON_PPM ||
                      UNITP == CON_VOLFR || UNITP == CON_VOLPROC || UNITP == CON_pVOLFR )
             {
                PhFract = Xincr / R1;
             }
             else {
                 PhFract = 0.0; // Other units are not allowed here (for the phase in MGP)
                 // - issue an error message!
             }
          }

          if( (x_gf >= 0 && k > x_gf) || (x_gf < 0 && x_aq >= 0 && k > x_aq) )    // Bugfix 1.11.21  DK
          {
              // Condensed phase other than aqueous or gas
              if( mtp->PsMode == RMT_MODE_S )
              {
                  PhFract = 0.0;  // For sequential reactors S mode, only aq and/or gas fluid phases are moving in MGP!
              }
              else if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL ||
                       UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM ||
                       UNITP == CON_MOLFR || UNITP == CON_MOLPROC || UNITP == CON_pMOLFR ||
                       UNITP == CON_WTFR || UNITP == CON_WTPROC || UNITP == CON_PPM ||
                       UNITP == CON_VOLFR || UNITP == CON_VOLPROC || UNITP == CON_pVOLFR ||
                       UNITP == CON_MOLAL || UNITP == CON_MMOLAL || UNITP == CON_pMOLAL ||
                       UNITP == CON_MOLAR || UNITP == CON_MMOLAR || UNITP == CON_pMOLAR ||
                       UNITP == CON_AQWFR || UNITP == CON_AQWPROC || UNITP == CON_AQPPM)
              {
                 PhFract = Xincr / R1;
              }
              else {
                  PhFract = 0.0; // Other units are not allowed here (for the phase in MGP)
                  // - issue an error message!
              }
/*
              else  // B or F mode
                  if( UNITP == QUAN_MOL || UNITP == QUAN_MMOL || UNITP == QUAN_MKMOL )
              {
                  // Mole amount of condensed phase used in MGP
                  // QUAN_MKMOL = 'Y',  QUAN_MMOL = 'h',  QUAN_MOL = 'M',  // NUMBER OF MOLES
                  PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if( UNITP == QUAN_KILO || UNITP == QUAN_GRAM || UNITP == QUAN_MGRAM )
              {
                  // Mass of condensed phase to be taken into MGP
                  // QUAN_MGRAM = 'y',  QUAN_GRAM = 'g',  QUAN_KILO = 'G', // MASS
                  PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if( UNITP == CON_MOLFR || UNITP == CON_MOLPROC || UNITP == CON_pMOLFR )
              {
                  // Mole fraction of condensed phase relative to its amount ( <= 1 ) used in MGP
                  // CON_MOLFR = 'n', CON_MOLPROC = 'N', CON_pMOLFR = 'f', // MOLE FRACTION
                  PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if(UNITP == CON_WTFR || UNITP == CON_WTPROC || UNITP == CON_PPM)
              {
                  // Mass fraction of condensed phase relative to mass of the system ( <= 1 ) used in MGP
                  // CON_WTFR  = 'w', CON_WTPROC =  '%', CON_PPM =    'P', // MASS FRACTION of aq phase used in MGP
                  PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if(UNITP == CON_VOLFR || UNITP == CON_VOLPROC || UNITP == CON_pVOLFR)
              {
                  // Volume fraction of condensed phase relative to volume of system ( <= 1 ) used in MGP
                  // CON_VOLFR = 'v', CON_VOLPROC = 'V', CON_pVOLFR = 'u', // VOLUME FRACTION
                  PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if( (x_aq == 0) && (UNITP == CON_MOLAL || UNITP == CON_MMOLAL || UNITP == CON_pMOLAL) )
              {
                 // Molality of condensed phase ( <= 1 ) relative to mass of water-solvent in aqueous phase used in MGP
                 // CON_MOLAL = 'm', CON_MMOLAL =  'i', CON_pMOLAL = 'p', // MOLALITY
                 PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if( (x_aq == 0) && (UNITP == CON_MOLAR || UNITP == CON_MMOLAR || UNITP == CON_pMOLAR) )
              {
                 // Molarity of condensed phase ( <= 1 ) relative to volume of aqueous phase used in MGP
                 // CON_MOLAR = 'L', CON_MMOLAR =  'j', CON_pMOLAR = 'q', // MOLARITY
                 PhFract = Xincr / R1; // node1_xPH( q, k );
              }
              else if( (x_aq == 0) && (UNITP == CON_AQWFR || UNITP == CON_AQWPROC || UNITP == CON_AQPPM) )
              {
                 // Concentration of condensed phase ( <= 1 ) per mass of aqueous phase used in MGP
                 // CON_AQWFR = 'C', CON_AQWPROC = 'c', CON_AQPPM =  'a', // CONCENTRATION
                 PhFract = Xincr / R1; // node1_xPH( q, k );
              }
*/
          }

          // Adding IC composition of phase k to MGP f and calculating bulk MGP IC partition coeffs
          if( fabs(PhFract) > 0.0 )
          {
//            double yqfi, n1bPHqki, n1bICqi, gqfi; // debugging
            if( k < na->pCSD()->nPSb )
            {    // phase - solutions
                 for(i=0; i<mtp->Nf; i++ )
                 {
                     if(node1_bPS( q, k, i ) < 1e-19)
                         continue;
                     // add PhFract part of amount of IC i in phase k to MGP f composition
                     y(q,f,i) += node1_bPS( q, k, i ) * PhFract;
                     // update partition coefficient of IC i in phase k taken to MGP f relative to bulk amount in system
                     g(q,f,i) = ( node1_bIC( q, i ) > 1e-19 ? y(q,f,i)/node1_bIC( q, i ) : 0.0);
                     if( g(q,f,i) > 1.0 )
                         g(q,f,i) = 1.0;
//   Alternative for debugging
//                    n1bPHqki = node1_bPS( q, k, i );    // amount of IC i in phase k
//                    if(n1bPHqki < 1e-19)
//                       continue;
//                    yqfi = y(q,f,i);
//                    yqfi += n1bPHqki * PhFract;  // add part of the phase transferred to MGP
//                    y(q,f,i) = yqfi;
//                    gqfi = g(q,f,i); n1bICqi = node1_bIC( q, i );
//                    gqfi = ( n1bICqi > 1e-19 ? y(q,f,i)/n1bICqi : 0.);
                   // g(q,f,i) = ( node1_bIC( q, i ) > 0. ? y(q,f,i)*PhFract/node1_bIC( q, i ): 0.);
//                     if( gqfi > 1.0 )
//                        gqfi = 1.0;
//                    g(q,f,i) = gqfi;
                 } // i
            }
            else
            {  // pure phase
                 for(i=0; i<mtp->Nf; i++ )
                 {
                     if(node1_bPH( q, k, i ) < 1e-19)
                         continue;
                     // add PhFract part of amount of IC i in phase k to MGP f composition
                    y(q,f,i) += node1_bPH( q, k, i ) * PhFract;
                    // update partition coefficient of IC i in phase k taken to MGP f relative to bulk amount in system
                    g(q,f,i) = ( node1_bIC( q, i ) > 1e-19 ? y(q,f,i)/node1_bIC( q, i ): 0.);
                    if( g(q,f,i) > 1.0 )
                        g(q,f,i) = 1.0;
//   Alternative for debugging
//                    n1bPHqki = node1_bPH( q, k, i );
//                    if(n1bPHqki < 1e-19)
//                       continue;
//                    yqfi = y(q,f,i);
//                    yqfi += n1bPHqki * PhFract;
//                    y(q,f,i) = yqfi;
//                    gqfi = g(q,f,i); n1bICqi = node1_bIC( q, i );
//                    gqfi = ( n1bICqi > 1e-19 ? y(q,f,i)/n1bICqi: 0.);
//                    if( gqfi > 1.0 )
//                       gqfi = 1.0;
//                    g(q,f,i) = gqfi;
                 }  // i
            }
            if( x_aq == 0 && k == 0 )
            {
                y(q,f,mtp->Nf-1) = 0.0;    // provisional - zeroing-off charge in MGP, may change in future
                g(q,f,mtp->Nf-1) = 0.0;
            }
          }
       } // k
    }// f
}

// Calculation of box masses, MGP compositions and MGP distribution coefficients in all boxes
//
void TGEM2MT::CalcMGPdata()
{
  long int q;
    // Calculation of current box reactive IC masses in kg
       for( q=0; q <mtp->nC; q++ )
       {
           mtp->BM[q] = BoxMasses( q );
       }
    // Calculation of MGP bulk compositions in boxes (in moles of ICs)
       for( q=0; q <mtp->nC; q++ )
       {
           ComposMGPinBox( q );
       }
    // Calculation of MGP masses in boxes (in kg)

 }

// Calculate new equilibrium states in the boxes for tcur = t
//  Ni   number of independent components
//  pr
//  tcur - current time
//  step - current time step
//
bool TGEM2MT::BoxEqStatesUpdate(  long int Ni, long int /*pr*/, double tcur, double step )
{
  bool iRet = true;
  FILE* diffile = nullptr;

    mtp->dTau = step*mtp->tf;     // This is a doubtful circular dependence
    mtp->cTau = tcur*mtp->tf;
    mtp->oTau = mtp->cTau;

 if( Ni >= 0)
 { // Update bulk compositions in boxes at current time point
   BoxesBCupdate();
 } 
 // Calculate new box equilibrium states at current time tcur
 if( mtp->PsSIA != S_ON )
     CalcIPM( NEED_GEM_AIA, 0, mtp->nC-1, diffile );
 else
     CalcIPM( NEED_GEM_SIA, 0, mtp->nC-1, diffile );
 
 if( Ni >= 0 )
 {
   // Here one has to compare old and new equilibrium phase assemblage
   // and pH/pe in all nodes and decide if the time step was Ok or it
   // should be decreased. If so then the nodes from C0 should be
   // copied to C1 (to be implemented)

   // Output of the results if step accepted
   if( mtp->PsMO != S_OFF  )
      PrintPoint( 1 );
 }

   
  // copy node array for T0 into node array for T1
  mtp->oTau = mtp->cTau;
  copyNodeArrays();

  CalcMGPdata();  // Calculation of masses, MGP compositions and part coeff tables

  return iRet;
}

// Initialization of box-flux transport calculations
//
void TGEM2MT::BoxFluxTransportStart()
{
    if( mtp->tf <= 0 ) // foolproof - time step reduction factor (decimal)
        mtp->tf = 1.;
    mtp->dTau = mtp->Tau[STEP_];
    mtp->oTau =  mtp->Tau[START_];
    // mtp->cTau = mtp->Tau[START_];
//    long int ntm_old = mtp->ntM;
//    long int ntm_new = (mtp->Tau[STOP_]-mtp->Tau[START_])/mtp->Tau[STEP_];
//    if(ntm_old > 0)
    mtp->ntM = (mtp->Tau[STOP_]-mtp->Tau[START_])/(mtp->Tau[STEP_])+1;
    mtp->ct = 0;
    mtp->qf = 0;

    if( mtp->gfc )
        delete[] mtp->gfc;
//    else
    mtp->gfc = new double[mtp->nC * mtp->nPG * mtp->Nf];

    if( mtp->yfb )
            delete[] mtp->yfb;
//    else
    mtp->yfb = new double[mtp->nC * mtp->nPG * mtp->Nf];

    if( mtp->tt )
        delete[] mtp->tt;
//    else
    mtp->tt = new double[mtp->nC * mtp->Nf][9];

    if( mtp->BM)
       delete[] mtp->BM;
//    else
    mtp->BM = new double[mtp->nC];

    if( mtp->BdM)
       delete[] mtp->BdM;
//    else
    mtp->BdM = new double[mtp->nC];

    if( mtp->BmgpM)
       delete[] mtp->BmgpM;
//    else
    mtp->BmgpM = new double[mtp->nC*mtp->nPG];

    if( mtp->FmgpJ)
       delete[] mtp->FmgpJ;
//    else
    mtp->FmgpJ = new double[mtp->nFD];

    long int q, f, i;
    // Cleaning work arrays
       for( q=0; q <mtp->nC; q++ )
       {
            mtp->BM[q] = 0.0;
            mtp->BdM[q] = 0.0;
            for(i=0; i<mtp->Nf; i++ )
            {
                 Mb(q,i) = 0.;
                 dMb(q,i) = 0.0;
            }
       }
       for( q=0; q <mtp->nC; q++ )
       {
           for(f=0; f<mtp->nPG; f++ )
           {
               bmgpm(q,f) = 0.0;
               for(i=0; i<mtp->Nf; i++ )
               {
                   y(q,f,i) = 0.0;
                   g(q,f,i) = 0.0;
               }
           }
       }

//    CalcMGPdata(); // Initial calculation of masses, MGP compositions and part coeff tables
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
    long int p, i, kk, x_aq=-1, x_gf=-1;
    bool iRet = false;
    clock_t outp_time = (clock_t)0;

    BoxFluxTransportStart();


// na->CopyWorkNodeFromArray( 0, mtp->nC,  na->pNodT1() );
// na->GEM_write_dbr( "node0000.dat",  0, false );

    // In this mode, only calculation of equilibria in nodes one after another!
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
                  dMBflux( kk+1, mtp->dMB );
                  mtp->FDLi[kk+1][1] = oldp;
                  // mtp->FDLf[kk+1][0] = oldo;
               }
               // Adding the incoming MGP to p-th box
               dMBflux( kk, mtp->dMB );
             }
         }   // kk

         // change bulk composition in the box p
         BoxComposUpdate( p );

//sprintf(buf, "node_%4.4d_wave_%4.4d_in_dbr.dat", p, mtp->ct );
//na->CopyWorkNodeFromArray( p, mtp->nC,  na->pNodT1() );
//na->GEM_write_dbr( buf,  0, false);

         // calculate equilibrium state in p-th box
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

       // copy node array T1 into node array T0
       copyNodeArrays();

        mtp->ct += 1;
        mtp->oTau = mtp->cTau;
        mtp->cTau += mtp->dTau;

        if( mtp->PsVTK != S_OFF )
            outp_time += PrintPoint( 0 );

   } while ( mtp->cTau < mtp->Tau[STOP_] && mtp->ct < mtp->ntM );

  }
  catch( TError& xcpt )
  {
      gems_logger->error(" {}  {}", xcpt.title, xcpt.mess);
       return 1;
  }

  return 0;
}

#undef dMb
#undef Mb

const double MAXSTEP = 1.7;  // as was initially set - 10 sec is too much for the integrator!
//Calculate generic box megasystem of 'B' or 'F' type with MGP fluxes
//   mode - SIA or AIA
//
bool TGEM2MT::CalcBoxFluxModel( char /*mode*/ )
{
  try
  {	 
    BoxFluxTransportStart();

  // Initial calculation of chemical equilibria in all nodes with AIA initial approximation
  BoxEqStatesUpdate( -1,  0, mtp->cTau/mtp->tf, mtp->dTau/mtp->tf );

  // time iteration part
   nfcn = nstep = naccept = nrejct = 0;
   double cfactor = (mtp->dTau/mtp->tf/MAXSTEP > 1.0? mtp->dTau/mtp->tf/MAXSTEP: 1.0);
   MaxIter = mtp->ntM*cfactor;   // Need to check the multiplicator for maximum number of iterations
   double new_dTau = 0.0;
   
  // Getting into a time integration loop
  new_dTau = INTEG( 1e-3, /*mtp->cdv,*/ mtp->dTau/mtp->tf, mtp->Tau[START_]/mtp->tf, mtp->Tau[STOP_]/mtp->tf );
  new_dTau *= mtp->tf;

  gems_logger->debug("GEM2MT B mode: new dTau = {} ;  old dTau = {}", new_dTau, mtp->dTau);
    
  }
  catch( TError& xcpt )
  {
       gems_logger->error(" {}  {}", xcpt.title, xcpt.mess);
       return 1;
  }

  return 0;
}

//--------------------------------------------------------------------
// ODE integration process implementation

//const long int NMAX = 800;
const long int KM = 8;
const double UROUND = 1.73e-18;
const double FAC1 = 2.e-2;
const double FAC2 = 4.0;
const double FAC3 = .9;
const double FAC4 = .8;
const double SAFE1 = .65;
const double SAFE2 = .94;
// const double MAXSTEP = 1.7;
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
//  over the time interval [t_begin, t_end] with a desired time step length
//     step and the precision eps (the step can be reduced inside implicitly also as mtp->dTau)
//  Uses the implicit method of rational extrapolation (Bulirsch-Stoer)
//
// returns current (possibly reduced) step value or negative value in case of error
// 
double
TGEM2MT::INTEG( double eps, double step, double t_begin, double t_end )
{
    double  t, h1, h, v1;
    long int  j,i,reject,last,kc=0,kopt,k;

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
    // BoxEqStatesUpdate( 0, k, t, h ); // 14/12/2007 ????? may be done before in calc
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
//        ErrorIf( nstep >= mtp->ntM /*MaxIter*/, "INTEG", "No convergence - too many iterations" ); // 14/12/2007 !!!!
        ErrorIf( nstep >= MaxIter, "INTEG", "No convergence - too many iterations" );
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
        // step = h;   // Doubtful - indirect dependence!
//        mtp->cTau = t;
        for( i=0; i < mtp->nC*mtp->Nf; i++ )
            x[i] = mtp->tt[i][0];
        Solut(  x, dx, t );
        naccept++;
        mtp->ct++;
        BoxEqStatesUpdate( naccept, kc, t, h );  // here current time and step is set in mtp data structure!
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
    *tv = t*mtp->tf;     // set current time via pointer
    step = h;    // may be reduced
    return step;
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


