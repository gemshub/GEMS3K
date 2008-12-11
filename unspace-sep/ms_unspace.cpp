//-------------------------------------------------------------------
// $Id: ms_unspace.cpp 792 2006-09-19 08:10:41Z gems $
//
// Standalone variant of UnSpace module service functions 
//
// Copyright (C) 2007,2008 S.Dmytrieva, D.Kulik
//
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>

#include "ms_unspace.h"
#include "io_arrays.h"
#include "node.h"
//#include <iomanip>

//extern bool _comment;


TUnSpace* TUnSpace::pu;

TUnSpace::TUnSpace()
{
    usp=&us[0];

    // default data
    memset(usp, 0, sizeof(UNSPACE) );
    setup_defaults();
    na = 0;
}  

// Setting defaults before reading TUnSpace structure from text file
// 
void TUnSpace::setup_defaults()
{
	usp->Gstat = 0;
	usp->Astat = 0;
	memcpy( &usp->PsUnInt, "%%AA+------", 11 );
	memcpy( &usp->Pa_f_pha, "------0A1+", 10 );
	usp->Pa_Crit = UNSP_CRIT_PA;
	memcpy( &usp->PvPOM, "+-+-", 4 );
    usp->quan_lev = 0.05;
    usp->pH_lo = 0.;
    usp->pH_up = 14.;
    usp->Eh_lo = -1.;
    usp->Eh_up = 1.;
    usp->IC_lo = 0.;
    usp->IC_up = 3.;

}

TUnSpace::~TUnSpace()
{
  Free();
  if( na )
   delete na;
}


// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up Node classes
int TUnSpace::TaskSystemInit( const char *chbr_in1 )
{
  // The NodeArray must be allocated here
  TNode::na = na = new TNode();

  // Here we read the files needed as input for initializing GEMIPM2K
  // This files anloaded from UnSpace print comand
  if( TNode::na->GEM_init( chbr_in1 ) )
  {
    cout << "Error: Chemical system definition files listed in " << endl
         << chbr_in1 << endl
         << " cannot be found or are corrupt" << endl;
     return 1;  // error occured during reading the files
  }
  TNode::na->pCNode()->NodeStatusCH = NEED_GEM_AIA; // activating GEM IPM for automatic initial
                                    // approximation
// re-calculating equilibrium by calling GEMIPM
  TNode::na->GEM_run( false );
  // setup some internal data 
  pmu = TProfil::pm->pmp;
  
  // !!! internal from Profile Must be checked
  //    RMULTS* mup;
  for( int ii=0; ii<pmu->FIs; ii++)
  {
    if( pmu->PHC[ii] == 'a' )
    	mup_Laq = pmu->L1[ii];
  }
  // mup_Laq = pmu->LO+1;     // LO  - of DC in aqueous phase
   mup_Pg = pmu->PG;      // PG  - of DC in gaseous phase
   mup_SF   =  pmu->SF;   // List of PHASE definition keys [0:Fi-1]             DB   
   mup_Ll = pmu->L1;      // L1 vector, shows a number of DC included to each phase [0:Fi-1] DB

  //    SYSTEM *syu; // syu = TProfil::pm->syp;
  // char  *syu_Dcl;  // DC selection switches { + * - } [0:mu.L-1]
    syu_B = pmu->B;  // Vector b of bulk chemical composition of the system, moles [0:mu.N-1]
     
  //    MTPARM *tpp;
    pmu->tpp_G = new double[pmu->L]; // Partial molar(molal) Gibbs energy g(TP) (always), J/mole 
    pmu->tpp_S = new double[pmu->L];    // Partial molar(molal) entropy s(TP), J/mole/K
    pmu->tpp_Vm = new double[pmu->L];   // Partial molar(molal) volume Vm(TP) (always), J/bar
    fillValue( pmu->tpp_G, 0., pmu->L );
    fillValue( pmu->tpp_S, 0., pmu->L );
    fillValue( pmu->tpp_Vm, 0., pmu->L );
    tpp_S = pmu->S0;  // overload

     
     pmu->Guns = new double[pmu->L];  //  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
     pmu->Vuns = new double[pmu->L];  //  mu.L work vector of uncertainty space increments to tp->Vm
     fillValue( pmu->Guns, 0., pmu->L );
     fillValue( pmu->Vuns, 0., pmu->L );

     
  return 0;
}

// Calculation of the UnSpace task
void
TUnSpace::CalcTask( const char *key )
{
    int nAdapt = 1;

    if( usp->Pa_Adapt > '1')
    {  nAdapt = (int)(usp->Pa_Adapt-'0');
       usp->Gstat = GS_INDEF; // for Adapt mode need   buildTestedArrays
                                  // for each cicle
    }

   for( int ii=0; ii<(usp->L* usp->N); ii++ )
     usp->A[ii] = pmu->A[ii];
    
    for( int ii=0; ii<nAdapt; ii++ )
    {
      if( usp->Gstat != GS_DONE )
      {
          usp->Gstat = GS_GOING;
          usp->Astat = AS_INDEF;
          buildTestedArrays();
      }
      if( usp->PsGen[0] == S_ON )
      {
        analyseArrays();
        if( usp->Pa_Adapt > '1')
           AdapG();                    // !!!! test Kostin break ob =0 or ob>Q*0.95
      }
    }

  if( usp->Pa_Adapt > '1')
    usp->Pa_Adapt = '1';
  usp->Gstat = GS_DONE;   // for Adapt mode need   buildTestedArrays
                                 // for each cicle
  usp->Astat = AS_DONE;
}

// read TUnSpace structure from file
int TUnSpace::ReadTask( const char *unsp_in1 )
{
 // read GEM2MT structure from file
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
   fstream ff(unsp_in1, ios::in );
   ErrorIf( !ff.good() , unsp_in1, "Fileopen error");
   from_text_file( ff );
   return 0;
  }
  catch(TError& err)
  {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
  }
  return 1;
}

// Write TUnSpace data and results to file
//
int TUnSpace::WriteTask( const char *unsp_out_file )
{
 // write UnSpace task setup data to file
  fstream f_log("ipmlog.txt", ios::out|ios::app );
  try
  {
   fstream ff(unsp_out_file, ios::out );
   ErrorIf( !ff.good() , unsp_out_file, "Fileopen error");
   to_text_file( ff, true);
   gstring filename = unsp_out_file;
           filename += ".res";
   fstream ff1( filename.c_str(), ios::out );
   ErrorIf( !ff1.good() , filename.c_str(), "Fileopen error");
   result_to_text_file( ff1, true );
   return 0;
  }
  catch(TError& err)
  {
      f_log << err.title.c_str() << "  : " << err.mess.c_str() << endl;
  }
  return 1;
}

//--------------------------------------------------------------------------------------------------
outField TUnSpace_static_fields[25] =  {
 { "PsGen", 1,0 },  // 7
 { "Pa_Adapt" , 1,0 },
 { "Pa_OF", 1,0 },
 { "Pa_Crit", 1,0 },
 { "Pa_Zcp", 1,0 },
 { "PvSi", 1,0 },
 { "Pa_f_pha", 0,0 }, // default '-' !! do not forget to set default in the constructor
 { "Pa_f_mol", 0,0 },
 { "Pa_f_fug", 0,0 },
 { "Pa_f_mfr", 0,0 },
 { "Pa_f_pH", 0,0 },
 { "Pa_f_Eh", 0,0 },
 { "Pa_f_IC", 0,0 },
 { "PvPOM", 0,0 }, // default '+'
 { "PvPOR", 0,0 },  // default '-'
 // Dimensionalities related to the UnSpace problem
 { "Q", 1,0 },
 { "quan_lev", 0,0 },  // float 0.05
 { "nPG", 1,0 },
 { "nGB", 1,0 },
 { "nGN", 1,0 },
 { "nGR", 1,0 },
 { "N", 1,0 },  // take these dims from TNode (actually no need to read here)
 { "L", 1,0 },
 { "Ls", 1,0 },
 { "Fi", 1,0 }
};


outField TUnSpace_dynamic_fields[46] =  {
		{ "SGp", 1, 0 },
		{ "PbD", 1, 0 },
		{ "OVB", 1, 0 },
		{ "OVR", 1, 0 },
		{ "OVN", 1, 0 },
		{ "pH_lo", 0, 0 },
		{ "pH_up", 0, 0 },
		{ "Eh_lo", 0, 0 },
		{ "Eh_up", 0, 0 },
		{ "IC_lo", 0, 0 },
		{ "IC_up", 0, 0 },
		{ "m_t_lo", 0, 0 },
		{ "m_t_up", 0, 0 },
		{ "fug_lo", 0, 0 },
		{ "fug_up", 0, 0 },
		{ "f_PhA", 0, 0 },
		{ "NgT", 0, 0 },
		{ "T", 0, 0 },
		{ "IntT", 0, 0 },
		{ "NgP", 0, 0 },
		{ "P", 0, 0 },
		{ "IntP", 0, 0 },
		{ "NgV", 0, 0 },
		{ "V", 0, 0 },
		{ "IntV", 0, 0 },
		{ "NgLg", 1, 0 },
		{ "Gs", 1, 0 },
		{ "IntLg", 1, 0 },
		{ "NgGam", 0, 0 },
		{ "GAMs", 0, 0 },
		{ "IntGam", 0, 0 },
		{ "NgLs", 0, 0 },
		{ "Ss", 0, 0 },
		{ "IntLs", 0, 0 },
		{ "NgLv", 0, 0 },
		{ "Vs", 0, 0 },
		{ "IntLv", 0, 0 },
		{ "NgNb", 0, 0 },
		{ "Bs", 0, 0 },
		{ "IntNb", 0, 0 },
		{ "UnICn", 1, 0 },
		{ "UgDCn", 1, 0 },
		{ "UaDCn", 1, 0 },
		{ "UnDCAn", 1, 0 },
		{ "ParNames", 1, 0 },
		{ "ncp", 0, 0 }
};

   
// Reading TUnSpace structure from text file
void TUnSpace::from_text_file(fstream& ff)
{
// static arrays
 TReadArrays  rdar( 25, TUnSpace_static_fields, ff);
 short nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
   case 0: rdar.readArray( "PsGen", usp->PsGen, 7, 1);
           break;
   case 1: rdar.readArray( "Pa_Adapt", &usp->Pa_Adapt, 1, 1);
           break;
   case 2: rdar.readArray( "Pa_OF", &usp->Pa_OF, 1, 1);
           break;
   case 3: rdar.readArray( "Pa_Crit", &usp->Pa_Crit, 1, 1);
           break;
   case 4: rdar.readArray( "Pa_Zcp", &usp->Pa_Zcp, 1, 1);
           break;
   case 5: rdar.readArray( "PvSi", &usp->PvSi, 1, 1);
           break;
   case 6: rdar.readArray( "Pa_f_pha", &usp->Pa_f_pha, 1, 1);
           break;
   case 7: rdar.readArray( "Pa_f_mol", &usp->Pa_f_mol, 1, 1);
           break;
   case 8: rdar.readArray( "Pa_f_fug", &usp->Pa_f_fug, 1, 1);
           break;
   case 9: rdar.readArray( "Pa_f_mfr", &usp->Pa_f_mfr, 1, 1);
           break;
   case 10: rdar.readArray( "Pa_f_pH", &usp->Pa_f_pH, 1, 1);
           break;
   case 11: rdar.readArray( "Pa_f_Eh", &usp->Pa_f_Eh, 1, 1);
           break;
   case 12: rdar.readArray( "Pa_f_IC", &usp->Pa_f_IC, 1, 1);
           break;
   case 13: rdar.readArray( "PvPOM", &usp->PvPOM, 1, 1);
           break;
   case 14: rdar.readArray( "PvPOR", &usp->PvPOR, 1, 1);
           break;
   case 15: rdar.readArray( "Q", &usp->Q, 1);
           break;
   case 16: rdar.readArray( "quan_lev", &usp->quan_lev, 1);
           break;
   case 17: rdar.readArray( "nPG", &usp->nPG, 1);
           break;
   case 18: rdar.readArray( "nGB", &usp->nGB, 1);
           break;
   case 19: rdar.readArray( "nGN", &usp->nGN, 1);
           break;
   case 20: rdar.readArray( "nGR", &usp->nGR, 1);
           break;
   case 21: rdar.readArray( "N", &usp->N, 1);
           break;
   case 22: rdar.readArray( "L", &usp->L, 1);
           break;
   case 23: rdar.readArray( "Ls", &usp->Ls, 1);
           break;
   case 24: rdar.readArray( "Fi", &usp->Fi, 1);
           break;
 }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be readed from TUnSpace structure";
    Error( "Error", ret);
  }

 usp->nG = usp->nGB +usp->nGR+usp->nGN; 	
  Alloc();

//dynamic data
 TReadArrays  rddar( 46, TUnSpace_dynamic_fields, ff);

// Set up flags
 if(usp->PsGen[0] != S_ON )
 {
    rddar.setNoAlws( 26 /*"Gs"*/);
    rddar.setNoAlws( 25 /*"NgLg"*/);
    rddar.setNoAlws( 27 /*"IntLg"*/);
 }
 if(usp->PsGen[0] == S_ON  && usp->Pa_f_mol == S_ON)
 {
    rddar.setAlws( 11 /*"m_t_lo"*/);
    rddar.setAlws( 12 /*"m_t_up"*/);
 }
 if( usp->PsGen[0] == S_ON  && usp->Ls && usp->Pa_f_fug== S_ON )
 {
    rddar.setAlws( 13 /*"fug_lo"*/);
    rddar.setAlws( 14 /*"fug_up"*/);
 }
 if(usp->PsGen[1]== S_ON)
 {
   rddar.setAlws( 32 /*"Ss"*/);
   rddar.setAlws( 31 /*"NgLs"*/);
   rddar.setAlws( 33 /*"IntLs"*/);
 }
 if(usp->PsGen[5]== S_ON)
 {
  rddar.setAlws( 35 /*"Vs"*/);
  rddar.setAlws( 34 /*"NgLv"*/);
  rddar.setAlws( 36 /*"IntLv"*/);
 }
 if(usp->PsGen[2]== S_ON)
 {
  rddar.setAlws( 38 /*"Bs"*/);
  rddar.setAlws( 37 /*"NgNb"*/);
  rddar.setAlws( 39 /*"IntNb"*/);
 } 
 if(usp->PsGen[6]== S_ON && usp->Ls )   
 {
  rddar.setAlws( 29 /*"GAMs"*/);
  rddar.setAlws( 28 /*"NgGam"*/);
  rddar.setAlws( 30 /*"IntGam"*/);
 } 
 if( usp->Pa_f_pha == S_ON)
	  rddar.setAlws( 15 /*"f_PhA"*/);

  nfild = rddar.findNext();
  while( nfild >=0 )
  {
   switch( nfild )
   {
   case 0: rddar.readArray(  "SGp", usp->SGp[0], usp->nG, MAXPHNAME);
            break;
   case 1: rddar.readArray(  "PbD", usp->PbD, usp->nG);
            break;
   case 2: rddar.readArray(  "OVB", usp->OVB, usp->nGB);
            break;
   case 3: rddar.readArray(  "OVR", usp->OVR, usp->nGR);
            break;
   case 4: rddar.readArray(  "OVN", usp->OVN, usp->nGN);
            break;
   case 5: rddar.readArray(  "pH_lo", &usp->pH_lo, 1);
             break;
   case 6: rddar.readArray(  "pH_up", &usp->pH_up, 1);
             break;
   case 7: rddar.readArray(  "Eh_lo", &usp->Eh_lo, 1);
             break;
   case 8: rddar.readArray(  "Eh_up", &usp->Eh_up, 1);
             break;
   case 9: rddar.readArray(  "IC_lo", &usp->IC_lo, 1);
             break;
   case 10: rddar.readArray(  "IC_up", &usp->IC_up, 1);
             break;
   case 11:if( !usp->m_t_lo ) Error( "Error", "Array m_t_lo not used in task"); 
	        rddar.readArray(  "m_t_lo", usp->m_t_lo, usp->N);
            break;
   case 12:if( !usp->m_t_up ) Error( "Error", "Array m_t_up not used in task"); 
            rddar.readArray(  "m_t_up", usp->m_t_up, usp->N);
            break;
   case 13: if( !usp->fug_lo ) Error( "Error", "Array fug_lo not used in task"); 
             rddar.readArray(  "fug_lo", usp->fug_lo, usp->Ls);
            break;
   case 14: if( !usp->fug_up ) Error( "Error", "Array fug_up not used in task"); 
            rddar.readArray(  "fug_up", usp->fug_up, usp->Ls);
            break;
   case 15: if( !usp->f_PhA ) Error( "Error", "Array f_PhA not used in task"); 
            rddar.readArray(  "f_PhA", usp->f_PhA, usp->N);
            break;
   case 16: rddar.readArray(  "NgT", &usp->NgT, 1);
             break;
   case 17: rddar.readArray(  "T", usp->T, 2);
            break;
   case 18: rddar.readArray(  "IntT", usp->IntT, 2);
            break;
   case 19: rddar.readArray(  "NgP", &usp->NgP, 1);
             break;
   case 20: rddar.readArray(  "P", usp->P, 2);
            break;
   case 21: rddar.readArray(  "IntP", usp->IntP, 2);
            break;
   case 22: rddar.readArray(  "NgV", &usp->NgV, 1);
            break;
   case 23: rddar.readArray(  "V", usp->V, 2);
            break;
   case 24: rddar.readArray(  "IntV", usp->IntV, 2);
            break;
   case 25: if( !usp->NgLg ) Error( "Error", "Array NgLg not used in task"); 
            rddar.readArray(  "NgLg", usp->NgLg, usp->L);
            break;
   case 26: if( !usp->Gs ) Error( "Error", "Array Gs not used in task"); 
            rddar.readArray(  "Gs", usp->Gs[0], usp->L*2);
            break;
   case 27: if( !usp->IntLg ) Error( "Error", "Array IntLg not used in task"); 
            rddar.readArray(  "IntLg", usp->IntLg[0], usp->L*2);
            break;
   case 28: if( !usp->NgGam ) Error( "Error", "Array NgGam not used in task"); 
            rddar.readArray(  "NgGam", usp->NgGam, usp->Ls);
            break;
   case 29: if( !usp->GAMs ) Error( "Error", "Array GAMs not used in task"); 
            rddar.readArray(  "GAMs", usp->GAMs[0], usp->Ls*2);
            break;
   case 30: if( !usp->IntGam ) Error( "Error", "Array IntGam not used in task"); 
            rddar.readArray(  "IntGam", usp->IntGam[0], usp->Ls*2);
            break;
   case 31: if( !usp->NgLs ) Error( "Error", "Array NgLs not used in task"); 
            rddar.readArray(  "NgLs", usp->NgLs, usp->L);
            break;
   case 32: if( !usp->Ss ) Error( "Error", "Array Ss not used in task"); 
            rddar.readArray(  "Ss", usp->Ss[0], usp->L*2);
            break;
   case 33: if( !usp->IntLs ) Error( "Error", "Array IntLs not used in task"); 
            rddar.readArray(  "IntLs", usp->IntLs[0], usp->L*2);
            break;
   case 34: if( !usp->NgLv ) Error( "Error", "Array NgLv not used in task"); 
            rddar.readArray(  "NgLv", usp->NgLv, usp->L);
            break;
   case 35: if( !usp->Vs ) Error( "Error", "Array Vs not used in task"); 
             rddar.readArray(  "Vs", usp->Vs[0], usp->L*2);
            break;
   case 36: if( !usp->IntLv ) Error( "Error", "Array IntLv not used in task"); 
             rddar.readArray(  "IntLv", usp->IntLv[0], usp->L*2);
            break;
   case 37: if( !usp->NgNb ) Error( "Error", "Array NgNb not used in task"); 
            rddar.readArray(  "NgNb", usp->NgNb, usp->N);
            break;
   case 38: if( !usp->Bs ) Error( "Error", "Array Bs not used in task"); 
            rddar.readArray(  "Bs", usp->Bs[0], usp->N*2);
            break;
   case 39:if( !usp->IntNb ) Error( "Error", "Array IntNb not used in task"); 
            rddar.readArray(  "IntNb", usp->IntNb[0], usp->N*2);
            break;
   case 40: rddar.readArray(  "UnICn", usp->UnICn[0], UNSP_SIZE1, NAME_SIZE);
            break;
   case 41: rddar.readArray(  "UgDCn", usp->UgDCn[0], UNSP_SIZE1, NAME_SIZE);
            break;
   case 42: rddar.readArray(  "UaDCn", usp->UaDCn[0], UNSP_SIZE1, NAME_SIZE);
            break;
   case 43: rddar.readArray(  "UnDCAn", usp->UnICn[0], UNSP_SIZE2, NAME_SIZE);
            break;
   case 44: rddar.readArray(  "ParNames", usp->ParNames[0], usp->nPG, PARNAME_SIZE);
            break;
   case 45: // for test
	   rddar.readArray(  "ncp", usp->ncp, usp->Q*usp->nG );
  }
     nfild = rddar.findNext();
 }
 
 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from TUnSpace structure";
    Error( "Error", ret);
  }

}

//realloc dynamic memory for work arrays (do not reading)
void TUnSpace::phase_lists_new()
{
// alloc memory for nPhA size
if( usp->nPhA > 0 )
  {
    usp->PhAndx = new short[usp->nPhA*usp->N];
    usp->PhNum = new short[usp->nPhA];
    usp->PhAID = new char[usp->nPhA][8];
    usp->PhAlst = new char[usp->nPhA][80];
    usp->PhAfreq = new float[ usp->nPhA];
  }
}

//realloc dynamic memory for work arrays (do not reading)
void TUnSpace::work_dyn_new()
{
//  work (not in record)
    usp->A = new float[ usp->L* usp->N];
    usp->sv = new short[ usp->Q];

    usp->Zcp = new double[ usp->Q];
    usp->Zmin = new double[ usp->Q];
    usp->Zmax = new double[ usp->Q];
    usp->ZmaxAbs = new double[ usp->Q];
    usp->Hom = new double[ usp->Q];
    usp->Prob = new double[ usp->Q];

    usp->pmr = 0;
    if( usp->PvPOR == S_ON )
    {  usp->POR = new float[ usp->Q];
       usp->pmr = usp->POR;
    }
    
    if( usp->PvPOM == S_ON )
    {
      usp->POM = new float[ usp->Q*usp->Q]; 
      usp->pmr = usp->POM;
    }

  usp->UnIC = new double[usp->N][UNSP_SIZE1];
  usp->UgDC = new double[usp->nPG][UNSP_SIZE1];
  usp->UaDC = new double[usp->Ls][UNSP_SIZE1];
  usp->UnDCA = new double[usp->nPG][UNSP_SIZE2];
}

//realloc dynamic memory for work arrays (reading)
void TUnSpace::nG_dyn_new()
{
	usp->nG = usp->nGB +usp->nGR+usp->nGN; 	
  if(usp->nG)
  {
    usp->PbD =  new short[ usp->nG];
    usp->SGp =  new char[ usp->nG][MAXPHNAME];
    usp->ncp =  new float [ usp->Q* usp->nG];
  }
  if(usp->nGB)
    usp->OVB = new float[ usp->nGB+1 ];
  if(usp->nGR)
    usp->OVR = new float[ usp->nGR+1];
  if(usp->nGN)
    usp->OVN = new float[ usp->nGN+1];

  if( usp->nPG )
    usp->ParNames = new char[usp->nPG][PARNAME_SIZE];
}

// realloc dynamic memory
void TUnSpace::Alloc()
{

    usp->UnICn = new char[UNSP_SIZE1][NAME_SIZE];
    usp->UgDCn = new char[UNSP_SIZE1][NAME_SIZE];
    usp->UaDCn = new char[UNSP_SIZE1][NAME_SIZE];
    usp->UnDCAn = new char[UNSP_SIZE2][NAME_SIZE];
  // allocation memory for arrays

  if(usp->PsGen[0] == S_ON )
  {
    usp->Gs = new float[usp->L][2];
    usp->NgLg = new short[usp->L];
    usp->IntLg = new float[usp->L][2];
    
    // do not read 
    usp->vG = new double[ usp->Q* usp->L];
    usp->vY = new double[ usp->Q* usp->L];
    usp->vYF = new double[ usp->Q* usp->Fi];
    usp->vGam = new double[ usp->Q* usp->L];
    usp->vMol = new double[ usp->Q* usp->N];
    usp->vU = new double[ usp->Q* usp->N];
    usp->vpH = new float[usp->Q][3];
    usp->vT = new float[usp->Q];
    usp->vP = new float[usp->Q];
    usp->vV = new float[usp->Q];

    usp->qQ = (short)(usp->quan_lev*usp->Q);
    if(usp->qQ<1)
        usp->qQ=1;
    usp->quanCx = new short[usp->Q][4];
    usp->quanCv = new float[usp->Q][4];
  }

  if(usp->PsGen[0] == S_ON  && usp->Pa_f_mol == S_ON)
  {
    usp->m_t_lo = new float[ usp->N];
    usp->m_t_up = new float[ usp->N];
  }

  // do not read 
  if( usp->PsGen[0] == S_ON  && usp->Ls )
       usp->vFug = new double[ usp->Q* usp->Ls];

  if( usp->PsGen[0] == S_ON  && usp->Ls && usp->Pa_f_fug== S_ON )
  {
    usp->fug_lo = new float[ usp->Ls];
    usp->fug_up = new float[ usp->Ls];
  }

  if(usp->PsGen[1]== S_ON)
  {
    usp->Ss = new float[usp->L][2];
    usp->NgLs = new short[usp->L];
    usp->IntLs = new float[usp->L][2];
    usp->vS = new double[ usp->Q* usp->L]; // do not read
  }
  
  if(usp->PsGen[5]== S_ON)
  {
    usp->Vs = new float[usp->L][2];
    usp->NgLv = new short[usp->L];
    usp->IntLv = new float[usp->L][2];
    usp->vmV = new double[ usp->Q* usp->L]; // do not read
  }

  if(usp->PsGen[2]== S_ON)
  {
    usp->NgNb = new short[usp->N];
    usp->IntNb = new float[usp->N][2];
    usp->Bs = new double[usp->N][2];
    usp->vB = new double[ usp->Q* usp->N];  // do not read
  } 

  if(usp->PsGen[6]== S_ON && usp->Ls )   
  {
    usp->NgGam = new short[usp->Ls];
    usp->IntGam = new float[usp->Ls][2];
    usp->GAMs = new float[usp->Ls][2];
    usp->vNidP = new double[ usp->Q* usp->Ls]; // do not read
  } 
  
  if(/*usp->PsGen[0??] == S_ON  &&*/ usp->Pa_f_pha == S_ON)
    usp->f_PhA = new short[ usp->N ];

   nG_dyn_new();               
   work_dyn_new();
   phase_lists_new();
}

void TUnSpace::Free()
{
  if( usp->PhAndx  )
    delete[] usp->PhAndx;
  if( usp->PhNum  )
    delete[] usp->PhNum;
  if( usp->PhAID  )
     delete[] usp->PhAID;
  if(usp->PhAlst  )
      delete[] usp->PhAlst;
  if( usp->PhAfreq  )
      delete[] usp->PhAfreq;

  //  work (not in record)
  if( usp->A  )
     delete[] usp->A;
  if( usp->sv  )
     delete[] usp->sv;
  if( usp->Zcp  )
    delete[] usp->Zcp;
  if( usp->Zmin  )
      delete[] usp->Zmin;
  if( usp->Zmax  )
    delete[] usp->Zmax;
  if( usp->ZmaxAbs  )
     delete[] usp->ZmaxAbs;
  if(usp->Hom  )
      delete[] usp->Hom;
  if( usp->Prob  )
      delete[] usp->Prob;
   if( usp->POR  )
      delete[] usp->POR;
   if( usp->POM  )
      delete[] usp->POM;
  if( usp->UnIC  )
    delete[] usp->UnIC;
  if( usp->UgDC  )
      delete[] usp->UgDC;
  if( usp->UaDC  )
    delete[] usp->UaDC;
  if( usp->UnDCA  )
     delete[] usp->UnDCA;

//realloc dynamic memory for work arrays (reading)
   if( usp->PbD )
      delete[] usp->PbD;
   if( usp->SGp  )
      delete[] usp->SGp;
   if( usp->ncp  )
     delete[] usp->ncp;
   if( usp->OVB  )
       delete[] usp->OVB;
   if( usp->OVR  )
     delete[] usp->OVR;
   if( usp->OVN  )
      delete[] usp->OVN;
   if(usp->ParNames  )
       delete[] usp->ParNames;

// allocation memory for arrays
   if( usp->UnICn )
      delete[] usp->UnICn;
   if( usp->UgDCn  )
      delete[] usp->UgDCn;
   if( usp->UaDCn  )
     delete[] usp->UaDCn;
   if( usp->UnDCAn  )
       delete[] usp->UnDCAn;
   
   if( usp->Gs )
      delete[] usp->Gs;
   if( usp->NgLg  )
      delete[] usp->NgLg;
   if( usp->IntLg  )
     delete[] usp->IntLg;
   if( usp->vG  )
       delete[] usp->vG;
  if( usp->vY)
     delete[] usp->vY;
  if( usp->vYF  )
     delete[] usp->vYF;
  if( usp->vGam  )
    delete[] usp->vGam;
  if( usp->vMol  )
      delete[] usp->vMol;
   if( usp->vU )
      delete[] usp->vU;
   if( usp->vpH  )
      delete[] usp->vpH;
   if( usp->vT  )
     delete[] usp->vT;
   if( usp->vP  )
       delete[] usp->vP;
   if( usp->vV )
      delete[] usp->vV;
   if( usp->quanCx  )
      delete[] usp->quanCx;
   if( usp->quanCv  )
     delete[] usp->quanCv;
   if( usp->m_t_lo  )
       delete[] usp->m_t_lo;
  if( usp->m_t_up)
     delete[] usp->m_t_up;
  if( usp->vFug  )
     delete[] usp->vFug;
  if( usp->fug_lo  )
    delete[] usp->fug_lo;
  if( usp->fug_up  )
      delete[] usp->fug_up;
   if( usp->Ss )
      delete[] usp->Ss;
   if( usp->NgLs  )
      delete[] usp->NgLs;
   if( usp->IntLs  )
     delete[] usp->IntLs;
   if( usp->vS  )
       delete[] usp->vS;
    if( usp->Vs  )
       delete[] usp->Vs;
    if( usp->NgLv  )
      delete[] usp->NgLv;
    if( usp->IntLv  )
        delete[] usp->IntLv;
   if( usp->vmV)
      delete[] usp->vmV;
   if( usp->NgNb  )
      delete[] usp->NgNb;
   if( usp->IntNb  )
     delete[] usp->IntNb;
   if( usp->Bs  )
       delete[] usp->Bs;
   if( usp->vB )
       delete[] usp->vB;
   if( usp->NgGam  )
       delete[] usp->NgGam;
   if( usp->IntGam  )
      delete[] usp->IntGam;
   if( usp->GAMs  )
        delete[] usp->GAMs;
   if( usp->vNidP  )
        delete[] usp->vNidP;
   if( usp->f_PhA  )
       delete[] usp->f_PhA;
}

//---------------------------------------------------------------------------

