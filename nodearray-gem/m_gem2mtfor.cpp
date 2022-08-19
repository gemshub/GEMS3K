//-------------------------------------------------------------------
// $Id: main.cpp 792 2006-09-19 08:10:41Z gems $
//
// Debugging version of a finite-difference 1D advection-diffusion
// mass transport model supplied by Dr. Frieder Enzmann (Uni Mainz)
// coupled with GEMIPM2K module for calculation of chemical equilibria
//
// Direct access to the TNodeArray class for storing all data for nodes
//
// Copyright (C) 2005,2007, 2021 S.Dmytriyeva, F.Enzmann, D.Kulik
//
//-------------------------------------------------------------------

#include "m_gem2mt.h"

#include "GEMS3K/io_template.h"
#include "GEMS3K/io_nlohmann.h"
#include "GEMS3K/io_simdjson.h"
#include "GEMS3K/io_keyvalue.h"

extern const char* _GEMIPM_version_stamp;

//=============================================================

io_formats::outField GEM2MT_static_fields[57] =  {
    // Allocation and setup flags
     { "PvPGD" , 1, 0, 0, "# PvPGD: Use mobile phase groups definitions (+ -)" },
     { "PvFDL" , 1, 0, 0, "# PvFDL: Use MGP flux definition list (+ -)" },
     { "PvSFL" , 1, 0, 0, "# PvSFL: Use source fluxes and elemental stoichiometries for them (+ -)" },
     { "PvGrid" , 1, 0, 0, "# PvGrid: Use array of grid point locations (+ -)" },
     { "PvDDc" , 1, 0, 0, "# PvDDc:  Use diffusion coefficients for DC - DDc vector (+ -)" },
     { "PvDIc" , 1, 0, 0, "# PvDIc:  Use diffusion coefficients for IC - DIc vector (+ -)" },
     { "PvnVTK" , 1, 0, 0, "# PvnVTK: Use selected fields to VTK format (+ -)" },
   // Controls on operation
     { "PsMode", 1, 0, 0, "# PsMode: Code of GEM2MT mode of operation { S F A D T W V }" },
     { "PsSIA" , 0, 0, 0, "# PsSIA: Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)" },
     { "PsSdat" , 0, 0, 0, "# PsSdat:Save DataCH and inital DataBR files as text json files (j|f), as key value files (t|o) or binary (b)" },
     { "PsSdef" , 0, 0, 0, "# PsSdef:Do not write data items that contain only default values (+ -)" },
     { "PsScom" , 0, 0, 0, "# PsScom:Write files with comments for all data entries ( text mode ) or as pretty JSON (+ -)" },
     { "PsMO" , 0, 0, 0, "# PsMO: Use non stop debug output for nodes (+ -)" },
     { "PsVTK" , 1, 0, 0, "# PsVTK: Use non stop debug output nodes to VTK format(+ -)" },
     { "PsMPh" , 1, 0, 0, "# PsMPh: Type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )" },
      // sizes
     { "nC", 1, 0, 0, "# nC:  Input number of local equilibrium cells (nodes)" },
     { "nIV", 1, 0, 0, "# nIV:  Number of initial variants of the chemical system, nIV <= nC" },
     { "nMGP", 0, 0, 0, "# nMGP:  Number of mobile groups of phases, nMGP >= 0" },
     { "nFD", 0, 0, 0, "# nFD: Number of MGP fluxes defined in the megasystem, nFD >= 0" },
     { "nSFD", 0, 0, 0, "# nSFD:  Number of IC source flux compositions defined in megasystem, nSFD >= 0" },
     { "nEl" , 1, 0, 0, "# nEl: Number of electrolytes for diffusion coefficients in mDEl" },
     { "nPTypes", 0, 0, 0, "# nPTypes:  Number of allocated particle types < 20" },
     { "nProps", 0, 0, 0, "# nProps:  Number of particle statistic properties (for monitoring) >= nPTypes" },
     { "Nsd",   1, 0, 0, "# Nsd:   Number of references to data sources" },
     { "bTau" , 1, 0, 0, "# bTau:   Time point for the simulation break (Tau[0] at start)" },
     { "ntM" , 0, 0, 0, "# ntM:  Maximum allowed number of time iteration steps" },
     { "nVTKfld" , 0, 0, 0, "# nVTKfld:   Number of selected fields to VTK format" },
    { "nPai",   1, 0, 0, "# nPai: Number of P points in MTP interpolation array in DataCH ( 1 to 10 )" },
    { "nTai",   1, 0, 0, "# nTai: Number of T points in MTP interpolation array in DataCH ( 1 to 20 )" },
    { "Lsf",   1, 0, 0, "# Lsf:  Number of DCs in phases-solutions in Multi (DATACH)" },
    { "Nf",   1, 0, 0, "# Nf:   nICb number of ICs in  (DATABR) for setting box-fluxes" },
    { "FIf",   1, 0, 0, "# FIf:  nPHb number of phases in (DATABR) for setting box-fluxes" },
    { "Tau" , 1, 0, 0, "# Tau:   Physical time iterator (start,end,step)" },
    { "sizeLc" , 1, 0, 0, "# sizeLc:  Spatial dimensions of the medium defines topology of nodes ( x y z )" },
    //  Input for compositions of initial systems
    { "InpSys" , 0 , 0, 0, "# InpSys: Masses (kg) and volumes (L) for initial systems: Ms (total mass, normalize)" },
    { "Vsysb" , 0 , 0, 0, "# Vsysb:  Vs (total volume of the object, for volume concentrations)" },
    { "Mwatb" , 0 , 0, 0, "# Mwatb:  M(H2O) (mass of water-solvent for molalities)" },
    { "Maqb" , 0 , 0, 0, "# Maqb:   Maq (mass of aqueous solution for ppm etc.)" },
    { "Vaqb" , 0 , 0, 0, "# Vaqb:   Vaq (volume of aqueous solution for molarities)" },
    { "Pgb" , 0 , 0, 0, "# Pgb:   Pg (pressure in gas, for partial pressures)" },
    { "Tmolb" , 0 , 0, 0, "# Tmolb:  MOL total mole amount for basis sub-system composition calculations" },
    { "WmCb" , 0 , 0, 0, "# WmCb:  mole fraction of the carrier DC (e.g. sorbent or solvent)" },
    { "Asur" , 0 , 0, 0, "# TauAsur:  Specific surface area of the sorbent (for adsorbed species)" },
    { "tf" , 1 , 0, 0, "# tf:  Advection/diffusion mass transport: time step reduction factor (usually 1)" },
    { "Vt" , 0 , 0, 0, "# Vt:  Initial total node volume (m^3)" },
    { "vp" , 1, 0, 0, "# vp:  Fluid advection velocity (m/sec)" },
    { "eps" , 0 , 0, 0, "# eps:  Initial node effective porosity (0 < eps < 1), usually 1" },
    { "Km" , 0 , 0, 0, "# Km:  Initial effective permeability, m2, usually 1" },
    { "al" , 1 , 0, 0, "# al:  Initial value of specific longitudinal dispersivity (m), usually 1e-3" },
    { "Dif" , 1 , 0, 0, "# Dif:  Initial general aqueous medium diffusivity (m2/sec), usually 1e-9" },
    { "nto" , 0 , 0, 0, "# nto:  Initial tortuosity factor, usually 1" },
    { "cdv" , 1 , 0, 0, "# cdv:   Cutoff for IC amount differences in the node between time steps (mol), usually 1e-9" },
    { "cez" , 1 , 0, 0, "# cez:   Cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12" },
    { "Name" , 0 , 0, 0, "# Name: Comment line for the full name of this GEM2MT task " },
    { "Note" , 0 , 0, 0, "# Note: Comment line for optional enhancement of the dtName line " },
    { "mtWrkS" , 0 , 0, 0, "# mtWrkS: internal" },
    { "mtWrkF" , 0 , 0, 0, "# mtWrkF: internal" }
};

io_formats::outField GEM2MT_dynamic_fields[26] =
{  // write/read dynamic (array) data to/from the text-format IPM file
   {  "SDref",  0, 0, 0, "# SDref: List of SDref keys to data sources " },
   {  "SDval",  0, 0, 0, "# SDval: List of short comments to SDref  record keys in the GSDref list " },
   {  "DiCp",   1, 0, 0, "# DiCp:  Array of indexes of initial system variants for distributing to nodes [nC]" },
   {  "FDLi",   0, 0, 0, "# FDLi: Source/Receive box index in the flux definition" },
   {  "xFlds",  0, 0, 0, "# xFlds: List of selected fields to VTK format" },
   {  "mDDc",  0, 0, 0, "# mDDc: [Ls] diffusion coefficients for DC" },
   {  "mDIc",  0, 0, 0, "# mDIc: [N] diffusion coefficients for IC" },
   {  "mDEl",  0, 0, 0, "# mDEl: [nE] diffusion coefficients for electrolyte salts" },
   {  "HydP",  0, 0, 0, "# HydP:  Initial hydraulic parameters in nodes: Vt, vp, eps, Km, al, Dif,  nto" },
    {  "BSF",  0, 0, 0, "# BSF: [nSFD][N] table of bulk compositions of source fluxes " },
    {  "MB",  0, 0, 0, "\n# MB: [nC]  column of current masses of boxes or reservoirs (in kg) " },
    {  "dMB",  0, 0, 0, "# dMB: nC][Nb]  Table of current derivatives dM for elements in reservoirs " },
    {  "FDLf",  0, 0, 0, "# FDLf:  [nFD][4] Part of the flux defnition list: flux order ,flux rate, MGP quantity" },
    {  "PGT",  0, 0, 0, "# PGT: Quantities of phases in MGP [Fi][nPG]" },
    {  "nam_i",  0, 0, 0, "\n# nam_i: [nIV][12] id names of initial systems" },
    {  "for_e ",  0, 0, 0, "# for_e: [nE][40] formulae for diffusing dissolved electrolytes" },
    {  "FDLid",  0, 0, 0, "# FDLid: [nFD] IDs of fluxes" },
    {  "FDLop",  0, 0, 0, "# FDLop: [nFD] Operation codes (letters)" },
    {  "FDLmp",  0, 0, 0, "# FDLmp: [nFD] ID of MGP to move in this flux " },
    {  "MGPid",  0, 0, 0, "# MGPid: [nPG] ID list of mobile phase groups" },
    {  "UMGP",  0, 0, 0, "# UMGP: [nFi] units for setting phase quantities in MGP (see PGT )" },
    {  "mGrid",  0, 0, 0, "# mGrid: Array of grid point locations, size is nC*3" },
    {  "NPmean",  0, 0, 0, "# NPmean: Array of initial mean particle type numbers per node [nPTypes]" },
    {  "nPmin",  0, 0, 0, "# nPmin: Minimum average total number of particles of each type per one node [nPTypes]" },
    {  "nPmax",  0, 0, 0, "# nPmax: Maximum average total number of particles of each type per one node [nPTypes]" },
    {  "ParTD",  0, 0, 0, "# ParTD: Array of particle type definitions at t0 or after interruption [nPTypes]" }
};

// reset mt counters
void TGEM2MT::mt_reset()
{
// setup  counters
//  mtp->cT = mtp->Tai[START_];
//  mtp->cP = mtp->Pai[START_];
  mtp->cV = 0.;
  mtp->cTau = mtp->Tau[START_];
  mtp->ctm = mtp->tmi[START_];
  mtp->cnv = mtp->NVi[START_];
  mtp->qc = 0;
  mtp->kv = 0;
  mtp->jt = 0;
#ifndef IPMGEMPLUGIN
  mtp->cT = mtp->PTVm[START_][1];
  mtp->cP = mtp->PTVm[START_][0];
#endif
  mtp->ct = 0;
}

//internal calc record structure
bool TGEM2MT::internalCalc()
{
     bool iRet = 0;
     calcFinished = false;

     if( mtp->PsMode == RMT_MODE_B ) // || mtp->PsMode == RMT_MODE_F  ) // Flux-box integrated model
     {
         iRet = CalcBoxFluxModel( NEED_GEM_SIA );
     }
     else if( mtp->PsMode == RMT_MODE_S )
     {
         iRet = CalcSeqReacModel( NEED_GEM_SIA );
     }
     else if( mtp->PsMode == RMT_MODE_A || mtp->PsMode == RMT_MODE_C || mtp->PsMode == RMT_MODE_W
           || mtp->PsMode == RMT_MODE_F )  // 1D RMT models or simple 1D flux-box pipe sequence w/o integration
     {
         iRet =  Trans1D( NEED_GEM_SIA );  // here A,W,D and also F modes
     }
     else
     {
         ;  // Wrong model code - error message to be issued
     }
    calcFinished = true;
    return iRet;
}



//==========================================================================================
//set default information
void TGEM2MT::set_def(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
        "E03GTrem: Attempt to access corrupted dynamic memory.");

//    TProfil *aPa= TProfil::pm;
    memcpy( &mtp->PunE, "jjbC", 4 );
    memcpy( &mtp->PvICi, "++-----------+S00--+--+-----", 28 );
    strcpy( mtp->name,  "`" );
    strcpy( mtp->notes, "`" );
    strcpy( mtp->xNames, "X" );
    strcpy( mtp->yNames, "Y" );
    //memcpy( mtp->xNames, TProfil::pm->pa.GDpcc[0], MAXAXISNAME );
    //memcpy( mtp->yNames, TProfil::pm->pa.GDpcc[1], MAXAXISNAME );
    memset( &mtp->nC, 0, sizeof(long int)*32 );
    memset( &mtp->Msysb, 0, sizeof(double)*20 );
    memset( mtp->size[0], 0, sizeof(float)*8 );
    memset( mtp->sizeLc, 0, sizeof(double)*3 );
    memset( mtp->sykey, 0, sizeof(char)*(EQ_RKLEN+10) );
    mtp->nC = 21;
    mtp->nIV =2;
    mtp->ntM =1000;
    mtp->cdv = 1e-9;
    mtp->cez = 1e-12;
    mtp->nYS =0;
    mtp->nYE =1;
    mtp->nPai =1;
    mtp->nTai =1;
    mtp->tmi[START_] = 1000;
    mtp->tmi[STOP_] = 1200;
    mtp->tmi[STEP_] = 1;
    mtp->NVi[START_] = 0;
    mtp->NVi[STOP_] = 0;
    mtp->NVi[STEP_] = 0;
    mtp->Pai[START_] = 1.;
    mtp->Pai[STOP_] = 1.;
    mtp->Pai[STEP_] = 0.;
    mtp->Pai[3] = .5;      //Ptol
    mtp->Tai[START_] = 25.;
    mtp->Tai[STOP_] = 25.;
    mtp->Tai[STEP_] = 0.;
    mtp->Tai[3] = 1.;     //Ttol
    mtp->Tau[START_] = 0.;
    mtp->Tau[STOP_] = 1000.;
    mtp->Tau[STEP_] = 1.;
// pointers
    mtp->lNam = nullptr;
    mtp->lNamE = nullptr;
    mtp->tExpr = nullptr;
    mtp->gExpr = nullptr;
    mtp->sdref = nullptr;
    mtp->sdval = nullptr;
    mtp->DiCp = nullptr;
    mtp->FDLi = nullptr;
    mtp->PTVm = nullptr;
    mtp->StaP = nullptr;
    mtp->xVTKfld = nullptr;
    mtp->xEt = nullptr;
    mtp->yEt = nullptr;
    mtp->Bn = nullptr;
    mtp->HydP = nullptr;
    mtp->qpi = nullptr;
    mtp->qpc = nullptr;
    mtp->xt = nullptr;
    mtp->yt = nullptr;
    mtp->CIb = nullptr;
    mtp->CAb = nullptr;
    mtp->FDLf = nullptr;
    mtp->PGT = nullptr;
    mtp->Tval = nullptr;
    mtp->Pval = nullptr;
    mtp->nam_i = nullptr;
    mtp->for_i = nullptr;
    mtp->stld = nullptr;
    mtp->CIclb = nullptr;
    mtp->AUcln = nullptr;
    mtp->FDLid = nullptr;
    mtp->FDLop = nullptr;
    mtp->FDLmp = nullptr;
    mtp->MGPid = nullptr;
    mtp->UMGP = nullptr;
    mtp->SBM = nullptr;
#ifndef IPMGEMPLUGIN
    plot = nullptr;
#endif
    mtp->BSF = nullptr;
    mtp->MB = nullptr;
    mtp->dMB = nullptr;
    mtp->DDc = nullptr;
    mtp->DIc = nullptr;
    mtp->DEl = nullptr;
    mtp->for_e = nullptr;
    mtp->xIC = nullptr;
    mtp->xDC = nullptr;
    mtp->xPH = nullptr;
    mtp->grid = nullptr;
    mtp->NPmean = nullptr;
    mtp->nPmin = nullptr;
    mtp->nPmax = nullptr;
    mtp->ParTD = nullptr;
// work
    mtp->BM = nullptr;
    mtp->BdM = nullptr;
    mtp->FmgpJ = nullptr;
    mtp->BmgpM = nullptr;
//
    mtp->An = nullptr;
    mtp->Ae = nullptr;
    mtp->gfc = nullptr;
    mtp->yfb = nullptr;
    mtp->tt = nullptr;
    mtp->etext = nullptr;
    mtp->tprn = nullptr;
    na = nullptr;
    pa_mt = nullptr;
}


/*set default information
void TGEM2MT::set_def(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
        "E03GTrem: Attempt to access corrupted dynamic memory.");

//    TProfil *aPa= TProfil::pm;
    memcpy( &mtp->PunE, "jjbC", 4 );
    memcpy( &mtp->PvICi, "++------------S00--+-++-----", 28 );
    strcpy( mtp->name,  "`" );
    strcpy( mtp->notes, "`" );
    strcpy( mtp->xNames, "X" );
    strcpy( mtp->yNames, "Y" );
    memset( &mtp->nC, 0, sizeof(long int)*32 );
    memset( &mtp->Msysb, 0, sizeof(double)*20 );
    memset( mtp->size[0], 0, sizeof(float)*8 );
    memset( mtp->sizeLc, 0, sizeof(double)*3 );
    memset( mtp->sykey, 0, sizeof(char)*(EQ_RKLEN+10) );
    mtp->nC = 21;
    mtp->nIV =2;
    mtp->ntM =1000;
    mtp->cdv = 1e-9;
    mtp->cez = 1e-12;
    mtp->nYS =0;
    mtp->nYE =1;
    mtp->nPai =1;
    mtp->nTai =1;
    mtp->tmi[START_] = 1000;
    mtp->tmi[STOP_] = 1200;
    mtp->tmi[STEP_] = 1;
    mtp->NVi[START_] = 0;
    mtp->NVi[STOP_] = 0;
    mtp->NVi[STEP_] = 0;
    mtp->Pai[START_] = 1.;
    mtp->Pai[STOP_] = 1.;
    mtp->Pai[STEP_] = 0.;
    mtp->Pai[3] = .5;      //Ptol
    mtp->Tai[START_] = 25.;
    mtp->Tai[STOP_] = 25.;
    mtp->Tai[STEP_] = 0.;
    mtp->Tai[3] = 1.;     //Ttol
    mtp->Tau[START_] = 0.;
    mtp->Tau[STOP_] = 1000.;
    mtp->Tau[STEP_] = 1.;
// pointers
    mtp->lNam = NULL;
    mtp->lNamE = NULL;
    mtp->tExpr = 0;
    mtp->gExpr = 0;
    mtp->sdref = 0;
    mtp->sdval = 0;
    mtp->DiCp = 0;
    mtp->FDLi = 0;
    mtp->PTVm = 0;
    mtp->StaP = 0;
    mtp->xVTKfld = 0;
    mtp->xEt = 0;
    mtp->yEt = 0;
    mtp->Bn = 0;
    mtp->HydP = 0;
    mtp->qpi = 0;
    mtp->qpc = 0;
    mtp->xt = 0;
    mtp->yt = 0;
    mtp->CIb = 0;
    mtp->CAb = 0;
    mtp->FDLf = 0;
    mtp->PGT = 0;
    mtp->Tval = 0;
    mtp->Pval = 0;
    mtp->nam_i = 0;
    mtp->for_i = 0;
    mtp->stld = 0;
    mtp->CIclb = 0;
    mtp->AUcln = 0;
    mtp->FDLid = 0;
    mtp->FDLop = 0;
    mtp->FDLmp = 0;
    mtp->MGPid = 0;
    mtp->UMGP = 0;
    mtp->SBM = 0;
    mtp->BSF = 0;
    mtp->MB = 0;
    mtp->dMB = 0;
    mtp->DDc = 0;
    mtp->DIc = 0;
    mtp->DEl = 0;
    mtp->for_e = 0;
    mtp->xIC = 0;
    mtp->xDC = 0;
    mtp->xPH = 0;
    mtp->grid = 0;
    mtp->NPmean = 0;
    mtp->nPmin = 0;
    mtp->nPmax = 0;
    mtp->ParTD = 0;
    mtp->BM = 0;
    mtp->BdM = 0;
    mtp->BmgpM = 0;
    mtp->BmgpJ = 0;

// work
    mtp->An = 0;
    mtp->Ae = 0;
    mtp->gfc = 0;
    mtp->yfb = 0;
    mtp->tt = 0;
    mtp->etext = 0;
    mtp->tprn = 0;
    na = 0;
    pa = 0;
}
*/

//==============================================================================

void TGEM2MT::checkAlws(io_formats::TRWArrays&  prar1, io_formats::TRWArrays&  prar) const
{
    // Set always for task
    if( mtp->PsMode == RMT_MODE_W  )
    {
        prar1.setAlws( f_nPTypes);
        prar1.setAlws( f_nProps);
        prar.setAlws( f__NPmean);
        prar.setAlws( f__nPmin);
        prar.setAlws( f__nPmax);
        prar.setAlws( f__ParTD);
    }
    if( mtp->PvGrid == S_ON )
        prar.setAlws( f__mGrid);

    if( mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B )
        prar.setAlws( f__HydP);

    if( mtp->PvFDL == S_ON )
      {
        prar1.setAlws( f_nFD );
        prar.setAlws( f__FDLi);
        prar.setAlws( f__FDLf);
        prar.setAlws( f__FDLid);
        prar.setAlws( f__FDLop);
        prar.setAlws( f__FDLmp);
      }

    if( mtp->PvPGD == S_ON )
      {
        prar1.setAlws( f_nMGP );
        prar.setAlws( f__PGT);
        prar.setAlws( f__MGPid);
        prar.setAlws( f__UMGP);
     }

    if( mtp->PvSFL == S_ON )
    {
      prar1.setAlws( f_nSFD );
      prar.setAlws( f__BSF);
    }

    /*if( mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF )
      {
        prar.setAlws( f__MB);
        prar.setAlws( f__dMB);
      }
    */
    if( mtp->PvnVTK == S_ON )
    {
      prar1.setAlws( f_nVTKfld );
      prar.setAlws( f__xFlds);
    }

    if( mtp->PvDDc == S_ON )
        prar.setAlws( f__mDDc);
    if( mtp->PvDIc == S_ON )
        prar.setAlws( f__mDIc);

    if( mtp->nEl > 0  )
      {
        prar.setAlws( f__mDEl);
        prar.setAlws( f__for_e);
     }
}

template<typename TIO>
void TGEM2MT::to_text_file( TIO& out_format, bool with_comments, bool brief_mode ) const
{
    bool _comment = with_comments;

    out_format.put_head( "", "gem2mt");
    io_formats::TPrintArrays<TIO>  prar1(57, GEM2MT_static_fields, out_format );
    io_formats::TPrintArrays<TIO>  prar(26, GEM2MT_dynamic_fields, out_format );

    // Set always for task
    checkAlws(prar1,prar);

    if( _comment )
    {
        prar1.writeComment( _comment, std::string( "# ") + _GEMIPM_version_stamp);;
        //        << "# File: " << path << endl;
        prar1.writeComment( _comment, "# Comments can be marked with # $ ; as the first character in the line\n");
    }

    prar1.writeArrayF(f_Name, mtp->name, 1, MAXFORMULA, _comment, brief_mode );
    prar1.writeArrayF(f_Note, mtp->notes, 1, MAXFORMULA, _comment, brief_mode );

    if( _comment )
        prar1.writeComment( _comment, "\n## (1) Allocation and setup flags");

    prar1.writeField(f_PvPGD, mtp->PvPGD, _comment, brief_mode  );
    prar1.writeField(f_PvFDL, mtp->PvFDL, _comment, brief_mode  );
    prar1.writeField(f_PvSFL, mtp->PvSFL, _comment, brief_mode  );
    prar1.writeField(f_PvGrid, mtp->PvGrid, _comment, brief_mode  );
    prar1.writeField(f_PvDDc, mtp->PvDDc, _comment, brief_mode  );
    prar1.writeField(f_PvDIc, mtp->PvDIc, _comment, brief_mode  );
    prar1.writeField(f_PvnVTK, mtp->PvnVTK, _comment, brief_mode  );

    if( _comment )
        prar1.writeComment( _comment, "\n## (2) Controls on operation");

    prar1.writeField(f_PsMode, mtp->PsMode, _comment, brief_mode  );
    prar1.writeField(f_PsSIA, mtp->PsSIA, _comment, brief_mode  );
    prar1.writeField(f_PsSdat, mtp->PsSdat, _comment, brief_mode  );
    prar1.writeField(f_PsSdef, mtp->PsSdef, _comment, brief_mode  );
    prar1.writeField(f_PsScom, mtp->PsScom, _comment, brief_mode  );
    prar1.writeField(f_PsMO, mtp->PsMO, _comment, brief_mode  );
    prar1.writeField(f_PsVTK, mtp->PsVTK, _comment, brief_mode  );
    prar1.writeField(f_PsMPh, mtp->PsMPh, _comment, brief_mode  );

    if( _comment )
        prar1.writeComment( _comment, "\n## (3) Dimensions for gem2mt (memory allocation)");

    prar1.writeField(f_nC, mtp->nC, _comment, brief_mode  );
    prar1.writeField(f_nIV, mtp->nIV, _comment, brief_mode  );
    prar1.writeField(f_nMGP, mtp->nPG, _comment, brief_mode  );
    prar1.writeField(f_nFD, mtp->nFD, _comment, brief_mode  );
    prar1.writeField(f_nSFD, mtp->nSFD, _comment, brief_mode  );
    prar1.writeField(f_nEl, mtp->nEl, _comment, brief_mode  );
    prar1.writeField(f_nPTypes, mtp->nPTypes, _comment, brief_mode  );
    prar1.writeField(f_nProps, mtp->nProps, _comment, brief_mode  );
    prar1.writeField(f_Nsd, mtp->Nsd, _comment, brief_mode  );
    prar1.writeField(f_bTau, mtp->bTau, _comment, brief_mode  );
    prar1.writeField(f_ntM, mtp->ntM, _comment, brief_mode  );
    prar1.writeField(f_nVTKfld, mtp->nVTKfld, _comment, brief_mode  );

    prar1.writeField(f_nPai, mtp->nPai, _comment, brief_mode  );
    prar1.writeField(f_nTai, mtp->nTai, _comment, brief_mode  );
    prar1.writeField(f_Lsf, mtp->Lsf, _comment, brief_mode  );
    prar1.writeField(f_Nf, mtp->Nf, _comment, brief_mode  );
    prar1.writeField(f_FIf, mtp->FIf, _comment, brief_mode  );

    prar1.writeArray(f_Tau, mtp->Tau, 3, 3, _comment, brief_mode  );
    prar1.writeArray(f_sizeLc, mtp->sizeLc, 3, 3, _comment, brief_mode  );

    if( _comment )
        prar1.writeComment( _comment, "\n## (4) Input for compositions of initial systems");

    prar1.writeField(f_InpSys, mtp->Msysb, _comment, brief_mode  );
    prar1.writeField(f_Vsysb, mtp->Vsysb, _comment, brief_mode  );
    prar1.writeField(f_Mwatb, mtp->Mwatb, _comment, brief_mode  );
    prar1.writeField(f_Maqb, mtp->Maqb, _comment, brief_mode  );
    prar1.writeField(f_Vaqb, mtp->Vaqb, _comment, brief_mode  );
    prar1.writeField(f_Pgb, mtp->Pgb, _comment, brief_mode  );
    prar1.writeField(f_Tmolb, mtp->Tmolb, _comment, brief_mode  );
    prar1.writeField(f_WmCb, mtp->WmCb, _comment, brief_mode  );
    prar1.writeField(f_Asur, mtp->Asur, _comment, brief_mode  );
    prar1.writeField(ff_tf, mtp->tf, _comment, brief_mode  );
    prar1.writeField(ff_Vt, mtp->vol_in, _comment, brief_mode  );
    prar1.writeField(ff_vp, mtp->fVel, _comment, brief_mode  );
    prar1.writeField(ff_eps, mtp->eps_in, _comment, brief_mode  );
    prar1.writeField(ff_Km, mtp->Km_in, _comment, brief_mode  );
    prar1.writeField(ff_al, mtp->al_in, _comment, brief_mode  );
    prar1.writeField(ff_Dif, mtp->Dif_in, _comment, brief_mode  );
    prar1.writeField(ff_nto, mtp->nto_in, _comment, brief_mode  );
    prar1.writeField(ff_cdv, mtp->cdv, _comment, brief_mode  );
    prar1.writeField(ff_cez, mtp->cez, _comment, brief_mode  );


    if( _comment )
        prar1.writeComment( _comment, "\n## (5) Internal stop point");

    prar1.writeArray(f_mtWrkS, &mtp->ctm, 12, 6, _comment, brief_mode  );
    prar1.writeArray(f_mtWrkF, &mtp->cT, 10, 6, _comment, brief_mode  );

    prar1.writeComment( true,  "\n<END_DIM>\n" );

    // dynamic arrays - must follow static data
    if( mtp->PsMode == RMT_MODE_W  )
    {
        if( _comment )
            prar.writeComment( _comment,  "\n## W random-walk advection-diffusion coupled RMT model");
        prar.writeArray(  f__NPmean, mtp->NPmean, mtp->nPTypes, -1L,_comment, brief_mode);
        prar.writeArray(  f__nPmin, mtp->nPmin, mtp->nPTypes, -1L,_comment, brief_mode);
        prar.writeArray(  f__nPmax, mtp->nPmax, mtp->nPTypes, -1L,_comment, brief_mode);
        prar.writeArray(  f__ParTD, &mtp->ParTD[0][0], mtp->nPTypes*6, 6L,_comment, brief_mode);
    }
    if( mtp->PvGrid == S_ON )
        prar.writeArray(  f__mGrid, &mtp->grid[0][0], mtp->nC*3, 3L,_comment, brief_mode);
    prar.writeArray(  f__DiCp, &mtp->DiCp[0][0], mtp->nC*2, 2L,_comment, brief_mode);

    if( mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B )
        prar.writeArray(  f__HydP, &mtp->HydP[0][0], mtp->nC*SIZE_HYDP, SIZE_HYDP,_comment, brief_mode);

    if( mtp->PvFDL == S_ON )
    {
        if( _comment )
            prar.writeComment( _comment, "\n## Use flux definition list");
        prar.writeArray(  f__FDLi, &mtp->FDLi[0][0],  mtp->nFD*2, 2L,_comment, brief_mode);
        prar.writeArray(  f__FDLf, &mtp->FDLf[0][0],  mtp->nFD*4, 4L,_comment, brief_mode);
        prar.writeArrayF(  f__FDLid, &mtp->FDLid[0][0], mtp->nFD, MAXSYMB,_comment, brief_mode);
        prar.writeArrayF(  f__FDLop, &mtp->FDLop[0][0], mtp->nFD, MAXSYMB,_comment, brief_mode);
        prar.writeArrayF(  f__FDLmp, &mtp->FDLmp[0][0], mtp->nFD, MAXSYMB,_comment, brief_mode);
    }

    if( mtp->PvPGD == S_ON )
    {
        if( _comment )
            prar.writeComment( _comment, "\n## Use phase groups definitions");
        prar.writeArray(  f__PGT, mtp->PGT,  mtp->FIf*mtp->nPG, mtp->nPG,_comment, brief_mode);
        prar.writeArrayF(  f__MGPid, &mtp->MGPid[0][0],  mtp->nPG, MAXSYMB,_comment, brief_mode);
        prar.writeArrayF(  f__UMGP, mtp->UMGP, mtp->FIf, 1L,_comment, brief_mode);
    }

    if( mtp->PvSFL == S_ON )
        prar.writeArray(  f__BSF, mtp->BSF,  mtp->nSFD*mtp->Nf, mtp->Nf,_comment, brief_mode);

    if( mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF )
    {
        prar.writeArray(  f__MB, mtp->MB,  mtp->nC*mtp->Nf, mtp->Nf,_comment, brief_mode);
        prar.writeArray(  f__dMB, mtp->dMB, mtp->nC*mtp->Nf,mtp->Nf,_comment, brief_mode);
    }

    if( mtp->PvnVTK == S_ON )
        prar.writeArray(  f__xFlds, &mtp->xVTKfld[0][0],  mtp->nVTKfld*2, 2L,_comment, brief_mode);

    if( mtp->PvDDc == S_ON )
        prar.writeArray(  f__mDDc, mtp->DDc, mtp->Lsf, -1L,_comment, brief_mode);
    if( mtp->PvDIc == S_ON )
        prar.writeArray(  f__mDIc, mtp->DIc, mtp->Nf, -1L,_comment, brief_mode);
    if( mtp->nEl > 0  )
    {
        prar.writeArray(  f__mDEl, mtp->DEl, mtp->nEl, -1L,_comment, brief_mode);
        prar.writeArrayF(  f__for_e, mtp->for_e[0],mtp->nEl, MAXFORMUNITDT,_comment, brief_mode);
    }

    prar.writeArrayF(  f__nam_i, mtp->nam_i[0],mtp->nIV, MAXIDNAME,_comment, brief_mode);

    if( mtp->Nsd > 0 )
    {
        prar.writeArrayF(  f_SDref, mtp->sdref[0],mtp->Nsd, V_SD_RKLEN,_comment, brief_mode);
        prar.writeArrayF(  f__SDval, mtp->sdval[0],mtp->Nsd, V_SD_RKLEN,_comment, brief_mode);
    }

    //!!!mtp->Tval  = new double[ mtp->nTai ];  // from DataCH
    //!!!mtp->Pval  = new double[ mtp->nPai ];

    out_format.dump( _comment );
}

// Reading dataCH structure from text file
template<typename TIO>
void TGEM2MT::from_text_file(TIO& in_format)
{
    io_formats::TReadArrays<TIO> rdar( 57, GEM2MT_static_fields, in_format);

    long int nfild = rdar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f_PvPGD: rdar.readArray( "PvPGD",  &mtp->PvPGD, 1, 1);
            break;
        case f_PvFDL: rdar.readArray( "PvFDL", &mtp->PvFDL, 1, 1);
            break;
        case f_PvSFL: rdar.readArray( "PvSFL", &mtp->PvSFL, 1, 1);
            break;
        case f_PvGrid: rdar.readArray( "PvGrid",  &mtp->PvGrid, 1, 1);
            break;
        case f_PvDDc: rdar.readArray( "PvDDc", &mtp->PvDDc, 1, 1);
            break;
        case f_PvDIc: rdar.readArray( "PvDIc", &mtp->PvDIc, 1, 1);
            break;
        case f_PvnVTK: rdar.readArray( "PvnVTK", &mtp->PvnVTK, 1, 1);
            break;
        case f_PsMode: rdar.readArray( "PsMode",  &mtp->PsMode, 1, 1);
            break;
        case f_PsSIA: rdar.readArray( "PsSIA",  &mtp->PsSIA, 1, 1);
            break;
        case f_PsSdat: rdar.readArray( "PsSdat",  &mtp->PsSdat, 1, 1);
            break;
        case f_PsSdef: rdar.readArray( "PsSdef",  &mtp->PsSdef, 1, 1);
            break;
        case f_PsScom: rdar.readArray( "PsScom",  &mtp->PsScom, 1, 1);
            break;
        case f_PsMO: rdar.readArray( "PsMO",  &mtp->PsMO, 1, 1);
            break;
        case f_PsVTK: rdar.readArray( "PsVTK",  &mtp->PsVTK, 1, 1);
            break;
        case f_PsMPh: rdar.readArray( "PsMPh",  &mtp->PsMPh, 1, 1);
            break;
        case f_nC: rdar.readArray( "nC",  &mtp->nC, 1);
            break;
        case f_nIV: rdar.readArray( "nIV",  &mtp->nIV, 1);
            break;
        case f_nMGP: rdar.readArray( "nMGP",  &mtp->nPG, 1);
            break;
        case f_nFD: rdar.readArray( "nFD",  &mtp->nFD, 1);
            break;
        case f_nSFD: rdar.readArray( "nSFD",  &mtp->nSFD, 1);
            break;
        case f_nEl: rdar.readArray( "nEl",  &mtp->nEl, 1);
            break;
        case f_nPTypes: rdar.readArray( "nPTypes",  &mtp->nPTypes, 1);
            break;
        case f_nProps: rdar.readArray( "nProps",  &mtp->nProps, 1);
            break;
        case f_Nsd: rdar.readArray( "Nsd",  &mtp->Nsd, 1);
            break;
        case f_bTau: rdar.readArray( "bTau",  &mtp->bTau, 1);
            break;
        case f_ntM: rdar.readArray( "ntM",  &mtp->ntM, 1);
            break;
        case f_nVTKfld: rdar.readArray( "nVTKfld",  &mtp->nVTKfld, 1);
            break;

        case f_nPai: rdar.readArray( "nPai",  &mtp->nPai, 1);
            break;
        case f_nTai: rdar.readArray( "nTai",  &mtp->nTai, 1);
            break;
        case f_Lsf: rdar.readArray( "Lsf",  &mtp->Lsf, 1);
            break;
        case f_Nf: rdar.readArray( "Nf",  &mtp->Nf, 1);
            break;
        case f_FIf: rdar.readArray( "FIf",  &mtp->FIf, 1);
            break;

        case f_Tau: rdar.readArray( "Tau",  mtp->Tau, 3);
            break;
        case f_sizeLc: rdar.readArray( "sizeLc",  mtp->sizeLc, 3);
            break;
        case f_InpSys: rdar.readArray( "InpSys",  &mtp->Msysb, 1);
            break;
        case f_Vsysb: rdar.readArray( "Vsysb",  &mtp->Vsysb, 1);
            break;
        case f_Mwatb: rdar.readArray( "Mwatb",  &mtp->Mwatb, 1);
            break;
        case f_Maqb: rdar.readArray( "Maqb",  &mtp->Maqb, 1);
            break;
        case f_Vaqb: rdar.readArray( "Vaqb",  &mtp->Vaqb, 1);
            break;
        case f_Pgb: rdar.readArray( "Pgb",  &mtp->Pgb, 1);
            break;
        case f_Tmolb: rdar.readArray( "Tmolb",  &mtp->Tmolb, 1);
            break;
        case f_WmCb: rdar.readArray( "WmCb",  &mtp->WmCb, 1);
            break;
        case f_Asur: rdar.readArray( "Asur",  &mtp->Asur, 1);
            break;
        case ff_tf: rdar.readArray( "tf",  &mtp->tf, 1);
            break;
        case ff_Vt: rdar.readArray( "Vt",  &mtp->vol_in, 1);
            break;
        case ff_vp: rdar.readArray( "vp",  &mtp->fVel, 1);
            break;
        case ff_eps: rdar.readArray( "eps",  &mtp->eps_in, 1);
            break;
        case ff_Km: rdar.readArray( "Km",  &mtp->Km_in, 1);
            break;
        case ff_al: rdar.readArray( "al",  &mtp->al_in, 1);
            break;
        case ff_Dif: rdar.readArray( "Dif",  &mtp->Dif_in, 1);
            break;
        case ff_nto: rdar.readArray( "nto",  &mtp->nto_in, 1);
            break;
        case ff_cdv: rdar.readArray( "cdv",  &mtp->cdv, 1);
            break;
        case ff_cez: rdar.readArray( "cez",  &mtp->cez, 1);
            break;

        case f_Name: rdar.readArray( "Name",  mtp->name, 1, MAXFORMULA);
            break;
        case f_Note: rdar.readArray( "Note",  mtp->notes, 1, MAXFORMULA);
            break;
        case f_mtWrkS: rdar.readArray( "mtWrkS",  &mtp->ctm, 12);
            break;
        case f_mtWrkF: rdar.readArray( "mtWrkF",  &mtp->cT, 10);
        }
        nfild = rdar.findNext();
    }

    //dynamic data
    io_formats::TReadArrays<TIO>   rddar(  26, GEM2MT_dynamic_fields, in_format);

    // set alwase flags for gems2mt
    checkAlws(rdar, rddar);

    // testing read
    auto ret = rdar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from gem2mt structure";
        Error( "Error", ret);
    }

    // realloc memory
#ifndef IPMGEMPLUGIN
    dyn_new(0);
#else
    mem_new(0);
#endif

    nfild = rddar.findNext();
    while( nfild >=0 )
    {
        switch( nfild )
        {
        case f__NPmean: rddar.readArray( "NPmean", mtp->NPmean, mtp->nPTypes);
            break;
        case f__nPmin: rddar.readArray( "nPmin", mtp->nPmin, mtp->nPTypes);
            break;
        case f__nPmax: rddar.readArray( "nPmax", mtp->nPmax, mtp->nPTypes);
            break;
        case f__ParTD: rddar.readArray( "ParTD", &mtp->ParTD[0][0], mtp->nPTypes*6);
            break;
        case f__mGrid: rddar.readArray( "mGrid", &mtp->grid[0][0], mtp->nC*3);
            break;
        case f__DiCp: rddar.readArray( "DiCp", &mtp->DiCp[0][0], mtp->nC*2);
            break;
        case f__HydP: rddar.readArray( "HydP", &mtp->HydP[0][0], mtp->nC*SIZE_HYDP);
            break;
        case f__FDLi: rddar.readArray( "FDLi", &mtp->FDLi[0][0], mtp->nFD*2);
            break;
        case f__FDLf: rddar.readArray( "FDLf", &mtp->FDLf[0][0], mtp->nFD*4);
            break;
        case f__FDLid: rddar.readArray( "FDLid", &mtp->FDLid[0][0], mtp->nFD, MAXSYMB);
            break;
        case f__FDLop: rddar.readArray( "FDLop", &mtp->FDLop[0][0], mtp->nFD, MAXSYMB);
            break;
        case f__FDLmp: rddar.readArray( "FDLmp", &mtp->FDLmp[0][0], mtp->nFD, MAXSYMB);
            break;
        case f__PGT: rddar.readArray( "PGT", mtp->PGT, mtp->FIf*mtp->nPG);
            break;
        case f__MGPid: rddar.readArray( "MGPid", &mtp->MGPid[0][0], mtp->nPG, MAXSYMB);
            break;
        case f__UMGP: rddar.readArray( "UMGP", mtp->UMGP, mtp->FIf, 1);
            break;
        case f__BSF: rddar.readArray( "BSF", mtp->BSF, mtp->nSFD*mtp->Nf);
            break;
        case f__MB: rddar.readArray( "MB", mtp->MB, mtp->nC*mtp->Nf);
            break;
        case f__dMB: rddar.readArray( "dMB", mtp->dMB, mtp->nC*mtp->Nf);
            break;
        case f__xFlds: rddar.readArray( "xFlds", &mtp->xVTKfld[0][0], mtp->nVTKfld*2);
            break;
        case f__mDDc: rddar.readArray( "mDDc", mtp->DDc, mtp->Lsf);
            break;
        case f__mDIc: rddar.readArray( "mDIc", mtp->DIc, mtp->Nf);
            break;
        case f__mDEl: rddar.readArray( "mDEl", mtp->DEl, mtp->nEl);
            break;
        case f__for_e: rddar.readArray( "for_e", mtp->for_e[0],mtp->nEl, MAXFORMUNITDT);
            break;
        case f__nam_i: rddar.readArray( "nam_i", mtp->nam_i[0],mtp->nIV, MAXIDNAME);
            break;
        case f_SDref: rddar.readArray( "SDref",mtp->sdref[0],mtp->Nsd, V_SD_RKLEN);
            break;
        case f__SDval: rddar.readArray( "SDval", mtp->sdval[0],mtp->Nsd, V_SD_RKLEN);
            break;
        }
        nfild = rddar.findNext();
    }

    // testing read
    ret = rddar.testRead();
    if( !ret.empty() )
    { ret += " - fields must be read from gem2mt structure";
        Error( "Error", ret);
    }

}

#ifdef USE_NLOHMANNJSON
template void  TGEM2MT::to_text_file<io_formats::NlohmannJsonWrite>( io_formats::NlohmannJsonWrite& out_format, bool with_comments, bool brief_mode ) const;
template void  TGEM2MT::from_text_file<io_formats::NlohmannJsonRead>( io_formats::NlohmannJsonRead& out_format );
#else
template void  TGEM2MT::to_text_file<io_formats::SimdJsonWrite>( io_formats::SimdJsonWrite& out_format, bool with_comments, bool brief_mode ) const;
template void  TGEM2MT::from_text_file<io_formats::SimdJsonRead>( io_formats::SimdJsonRead& out_format );
#endif
template void  TGEM2MT::to_text_file<io_formats::KeyValueWrite>( io_formats::KeyValueWrite& out_format, bool with_comments, bool brief_mode ) const;
template void  TGEM2MT::from_text_file<io_formats::KeyValueRead>( io_formats::KeyValueRead& out_format );


// --------------------- end of m_gem2mtbox.cpp ---------------------------

