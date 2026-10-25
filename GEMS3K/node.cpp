//-------------------------------------------------------------------
// $Id$
//
/// \file node.cpp
/// Implementation of TNode class functionality including initialization
/// and execution of the GEM IPM 3 kernel
/// Works with DATACH and DATABR structures
//
// Copyright (c) 2005-2012 S.Dmytriyeva, D.Kulik, G.Kosakowski, F.Hingerl
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include "gems3k_timing.h"
#include "node.h"
#include "datach_api.h"
#include "num_methods.h"
#include "kinetics.h"
#include "v_service.h"
#include "jsonconfig.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include "lazy_file_sink.h"
#ifdef USE_OPTIMA_SOLVER
#include <Optima/Exception.hpp>
#endif


// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> TNode::node_logger = spdlog::stdout_color_mt("tnode");
// Thread-safe logger to file !!!!  to do change
std::shared_ptr<spdlog::logger> TNode::ipmlog_file;// = spdlog::stdout_color_mt("ipmlog");

// Conversion factors
const double bar_to_Pa = 1e5,
m3_to_cm3 = 1e6,
kg_to_g = 1e3;

void TNode::write_ThermoFun_format_stream(std::iostream &stream, bool compact) {
    stream << thermo_json_string;
}

std::once_flag once_flag;

// Constructor of the class instance in memory for standalone GEMS3K or coupled program
TNode::TNode()
{
    std::call_once(once_flag, [](){
        // Thread-safe one-time operation
    if(!ipmlog_file) {
        // ipmlog.txt is created when the first message is written, not here
        ipmlog_file = std::make_shared<spdlog::logger>("ipmlog",
                          gems3k_detail::make_lazy_file_sink(GemsSettings::with_directory("ipmlog.txt"),
                                                             GemsSettings::log_file_size,
                                                             GemsSettings::log_file_count));
        spdlog::register_logger(ipmlog_file);
    }
#ifdef USE_OPTIMA_SOLVER
    // Optima prints its warnings to std::cout unless told otherwise: send them to the "tnode" logger, so that
    // update_loggers(use_stdout, logfile, level) decides whether they reach the console, the log file, or nowhere
    Optima::setWarningHandler([](const std::string& message) {
        if(node_logger) node_logger->warn("Optima: {}", message);
    });
#endif
    });
    CSD = NULL;
    CNode = NULL;
    allocMemory();
    dbr_file_name = "dbr_file_name";
    load_thermodynamic_data = false;
}

TNode::TNode( const TNode& otherNode )
{
    CSD = 0;
    CNode = 0;

    allocMemory();
    dbr_file_name = otherNode.dbr_file_name;
    // copy data from otherNode
    datach_copy( otherNode.CSD );
    databr_copy( otherNode.CNode );

    multi_ptr()->copyMULTI( *otherNode.multi_base );
    // the targets go with the copy: the registered conditions and which of them came from the DBR
#ifdef USE_OPTIMA_SOLVER
    multi_ptr()->RestoreControlConditions( otherNode.multi_base->GetControlConditions() );
#endif
    dbrTarget_pH = otherNode.dbrTarget_pH;
    dbrTarget_Eh = otherNode.dbrTarget_Eh;
    dbrTargetWarned = otherNode.dbrTargetWarned;
    apiTarget_pH = otherNode.apiTarget_pH;
    apiTarget_Eh = otherNode.apiTarget_Eh;

    // copy intervals for minimizatiom
    pmm->Pai[0] = CSD->Pval[0]/bar_to_Pa;
    pmm->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
    pmm->Pai[2] = getStep( pmm->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
    pmm->Pai[3] = CSD->Ptol/bar_to_Pa;

    pmm->Tai[0] = CSD->TKval[0]-C_to_K;
    pmm->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
    pmm->Tai[2] = getStep( pmm->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
    pmm->Tai[3] = CSD->Ttol;

    pmm->Fdev1[0] = 0.;
    pmm->Fdev1[1] = 1e-6;   // 24/05/2010 must be copy from GEMS3 structure
    pmm->Fdev2[0] = 0.;
    pmm->Fdev2[1] = 1e-6;

    node_logger->debug("copy constructor...");
}

TNode::~TNode()
{
    clear_ThermoEngine();
    freeMemory();
}

void TNode::allocMemory()
{
    // memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;
    // mem_set( CSD, 0, sizeof(DATACH) );
    dbr_dch_api::datach_reset(CSD);
    // mem_set( CNode, 0, sizeof(DATABR) );
    dbr_dch_api::databr_reset(CNode, 2);

    // allocation internal structures
    internal_multi.reset(new TMultiBase(this));
    multi_base = internal_multi.get();
    multi_base->set_def();
    pmm = multi_base->GetPM();
    atp.reset( new TActivity( CSD, CNode, this ) );
    //    atp->set_def();
    kip.reset( new TKinetics( CSD, CNode, this ) );
    kip->set_def();
}

void TNode::freeMemory()
{
    dbr_dch_api::datach_free(CSD);
    // CSD = 0;
    delete CSD;
    CNode = databr_free(CNode);
}

std::string TNode::system_id() const
{
    return char_array_to_string(multi_ptr()->GetPM()->stkey, EQ_RKLEN);
}

// Checks if given temperature TK and pressure P fit within the interpolation
// intervals of the DATACH lookup arrays (returns true) or not (returns false)
bool  TNode::check_TP( double TK, double P ) const
{
    auto error_mesg = dbr_dch_api::check_TP(CSD, TK, P);
    if( !error_mesg.empty() ) {
        ipmlog_file->info("In node {}: {}", CNode->NodeHandle, error_mesg);
        return false;
    }
    return true;
}

void TNode::Get_sMod(int ndx, std::string &sMod)
{
    sMod = char_array_to_string(pmm->sMod[ndx], 8);
}

//-------------------------------------------------------------------------------------------------------------------------------
// (2) Main call for GEM IPM calculations using the input bulk composition, temperature, pressure
//   and metastability constraints provided in the work instance of DATABR structure.
//   Actual calculation will be performed only when dBR->NodeStatusCH == NEED_GEM_SIA (5) or dBR->NodeStatusCH = NEED_GEM_AIA (1).
//   By other values of NodeStatusCH, no calculation will be performed and the status will remain unchanged.
//  In "smart initial approximation" (SIA) mode, the program can automatically switch into the "automatic initial
//  approximation" (AIA) mode and return  OK_GEM_AIA instead of OK_GEM_SIA.
//  Parameter:
//   uPrimalSol  flag to define the mode of GEM smart initial approximation
//               (only if dBR->NodeStatusCH = NEED_GEM_SIA has been set before GEM_run() call).
//               false  (0) -  use speciation and activity coefficients from previous GEM_run() calculation
//               true  (1)  -  use speciation provided in the DATABR memory structure (e.g. after reading the DBR file)
//  Return values:    NodeStatusCH  (the same as set in dBR->NodeStatusCH). Possible values (see "databr.h" file for the full list)
// Defined here, not in ipm_optima.cpp, so a build without USE_OPTIMA_SOLVER still links.
thread_local bool g_optimaLineSearchRetry = false;
thread_local bool g_optimaLSEscapeOff = false;
thread_local bool g_optimaEscRetryNext = false;   // a failed attempt is re-run with the stall escape off (escretry)
thread_local bool g_optimaCgSeedArmed = false;
thread_local bool g_optimaTargetWarned = false;   // the unusual-target warning is shown once per GEM_run()
thread_local long int g_optimaIterBudget = -1;     // Optima iterations left for this call; -1 = no limit
thread_local bool g_optimaFixAccept = false;     // Optima takes the full step, clipped to the bounds
thread_local OptimaQuickTryWatch g_optimaQuickTry;   // the quick first try's frozen and divergence checks
extern thread_local bool g_aqueousAbsent;
extern thread_local bool g_ipmLevelUsed;   // pa_IpmLevelStart (ipm_main.cpp)
extern thread_local bool g_ipmNoProgArmed;   // pa_IpmNoProgress (ipm_main.cpp)
extern thread_local bool g_ipmNoProgFired;
// Arms pa_IpmNoProgress for one native solve.
struct NoProgArm
{
    explicit NoProgArm( bool on ) { g_ipmNoProgArmed = on; g_ipmNoProgFired = false; }
    ~NoProgArm() { g_ipmNoProgArmed = false; }
};
extern thread_local bool g_ipmLevelOff;
thread_local bool g_optimaPreSolveEmpty = false; // the reduced pre-solve found nothing earlier in this GEM_run()
thread_local long int g_lvlFullCalls = 0;   // pa_OptimaColdLevel / pa_OptimaLevelStart: cold Optima calls past the quick try in this GEM_run()
thread_local bool g_optimaCapCall = false;   // inside a capped Optima call (pa_OptimaCallCap)
thread_local long int g_optimaCapAllowance = 0;   // extra iterations for the repair step once the cap is used up; -1 = used
thread_local bool g_optimaDryRedoCall = false;   // inside the dry-retry redo (pa_OptimaDryRetry)
thread_local bool g_optimaDryFailed = false;     // the dry stop fired and the water-first retry failed
thread_local bool g_optimaDryOff = false;        // the dry stop and its retry are off for this call
bool TNode::SyncMultiWithDataBr()
{
    bool same = true;
    for( long int i = 0; i < CSD->nDCb && same; i++ )
        same = ( CNode->xDC[i] == pmm->X[ CSD->xdc[i] ] );
    for( long int i = 0; i < CSD->nICb && same; i++ )
        same = ( CNode->uIC[i] == pmm->U[ CSD->xic[i] ] );
    if( same )
        return false;
    unpackDataBr( true );
    for( long int i = 0; i < CSD->nICb; i++ )
        pmm->U[ CSD->xic[i] ] = CNode->uIC[i];
    multi_ptr()->RefreshFromAmounts();
    return true;
}

void TMultiBase::BeginCalculation()
{
    g_optimaTargetWarned = false;
    g_optimaPreSolveEmpty = false;
    g_lvlFullCalls = 0;
}

// Warn once when the system holds more than 1e6 mol of elements, more than the solvers are tested for.
void TNode::WarnIfSystemTooLarge()
{
    const double bigLimit = 1e6;
    double total = 0., largest = 0.;
    long int iLargest = -1;
    for( long int i = 0; i < CSD->nICb; i++ )
    {
        if( CSD->ccIC[ IC_xDB_to_xCH( i ) ] == IC_CHARGE ) continue;
        total += CNode->bIC[i];
        if( CNode->bIC[i] > largest ) { largest = CNode->bIC[i]; iLargest = i; }
    }
    if( total > bigLimit )
    {
        native_trace_decide( "bigsystem total=%.3e limit=%.3e", total, bigLimit );
        static thread_local bool bigWarned = false;
        if( !bigWarned && iLargest >= 0 )
            node_logger->warn( "GEM_run(): the system holds {:.3g} mol of elements ({:.3g} mol of {}), more than the {:.3g} mol "
                               "the solvers are tested for; the answer may be inaccurate", total, largest,
                               CSD->ICNL[ IC_xDB_to_xCH( iLargest ) ], bigLimit );
        bigWarned = true;
    }
}

// The layers of one call, outermost first: the call cap, the dry redo, then the quick first try and the
// solve with its retries. Each layer calls the next one.
long int TNode::GEM_run( bool uPrimalSol )
{
    WarnIfSystemTooLarge();
    const long int requested = CNode->NodeStatusCH;
    ApplyDbrTargets( requested );
    restoreEnteredComposition( requested );
    const std::vector<double> bBefore( CNode->bIC, CNode->bIC + CSD->nICb );
    long int st;
    {
        gems3k_timing::Scope t( gems3k_timing::Call );
        st = GEM_run_rescue( uPrimalSol );
    }
    rememberTitration( requested, bBefore );
    gems3k_timing::flush( requested, st, CNode->IterDone );
    return st;
}

bool TNode::modeUsesTargets( long int mode )
{
#ifdef USE_OPTIMA_SOLVER
    return mode == NEED_GEM_AOP || mode == NEED_GEM_SOP || mode == NEED_GEM_HOP || mode == NEED_GEM_SHP;
#else
    (void)mode;
    return false;
#endif
}

void TNode::restoreEnteredComposition( long int mode )
{
    if( titrSolved_.empty() )
        return;
    const size_t n = titrSolved_.size();
    bool untouched = n == (size_t)CSD->nICb;
    for( size_t i = 0; i < n && untouched; i++ )
        untouched = CNode->bIC[i] == titrSolved_[i];
    if( !untouched )
    {   // the caller set a new composition: nothing to put back
        titrEntered_.clear();
        titrSolved_.clear();
        return;
    }
#ifdef USE_OPTIMA_SOLVER
    if( modeUsesTargets( mode ) && !multi_ptr()->GetControlConditions().empty() )
        return;   // the target is applied again from this composition
#endif
    for( size_t i = 0; i < n; i++ )
        CNode->bIC[i] = titrEntered_[i];
    native_trace_decide( "titration-restore mode=%ld", mode );
    node_logger->info( "The last solve used a pH or Eh target and changed the bulk composition; mode {} does not use "
                       "a target, so the composition that was entered is solved", mode );
    titrEntered_.clear();
    titrSolved_.clear();
}

void TNode::rememberTitration( long int mode, const std::vector<double>& before )
{
#ifdef USE_OPTIMA_SOLVER
    if( !modeUsesTargets( mode ) || multi_ptr()->GetControlConditions().empty() )
        return;
    const size_t n = before.size();
    // a call that starts from the composition an earlier target solve wrote keeps that solve's entered composition
    const bool continuing = titrSolved_.size() == n && std::equal( before.begin(), before.end(), titrSolved_.begin() );
    if( !continuing )
        titrEntered_ = before;
    titrSolved_.assign( CNode->bIC, CNode->bIC + n );
#else
    (void)mode; (void)before;
#endif
}

void TNode::ApplyDbrTargets( long int mode )
{
    auto given = []( double v, const char* name ) {
        if( IsDoubleEmpty( v ) ) return false;
        if( !std::isfinite( v ) )
        {
            node_logger->warn( "The DBR {} target is not a finite number; it is not used", name );
            return false;
        }
        return true; };
    const bool pHgiven = given( CNode->pHtarget, "pH" ), Ehgiven = given( CNode->Ehtarget, "Eh" );
#ifdef USE_OPTIMA_SOLVER
    const bool used = mode == NEED_GEM_AOP || mode == NEED_GEM_SOP || mode == NEED_GEM_HOP || mode == NEED_GEM_SHP;
#else
    const bool used = false;   // without Optima no mode uses targets
    (void)mode;
#endif
    if( ( pHgiven || Ehgiven ) && !used && !dbrTargetWarned )
    {
        node_logger->warn( "The DBR has a pH or Eh target; only the Optima modes (AOP, SOP, HOP, SHP) use it" );
        dbrTargetWarned = true;
    }
    // a target is registered only for a mode that uses it, so the other modes run as without a target
    const bool pH = used && pHgiven, Eh = used && Ehgiven;
#ifdef USE_OPTIMA_SOLVER
    TMultiBase* m = multi_ptr();
    // a DBR target that this call does not use is removed; the target of the same name that the program
    // set (Set_pH_target(), Set_Eh_target()) and the DBR target replaced comes back
    auto expire = [m]( const char* name, bool was, bool now, std::vector<EqControlCondition>& saved ) {
        if( !was || now )
            return;
        std::vector<EqControlCondition> keep;
        for( const auto& c : m->GetControlConditions() )
            if( c.name != name ) keep.push_back( c );
        for( const auto& c : saved ) keep.push_back( c );
        saved.clear();
        m->RestoreControlConditions( keep );
    };
    expire( "pH", dbrTarget_pH, pH, apiTarget_pH );
    expire( "Eh", dbrTarget_Eh, Eh, apiTarget_Eh );
    // keep the program's target of the same name before a DBR target replaces it
    auto keepApi = [m]( const char* name, bool was, bool now, std::vector<EqControlCondition>& saved ) {
        if( was || !now )
            return;
        saved.clear();
        for( const auto& c : m->GetControlConditions() )
            if( c.name == name ) saved.push_back( c );
    };
    keepApi( "pH", dbrTarget_pH, pH, apiTarget_pH );
    keepApi( "Eh", dbrTarget_Eh, Eh, apiTarget_Eh );
    dbrTarget_pH = pH;
    dbrTarget_Eh = Eh;
    if( pH || Eh )
    {
        // the target setters find the H, charge and H+ indices through the node
        TNode* const node0 = m->node1;
        m->node1 = this;
        try {
            if( pH ) m->SetControlCondition_pH( CNode->pHtarget );
            if( Eh ) m->SetControlCondition_Eh( CNode->Ehtarget );
        }
        catch( ... ) { m->node1 = node0; throw; }
        m->node1 = node0;
    }
#endif
}

// pa_FailRescue: a call that ends with an error or a bad result is tried again in steps, stopping at the
// first that passes: (1) a cold call solved warm from the node's previous answer, (2) nudged compositions,
// (3) a walk from the node's previous answer to the new composition, (4) looser settings, (5) the other
// solver from the best point so far. Each passing answer comes from a solve at the exact composition with
// the normal settings. With pa_FailRescue = 2, if none passes, the attempt that best balances the
// composition (then lowest G) is returned as a bad result; if
// none balances (or the DATABR lacks species or elements, so balance cannot be shown), the original failure is returned. A kinetic node is finished with this step's kinetic
// limits held fixed, so the kinetics are not applied twice. pa_RescueBudget: after the call has used that many
// iterations, steps 1 to 4 start no new attempt; an attempt that passed still gets its finish, and step 5, the
// other solver, still runs once. DECIDE rescue, rescue-result, rescue-budget.
long int TNode::GEM_run_rescue( bool uPrimalSol )
{
    rescueStep = 0;
    BASE_PARAM* bp = multi_ptr()->base_param();
    const long int requested = CNode->NodeStatusCH;
    const long int level = bp->FailRescue;
    const bool sixModes = requested == NEED_GEM_AIA || requested == NEED_GEM_SIA || requested == NEED_GEM_AOP
                       || requested == NEED_GEM_SOP || requested == NEED_GEM_HOP || requested == NEED_GEM_SHP;
    if( level <= 0 || !sixModes )
        return GEM_run_capped( uPrimalSol );

    auto snap = [this]() {
        DATABR* d = new DATABR;
        dbr_dch_api::databr_reset( d, 1 );
        dbr_dch_api::databr_realloc( CSD, d );
        { DATABR* live = CNode; CNode = d; databr_copy( live ); CNode = live; }
        return d; };
    auto isOk = []( long int s ) { return s == OK_GEM_AIA || s == OK_GEM_SIA || s == OK_GEM_AOP || s == OK_GEM_SOP
                                          || s == OK_GEM_HOP || s == OK_GEM_SHP; };
    auto isBad = []( long int s ) { return s == BAD_GEM_AIA || s == BAD_GEM_SIA || s == BAD_GEM_AOP || s == BAD_GEM_SOP
                                           || s == BAD_GEM_HOP || s == BAD_GEM_SHP; };

    DATABR* D0 = snap();   // the call's inputs
    const long int ITau0 = pmm->ITau, pKMM0 = pmm->pKMM;   // kinetic state before the first attempt
    long int st = GEM_run_capped( uPrimalSol );
    // a warm call may report its outcome under the cold code of its solver, so any bad or error code counts
    const bool isErr = st == ERR_GEM_AIA || st == ERR_GEM_SIA || st == ERR_GEM_AOP || st == ERR_GEM_SOP
                    || st == ERR_GEM_HOP || st == ERR_GEM_SHP;
    if( !isBad( st ) && !isErr )
    {
        databr_free( D0 );
        return st;
    }

    const long int nIC = CSD->nICb;
    const std::vector<double> b0( CNode->bIC, CNode->bIC + nIC );
    const bool kinetics = D0->dt > 0.;
    const double dt0 = D0->dt;
    const long int warmMode = ( requested == NEED_GEM_AIA || requested == NEED_GEM_SIA ) ? NEED_GEM_SIA
                            : ( requested == NEED_GEM_AOP || requested == NEED_GEM_SOP ) ? NEED_GEM_SOP : NEED_GEM_SHP;
    DATABR* failed = snap();
    const std::string failedError = ipmlog_error;
    // pH/Eh control results live outside the DATABR, so they are saved with every answer that may be returned
    const std::vector<EqControlCondition> failedCC = multi_ptr()->GetControlConditions();
    long int itf = pmm->ITF, itg = pmm->ITG;
    long int itTotal = isErr ? pmm->ITF + pmm->ITG : CNode->IterDone;
    double seconds = CalcTime;
    long int nAttempts = 0;

    struct Entry { DATABR* d; double rel; double G; long int step; std::vector<EqControlCondition> cc; bool targetsMet; };
    std::vector<Entry> pool;
    auto addToPool = [&]( long int s, long int step ) {
        if( !isOk( s ) && !isBad( s ) )
            return;   // an error leaves no packed answer
        double rel = 0.;
        NudgedAnswerBalances( CNode, b0, rel );
        // under pH/Eh control an attempt counts only if the answer itself meets every target (a hybrid's BAD answer is
        // native's, while the control record is its Optima leg's, so the record alone cannot tell)
        bool met = true;
        for( const auto& c : multi_ptr()->GetControlConditions() )
        {
            const double tol = c.tolerance > 0. ? c.tolerance : ( c.defaultToleranceFn ? c.defaultToleranceFn() : 1e-3 );
            if( c.name == "pH" ) met = met && fabs( CNode->pH - c.target ) <= tol;
            else if( c.name == "Eh" ) met = met && fabs( CNode->Eh - c.target ) <= tol;
            else met = met && c.targetMet;
        }
        pool.push_back( { snap(), rel, CNode->Gs, step, multi_ptr()->GetControlConditions(), met } ); };
    if( isBad( st ) )
        addToPool( st, 0 );

    // One solve of the current CNode; a kinetic node with frozenBox is solved without the kinetics step,
    // inside the kinetic limits already in CNode.
    auto run = [&]( long int mode, bool primal, bool frozenBox ) {
        const long int ITau = pmm->ITau, pKMM = pmm->pKMM;
        const double kTau = pmm->kTau, kdT = pmm->kdT;
        if( kinetics && frozenBox ) CNode->dt = 0.;
        if( kinetics && !frozenBox ) { pmm->ITau = ITau0; pmm->pKMM = pKMM0; }   // repeat this time step, not the next one
        CNode->NodeStatusCH = mode;
        CNode->IterDone = 0;   // a failed solve does not write it, so a loaded count is not counted again
        const long int s = GEM_run_capped( primal );
        if( kinetics && frozenBox )
        {
            CNode->dt = dt0;
            pmm->ITau = ITau; pmm->pKMM = pKMM; pmm->kTau = kTau; pmm->kdT = kdT;
        }
        itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
        itTotal += CNode->IterDone > 0 ? CNode->IterDone : pmm->ITF + pmm->ITG;
        nAttempts++;
        return s; };
    auto setB = [&]( const std::vector<double>& b ) {
        for( long int i = 0; i < nIC; i++ ) CNode->bIC[i] = b[(size_t)i]; };
    // The answer now in CNode, finished warm at the exact composition with the normal settings.
    auto finish = [&]( long int step ) {
        setB( b0 );
        const long int s = run( warmMode, true, true );
        addToPool( s, step );
        return isOk( s ); };

    native_trace_quiet( true );
    bool passed = false;
    // pa_RescueBudget: true once the call has used its iterations; no new attempt of steps 1-4 starts then
    bool budgetUsed = false;
    auto overBudget = [&]() {
        if( !budgetUsed && bp->RescueBudget > 0 && itTotal >= bp->RescueBudget )
        {
            budgetUsed = true;
            native_trace_decide( "rescue-budget it=%ld budget=%ld attempts=%ld", itTotal, bp->RescueBudget, nAttempts );
        }
        return budgetUsed; };

    // (1) a cold call with a previous answer in its inputs is solved warm from it, as a transport code would
    double xPrev = 0.;
    for( long int jj = 0; jj < CSD->nDCb; jj++ ) xPrev += D0->xDC[jj];
    const bool coldReq = requested == NEED_GEM_AIA || requested == NEED_GEM_AOP || requested == NEED_GEM_HOP;
    if( coldReq && xPrev > 0. && !overBudget() )
    {
        databr_copy( D0 );
        const long int s = run( warmMode, true, false );
        addToPool( s, 1 );
        if( isOk( s ) ) { passed = true; rescueStep = 1; }
        native_trace_decide( "rescue step=1 mode=%ld status=%ld it=%ld", warmMode, s, itTotal );
    }

    // (2) nudged compositions, as pa_ColdRetryNudges: bIC[i] x (1 + s_i k 1e-15), k = +1, -1, +2, -2, ...
    for( long int n = 0; n < bp->ColdRetryNudges && !passed && !overBudget(); n++ )
    {
        const double k = ( n % 2 == 0 ? 1. : -1. ) * (double)( n / 2 + 1 );
        databr_copy( D0 );
        for( long int i = 0; i < nIC; i++ )
            CNode->bIC[i] = b0[(size_t)i] * ( 1. + ( ( i % 2 ) ? -1. : 1. ) * k * 1e-15 );
        const long int s = run( requested, uPrimalSol, false );
        long int f = -1;
        if( isOk( s ) )
        {
            addToPool( s, 2 );
            passed = finish( 2 );
            f = CNode->NodeStatusCH;
            if( passed ) rescueStep = 2;
        }
        native_trace_decide( "rescue step=2 attempt=%ld nudge=%+g status=%ld finish=%ld it=%ld", n + 1, k, s, f, itTotal );
    }

    // (3) walk from the bulk the previous answer balances to the new one, in equal warm steps; a failed step
    // is halved, at most three times. A kinetic node walks inside the kinetic limits of the first attempt
    // that produced them, and skips the walk if none did.
    if( !passed && !overBudget() && bp->RescueWalkSteps > 0 && CSD->nDCb == CSD->nDC && CSD->nICb == CSD->nIC )
    {
        std::vector<double> bPrev( (size_t)nIC, 0. );
        double xSum = 0.;
        for( long int jj = 0; jj < CSD->nDCb; jj++ ) xSum += D0->xDC[jj];
        for( long int i = 0; i < nIC; i++ )
            for( long int jj = 0; jj < CSD->nDCb; jj++ )
                bPrev[(size_t)i] += CSD->A[ CSD->xdc[jj] * CSD->nIC + CSD->xic[i] ] * D0->xDC[jj];
        double dist = 0.;
        for( long int i = 0; i < nIC; i++ )
            if( CSD->ccIC[ CSD->xic[i] ] != IC_CHARGE )
                dist = std::max( dist, fabs( b0[(size_t)i] - bPrev[(size_t)i] ) / std::max( fabs( b0[(size_t)i] ), 1e-300 ) );
        const DATABR* boxFrom = kinetics ? ( pool.empty() ? nullptr : pool.back().d ) : D0;
        if( xSum > 0. && dist > 1e-12 && boxFrom )
        {
            DATABR* good = snap();
            databr_copy( D0 );
            for( long int jj = 0; jj < CSD->nDCb; jj++ ) { CNode->dul[jj] = boxFrom->dul[jj]; CNode->dll[jj] = boxFrom->dll[jj]; }
            { DATABR* live = CNode; CNode = good; databr_copy( live ); CNode = live; }
            double tGood = 0., h = 1. / (double)bp->RescueWalkSteps;
            long int halvings = 0, steps = 0, s = -1;
            while( tGood < 1. && !overBudget() )
            {
                const double t = std::min( 1., tGood + h );
                databr_copy( good );
                std::vector<double> b( (size_t)nIC );
                for( long int i = 0; i < nIC; i++ )
                    b[(size_t)i] = ( t >= 1. ) ? b0[(size_t)i] : bPrev[(size_t)i] + t * ( b0[(size_t)i] - bPrev[(size_t)i] );
                setB( b );
                s = run( warmMode, true, true );
                steps++;
                if( isOk( s ) )
                {
                    tGood = t;
                    { DATABR* live = CNode; CNode = good; databr_copy( live ); CNode = live; }
                }
                else if( t >= 1. && isBad( s ) )
                {
                    addToPool( s, 3 );
                    break;
                }
                else if( ++halvings > 3 )
                    break;
                else
                    h *= 0.5;
            }
            if( tGood >= 1. ) { passed = true; rescueStep = 3; }
            databr_free( good );
            native_trace_decide( "rescue step=3 steps=%ld halvings=%ld reached=%.4f dist=%.3e status=%ld it=%ld",
                                 steps, halvings, tGood, dist, s, itTotal );
        }
        else
            native_trace_decide( "rescue step=3 skipped prev=%d dist=%.3e box=%d", xSum > 0. ? 1 : 0, dist, boxFrom ? 1 : 0 );
    }

    // A kinetic node where no attempt so far produced any result: its kinetic limits are likely not reachable
    // from this composition, which neither looser settings nor the other solver change, so steps 4 and 5 are skipped.
    const bool noResultKinetic = kinetics && pool.empty();
    if( !passed && noResultKinetic )
        native_trace_decide( "rescue skip=4,5 reason=kinetic-no-result attempts=%ld", nAttempts );

    // (4) one solve with looser settings, then the finish with the normal ones.
    if( !passed && !noResultKinetic && !overBudget() )
    {
        databr_copy( D0 );
        TMultiBase::LocalSettings ownSettings( *multi_ptr() );   // the other nodes keep their settings
        BASE_PARAM* bp4 = multi_ptr()->base_param();
        const BASE_PARAM saved = *bp4;
        bp4->DK *= 10.;
        bp4->OptimaTol = 1e-6;
        bp4->OptimaDcFloor = 1e-12;
        bp4->LogBarrierTau = 1e-12;
        const long int s = run( requested, uPrimalSol, false );
        *bp4 = saved;
        long int f = -1;
        if( isOk( s ) || isBad( s ) )
        {
            addToPool( s, 4 );   // kept as a candidate in case its finish gives no result
            passed = finish( 4 );
            f = CNode->NodeStatusCH;
            if( passed ) rescueStep = 4;
        }
        native_trace_decide( "rescue step=4 status=%ld finish=%ld it=%ld", s, f, itTotal );
    }

    // (5) the other solver, warm at the exact composition from the best balanced answer so far (or from the inputs).
    auto bestIndex = [&]() {
        long int best = -1;
        // a DATABR without every species and element cannot show that an attempt balances, so none is offered
        if( CSD->nDCb != CSD->nDC || CSD->nICb != CSD->nIC )
            return best;
        for( size_t e = 0; e < pool.size(); e++ )
        {
            if( pool[e].rel > 1. || !pool[e].targetsMet ) continue;
            if( best < 0 || pool[e].G < pool[(size_t)best].G ) best = (long int)e;
        }
        return best; };
    if( !passed && !noResultKinetic )
    {
        std::vector<long int> others;
        // native does not apply pH/Eh control conditions, so a controlled call never goes to it
        const bool controlled = !multi_ptr()->GetControlConditions().empty();
        if( warmMode == NEED_GEM_SIA ) others.push_back( NEED_GEM_SOP );
        else if( warmMode == NEED_GEM_SOP ) { if( !controlled ) others.push_back( NEED_GEM_SIA ); }
        else { if( !controlled ) others.push_back( NEED_GEM_SIA ); others.push_back( NEED_GEM_SOP ); }
        // step 5 runs even when the budget is used: it is one solve per solver, and often the one that passes
        for( long int code : others )
        {
            if( passed ) break;
            const long int other = code;
            const long int bi = bestIndex();
            long int s;
            if( bi >= 0 )
            {
                databr_copy( pool[(size_t)bi].d );
                setB( b0 );
                s = run( other, true, true );
            }
            else
            {
                databr_copy( D0 );
                s = run( other, uPrimalSol, false );
            }
            addToPool( s, 5 );
            if( isOk( s ) ) { passed = true; rescueStep = 5; }
            native_trace_decide( "rescue step=5 solver=%ld from=%s status=%ld it=%ld", other, bi >= 0 ? "best" : "inputs", s, itTotal );
        }
    }
    native_trace_quiet( false );

    long int result;
    if( passed )
    {
        result = requested + 1;
        CNode->NodeStatusCH = result;
        native_trace_run_result( *pmm, requested, result, multi_base );
    }
    else
    {
        const long int bi = ( level >= 2 ) ? bestIndex() : -1;
        if( bi >= 0 )
        {
            databr_copy( pool[(size_t)bi].d );
            multi_ptr()->RestoreControlConditions( pool[(size_t)bi].cc );
            setB( b0 );
            CNode->dt = dt0;
            result = requested + 2;
            CNode->NodeStatusCH = result;
            rescueStep = 6;
            node_logger->warn( "GEM_run(): mode {} failed ({}); every rescue step failed, returned the best balanced attempt "
                               "(mass balance {:.3g} of the tolerance) as a bad result", requested, st, pool[(size_t)bi].rel );
        }
        else
        {
            databr_copy( failed );
            multi_ptr()->RestoreControlConditions( failedCC );
            ipmlog_error = failedError;
            result = st;
        }
    }
    native_trace_decide( "rescue-result requested=%ld first=%ld step=%ld status=%ld attempts=%ld pool=%ld it=%ld budget=%d",
                         requested, st, rescueStep, result, nAttempts, (long)pool.size(), itTotal, budgetUsed ? 1 : 0 );
    node_logger->log( passed ? spdlog::level::debug : spdlog::level::info, "GEM_run(): mode {} failed ({}) - rescue: {}", requested, st,
                       passed ? "recovered" : ( rescueStep == 6 ? "best attempt returned" : "still failed" ) );
    if( budgetUsed && !passed )
        node_logger->info( "GEM_run(): the rescue stopped at its iteration limit (pa_RescueBudget = {})", bp->RescueBudget );
    for( auto& e : pool ) databr_free( e.d );
    databr_free( failed );
    databr_free( D0 );
    pmm->ITF = itf; pmm->ITG = itg; NumIterFIA = itf; NumIterIPM = itg;
    CNode->IterDone = itTotal; CalcTime = seconds;
    return result;
}

// pa_OptimaCallCap: an AOP or SOP call gets a whole-call iteration cap. A call that fails with
// the cap used up is re-solved at up to pa_OptimaCapNudges 1e-15 nudges of the bulk
// composition, each with a fresh cap; a nudged answer is finished warm at the exact
// composition. If no nudge solves, the call is run once more without the cap, so the cap
// never loses an answer. DECIDE callcap, capretry, caplast.
long int TNode::GEM_run_capped( bool uPrimalSol )
{
    const BASE_PARAM* capPa = multi_ptr()->base_param();
    const long int cap = optima_call_cap( capPa->OptimaCallCap, pmm->L );
    if( cap > 0 && !g_optimaCapCall && CNode->dt <= 0.
        && ( CNode->NodeStatusCH == NEED_GEM_AOP || CNode->NodeStatusCH == NEED_GEM_SOP ) )
    {
        const long int req = CNode->NodeStatusCH;
        const long int capNudges = capPa->OptimaCapNudges;
        g_optimaCapCall = true;
        g_optimaIterBudget = cap;
        g_optimaCapAllowance = cap / 2;
        long int st = GEM_run_dry_redo( uPrimalSol );
        long int total = CNode->IterDone;
        const bool capHit = ( g_optimaIterBudget == 0 || g_optimaCapAllowance < 0 );
        const bool okFirst = ( st == OK_GEM_AOP || st == OK_GEM_SOP );
        native_trace_decide( "callcap cap=%ld status=%ld it=%ld hit=%d", cap, (long)st, (long)total, capHit ? 1 : 0 );
        if( !okFirst && capHit )
        {
            const long int nIC = CSD->nICb;
            DATABR* failed = new DATABR;
            dbr_dch_api::databr_reset( failed, 1 );
            dbr_dch_api::databr_realloc( CSD, failed );
            { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
            const std::string failedError = ipmlog_error;
            const std::vector<double> b0( CNode->bIC, CNode->bIC + nIC );
            long int result = st;
            bool recovered = false;
            for( long int n = 0; n < capNudges && !recovered; n++ )
            {
                const double k = ( n % 2 == 0 ? 1. : -1. ) * (double)( n / 2 + 1 );
                databr_copy( failed );
                for( long int i = 0; i < nIC; i++ )
                    CNode->bIC[i] = b0[(size_t)i] * ( 1. + ( ( i % 2 ) ? -1. : 1. ) * k * 1e-15 );
                CNode->NodeStatusCH = NEED_GEM_AOP;
                g_optimaIterBudget = cap;
                g_optimaCapAllowance = cap / 2;
                const long int cold = GEM_run_dry_redo( uPrimalSol );
                total += CNode->IterDone;
                long int warm = -1;
                if( cold == OK_GEM_AOP )
                {
                    DATABR* nudgedOk = new DATABR;
                    dbr_dch_api::databr_reset( nudgedOk, 1 );
                    dbr_dch_api::databr_realloc( CSD, nudgedOk );
                    { DATABR* live = CNode; CNode = nudgedOk; databr_copy( live ); CNode = live; }
                    for( long int i = 0; i < nIC; i++ )
                        CNode->bIC[i] = b0[(size_t)i];
                    CNode->NodeStatusCH = NEED_GEM_SOP;
                    g_optimaIterBudget = cap;
                    g_optimaCapAllowance = cap / 2;
                    warm = GEM_run_dry_redo( uPrimalSol );
                    total += CNode->IterDone;
                    bool keep = true;
                    if( warm != OK_GEM_SOP && warm != OK_GEM_AOP )
                    {   // keep the nudged answer, reported at the exact composition, if it balances that composition
                        double worstRel = 0.;
                        keep = NudgedAnswerBalances( nudgedOk, b0, worstRel );
                        if( keep )
                        {
                            databr_copy( nudgedOk );
                            for( long int i = 0; i < nIC; i++ )
                                CNode->bIC[i] = b0[(size_t)i];
                            native_trace_decide( "capretry-keepnudged nudge=%+g mb_rel=%.3e", k, worstRel );
                        }
                        else
                            native_trace_decide( "capretry-keepnudged-refused nudge=%+g mb_rel=%.3e", k, worstRel );
                    }
                    if( keep )
                    {
                        result = ( req == NEED_GEM_AOP ) ? OK_GEM_AOP : OK_GEM_SOP;
                        recovered = true;
                    }
                    databr_free( nudgedOk );
                }
                native_trace_decide( "capretry attempt=%ld nudge=%+g aop=%ld sop=%ld it=%ld",
                                     (long)( n + 1 ), k, (long)cold, (long)warm, (long)total );
            }
            if( !recovered )
            {   // last try: the original call without the cap
                databr_copy( failed );
                for( long int i = 0; i < nIC; i++ )
                    CNode->bIC[i] = b0[(size_t)i];
                CNode->NodeStatusCH = req;
                g_optimaIterBudget = -1;
                g_optimaCapAllowance = 0;
                const long int last = GEM_run_dry_redo( uPrimalSol );
                total += CNode->IterDone;
                native_trace_decide( "caplast status=%ld it=%ld", (long)last, (long)total );
                if( last == OK_GEM_AOP || last == OK_GEM_SOP )
                {
                    result = last;
                    recovered = true;
                }
            }
            if( !recovered )
            {
                databr_copy( failed );
                ipmlog_error = failedError;
                result = st;
            }
            CNode->NodeStatusCH = result;
            st = result;
            databr_free( failed );
            node_logger->info( "GEM_run(): Optima call reached its iteration cap ({}) - retried at nudged compositions"
                               " and without the cap: {}", cap, recovered ? "recovered" : "still failed" );
        }
        CNode->IterDone = total;
        g_optimaIterBudget = -1;
        g_optimaCapAllowance = 0;
        g_optimaCapCall = false;
        return st;
    }
    return GEM_run_dry_redo( uPrimalSol );
}

// pa_OptimaDryRetry: if the dry stop fired and the water-first retry failed, the low water may be
// the real answer, so the call is run again with the dry stop and its retry off. DECIDE dryredo.
long int TNode::GEM_run_dry_redo( bool uPrimalSol )
{
    if( multi_ptr()->base_param()->OptimaDryRetry && !g_optimaDryRedoCall && !g_optimaDryOff && CNode->dt <= 0.
        && ( CNode->NodeStatusCH == NEED_GEM_AOP || CNode->NodeStatusCH == NEED_GEM_SOP ) )
    {
        const long int req = CNode->NodeStatusCH;
        g_optimaDryRedoCall = true;
        g_optimaDryFailed = false;
        const long int budget0 = g_optimaIterBudget;
        const long int allowance0 = g_optimaCapAllowance;
        long int st = GEM_run_tries( uPrimalSol );
        long int total = CNode->IterDone;
        if( st != OK_GEM_AOP && st != OK_GEM_SOP && g_optimaDryFailed )
        {
            g_optimaIterBudget = budget0 < 0 ? -1 : std::max( 0L, budget0 - total );
            g_optimaCapAllowance = allowance0;
            g_optimaDryOff = true;
            CNode->NodeStatusCH = req;
            st = GEM_run_tries( uPrimalSol );
            total += CNode->IterDone;
            g_optimaDryOff = false;
            native_trace_decide( "dryredo status=%ld it=%ld", (long)st, (long)total );
        }
        CNode->IterDone = total;
        g_optimaIterBudget = budget0 < 0 ? -1 : std::max( 0L, budget0 - total );   // undo the forced stop
        g_optimaDryRedoCall = false;
        return st;
    }
    return GEM_run_tries( uPrimalSol );
}

// The quick first try, then the solve with its retries, then the notes on the answer.
long int TNode::GEM_run_tries( bool uPrimalSol )
{
    const long int requested = CNode->NodeStatusCH;
    const bool kinetics = CNode->dt > 0.;
    TMultiBase::BeginCalculation();
    // pa_OptimaFastFirst > 0: a cold Optima call first tries a quick solve - the full step clipped to the
    // bounds, plain settings and a small iteration cap. Anything but OK sends the call back to its inputs
    // for the normal solve, and both tries' iterations are counted. DECIDE fastfirst.
    long int fastF = 0, fastG = 0;
    double fastSec = 0.;
    bool fastFailed = false;
    const long int fastCap = multi_ptr()->base_param()->OptimaFastFirst;
#ifdef USE_OPTIMA_SOLVER   // the quick try is an Optima solve
    if( requested == NEED_GEM_AOP && !kinetics && fastCap > 0 )
    {
        DATABR* fastBefore = new DATABR;
        dbr_dch_api::databr_reset( fastBefore, 1 );
        dbr_dch_api::databr_realloc( CSD, fastBefore );
        { DATABR* live = CNode; CNode = fastBefore; databr_copy( live ); CNode = live; }
        TMultiBase::LocalSettings ownSettings( *multi_ptr() );   // the other nodes keep their settings
        BASE_PARAM* bp = multi_ptr()->base_param();
        const BASE_PARAM saved = *bp;
        bp->DG = 1e-5;                  // no rescaling
        bp->OptimaDcFloor = 1e-16;
        bp->LogBarrierTau = 1e-16;
        bp->OptimaLineSearch = 0.;
        bp->OptimaStallWindow = 0;
        bp->OptimaDimReduce = -1;
        bp->OptimaFDHessianDelay = 0;
        const long int budget0 = g_optimaIterBudget;
        g_optimaQuickTry = OptimaQuickTryWatch();
        g_optimaIterBudget = budget0 < 0 ? fastCap : std::min( fastCap, budget0 );
        const long int budgetSet = g_optimaIterBudget;
        const long int allowance0 = g_optimaCapAllowance;   // the repair allowance is not for this attempt
        g_optimaCapAllowance = 0;
        g_optimaFixAccept = true;
        g_errorsHandledByCaller++;   // a failed quick try is followed by the full solve
        const long int fast = GEM_run_single( uPrimalSol );
        g_errorsHandledByCaller--;
        g_optimaFixAccept = false;
        g_optimaCapAllowance = allowance0;
        *bp = saved;
        const long int used = budgetSet - g_optimaIterBudget;
        // Only a quick try that finished within its cap counts; one cut at the cap goes to the normal solve,
        // unless it ended OK and no missing phase would form (pa_OptimaFastAcceptCut).
        bool fastCut = ( used >= budgetSet );
        if( fastCut && fast == OK_GEM_AOP && saved.OptimaFastAcceptCut > 0 )
        {
            long int wk = -1;
            const double worst = multi_ptr()->AbsentPhaseWorstDrive( wk );
            // A mixed phase with a twin copy (a possible miscibility gap): two copies present, or an absent copy with
            // any drive to form, is left to the full solve.
            const std::vector<char> twin = multi_ptr()->PhaseHessianTwins();
            long int twinsPresent = 0;
            double tot = 0., twinDrive = 0.;
            for( long int k = 0; k < pmm->FI; k++ ) tot += pmm->XF[k];
            for( long int k = 0, jb = 0; k < pmm->FIs && k < (long int)twin.size(); jb += pmm->L1[k], k++ )
            {
                if( !twin[(size_t)k] ) continue;
                if( pmm->XF[k] > 1e-9 * tot ) { twinsPresent++; continue; }
                extern thread_local bool g_tpdStrongSearch;
                std::vector<double> yb;
                g_tpdStrongSearch = true;
                twinDrive = std::min( twinDrive, multi_ptr()->NativeTpdPhase( k, jb, yb ) );
                g_tpdStrongSearch = false;
            }
            const bool keep = !( worst < -1e-6 ) && twinsPresent < 2 && !( twinsPresent > 0 && twinDrive < -1e-10 );
            native_trace_decide( "fastcutaccept keep=%d worst=%.3e phase=%ld twins=%ld twindrive=%.3e", keep ? 1 : 0, worst, wk,
                                 twinsPresent, twinDrive );
            if( keep ) fastCut = false;
        }
        native_trace_decide( "fastfirst cap=%ld status=%ld it=%ld used=%ld cut=%d", fastCap, (long)fast,
                             (long)( pmm->ITF + pmm->ITG ), used, fastCut ? 1 : 0 );
        if( fast == OK_GEM_AOP && !fastCut )
        {
            g_optimaIterBudget = budget0;
            databr_free( fastBefore );
            return fast;
        }
        fastF = pmm->ITF; fastG = pmm->ITG; fastSec = CalcTime; fastFailed = true;
        g_optimaIterBudget = budget0 < 0 ? -1 : std::max( 0L, budget0 - used );
        databr_copy( fastBefore );
        databr_free( fastBefore );
        CNode->NodeStatusCH = requested;
        clearipmLogError();
    }
#endif
    // pa_OptimaLineSearch < 0: an Optima call that fails without the line search is re-run once
    // with it (g_optimaLineSearchRetry). The DATABR is snapshotted before the first attempt; a
    // retry that does not return OK hands the first failure back unchanged. DECIDE lsretry.
    // pa_OptimaLSStallEscape > 0 with the line search on: a failed call is re-run once with the
    // escape off (g_optimaLSEscapeOff), same snapshot/restore. DECIDE escretry.
    const bool optimaReq = !kinetics &&
        ( requested == NEED_GEM_AOP || requested == NEED_GEM_SOP || requested == NEED_GEM_HOP || requested == NEED_GEM_SHP );
    const bool lsRetryArmed = optimaReq && multi_ptr()->base_param()->OptimaLineSearch < 0.;
    // With the line search on, the objective memo is always active, so the legacy retry -
    // escape off and memo off - is armed whenever the line search is on.
    // With pa_OptimaColdRetry = 2 a failed warm call goes straight to the cold retry, so it gets no second warm try.
    const bool warmGoesCold = ( requested == NEED_GEM_SOP || requested == NEED_GEM_SHP )
                              && multi_ptr()->base_param()->OptimaColdRetry == 2;
    const bool escRetryArmed = optimaReq && !lsRetryArmed && !warmGoesCold && multi_ptr()->base_param()->OptimaLineSearch > 0.;
    DATABR* lsBefore = nullptr;
    if( lsRetryArmed || escRetryArmed )
    {
        lsBefore = new DATABR;
        dbr_dch_api::databr_reset( lsBefore, 1 );
        dbr_dch_api::databr_realloc( CSD, lsBefore );
        { DATABR* live = CNode; CNode = lsBefore; databr_copy( live ); CNode = live; }
    }
    // pa_IpmLevelStart: keep the inputs of a call with a native part, in case its start from the leveled solution fails.
    DATABR* lvlBefore = nullptr;
    g_ipmLevelUsed = false;
    if( !kinetics && ( requested == NEED_GEM_AIA || requested == NEED_GEM_SIA || requested == NEED_GEM_HOP || requested == NEED_GEM_SHP )
        && multi_ptr()->base_param()->IpmLevelStart > 0 )
    {
        lvlBefore = new DATABR;
        dbr_dch_api::databr_reset( lvlBefore, 1 );
        dbr_dch_api::databr_realloc( CSD, lvlBefore );
        { DATABR* live = CNode; CNode = lvlBefore; databr_copy( live ); CNode = live; }
    }
    g_optimaEscRetryNext = escRetryArmed;   // the first attempt leaves its last-resort re-solve to that retry
    long int status = GEM_run_single( uPrimalSol );
    g_optimaEscRetryNext = false;
    if( warmGoesCold && optimaReq && !lsRetryArmed && multi_ptr()->base_param()->OptimaLineSearch > 0.
        && status != OK_GEM_SOP && status != OK_GEM_SHP )
        native_trace_decide( "escretry skipped=warm requested=%ld first=%ld", (long)requested, (long)status );
    auto isOkStatus =[]( long int st ) { return st == OK_GEM_AOP || st == OK_GEM_SOP || st == OK_GEM_HOP || st == OK_GEM_SHP
                                                || st == OK_GEM_AIA || st == OK_GEM_SIA; };
    if( lsBefore )
    {
        const bool okNow = status == OK_GEM_AOP || status == OK_GEM_SOP || status == OK_GEM_HOP || status == OK_GEM_SHP
                        || status == OK_GEM_AIA || status == OK_GEM_SIA;
        if( !okNow )
        {
            DATABR* failed = new DATABR;
            dbr_dch_api::databr_reset( failed, 1 );
            dbr_dch_api::databr_realloc( CSD, failed );
            { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
            const std::string failedError = ipmlog_error;
            long int itf = pmm->ITF, itg = pmm->ITG; double seconds = CalcTime;
            databr_copy( lsBefore );                  // the first attempt's INPUTS
            CNode->NodeStatusCH = requested;
            if( lsRetryArmed ) g_optimaLineSearchRetry = true; else g_optimaLSEscapeOff = true;
            const long int second = GEM_run_single( uPrimalSol );
            g_optimaLineSearchRetry = false; g_optimaLSEscapeOff = false;
            itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
            native_trace_decide( "%s requested=%ld first=%ld second=%ld", lsRetryArmed ? "lsretry" : "escretry",
                                 (long)requested, (long)status, (long)second );
            node_logger->info( "GEM_run(): mode {} failed ({}) - re-solved once from the same inputs with {}: {}",
                               (long)requested, (long)status,
                               lsRetryArmed ? "the line-search retry" : "the stall escape off",
                               isOkStatus( second ) ? "recovered" : "still failed" );
            const bool ok2 = second == OK_GEM_AOP || second == OK_GEM_SOP || second == OK_GEM_HOP || second == OK_GEM_SHP
                          || second == OK_GEM_AIA || second == OK_GEM_SIA;
            if( !ok2 ) { databr_copy( failed ); ipmlog_error = failedError; }
            databr_free( failed );
            pmm->ITF = itf; pmm->ITG = itg; NumIterFIA = itf; NumIterIPM = itg;
            CNode->IterDone = itf + itg; CalcTime = seconds;
            status = CNode->NodeStatusCH;
        }
        databr_free( lsBefore );
    }
    // pa_IpmLevelStart: a call whose native part started from the leveled solution and did not end OK is solved
    // again once from the ordinary start; if that fails too, the first result is kept.
    if( lvlBefore )
    {
        if( g_ipmLevelUsed && !isOkStatus( status ) )
        {
            DATABR* failed = new DATABR;
            dbr_dch_api::databr_reset( failed, 1 );
            dbr_dch_api::databr_realloc( CSD, failed );
            { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
            const std::string failedError = ipmlog_error;
            long int itf = pmm->ITF, itg = pmm->ITG; double seconds = CalcTime;
            databr_copy( lvlBefore );
            CNode->NodeStatusCH = requested;
            clearipmLogError();
            g_ipmLevelOff = true;
            const long int second = GEM_run_single( uPrimalSol );
            g_ipmLevelOff = false;
            itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
            native_trace_decide( "natlvlretry requested=%ld first=%ld second=%ld", (long)requested, (long)status, (long)second );
            node_logger->info( "GEM_run(): mode {} failed ({}) after the leveled start - solved again from the ordinary start: {}",
                               (long)requested, (long)status, isOkStatus( second ) ? "recovered" : "still failed" );
            if( !isOkStatus( second ) ) { databr_copy( failed ); ipmlog_error = failedError; }
            databr_free( failed );
            pmm->ITF = itf; pmm->ITG = itg; NumIterFIA = itf; NumIterIPM = itg;
            CNode->IterDone = itf + itg; CalcTime = seconds;
            status = CNode->NodeStatusCH;
        }
        databr_free( lvlBefore );
    }
    if( requested == NEED_GEM_AIA && status == ERR_GEM_AIA && !kinetics && multi_ptr()->base_param()->FailRescue <= 0 )
    {
        const long int nudges = multi_ptr()->base_param()->ColdRetryNudges;
        if( nudges > 0 )
            status = GEM_run_cold_retry( nudges );
    }
    // pa_OptimaColdRetry > 0: a warm Optima call (SOP, or SHP whose Optima leg failed) that did
    // not return OK is re-solved cold as AOP from the same inputs. With pa_OptimaCgSeed > 0 a
    // failed cold call is tried once more from the column-generation seed. A retry that is not OK
    // hands the failed call back exactly as it was.
    const bool optimaColdRetry = ( multi_ptr()->base_param()->OptimaColdRetry > 0 );   // pa_OptimaColdRetry
    if( optimaColdRetry && !kinetics && ( requested == NEED_GEM_SOP || requested == NEED_GEM_SHP )
        && status != OK_GEM_SOP && status != OK_GEM_SHP && status != OK_GEM_AOP
        && status != OK_GEM_SIA && status != OK_GEM_AIA )
        status = GEM_run_optima_cold_retry( requested );
    // pa_OptimaCgSeed is used only as a second attempt of a failed cold AOP call (see
    // g_optimaCgSeedArmed). A retry that is not OK hands the first failure back unchanged.
    if( multi_ptr()->base_param()->OptimaCgSeed > 0. && !kinetics && requested == NEED_GEM_AOP && status != OK_GEM_AOP )
    {
        DATABR* failed = new DATABR;
        dbr_dch_api::databr_reset( failed, 1 );
        dbr_dch_api::databr_realloc( CSD, failed );
        { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
        const std::string failedError = ipmlog_error;
        long int itf = pmm->ITF, itg = pmm->ITG; double seconds = CalcTime;
        const long int second = GEM_run_aop_cgseed( itf, itg, seconds );
        if( second != OK_GEM_AOP ) { databr_copy( failed ); ipmlog_error = failedError; }
        databr_free( failed );
        pmm->ITF = itf; pmm->ITG = itg; NumIterFIA = itf; NumIterIPM = itg;
        CNode->IterDone = itf + itg; CalcTime = seconds;
        status = CNode->NodeStatusCH;
    }
    if( fastFailed )
    {
        pmm->ITF += fastF; pmm->ITG += fastG;
        NumIterFIA = pmm->ITF; NumIterIPM = pmm->ITG;
        CNode->IterDone = pmm->ITF + pmm->ITG;
        CalcTime += fastSec;
    }
#ifdef USE_OPTIMA_SOLVER
    if( isOkStatus( status ) )
        multi_ptr()->WarnCollapsedTwins();
#endif
    // The answer has an aqueous phase with no water: pH, pe and ionic strength are reported as 0.
    if( isOkStatus( status ) && g_aqueousAbsent )
    {
        native_trace_decide( "aqgone status=%ld", (long)status );
        static thread_local bool warned = false;
        if( !warned )
            node_logger->warn( "GEM_run(): the aqueous phase has no water left, so pH, pe and ionic strength are reported as 0" );
        warned = true;
    }
    return status;
}

// One GEM_run_single() with the legacy retry: with the line search on, a failed call is re-run
// once from the same inputs with the stall escape, the objective memo and reject_if_worse off
// (g_optimaLSEscapeOff); a retry that is not OK hands the first failure back. Used by the later
// retry steps (the warm->cold leg, the column-generation seed). Accumulates into the caller's
// totals.
long int TNode::GEM_run_single_legacy_retry( bool uPrimalSol, long int& itf, long int& itg, double& seconds, const char* step )
{
    const long int requested = CNode->NodeStatusCH;
    auto isOk = []( long int st ) { return st == OK_GEM_AOP || st == OK_GEM_SOP || st == OK_GEM_HOP || st == OK_GEM_SHP
                                           || st == OK_GEM_AIA || st == OK_GEM_SIA; };
    const bool armed = multi_ptr()->base_param()->OptimaLineSearch > 0. && !g_optimaLSEscapeOff;
    DATABR* before = nullptr;
    if( armed )
    {
        before = new DATABR;
        dbr_dch_api::databr_reset( before, 1 );
        dbr_dch_api::databr_realloc( CSD, before );
        { DATABR* live = CNode; CNode = before; databr_copy( live ); CNode = live; }
    }
    long int st = GEM_run_single( uPrimalSol );
    itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
    if( armed && !isOk( st ) )
    {
        DATABR* failed = new DATABR;
        dbr_dch_api::databr_reset( failed, 1 );
        dbr_dch_api::databr_realloc( CSD, failed );
        { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
        const std::string failedError = ipmlog_error;
        databr_copy( before );
        CNode->NodeStatusCH = requested;
        g_optimaLSEscapeOff = true;
        const long int second = GEM_run_single( uPrimalSol );
        g_optimaLSEscapeOff = false;
        itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
        native_trace_decide( "escretry step=%s requested=%ld first=%ld second=%ld", step, (long)requested, (long)st, (long)second );
        node_logger->info( "GEM_run(): {} attempt failed ({}) - re-solved with the stall escape off: {}",
                           step, (long)st, isOk( second ) ? "recovered" : "still failed" );
        if( !isOk( second ) ) { databr_copy( failed ); ipmlog_error = failedError; }
        databr_free( failed );
        st = CNode->NodeStatusCH;
    }
    if( before ) databr_free( before );
    return st;
}

// Cold AOP with the column-generation seed armed (one call); accumulates iterations/time into the caller's totals.
long int TNode::GEM_run_aop_cgseed( long int& itf, long int& itg, double& seconds )
{
    CNode->NodeStatusCH = NEED_GEM_AOP;
    g_optimaCgSeedArmed = true;
    const long int st = GEM_run_single_legacy_retry( false, itf, itg, seconds, "cgseed" );
    g_optimaCgSeedArmed = false;
    native_trace_decide( "cgseedretry aop=%ld", (long)st );
    node_logger->info( "GEM_run(): cold AOP failed - retried from the column-generation seed: {}",
                       st == OK_GEM_AOP ? "recovered" : "still failed" );
    return st;
}

// Cold AOP re-solve of a failed warm Optima call (pa_OptimaColdRetry). The inputs (bIC, T, P,
// dll/dul) are still in CNode, since packDataBr() writes only outputs. An OK_GEM_AOP result is
// reported under the cold code (as SIA's switch to OK_GEM_AIA); anything else restores the
// failed call's DATABR, status and error message.
long int TNode::GEM_run_optima_cold_retry( long int requested )
{
    DATABR* failed = new DATABR;
    dbr_dch_api::databr_reset( failed, 1 );
    dbr_dch_api::databr_realloc( CSD, failed );
    { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
    const std::string failedError = ipmlog_error;
    const long int failedStatus = failed->NodeStatusCH;
    long int itf = pmm->ITF, itg = pmm->ITG;
    double seconds = CalcTime;

    CNode->NodeStatusCH = NEED_GEM_AOP;
    long int cold = GEM_run_single_legacy_retry( false, itf, itg, seconds, "coldleg" );
    if( cold != OK_GEM_AOP && multi_ptr()->base_param()->OptimaCgSeed > 0. )   // the column-generation seed as a second attempt here too
        cold = GEM_run_aop_cgseed( itf, itg, seconds );
    native_trace_decide( "optimacoldretry requested=%ld failed=%ld aop=%ld", (long)requested,
                         (long)failedStatus, (long)cold );
    node_logger->info( "GEM_run(): warm Optima call (mode {}) failed ({}) - re-solved cold as AOP: {}",
                       (long)requested, (long)failedStatus, cold == OK_GEM_AOP ? "recovered" : "still failed" );
    if( cold != OK_GEM_AOP )
    {
        databr_copy( failed );
        ipmlog_error = failedError;
    }
    databr_free( failed );
    pmm->ITF = itf;
    pmm->ITG = itg;
    NumIterFIA = itf;
    NumIterIPM = itg;
    CNode->IterDone = itf + itg;
    CalcTime = seconds;
    return CNode->NodeStatusCH;
}

// Whether the amounts of a saved answer balance the bulk composition b (one value per IC of the DATABR) under the
// per-element test of the solve itself (relative to DHBM, or also absolute with pa_DT); the charge row is left out.
// False when not every species or element is in the DATABR: the balance cannot be computed, so the answer is not kept.
bool TNode::NudgedAnswerBalances( const DATABR* d, const std::vector<double>& b, double& worstRel ) const
{
    worstRel = 0.;
    if( CSD->nDCb != CSD->nDC || CSD->nICb != CSD->nIC || (long int)b.size() < CSD->nICb )
        return false;
    const double dtExp = fabs( (double)multi_ptr()->base_param()->DT );
    const double absCut = ( dtExp < 2. ) ? pmm->DHBM : pow( 10., -dtExp );
    bool ok = true;
    for( long int i = 0; i < CSD->nICb; i++ )
    {
        const long int ic = CSD->xic[i];
        if( CSD->ccIC[ic] == IC_CHARGE )
            continue;
        double c = b[(size_t)i];
        for( long int jj = 0; jj < CSD->nDCb; jj++ )
            c -= CSD->A[ CSD->xdc[jj] * CSD->nIC + ic ] * d->xDC[jj];
        const double a = fabs( c ), bar = b[(size_t)i] * pmm->DHBM;
        const double r = bar > 0. ? a / bar : ( a > 0. ? 1e300 : 0. );
        worstRel = std::max( worstRel, r );
        if( !multi_ptr()->base_param()->DT ? ( r > 1. ) : ( r > 1. && a > absCut ) )
            ok = false;
    }
    return ok;
}

// Recovery of a failed cold native call (pa_ColdRetryNudges): a converged solve one 1e-15 nudge
// away is a warm start from which SIA reaches the equilibrium at the exact composition. Every
// nudge is undone before the warm solve. If that warm solve fails, the nudged answer is returned
// instead. A call that no nudge recovers is handed back exactly as it failed.
long int TNode::GEM_run_cold_retry( long int maxNudges )
{
    const long int nIC = CSD->nICb;
    DATABR* failed = new DATABR;
    dbr_dch_api::databr_reset( failed, 1 );
    dbr_dch_api::databr_realloc( CSD, failed );
    { DATABR* live = CNode; CNode = failed; databr_copy( live ); CNode = live; }
    const std::string failedError = ipmlog_error;
    const std::vector<double> b0( CNode->bIC, CNode->bIC + nIC );
    long int itf = pmm->ITF, itg = pmm->ITG;
    double seconds = CalcTime;

    native_trace_quiet( true );
    long int result = ERR_GEM_AIA;
    // The first nudged cold answer, kept in case the warm re-solve at the exact composition fails.
    DATABR* nudgedOk = nullptr;
    double nudgedK = 0.;
    for( long int n = 0; n < maxNudges && result == ERR_GEM_AIA; n++ )
    {
        const double k = ( n % 2 == 0 ? 1. : -1. ) * (double)( n / 2 + 1 );
        for( long int i = 0; i < nIC; i++ )
            CNode->bIC[i] = b0[(size_t)i] * ( 1. + ( ( i % 2 ) ? -1. : 1. ) * k * 1e-15 );
        CNode->NodeStatusCH = NEED_GEM_AIA;
        const long int cold = GEM_run_single( false );
        itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
        long int warm = -1;
        if( cold == OK_GEM_AIA )
        {
            if( !nudgedOk )
            {
                nudgedOk = new DATABR;
                dbr_dch_api::databr_reset( nudgedOk, 1 );
                dbr_dch_api::databr_realloc( CSD, nudgedOk );
                { DATABR* live = CNode; CNode = nudgedOk; databr_copy( live ); CNode = live; }
                nudgedK = k;
            }
            for( long int i = 0; i < nIC; i++ )
                CNode->bIC[i] = b0[(size_t)i];
            CNode->NodeStatusCH = NEED_GEM_SIA;
            warm = GEM_run_single( false );
            itf += pmm->ITF; itg += pmm->ITG; seconds += CalcTime;
            if( warm == OK_GEM_SIA || warm == OK_GEM_AIA )
                result = OK_GEM_AIA;
            else if( warm == BAD_GEM_SIA || warm == BAD_GEM_AIA )
                result = BAD_GEM_AIA;
            else if( nudgedOk )
            {   // the nudged answer is off the exact composition by only |k|*1e-15 of each amount
                double worstRel = 0.;
                if( NudgedAnswerBalances( nudgedOk, b0, worstRel ) )
                {
                    databr_copy( nudgedOk );
                    for( long int i = 0; i < nIC; i++ )
                        CNode->bIC[i] = b0[(size_t)i];
                    result = OK_GEM_AIA;
                    native_trace_decide( "coldretry-keepnudged nudge=%+g mb_rel=%.3e", nudgedK, worstRel );
                    node_logger->info( "GEM_run(): the warm check at the exact composition failed; kept the answer "
                                       "solved at a composition nudged by {:g}e-15 of each amount", nudgedK );
                }
                else
                {
                    native_trace_decide( "coldretry-keepnudged-refused nudge=%+g mb_rel=%.3e", nudgedK, worstRel );
                    node_logger->info( "GEM_run(): the answer solved at a nudged composition does not balance the exact "
                                       "composition, so it is not kept" );
                }
            }
        }
        native_trace_decide( "coldretry attempt=%ld nudge=%+g aia=%ld sia=%ld itf=%ld itg=%ld",
                             (long)( n + 1 ), k, (long)cold, (long)warm, (long)itf, (long)itg );
    }
    native_trace_quiet( false );
    if( nudgedOk )
        databr_free( nudgedOk );
    node_logger->info( "GEM_run(): cold native call failed - retried from nudged compositions (up to {}): {}",
                       (long)maxNudges, result == ERR_GEM_AIA ? "still failed" : "recovered" );

    if( result == ERR_GEM_AIA )
    {
        databr_copy( failed );            // the failed call as it was, status ERR_GEM_AIA
        ipmlog_error = failedError;
    }
    else
    {
        CNode->NodeStatusCH = result;
        native_trace_run_result( *pmm, NEED_GEM_AIA, result, multi_base );
    }
    databr_free( failed );
    pmm->ITF = itf;
    pmm->ITG = itg;
    NumIterFIA = itf;
    NumIterIPM = itg;
    CNode->IterDone = itf + itg;
    CalcTime = seconds;
    return CNode->NodeStatusCH;
}

long int TNode::GEM_run_single( bool uPrimalSol )
{
    CalcTime = 0.0;
    PrecLoops = 0; NumIterFIA = 0; NumIterIPM = 0;
    clearipmLogError();

    try
    {
        ipmlog_file->debug(" GEM_run() begin Mode= {}", CNode->NodeStatusCH);

#ifndef NDEBUG
        // Bulk composition as received, before any internal processing (debug level).
        if (node_logger->should_log(spdlog::level::debug)) {
            for (long int i = 0; i < CSD->nICb; i++)
                node_logger->debug("bIC[{}] {} = {:.6e}", i, CSD->ICNL[CSD->xic[i]], CNode->bIC[i]);
        }
        // Warn about an IC at or below the numerical floor.
        for (long int i = 0; i < CSD->nICb; i++) {
            if (CNode->bIC[i] > 0. && CNode->bIC[i] <= pmm->DcMinM)
                node_logger->warn("bIC[{}] {} = {:.6e} is at/below the numerical floor (DcMinM={:.3e})",
                                   i, CSD->ICNL[CSD->xic[i]], CNode->bIC[i], pmm->DcMinM);
        }
#endif

        // Checking T and P  for interpolation intervals
        check_TP( CNode->TK, CNode->P);
        // Unpacking work DATABR structure into MULTI (GEM IPM structure): uses DATACH
        // setting up up PIA or AIA mode
        if( CNode->NodeStatusCH == NEED_GEM_SIA || CNode->NodeStatusCH == NEED_GEM_SOP )
        {
            pmm->pNP = 1;
            unpackDataBr( uPrimalSol );
        }
        else if( CNode->NodeStatusCH == NEED_GEM_AIA || CNode->NodeStatusCH == NEED_GEM_AOP
                 || CNode->NodeStatusCH == NEED_GEM_ROP || CNode->NodeStatusCH == NEED_GEM_HOP
                 || CNode->NodeStatusCH == NEED_GEM_SHP )
        {
            pmm->pNP = 0; // As default setting AIA/AOP/ROP/HOP/SHP mode - ROP is a single mode
                          // (no warm pair), HOP's first leg is the cold native solve, and
                          // SHP's warm native leg sets pmm->pNP = 1 itself inside
                          // CalculateEquilibriumStateHOP().
            if (CNode->dt > 0.)
                uPrimalSol = true;
            unpackDataBr( uPrimalSol );
        }
        else
            return CNode->NodeStatusCH;

        // Complete run configuration into the trace, once per call and for every mode: after
        // unpackDataBr() (pmm->B[] holds the caller's bulk vector) and before the dispatch
        // (NodeStatusCH is still the request). No-op unless GEMS3K_NATIVE_TRACE_FILE is set.
        native_trace_run_header( *pmm, multi_ptr()->base_param(), CNode->NodeStatusCH );
        // The mode as requested, kept for the result record below (the dispatch and
        // packDataBr() overwrite NodeStatusCH with the outcome).
        const long int traceRequestedMode = CNode->NodeStatusCH;

        // added 18.12.14 DK : setting chemical kinetics time counter and variables
        node_logger->debug("GEM_run dTime:{}, TimeStep: {} Time: {}", CNode->dt, CNode->NodeStatusFMT, CNode->Tm);
        if( CNode->dt <= 0. )
        {  // no kinetics to consider
            pmm->kTau = 0.;
            pmm->kdT = 0.;
            pmm->ITau = -1;
            pmm->pKMM = 2;  // no need to allocate TKinMet instances
        }
        else {   // considering kinetics
            pmm->kTau = CNode->Tm;
            pmm->kdT = CNode->dt;
            if( pmm->ITau < 0 && CNode->Tm/CNode->dt < 1e-9 )
            {   // we need to initialize TKinMet
                pmm->pKMM = -1;
                pmm->ITau = -1;
            }
            else  // TKinMet exists, simulation continues
                pmm->pKMM = 1; // pmm->ITau = CNode->Tm/CNode->dt;
        }

        // Invalidate the previous call's per-leg split before the dispatch, so only a mode that
        // actually ran two legs reports one.
        multi_ptr()->hop_split = TMultiBase::HopLegSplit();

        // Optima modes dispatch to ipm_optima.cpp. Without USE_OPTIMA_SOLVER they fall back to the
        // equivalent native AIA/SIA solve rather than failing.
        if( CNode->NodeStatusCH == NEED_GEM_HOP || CNode->NodeStatusCH == NEED_GEM_SHP )
        {
#ifdef USE_OPTIMA_SOLVER
            // Hybrid: native selects the phases, then a warm Optima leg finishes (restoring
            // native's answer as BAD_GEM_HOP if it fails) - TMultiBase::CalculateEquilibriumStateHOP().
            // SHP is HOP with a warm native leg (with a cold fallback inside).
            // A native leg stopped by pa_IpmNoProgress hands its last state to the Optima leg.
            NoProgArm arm( CNode->dt <= 0. && multi_ptr()->GetControlConditions().empty() );
            CalcTime = multi_ptr()->CalculateEquilibriumStateHOP( NumIterFIA, NumIterIPM,
                                       CNode->NodeStatusCH == NEED_GEM_SHP );
#else
            node_logger->warn("GEM_run(): NEED_GEM_HOP/SHP requested but GEMS3K was not built with "
                               "USE_OPTIMA_SOLVER - falling back to the equivalent native {} solve",
                               CNode->NodeStatusCH == NEED_GEM_SHP ? "SIA" : "AIA");
            CNode->NodeStatusCH = ( CNode->NodeStatusCH == NEED_GEM_SHP ) ? NEED_GEM_SIA : NEED_GEM_AIA;
            pmm->pNP = ( CNode->NodeStatusCH == NEED_GEM_SIA ) ? 1 : 0;
            CalcTime = multi_ptr()->CalculateEquilibriumState( NumIterFIA, NumIterIPM );
#endif
        }
        else if( CNode->NodeStatusCH == NEED_GEM_AOP || CNode->NodeStatusCH == NEED_GEM_SOP
            || CNode->NodeStatusCH == NEED_GEM_ROP )
        {
#ifdef USE_OPTIMA_SOLVER
            CalcTime = multi_ptr()->CalculateEquilibriumStateOptima( NumIterFIA, NumIterIPM,
                                       CNode->NodeStatusCH == NEED_GEM_ROP );
#else
            node_logger->warn("GEM_run(): NEED_GEM_AOP/SOP/ROP requested but GEMS3K was not built with "
                               "USE_OPTIMA_SOLVER - falling back to the equivalent native {} solve",
                               CNode->NodeStatusCH == NEED_GEM_SOP ? "SIA" : "AIA");
            CNode->NodeStatusCH = ( CNode->NodeStatusCH == NEED_GEM_SOP ) ? NEED_GEM_SIA : NEED_GEM_AIA;
            CalcTime = multi_ptr()->CalculateEquilibriumState( NumIterFIA, NumIterIPM );
#endif
        }
        else
        {
#ifdef USE_OPTIMA_SOLVER
            NoProgArm arm( CNode->dt <= 0. && multi_ptr()->GetControlConditions().empty()
                           && multi_ptr()->base_param()->NativeOptimaFinish > 0 );
            try
            {
                // an error here may be recovered below: it is logged there only if it is not
                struct HandledErrors { HandledErrors() { g_errorsHandledByCaller++; } ~HandledErrors() { g_errorsHandledByCaller--; } } handled;
                CalcTime = multi_ptr()->CalculateEquilibriumState( /*RefineLoops,*/ NumIterFIA, NumIterIPM  );
            }
            catch( TError& nerr )
            {
                // Native stopped at its iteration limit, or its linear system broke down: let Optima find the
                // answer; otherwise report the failure.
                double optSec = 0.;
                const bool limit = nerr.title.rfind( "E06IPM", 0 ) == 0;
                const bool broken = nerr.title.rfind( "E07IPM", 0 ) == 0;
                if( !( ( limit || broken ) && CNode->dt <= 0.
                       && ( CNode->NodeStatusCH == NEED_GEM_AIA || CNode->NodeStatusCH == NEED_GEM_SIA )
                       && multi_ptr()->NativeOptimaRescue( NumIterFIA, NumIterIPM, optSec, broken ) ) )
                {
                    // pa_IpmNoProgress stopped native and the finish failed: solve again cold without the stop.
                    if( !( limit && g_ipmNoProgFired ) )
                    {
                        LogError( nerr.title, nerr.mess );
                        throw;
                    }
                    g_ipmNoProgArmed = false;
                    pmm->pNP = 0;
                    const long int fia0 = NumIterFIA, ipm0 = NumIterIPM;
                    try
                    {
                        CalcTime = multi_ptr()->CalculateEquilibriumState( NumIterFIA, NumIterIPM );
                        NumIterFIA += fia0; NumIterIPM += ipm0;
                        native_trace_decide( "ipmnoprog-retry ok=1" );
                    }
                    catch( TError& )
                    {
                        native_trace_decide( "ipmnoprog-retry ok=0" );
                        LogError( nerr.title, nerr.mess );
                        throw nerr;
                    }
                }
                else
                    CalcTime = multi_ptr()->GetPM()->t_elap_sec + optSec;
            }
#else
            CalcTime = multi_ptr()->CalculateEquilibriumState( /*RefineLoops,*/ NumIterFIA, NumIterIPM  );
#endif
#ifdef USE_OPTIMA_SOLVER
            // A trace solid phase that still tends to grow: let Optima finish, keep only a better answer.
            if( ( CNode->NodeStatusCH == NEED_GEM_AIA || CNode->NodeStatusCH == NEED_GEM_SIA ) && CNode->dt <= 0. )
                CalcTime += multi_ptr()->NativeOptimaFinish( NumIterFIA, NumIterIPM );
#endif
        }

        // Extracting and packing GEM IPM results into work DATABR structure
        packDataBr();
        CNode->IterDone = NumIterFIA+NumIterIPM;
        //**************************************************************
        // only for testing output results for files
        //    GEM_write_dbr( "calculated_dbr.dat",  2 /*json*/ );
        //    GEM_print_ipm( "calc_multi.ipm" );
        // *********************************************************

        // test error result GEM IPM calculation of equilibrium state in MULTI
        long int erCode =  multi_ptr()->testMulti();

        if( erCode )
        {
            if( CNode->NodeStatusCH  == NEED_GEM_AIA )
                CNode->NodeStatusCH = BAD_GEM_AIA;
            else if( CNode->NodeStatusCH == NEED_GEM_AOP )
                CNode->NodeStatusCH = BAD_GEM_AOP;
            else if( CNode->NodeStatusCH == NEED_GEM_SOP )
                CNode->NodeStatusCH = BAD_GEM_SOP;
            else if( CNode->NodeStatusCH == NEED_GEM_ROP )
                CNode->NodeStatusCH = BAD_GEM_ROP;
            else if( CNode->NodeStatusCH == NEED_GEM_HOP )
                CNode->NodeStatusCH = BAD_GEM_HOP;
            else if( CNode->NodeStatusCH == NEED_GEM_SHP )
                CNode->NodeStatusCH = BAD_GEM_SHP;
            else
                CNode->NodeStatusCH = BAD_GEM_SIA;

            // internal multi error
            ipmlog_error = pmm->errorCode+ std::string(": ") +  pmm->errorBuf;
            // A thrown error was already logged by Error(); this status path throws nothing.
            // A failed quick try is followed by the full solve, so it is logged at debug, as in Error().
            if( g_errorsHandledByCaller > 0 )
                node_logger->debug( "{}", ipmlog_error );
            else
                node_logger->error( "{}", ipmlog_error );
        }
        else
        {
            if( CNode->NodeStatusCH  == NEED_GEM_AIA )
                CNode->NodeStatusCH = OK_GEM_AIA;
            else if( CNode->NodeStatusCH == NEED_GEM_AOP )
                CNode->NodeStatusCH = OK_GEM_AOP;
            else if( CNode->NodeStatusCH == NEED_GEM_SOP )
                CNode->NodeStatusCH = OK_GEM_SOP;
            else if( CNode->NodeStatusCH == NEED_GEM_ROP )
                CNode->NodeStatusCH = OK_GEM_ROP;
            else if( CNode->NodeStatusCH == NEED_GEM_HOP )
                CNode->NodeStatusCH = OK_GEM_HOP;
            else if( CNode->NodeStatusCH == NEED_GEM_SHP )
                CNode->NodeStatusCH = OK_GEM_SHP;
            else
                CNode->NodeStatusCH = OK_GEM_SIA;
        }

        // The regime this call reached, into the trace (after the status has settled). Not
        // emitted from the catch blocks: a thrown call's reported values are stale. No-op unless
        // GEMS3K_NATIVE_TRACE_FILE is set.
        native_trace_run_result( *pmm, traceRequestedMode, CNode->NodeStatusCH, multi_base );

        return CNode->NodeStatusCH;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
        if( CNode->NodeStatusCH  == NEED_GEM_AIA )
            CNode->NodeStatusCH = ERR_GEM_AIA;
        else if( CNode->NodeStatusCH == NEED_GEM_AOP )
            CNode->NodeStatusCH = ERR_GEM_AOP;
        else if( CNode->NodeStatusCH == NEED_GEM_SOP )
            CNode->NodeStatusCH = ERR_GEM_SOP;
        else if( CNode->NodeStatusCH == NEED_GEM_ROP )
            CNode->NodeStatusCH = ERR_GEM_ROP;
        else if( CNode->NodeStatusCH == NEED_GEM_HOP )
            CNode->NodeStatusCH = ERR_GEM_HOP;
        else if( CNode->NodeStatusCH == NEED_GEM_SHP )
            CNode->NodeStatusCH = ERR_GEM_SHP;
        else
            CNode->NodeStatusCH = ERR_GEM_SIA;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
        CNode->NodeStatusCH = T_ERROR_GEM;
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        CNode->NodeStatusCH = T_ERROR_GEM;
    }
    ipmlog_file->error("Error Node:{}  time:{}  dt:{}", CNode->NodeHandle, CNode->Tm, CNode->dt);
    if( multi_ptr()->base_param()->PSM >= 2  )
        ipmlog_file->error("{}", ipmlog_error);

    return CNode->NodeStatusCH;
}

// Returns GEMIPM2 calculation time in seconds elapsed during the last call of GEM_run() -
// can be used for monitoring the performance of calculations.
// Return value:  double number, may contain 0.0 if the calculation time is less than the
//                internal time resolution of C/C++ function
double TNode::GEM_CalcTime( long int& NumK2, long int& NumIterFIA1, long int& NumIterIPM1)
{
    NumK2 = pmm->K2;
    NumIterFIA1 = pmm->ITF;
    NumIterIPM1 = pmm->ITG;
    return CalcTime;
}

// To obtain the number of GEM IPM2 iterations performed during the last call of GEM_run() e.g. for monitoring the
// performance of GEMS3K in AIA or SIA modes, or for problem diagnostics.
// Parameters:  long int variables per reference (must be allocated before calling GEM_Iterations(), previous values will be lost. See Return values.
// Return values:
//   Function         Total number of EFD + IPM iterations from the last call to GEM_run()
//   PrecLoops        Number of performed IPM-2 precision refinement loops
//   NumIterFIA       Total number of performed MBR() iterations to obtain a feasible initial approximation for the IPM algorithm.
//   NumIterIPM       Total number of performed IPM main descent algorithm iterations.
long int TNode::GEM_Iterations( long int& PrecLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
    PrecLoops_ = PrecLoops;
    NumIterFIA_ = NumIterFIA;
    NumIterIPM_ = NumIterIPM;
    return NumIterFIA+NumIterIPM;
}

// The last two-leg (HOP/SHP) call's cost, native leg and Optima leg apart - see the
// declaration in node.h and TMultiBase::HopLegSplit.
bool TNode::GEM_IterationsHOP( long int& NumIterFIANative, long int& NumIterIPMNative,
                               long int& NumIterFIAOptima, long int& NumIterIPMOptima,
                               double& TimeNative, double& TimeOptima,
                               bool* Failed ) const
{
    const TMultiBase::HopLegSplit& s = multi_ptr()->hop_split;
    if( !s.valid )
        return false;               // the last solve was not two-leg - leave the outputs alone
    NumIterFIANative = s.fiaNative;
    NumIterIPMNative = s.ipmNative;
    NumIterFIAOptima = s.fiaOptima;
    NumIterIPMOptima = s.ipmOptima;
    TimeNative = s.timeNative;
    TimeOptima = s.timeOptima;
    if( Failed ) *Failed = s.failed;
    return true;
}

// Extracting and packing GEM IPM results into work DATABR structure
void TNode::packDataBr()
{
    long int ii;

    // set default data to DataBr
    //   CNode->NodeStatusCH = NEED_GEM_AIA;
    // pmm->pNP (0 = cold, 1 = warm) selects AIA/SIA vs. AOP/SOP; which pair is read from the
    // still unmodified NodeStatusCH.
    if( CNode->NodeStatusCH == NEED_GEM_ROP || CNode->NodeStatusCH == NEED_GEM_HOP
        || CNode->NodeStatusCH == NEED_GEM_SHP )
    {
        ; // ROP is a single mode, and HOP/SHP are told apart by the caller, not by pmm->pNP
          // (both legs set it themselves). They must be listed here, or the generic branch
          // below would rewrite the node to NEED_GEM_SIA.
    }
    else if( CNode->NodeStatusCH == NEED_GEM_AOP || CNode->NodeStatusCH == NEED_GEM_SOP )
    {
        if( pmm->pNP == 0 )
            CNode->NodeStatusCH = NEED_GEM_AOP;
        else
            CNode->NodeStatusCH = NEED_GEM_SOP;
    }
    else
    {
        if( pmm->pNP == 0 )
            CNode->NodeStatusCH = NEED_GEM_AIA;
        else
            CNode->NodeStatusCH = NEED_GEM_SIA;
    }

    CNode->TK = pmm->TCc+C_to_K; //25
    CNode->P = pmm->Pc*bar_to_Pa; //1
    //   CNode->IterDone = pmm->IT;
    CNode->IterDone = /*pmm->ITF+*/pmm->IT;   // Now complete number of FIA and IPM iterations
    // values
    CNode->Vs = pmm->VXc*1.e-6; // from cm3 to m3
    // pm.FX is seeded each call with kTotalGibbsEnergyUnset and both solver paths overwrite it,
    // so this branch should be unreachable: it reports an answer whose Gibbs energy was never
    // computed instead of publishing the marker. Checked with == because the marker is assigned
    // verbatim.
    if( pmm->FX == kTotalGibbsEnergyUnset )
        node_logger->error( "packDataBr(): total Gibbs energy is still kTotalGibbsEnergyUnset - the "
                           "solver path that produced NodeStatusCH={} did not set pm.FX, so "
                           "DATABR.Gs, GEM_to_MT()'s p_Gs and any exported -dbr carry a sentinel and "
                           "not an energy. Use TNode::Get_GibbsEnergy() meanwhile; this is a bug.",
                           CNode->NodeStatusCH );
    CNode->Gs = pmm->FX;
    CNode->Hs = pmm->HXc;
    CNode->IC = pmm->IC;
    CNode->pH = pmm->pH;
    CNode->pe = pmm->pe;
    //  CNode->Eh = pmm->FitVar[3];  Bugfix 19.12.2006  KD
    CNode->Eh = pmm->Eh;
    CNode->Ms = pmm->MBX;

    // arrays
    // An absent mixed phase gets a real saturation index in place of the cold-start placeholder.
    if( multi_ptr() )
    {
        try { multi_ptr()->ReportAbsentPhaseIndexes(); }
        catch( ... ) {}
    }
    for( ii=0; ii<CSD->nPHb; ii++ )
    {  CNode->xPH[ii] = pmm->XF[ CSD->xph[ii] ];
        if( CSD->nAalp >0 )
            CNode->aPH[ii] = pmm->Aalp[ CSD->xph[ii] ]*kg_to_g;
        CNode->omPH[ii] = pmm->Falp[ CSD->xph[ii] ];
    }
    floorParkedAmount.assign( (size_t)CSD->nPHb, -1. );
    for( ii=0; ii<CSD->nPHb; ii++ )
        if( multi_ptr() && multi_ptr()->PhaseFloorParked( CSD->xph[ii] ) )
            floorParkedAmount[(size_t)ii] = CNode->xPH[ii];
    for( ii=0; ii<CSD->nPSb; ii++ )
    {   CNode->vPS[ii] = pmm->FVOL[ CSD->xph[ii] ]/m3_to_cm3;
        CNode->mPS[ii] = pmm->FWGT[ CSD->xph[ii] ]/kg_to_g;
        CNode->xPA[ii] = pmm->XFA[ CSD->xph[ii] ];
        CNode->amru[ii] = pmm->PUL[ CSD->xph[ii] ];
        CNode->amrl[ii] = pmm->PLL[ CSD->xph[ii] ];
    }
    for( ii=0; ii<CSD->nPSb; ii++ )
        for(long int jj=0; jj<CSD->nICb; jj++ )
        { long int   new_ndx= (ii*CSD->nICb)+jj,
                    mul_ndx = ( CSD->xph[ii]*CSD->nIC )+ CSD->xic[jj];
            CNode->bPS[new_ndx] = pmm->BF[ mul_ndx ];
        }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->xDC[ii] = pmm->X[ CSD->xdc[ii] ];
        CNode->gam[ii] = pmm->Gamma[ CSD->xdc[ii] ];
        CNode->dul[ii] = pmm->DUL[ CSD->xdc[ii] ];// always for GEM2MT init
        CNode->dll[ii] = pmm->DLL[ CSD->xdc[ii] ];// always for GEM2MT init
    }
    for( ii=0; ii<CSD->nICb; ii++ )
    {  CNode->bIC[ii] = pmm->B[ CSD->xic[ii] ]; // always for GEM2MT  init
        CNode->rMB[ii] = pmm->C[ CSD->xic[ii] ];
        CNode->uIC[ii] = pmm->U[ CSD->xic[ii] ];
        CNode->bSP[ii] = pmm->BFC[ CSD->xic[ii] ];
    }
}

// Unpacking work DATABR structure into MULTI
//(GEM IPM work structure): uses DATACH
//  if uPrimalSol is true then the primal solution (vectors x, gamma, IC etc.)
//  will be unpacked - as an option for PIA mode with previous GEM solution from
//  the same node.
//  If uPrimalSol = false then the primal solution data will not be unpacked
//  into the MULTI structure (AIA mode or SIA mode with primal solution retained
//    in the MULTI structure from previous IPM calculation)
void TNode::unpackDataBr( bool uPrimalSol )
{
    long int ii;

    CheckMtparam(); // T or P change detection - moved to here from InitalizeGEM_IPM_Data() 11.10.2012
    pmm->kTau = CNode->Tm;  // added 18.12.14 DK
    pmm->kdT = CNode->dt;   // added 18.12.14 DK

    pmm->TCc = CNode->TK-C_to_K;
    pmm->Tc = CNode->TK;
    pmm->Pc  = CNode->P/bar_to_Pa;
    pmm->VXc = CNode->Vs/1.e-6; // from cm3 to m3
    // Obligatory arrays - always unpacked!
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        pmm->DUL[ CSD->xdc[ii] ] = CNode->dul[ii];
        pmm->DLL[ CSD->xdc[ii] ] = CNode->dll[ii];
        if( pmm->DUL[ CSD->xdc[ii] ] < pmm->DLL[ CSD->xdc[ii] ] )
        {
            Error("unpackDataBr", std::string("Upper kinetic restriction less than the lower one for DC&RC")
                  +char_array_to_string( pmm->SM[CSD->xdc[ii]], MAXDCNAME));
        }
    }
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        pmm->B[ CSD->xic[ii] ] = CNode->bIC[ii];
        if( ii < CSD->nICb-1 && pmm->B[ CSD->xic[ii] ] < multi_ptr()->base_param()->DB )
        {
            Error("unpackDataBr", std::string("Bulk mole amount of IC ")+
                  char_array_to_string(pmm->SB[CSD->xic[ii]], 6)+" is "+
                    std::to_string(pmm->B[ CSD->xic[ii] ])+" - out of range" );
        }
    }
    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        if( CSD->nAalp > 0 )
            pmm->Aalp[ CSD->xph[ii] ] = CNode->aPH[ii]/kg_to_g;
        pmm->Falp[ CSD->xph[ii] ] = CNode->omPH[ii];
    }

    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        pmm->XF[ CSD->xph[ii] ] = CNode->xPH[ii];
        //pmm->YF[ CSD->xph[ii] ] = CNode->xPH[ii];
    }

    if( !uPrimalSol )
    {
        //  Using primal solution retained in the MULTI structure instead -
        ; // the primal solution data from the DATABR structure are not unpacked

        //   pmm->IT = 0;
    }
    else {   // Unpacking primal solution provided in the node DATABR structure
        pmm->IT = CNode->IterDone; // ?  pmm->ITF+pmm->IT;
        //pmm->IT = 0;
        pmm->MBX = CNode->Ms;
        pmm->IC = CNode->IC;
        pmm->Eh = CNode->Eh;
        for( ii=0; ii<CSD->nDCb; ii++ )
            /*    pmm->X[ CSD->xdc[ii] ] = */
            pmm->Y[ CSD->xdc[ii] ] = CNode->xDC[ii];

        for( ii=0; ii<CSD->nPSb; ii++ )
        {
            //      pmm->FVOL[ CSD->xph[ii] ] = CNode->vPS[ii]*m3_to_cm3;
            //      pmm->FWGT[ CSD->xph[ii] ] = CNode->mPS[ii]*kg_to_g;
            pmm->PUL[ CSD->xph[ii] ] = CNode->amru[ii];
            pmm->PLL[ CSD->xph[ii] ] = CNode->amrl[ii];
        }

        for( ii=0; ii<CSD->nPHb; ii++ )
        {
            pmm->XF[ CSD->xph[ii] ] =
                    pmm->YF[ CSD->xph[ii] ] = this->Ph_Mole(ii);
            pmm->FVOL[ CSD->xph[ii] ] = this->Ph_Volume(ii)*m3_to_cm3;
            pmm->FWGT[ CSD->xph[ii] ] = this->Ph_Mass(ii)*kg_to_g;
        }

        for( long int k=0; k<CSD->nPSb; k++ )
            for(long int i=0; i<CSD->nICb; i++ )
            { long int dbr_ndx= (k*CSD->nICb)+i,
                        mul_ndx = ( CSD->xph[k]*CSD->nIC )+ CSD->xic[i];
                pmm->BF[ mul_ndx ] = CNode->bPS[dbr_ndx];
            }

        for( ii=0; ii<CSD->nPSb; ii++ )
            pmm->XFA[ CSD->xph[ii] ] = pmm->YFA[ CSD->xph[ii] ] = CNode->xPA[ii];

        for( ii=0; ii<CSD->nICb; ii++ )
            pmm->C[ CSD->xic[ii] ] = CNode->rMB[ii];
        for( ii=0; ii<CSD->nICb; ii++ )
            pmm->U[ CSD->xic[ii] ] = CNode->uIC[ii];

        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            pmm->Gamma[ CSD->xdc[ii] ] = CNode->gam[ii];
        }

        long int jb, je = 0;
        for( long int k=0; k<pmm->FIs; k++ )
        { // loop on solution phases
            jb = je;
            je += pmm->L1[ k ];
            // Load activity coeffs for phases-solutions
            for( ii=jb; ii<je; ii++ )
            {
                pmm->lnGam[ii] =  multi_ptr()->PhaseSpecificGamma( ii, jb, je, k, 1L );
            } // ii
        }
    }
    //  End
}


// (9) Optional, for passing the current mass transport time and time step into the work
// DATABR structure (for using it in TKinMet, or tracing/debugging, or in writing DBR files for nodes)
// This call should be used instead of obsolete GEM_set_MT() (provided below for compatibility with older codes)
//
void TNode::GEM_from_MT_time(
        double p_Tm,      // actual total simulation time, s          +       -      -
        double p_dt       // actual time step, s                      +       -      -
        )
{
    CNode->Tm = p_Tm;
    CNode->dt = p_dt;
}

void TNode::GEM_set_MT( // misleading name of the method - use GEM_from_MT_time() instead, see above
                        double p_Tm,      // actual total simulation time, s          +       -      -
                        double p_dt       // actual time step, s                      +       -      -
                        )
{
    CNode->Tm = p_Tm;
    CNode->dt = p_dt;
}

// (7a) Optional, to check if the time step in the work DATABR structure was o.k. for TKinMet calculations,
//  compared with the time step p_dt given before the GEM calculation. Checks the criteria for the validity
//  of time step. If time step was acceptable by a TKinMet model used, returns the actual time step after
//  copying (changed) AMRs into p_dul and p_dll vectors, as well as (changed) specific surface areas of
//  some (kinetically controlled) phases. Otherwise, returns a (smaller) suggested time step, while the
//  p_dul, p_pll, and p_asPH vectors remain unchanged.
//  Returns 0 or a negative number (unchanged p_dul and p_dll), if TKinMet calculations failed.
//
double GEM_to_MT_time(
        double p_dt,       ///< Actual time step, s                                     -       -     (+)   (+)
        double *p_dul,    ///< Upper AMR restrictions to amounts of DC [nDCb]          -       -      +     -
        double *p_dll     ///< Lower AMR restrictions to amounts of DC [nDCb]          -       -      +     -
        );

// (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
//  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
//   chemical data for all nodes used in the mass-transport model.
void TNode::GEM_restore_MT(
        long int  &p_NodeHandle,   // Node identification handle
        long int  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
        //                                    GEM input output  FMT control
        double &p_TK,      // Temperature T, Kelvin                       +       -      -
        double &p_P,      // Pressure P,  Pa                              +       -      -
        double &p_Vs,     // Volume V of reactive subsystem,  m3         (+)      -      +
        double &p_Ms,     // Mass of reactive subsystem, kg               -       -      +
        double *p_bIC,    // Bulk mole amounts of IC  [nICb]              +       -      -
        double *p_dul,    // Upper restrictions to amounts of DC [nDCb]   +       -      -
        double *p_dll,    // Lower restrictions to amounts of DC [nDCb]   +       -      -
        double *p_asPH    // Specific surface areas of phases,m2/kg[nPHb] +       -      -
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_TK = CNode->TK;
    p_P = CNode->P;
    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
        p_bIC[ii] = CNode->bIC[ii];
    for( ii=0; ii<CSD->nDCb; ii++ )
    {  p_dul[ii] = CNode->dul[ii];
        p_dll[ii] = CNode->dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            p_asPH[ii] = CNode->aPH[ii];
}

// (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
//  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
//   chemical data for all nodes used in the mass-transport model.
void TNode::GEM_restore_MT(
        long int  &p_NodeHandle,   // Node identification handle
        long int  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
        //                                    GEM input output  FMT control
        double &p_TK,      // Temperature T, Kelvin                       +       -      -
        double &p_P,      // Pressure P,  Pa                              +       -      -
        double &p_Vs,     // Volume V of reactive subsystem,  m3         (+)      -      +
        double &p_Ms,     // Mass of reactive subsystem, kg               -       -      +
        double *p_bIC,    // Bulk mole amounts of IC  [nICb]              +       -      -
        double *p_dul,    // Upper restrictions to amounts of DC [nDCb]   +       -      -
        double *p_dll,    // Lower restrictions to amounts of DC [nDCb]   +       -      -
        double *p_asPH,   // Specific surface areas of phases,m2/kg[nPHb] +       -      -
        double *p_omPH,  // Stability indices of phases,log10 scale [nPHb] (+)     +      -
        double *p_amru,   // Upper AMRs to masses of sol. phases [nPSb]   +       -      -
        double *p_amrl    // Lower AMRs to masses of sol. phases [nPSb]   +       -      -
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_TK = CNode->TK;
    p_P = CNode->P;
    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
        p_bIC[ii] = CNode->bIC[ii];
    for( ii=0; ii<CSD->nDCb; ii++ )
    {  p_dul[ii] = CNode->dul[ii];
        p_dll[ii] = CNode->dll[ii];
    }
    if( CSD->nAalp >0 )
    {  for( ii=0; ii<CSD->nPHb; ii++ )
            p_asPH[ii] = CNode->aPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        p_omPH[ii] = CNode->omPH[ii];
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        p_amru[ii] = CNode->amru[ii];
        p_amrl[ii] = CNode->amrl[ii];
    }
}

// (7)  Retrieves the GEMIPM2 chemical speciation calculation results from the work DATABR structure instance
//   into memory provided by the mass transport part. Dimensions and order of elements in the arrays must correspond
//   to those in currently existing DATACH memory structure.
void TNode::GEM_to_MT(
        long int &p_NodeHandle,    // Node identification handle
        long int &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
        long int &p_IterDone,      // Number of iterations performed in the last GEM IPM calculation
        //                                                  GEM input output  FMT control
        // Chemical scalar variables
        double &p_Vs,    // Total volume V of reactive subsystem at given P,T, m3    -      -      +     +
        double &p_Ms,    // Total mass of the reactive subsystem, kg                 -      -      +     +
        double &p_Gs,    // Total Gibbs energy of the reactive subsystem, J          -      -      +     +
        double &p_Hs,    // Total enthalpy of reactive subsystem, J (reserved)       -      -      +     +
        double &p_IC,    // Effective aqueous ionic strength, molal                  -      -      +     +
        double &p_pH,    // pH of aqueous solution                                   -      -      +     +
        double &p_pe,    // pe of aqueous solution                                   -      -      +     +
        double &p_Eh,    // Eh of aqueous solution, V                                -      -      +     +
        // Dynamic data - dimensions see in DATACH.H structure
        double  *p_rMB,  // Mole balance residuals for Independent Components [nICb] -      -       +     +
        double  *p_uIC,  // Dual solution: IC chemical potentials, mol/mol [nICb]    -      -       +     +
        double  *p_xDC,  // Primal solution: DC mole amounts  [nDCb]                 -      -       +     +
        double  *p_gam,  // External activity coefficients of DC [nDCb]              -      -       +     +
        double  *p_xPH,  // Total mole amounts of all phases [nPHb]                  -      -       +     +
        double  *p_vPS,  // Total volumes of multicomponent phases, m3   [nPSb]      -      -       +     +
        double  *p_mPS,  // Total mass of multicomponent phase (carrier),kg [nPSb]   -      -       +     +
        double  *p_bPS,  // Bulk compositions of phases  [nPSb][nICb]                -      -       +     +
        double  *p_xPA,  //Amount of carrier in a multicomponent asymmetric phase[nPSb]-    -       +     +
        double  *p_aPH,  //surface area for phases, m2                           [nPHb]-    -       +     +
        double  *p_bSP   //Bulk composition of all solids, moles [nICb]                -    -       +     +
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_IterDone = CNode->IterDone;

    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    p_Gs = CNode->Gs;
    p_Hs = CNode->Hs;
    p_IC = CNode->IC;
    p_pH = CNode->pH;
    p_pe = CNode->pe;
    p_Eh = CNode->Eh;

    for( ii=0; ii<CSD->nICb; ii++ )
    {
        p_rMB[ii] = CNode->rMB[ii];
        p_uIC[ii] = CNode->uIC[ii];
        p_bSP[ii] = CNode->bSP[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        p_xDC[ii] = CNode->xDC[ii];
        p_gam[ii] = CNode->gam[ii];
    }
    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        p_xPH[ii] = CNode->xPH[ii];
        p_aPH[ii] = CNode->aPH[ii]*Ph_Mass( ii );  // correction 9.10.2013 DK
    }
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        p_vPS[ii] = CNode->vPS[ii];
        p_mPS[ii] = CNode->mPS[ii];
        p_xPA[ii] = CNode->xPA[ii];
    }
    for( ii=0; ii<CSD->nPSb*CSD->nICb; ii++ )
        p_bPS[ii] = CNode->bPS[ii];
}

// (8) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH   // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
}

//(8a) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//This overloaded variant uses the xDC speciation vector for setting the
// new bulk chemical composition to be used in the next GEM_run() calculation.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,  // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
        double *p_xDC    // Mole amounts of DCs [nDCb] - will be convoluted
        // and added to the bIC GEM input vector (if full speciation
        // and not just increments then p_bIC vector must be zeroed off -
        // it will be calculated from p_xDC and stoichiometry matrix A
        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-simplex IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is ordered, leave if simplex is ordered (KD)

    // Optional part - convolution of xDC vector into bIC vector
    if( p_xDC )
    {  long int jj;
        // Correction of bIC vector by convoluting the amounts of DCs
        for( jj=0; jj<CSD->nDCb; jj++ )
            if( p_xDC[jj] > 0.0 )
                for( ii=0; ii<CSD->nICb; ii++ )
                    CNode->bIC[ii] += p_xDC[jj] * DCaJI( jj, ii );
    }
}

//(8b) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
//  structure contains a full list of Dependent Components used in GEM IPM2 calculations.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,   // Specific surface areas of phases, m2/kg [nPHb]  +       -      -
        double *p_xDC,  // Mole amounts of DCs [nDCb] - old primal soln.     +       -      -
        double *p_gam   // DC activity coefficients [nDCb] - old primal s.   +       -      -
        )
{
    long int ii;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];

    // Optional part - copying old primal solution from p_xDC and p_gam vectors
    if( p_xDC && p_gam )
    {
        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            CNode->xDC[ii] = p_xDC[ii];
            CNode->gam[ii] = p_gam[ii];
        }
    }
    else if( CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal solution provided!

    //  Discuss the policy!
    //   if( p_xDC )
    //   {  long int jj;
    //      // Correction of bIC vector by convoluting the amounts of DCs
    //      for( jj=0; jj<CSD->nDCb; jj++ )
    //        if( p_xDC[jj] )
    //          for( ii=0; ii<CSD->nICb; ii++ )
    //            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
    //   }

}


// (8c) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(long int  p_NodeHandle,   // Node identification handle
                        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                        //                                              GEM input output  FMT control
                        double p_TK,     // Temperature T, Kelvin                            +       -      -
                        double p_P,      // Pressure P, Pa                                   +       -      -
                        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
                        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_dll    // Lower restrictions to amounts of DC [nDCb]       +       -      -
                        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
}

// (8d) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(long int  p_NodeHandle,   // Node identification handle
                        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                        //                                              GEM input output  FMT control
                        double p_TK,     // Temperature T, Kelvin                            +       -      -
                        double p_P,      // Pressure P, Pa                                   +       -      -
                        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
                        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_asPH,  // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
                        double *p_omPH,   // Stability indices of phases,log10 scale [nPHb]  (+)      +      -
                        double *p_amru,   // Upper AMR to masses of sol. phases [nPSb]        +       -      -
                        double *p_amrl    // Lower AMR to masses of sol. phases [nPSb]        +       -      -
                        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
    if( CSD->nAalp >0 )
    { for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        CNode->omPH[ii] = p_omPH[ii];
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        CNode->amru[ii] = p_amru[ii];
        CNode->amrl[ii] = p_amrl[ii];
    }
}

//(8e) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
//  structure contains a full list of Dependent Components used in GEM IPM2 calculations.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,   // Specific surface areas of phases, m2/kg [nPHb]  +       -      -
        double *p_omPH,   // Stability indices of phases,log10 scale [nPHb]  (+)      +      -
        double *p_amru,   //< Upper AMR to masses of sol. phases [nPSb]      +       -      -
        double *p_amrl,   //< Lower AMR to masses of sol. phases [nPSb]      +       -      -
        double *p_xDC,  // Mole amounts of DCs [nDCb] - old primal soln.     +       -      -
        double *p_gam   // DC activity coefficients [nDCb] - old primal s.   +       -      -
        )
{
    long int ii;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
    { for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        CNode->omPH[ii] = p_omPH[ii];
    // Optional part - copying old primal solution from p_xDC and p_gam vectors
    if( p_xDC && p_gam )
    {
        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            CNode->xDC[ii] = p_xDC[ii];
            CNode->gam[ii] = p_gam[ii];
        }
    }
    else if( CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal solution provided!

    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        CNode->amru[ii] = p_amru[ii];
        CNode->amrl[ii] = p_amrl[ii];
    }
    //  Discuss the policy!
    //   if( p_xDC )
    //   {  long int jj;
    //      // Correction of bIC vector by convoluting the amounts of DCs
    //      for( jj=0; jj<CSD->nDCb; jj++ )
    //        if( p_xDC[jj] )
    //          for( ii=0; ii<CSD->nICb; ii++ )
    //            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
    //   }

}

//-----------------------End of node.cpp--------------------------


