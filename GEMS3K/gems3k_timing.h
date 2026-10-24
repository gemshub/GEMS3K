//-------------------------------------------------------------------
/// \file gems3k_timing.h
/// Optional time split of one GEM_run() call into solver parts.
/// Active only when GEMS3K_TIMING_FILE is set; each call then appends one line to that file.
/// Times are inclusive: a part called inside another is counted in both.
//-------------------------------------------------------------------
#ifndef GEMS3K_TIMING_H
#define GEMS3K_TIMING_H

#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace gems3k_timing {

enum Part {
    Call,          // whole GEM_run()
    Native,        // native IPM + MBR solve
    OptimaCall,    // one Optima-mode solve
    OSolve,        // inside Optima's solve()
    Obj,           // objective callback
    ObjFxx,        // objective calls that also build the Hessian
    FdMix,         // finite-difference Hessian columns of mixed-phase species
    FdPure,        // finite-difference Hessian columns of pure-phase species
    FdBlock,       // phase curvature blocks (pa_PhaseHessianFloor)
    ObjRestore,    // the activity call that restores the base state after the columns
    ActCoef,       // CalculateActivityCoefficients, anywhere
    Psf,           // PotentialSpaceFinish
    Tpd,           // NativeTpdPhase
    Wpsv,          // WorstPhaseStabilityViolation
    PreSolve,      // OptimaReducedPreSolve
    Level,         // leveled start (column generation)
    NParts
};

inline const char* name( int p )
{
    static const char* n[NParts] = { "call", "native", "optima", "osolve", "obj", "objfxx", "fdmix", "fdpure",
                                     "fdblock", "objrestore", "act", "psf", "tpd", "wpsv", "presolve", "level" };
    return n[p];
}

inline bool enabled()
{
    static const bool on = ( std::getenv( "GEMS3K_TIMING_FILE" ) != nullptr );
    return on;
}

struct Acc { double sec[NParts] = {}; long calls[NParts] = {}; };
inline thread_local Acc acc;

/// Adds the time of its own lifetime to one part.
class Scope
{
public:
    explicit Scope( Part p ) : part( p ), on( enabled() )
    {
        if( on ) t0 = std::chrono::steady_clock::now();
    }
    ~Scope()
    {
        if( !on ) return;
        acc.sec[part] += std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
        acc.calls[part]++;
    }
    Scope( const Scope& ) = delete;
    Scope& operator=( const Scope& ) = delete;
private:
    Part part;
    bool on;
    std::chrono::steady_clock::time_point t0;
};

/// Writes one line for the finished call and clears the sums.
inline void flush( long requested, long status, long iterations )
{
    if( !enabled() ) return;
    if( FILE* f = std::fopen( std::getenv( "GEMS3K_TIMING_FILE" ), "a" ) )
    {
        std::fprintf( f, "TIMING req=%ld status=%ld it=%ld", requested, status, iterations );
        for( int p = 0; p < NParts; p++ )
            if( acc.calls[p] )
                std::fprintf( f, " %s=%.6f/%ld", name( p ), acc.sec[p], acc.calls[p] );
        std::fprintf( f, "\n" );
        std::fclose( f );
    }
    acc = Acc();
}

} // namespace gems3k_timing

#endif // GEMS3K_TIMING_H
