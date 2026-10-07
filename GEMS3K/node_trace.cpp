//-------------------------------------------------------------------
/// \file node_trace.cpp
/// TNode::GEM_trace_regimes() - dilute-regime check for trace elements (see node.h).
//-------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <memory>
#include "node.h"
#include "datach_api.h"

std::vector<TNode::TraceRegime> TNode::GEM_trace_regimes( const std::vector<double>& factors, double traceRel,
                                                          double tol, long int mode,
                                                          const std::vector<std::string>& ofInterest )
{
    std::vector<TraceRegime> out;
    const long int nIC = CSD->nICb, nPH = CSD->nPHb, nICch = CSD->nIC;
    // no explicit list: the elements marked of interest on this node, if any
    const std::vector<std::string>& interest = ofInterest.empty() ? multi_base->elementsOfInterest : ofInterest;
    if( nIC < 1 || nPH < 1 ) return out;

    // ---- which ICs are trace (or of interest)
    std::vector<double> b0( CNode->bIC, CNode->bIC + nIC );
    auto isCharge = [&]( long int i ) { return CSD->ccIC[ IC_xDB_to_xCH( i ) ] == IC_CHARGE; };
    double total = 0.;
    for( long int i = 0; i < nIC; i++ ) if( !isCharge( i ) && b0[(size_t)i] > 0. ) total += b0[(size_t)i];
    std::vector<long int> trace;
    for( long int i = 0; i < nIC; i++ )
    {
        if( isCharge( i ) || !( b0[(size_t)i] > 0. ) ) continue;
        const std::string name = CSD->ICNL[ IC_xDB_to_xCH( i ) ];
        std::string trimmed = name; trimmed.erase( trimmed.find_last_not_of( " \t" ) + 1 );
        if( !interest.empty() )
        { if( std::find( interest.begin(), interest.end(), trimmed ) != interest.end() ) trace.push_back( i ); }
        else if( b0[(size_t)i] <= traceRel * total )
            trace.push_back( i );
    }
    if( trace.empty() ) return out;

    // ---- exact backup of the node (inputs and results): databr_copy() copies INTO CNode
    DATABR* backup = new DATABR;
    dbr_dch_api::databr_reset( backup, 1 );
    dbr_dch_api::databr_realloc( CSD, backup );
    { DATABR* live = CNode; CNode = backup; databr_copy( live ); CNode = live; }
    // ... and of MULTI: a warm GEM_run(false) starts from MULTI's retained primal, so restoring
    // DATABR alone would leave the next warm call starting from the last check solve. The
    // snapshot is a fresh TMultiBase filled by copyMULTI(); the restore copies values back
    // without reallocating (copyMULTIData(.., false)), since the TSolMod objects point into the
    // live arrays.
    std::unique_ptr<TMultiBase> savedMulti( new TMultiBase( this ) );
    savedMulti->set_def();
    savedMulti->copyMULTI( *multi_base );
    auto restoreNode = [&]() {
        databr_copy( backup );
        databr_free( backup );
        backup = nullptr;
        multi_base->copyMULTIData( *savedMulti, false );
    };

    std::vector<double> runFactors{ 1. };
    runFactors.insert( runFactors.end(), factors.begin(), factors.end() );
    struct Run { bool ok = false; std::vector<double> perPhase, xPH; std::vector<char> present; };  // perPhase[k*nIC+i]
    std::vector<Run> runs( runFactors.size() );
    std::vector<double> bc( (size_t)nIC );
    try {
    for( size_t r = 0; r < runFactors.size(); r++ )
    {
        databr_copy( backup );
        for( long int i : trace ) CNode->bIC[i] = b0[(size_t)i] * runFactors[r];
        CNode->NodeStatusCH = mode;
        const long int code = GEM_run( false );
        runs[r].ok = ( code == mode + 1 || code == mode + 2 );
        if( !runs[r].ok ) continue;
        runs[r].perPhase.assign( (size_t)(nPH * nIC), 0. );
        runs[r].present.assign( (size_t)nPH, 0 );
        runs[r].xPH.assign( CNode->xPH, CNode->xPH + nPH );
        for( long int k = 0; k < nPH; k++ )
        {
            Ph_BC( k, bc.data() );
            for( long int i = 0; i < nIC; i++ ) runs[r].perPhase[(size_t)(k*nIC + i)] = bc[(size_t)i];
            runs[r].present[(size_t)k] = CNode->xPH[k] > 0.;
        }
    }
    } catch( ... ) { restoreNode(); throw; }
    // phase amounts at the given amount (the factor-1 re-solve), for the stranded test
    double phaseTotal = 0.;
    const std::vector<double> xPH0 = runs[0].ok ? runs[0].xPH : std::vector<double>( (size_t)nPH, 0. );
    for( long int k = 0; k < nPH; k++ ) if( xPH0[(size_t)k] > 0. ) phaseTotal += xPH0[(size_t)k];

    // ---- species -> phase (DATACH order) and single-species phases
    std::vector<long int> phaseOfDC( (size_t)CSD->nDC, -1 );
    for( long int k = 0, j = 0; k < CSD->nPH; j += CSD->nDCinPH[k], k++ )
        for( long int jj = j; jj < j + CSD->nDCinPH[k]; jj++ ) phaseOfDC[(size_t)jj] = k;
    auto phaseName = [&]( long int kDB ) {
        std::string s = CSD->PHNL[ Ph_xDB_to_xCH( kDB ) ];
        s.erase( s.find_last_not_of( " \t" ) + 1 ); return s;
    };

    for( long int i : trace )
    {
        TraceRegime tr;
        tr.xIC = i;
        tr.name = CSD->ICNL[ IC_xDB_to_xCH( i ) ];
        tr.name.erase( tr.name.find_last_not_of( " \t" ) + 1 );
        tr.amount = b0[(size_t)i];
        if( !runs[0].ok ) { tr.verdict = "FAILED"; out.push_back( tr ); continue; }

        auto fractions = [&]( const Run& R ) {
            std::vector<double> f( (size_t)nPH, 0. ); double t = 0.;
            for( long int k = 0; k < nPH; k++ ) t += R.perPhase[(size_t)(k*nIC + i)];
            if( t > 0. ) for( long int k = 0; k < nPH; k++ ) f[(size_t)k] = R.perPhase[(size_t)(k*nIC + i)] / t;
            return f;
        };
        const std::vector<double> f0 = fractions( runs[0] );
        for( long int k = 0; k < nPH; k++ )
            if( f0[(size_t)k] > 1e-4 ) tr.phases.push_back( { phaseName( k ), f0[(size_t)k] } );
        std::sort( tr.phases.begin(), tr.phases.end(), []( auto& a, auto& b ) { return a.second > b.second; } );

        // carrier phases (DBR indices) of this IC
        const long int ich = IC_xDB_to_xCH( i );
        std::vector<char> carrier( (size_t)nPH, 0 );
        for( long int j = 0; j < CSD->nDC; j++ )
            if( CSD->A[ ich + j*nICch ] != 0. )
            {
                const long int kDB = Ph_xCH_to_xDB( phaseOfDC[(size_t)j] );
                if( kDB >= 0 && kDB < nPH ) carrier[(size_t)kDB] = 1;
            }

        bool failed = false;
        for( size_t r = 1; r < runs.size(); r++ )
        {
            if( !runs[r].ok ) { failed = true; continue; }
            const std::vector<double> f = fractions( runs[r] );
            for( long int k = 0; k < nPH; k++ )
                tr.maxFractionChange = std::max( tr.maxFractionChange, std::fabs( f[(size_t)k] - f0[(size_t)k] ) );
        }
        auto singleSpecies = [&]( long int kDB ) { return CSD->nDCinPH[ Ph_xDB_to_xCH( kDB ) ] == 1; };
        if( failed )
            tr.verdict = "FAILED";
        else if( tr.maxFractionChange <= tol )
            tr.verdict = "LINEAR";
        else
        {
            for( long int k = 0; k < nPH; k++ )
            {
                if( !carrier[(size_t)k] || !singleSpecies( k ) ) continue;
                bool seenIn = false, seenOut = false;
                for( const Run& R : runs ) if( R.ok ) ( R.present[(size_t)k] ? seenIn : seenOut ) = true;
                if( seenIn && seenOut ) tr.boundaryPhases.push_back( phaseName( k ) );
            }
            bool pureHost = false;
            for( long int k = 0; k < nPH; k++ )
                if( carrier[(size_t)k] && singleSpecies( k ) && runs[0].present[(size_t)k] ) pureHost = true;
            double dmin = 1e300, dmax = 0.;
            for( const Run& R : runs )
            {
                if( !R.ok ) continue;
                double d = 0.;
                for( long int k = 0; k < nPH; k++ ) if( !singleSpecies( k ) ) d += R.perPhase[(size_t)(k*nIC + i)];
                dmin = std::min( dmin, d ); dmax = std::max( dmax, d );
            }
            double fmin = 1e300, fmax = 0.;
            for( double f : runFactors ) { fmin = std::min( fmin, f ); fmax = std::max( fmax, f ); }
            const bool flat = dmin > 0. && dmax / dmin < 1.5 && fmax / fmin >= 10.;
            if( !tr.boundaryPhases.empty() ) tr.verdict = "BOUNDARY";
            else if( pureHost && flat )      tr.verdict = "SATURATED";
            else                             tr.verdict = "NONLINEAR";
        }
        // stranded: every carrier in ONE multi-species phase that is <= 1e-6 of the system
        long int only = -1; bool one = true;
        for( long int k = 0; k < nPH && one; k++ )
            if( carrier[(size_t)k] ) { if( only < 0 ) only = k; else one = false; }
        if( one && only >= 0 && !singleSpecies( only ) && phaseTotal > 0. && xPH0[(size_t)only] <= 1e-6 * phaseTotal )
            tr.verdict = ( tr.verdict == "LINEAR" ) ? std::string( "STRANDED" ) : "STRANDED+" + tr.verdict;
        out.push_back( tr );
    }

    // ---- restore the node exactly (DATABR and MULTI)
    restoreNode();
    return out;
}

long int TNode::GEM_set_elements_of_interest( const std::vector<std::string>& names )
{
    std::vector<std::string> kept;
    for( std::string n : names )
    {
        n.erase( 0, n.find_first_not_of( " \t" ) );
        n.erase( n.find_last_not_of( " \t" ) + 1 );
        bool found = false;
        for( long int i = 0; i < CSD->nIC && !found; i++ )
        {
            std::string ic = CSD->ICNL[i];
            ic.erase( ic.find_last_not_of( " \t" ) + 1 );
            found = ( ic == n );
        }
        if( !found )
            node_logger->warn( "GEM_set_elements_of_interest: '{}' is not an independent component of this system - skipped", n );
        else if( std::find( kept.begin(), kept.end(), n ) == kept.end() )
            kept.push_back( n );
    }
    multi_base->elementsOfInterest = kept;
    return (long int)kept.size();
}

const std::vector<std::string>& TNode::GEM_elements_of_interest() const
{
    return multi_base->elementsOfInterest;
}
