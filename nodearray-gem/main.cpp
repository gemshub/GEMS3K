//-------------------------------------------------------------------
// $Id$
//
// Debugging version of a finite-difference 1D advection-diffusion
// mass transport model supplied by Dr. Frieder Enzmann (Uni Mainz)
// coupled with GEMIPM2K module for calculation of chemical equilibria
//
// Direct access to the TNodeArray class for storing all data for nodes
//
// Copyright (C) 2005,2007 S.Dmytriyeva, F.Enzmann, D.Kulik
//
//-------------------------------------------------------------------

#include "m_gem2mt.h"

//---------------------------------------------------------------------------
// Test of 1D advection (finite difference method provided by Dr. F.Enzmann,
// Uni Mainz) coupled with GEMIPM2K kernel (PSI) using the TNodeArray class.
// Finite difference calculations split over independent components
// (through bulk composition of aqueous phase).
// Experiments with smoothing terms on assigning differences to bulk composition
// of nodes

// "TestVTK/GEM2MT-task.json" "TestVTK/CaWBoundC-dat.lst" "TestVTK/CaWBoundC-dbr.lst" "TestVTK/VTK"

int main( int argc, char* argv[] )
{

#ifndef USE_NLOHMANNJSON
    std::string gem2mt_in1 = "TestAD/CalColumnAD.dat";
    std::string ipm_lst = "TestAD/CalcColumn-dat.lst";
    std::string dbr_lst = "TestAD/CalcColumn-dbr.lst";
    std::string vtk_fname = "TestAD/VTK";
#else
    std::string gem2mt_in1 = "TestAD1/CalcColumnAD.json";
    std::string ipm_lst = "TestAD1/CalcColumn-dat.lst";
    std::string dbr_lst = "TestAD1/CalcColumn-dbr.lst";
    std::string vtk_fname = "TestAD1/VTK";
#endif


    // from argv
    if (argc >= 2 )
        gem2mt_in1 = argv[1];
    if (argc >= 3 )
        ipm_lst = argv[2];
    if (argc >= 4 )
        dbr_lst = argv[3];
    if (argc >= 5 )
        vtk_fname = argv[4];

    gems3k_update_loggers( true, "gems3k_logger.log", spdlog::level::info);
    auto ar_logger_gems3k = spdlog::get("ipm");
    ar_logger_gems3k->set_level(spdlog::level::warn);

    try{

        if(gem2mt_in1.empty() || ipm_lst.empty() || dbr_lst.empty() )
            Error( "Start task", "No inital files");

        // The NodeArray must be allocated here
        TGEM2MT::pm = new TGEM2MT( 0 );

        // Here we read the GEM2MT structure, prepared from GEMS or by hand
        if( TGEM2MT::pm->ReadTask( gem2mt_in1.c_str(), vtk_fname.c_str() ))
            return 1;  // error reading files

        // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
        if( TGEM2MT::pm->MassTransInit( ipm_lst.c_str(), dbr_lst.c_str() ) )
            return 1;  // error reading files

        //   TGEM2MT::pm-> WriteTask( "gem2mt_out.dat" );

        // here we call the mass-transport finite-difference coupled routine
        TGEM2MT::pm->RecCalc();
    }
    catch(TError& err)
    {
        TNode::ipmlog_file->error("Error {} : {}", err.title, err.mess);
        return 1;
    }

    return 0;
}


//---------------------------------------------------------------------------

