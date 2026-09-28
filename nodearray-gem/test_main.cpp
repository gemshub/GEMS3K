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
#include "GEMS3K/jsonconfig.h"

static int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_A(const std::string& ipm_lst, const std::string& dbr_lst);
static int task_C(const std::string& ipm_lst, const std::string& dbr_lst);
static int task_W(const std::string& ipm_lst, const std::string& dbr_lst);
static int task_F(const std::string& ipm_lst, const std::string& dbr_lst);
static int task_S(const std::string& ipm_lst, const std::string& dbr_lst);
static int task_B(const std::string& ipm_lst, const std::string& dbr_lst);

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
#else
    std::string gem2mt_in1 = "TestAD1/CalcColumnAD.json";
    std::string ipm_lst = "TestAD1/CalcColumn-dat.lst";
    std::string dbr_lst = "TestAD1/CalcColumn-dbr.lst";
#endif


    // from argv
    if (argc >= 2 )
        gem2mt_in1 = argv[1];
    if (argc >= 3 )
        ipm_lst = argv[2];
    if (argc >= 4 )
        dbr_lst = argv[3];

    GemsSettings().gems3k_update_loggers( true, "gems3k_logger.log", spdlog::level::info);


    try{
        //return task_from_file(gem2mt_in1, ipm_lst, dbr_lst);
        //return task_A(ipm_lst, dbr_lst);
        //return task_C(ipm_lst, dbr_lst);
        //return task_W(ipm_lst, dbr_lst);
        return task_F("TestF/CalcColumn-dat.lst", "TestF/CalcColumn-dbr.lst");
        //return task_S("TestS/CalcColumn-dat.lst", "TestS/CalcColumn-dbr.lst");
        //return task_B("TestB/CalcColumn-dat.lst", "TestB/CalcColumn-dbr.lst");
    }
    catch(TError& err) {
        TNode::ipmlog_file->error("Error {} : {}", err.title, err.mess);
        return 1;
    }

    return 0;
}

//---------------------------------------------------------------------------
// Test of 1D advection (finite difference method provided by Dr. F.Enzmann,
// Uni Mainz) coupled with GEMIPM2K kernel (PSI) using the TNodeArray class.
// Finite difference calculations split over independent components
// (through bulk composition of aqueous phase).
// Experiments with smoothing terms on assigning differences to bulk composition
// of nodes
// "TestVTK/GEM2MT-task.json" "TestVTK/CaWBoundC-dat.lst" "TestVTK/CaWBoundC-dbr.lst"
int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(gem2mt_file.empty() || ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT(0) );
    TGEM2MT::pm = mt_task.get();

    // Here we read the GEM2MT structure, prepared from GEMS or by hand
    if(TGEM2MT::pm->ReadTask(gem2mt_file, "")) {
        return 1;  // error reading files
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}


int task_A(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('A', 201) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_ON);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 500000, 5000);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task->setSpatialDimensions(0.2, 0, 0);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task->setTimeStepReductionFactor(1.);
    // Set initial total node volume (m^3)
    mt_task->setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task->setFluidAdvectionVelocity(1e-7);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(1.);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1.);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task->setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task->setInitialDiffusivity(0.);
    // Set initial tortuosity factor, usually 1
    mt_task->setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(1e-11);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}

// GUI record 'CalDolCol2:G:CalcColumn:0:0:1:25:0:1D-Cc-MgCl2-tmod:C:'
int task_C(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('C', 101) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_OFF);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 30000, 3);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task->setSpatialDimensions(0.2, 0, 0);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task->setTimeStepReductionFactor(2.);
    // Set initial total node volume (m^3)
    mt_task->setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task->setFluidAdvectionVelocity(2e-6);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(0.5);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1e-12);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task->setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task->setInitialDiffusivity(2e-09);
    // Set initial tortuosity factor, usually 1
    mt_task->setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(1e-11);

    // If need update input node distributing
    // Change index of initial system variant for the first node
    //mt_task->setDistributing(0, 0);
    // Change the type for for the first node
    //mt_task->setNodeType(0, 3);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // Change/define some other gem2mt arrays
    // Set list of selected fields and indexes to VTK format
    mt_task->setVTKfields({{33,0},{32,0},{38,1},{38,4}});

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}

// GUI record 'CalDolCol2:G:CalcColumn:0:0:1:25:0:1D-DifMgCl2-center:W:'
int task_W(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('W', 101) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_ON);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set number of allocated particle types < 20
    mt_task->setNumberParticleTypes(1);
    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 1000000, 1000);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task->setSpatialDimensions(0.2, 0, 0);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task->setTimeStepReductionFactor(5.);
    // Set initial total node volume (m^3)
    mt_task->setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task->setFluidAdvectionVelocity(2e-6);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(0.5);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1e-12);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task->setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task->setInitialDiffusivity(2e-09);
    // Set initial tortuosity factor, usually 1
    mt_task->setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(1e-11);

    // If need update input node distributing
    // Change index of initial system variant for the center node
    mt_task->setDistributing(50, 0);
    // Change the type for for the center node
    mt_task->setNodeType(50, 3);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // Change/define some other gem2mt arrays

    // Set list of selected fields and indexes to VTK format
    mt_task->setVTKfields({{33,0},{32,0},{38,1},{38,4}});
    // Set of particle statistic property
    mt_task->setParticle(0, 1000, 500, 1500, {0, 11, 21, 0, 0, 0});

    //If need update initial hydraulic parameters in node: Vt, vp, eps, Km, al, Dif,  nto
    //mt_task->setHydraulicParameters(1, 0.0015, 2.1e-06, 0.6, 1e-11, 0.0015, 3e-09, 1);

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}


//  CalDolCol2:G:CalcColumn:0:0:1:25:0:Test2:F:
int task_F(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('F', 51) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test F mode (flow-through reactors)");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_OFF);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set number of mobile groups of phases, nMGP >= 0
    mt_task->setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem, nFD >= 0
    mt_task->setNumberMGPfluxes(51);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 1200, 1);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task->setSpatialDimensions(0, 0, 0);

    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(1e-9);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1e-12);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(0);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(0);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // Change/define some other gem2mt arrays

    for(long int ii=0; ii<mt_task->nPhaseGroups(); ++ii) {
        // Set ID of mobile phase group
        mt_task->setPhaseGroupsID(ii, "phg"+std::to_string(ii+1));
        for(long int k=0; k<mt_task->nPhases(); ++k) {
            // Set quantities of phases in MGP
            mt_task->setPhaseGroupsQuantities(ii, k, (k==0 ? 1.: 0.));
        }
    }
    for(long int k=0; k<mt_task->nPhases(); ++k) {
        // Set units for setting phase quantities in MGP (see PGT )
        mt_task->setUnitsPhaseQuantities(k, 'M');
    }


    for(long int ii=0; ii<mt_task->nMGPfluxes(); ++ii) {
        // Set Source/Receive box index in the flux definition
        mt_task->setFluxSourceReceive(ii, ii, (ii<mt_task->nMGPfluxes()-1 ? ii+1: -1));
        // Set the flux defnition: flux order, flux rate, MGP quantity
        mt_task->setFluxSourceReceive(ii, 1., 0.5, 0, 0);
        //  Set the ID of MGP to move in this flux
        mt_task->setFluxMGPid(ii, "phg1");
        //  Set IDs of fluxes
        mt_task->setFluxIDs(ii, "qj");
     }

    // Set list of selected fields and indexes to VTK format
    mt_task->setVTKfields({{41,2},{41,6},{41,25},{41,26},{41,17},{41,13}});
    mt_task->setOutVTK(false);  // allocated but no write

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}


// CalDolCol2:G:CalcColumn:0:0:1:25:0:Test3:S:
int task_S(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('S', 51) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test S mode (initial 2.5 mmol/L calcite)  1 to 96 C");
    mt_task->setComment("`");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_OFF);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set number of mobile groups of phases, nMGP >= 0
    mt_task->setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem, nFD >= 0
    mt_task->setNumberMGPfluxes(51);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 30, 0.3);

    // ?? only to compare Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(1e-9);
    // ?? only to compare Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1e-12);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(0);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(0);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // Change/define some other gem2mt arrays
    // Set units for aq phase quantities
    mt_task->setUnitsPhaseQuantities(0, 'n');
    for(long int ii=0; ii<mt_task->nMGPfluxes(); ++ii) {
        // Set the flux defnition: flux order, flux rate, MGP quantity
        mt_task->setFluxSourceReceive(ii, 0., 1., 0, 0);
    }

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}


// CalDolCol2:G:CalcColumn:0:0:1:25:0:Test4:B:
int task_B(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('B', 50) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test B mode (flow-through reactors), T = 1 to 99 C");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA(S_OFF);
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set number of mobile groups of phases, nMGP >= 0
    mt_task->setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem, nFD >= 0
    mt_task->setNumberMGPfluxes(50);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 10000, 1);

    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(1e-9);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1e-12);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(0);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(0);

    // Change from default: Other nodes normal, the last two set as constant source and sink
    //mtp->DiCp[mtp->qc][1] = (mtp->qc<mtp->nC-2 ? mtp->DiCp[mtp->qc][1]: 3);
    mt_task->setNodeType(48, 0);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // Change/define some other gem2mt arrays

    // Set ID of the first mobile phase group
    mt_task->setPhaseGroupsID(0, "Pg1");
    for(long int ii=0; ii<mt_task->nMGPfluxes(); ++ii) {
        //  Set the ID of MGP to move in this flux
        mt_task->setFluxMGPid(ii, "Pg1");
    }
    // Change from default: flux order 1 (proportional to source MPG mass)
    mt_task->setFluxSourceReceive(0, 1., 0.1, 0, 0);

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}
//---------------------------------------------------------------------------

