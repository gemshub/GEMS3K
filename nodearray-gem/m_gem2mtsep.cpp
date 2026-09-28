//-------------------------------------------------------------------
// $Id: main.cpp 792 2006-09-19 08:10:41Z gems $
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
#include "GEMS3K/nodearray.h"
#include "GEMS3K/io_keyvalue.h"
#include "GEMS3K/io_simdjson.h"
#ifdef USE_NLOHMANNJSON
#include "GEMS3K/io_nlohmann.h"
#endif

TGEM2MT* TGEM2MT::pm;

TGEM2MT::TGEM2MT( size_t /*nrt*/ )
{
    mtp=&mt[0];
    set_def(0);
    ////mtp->PvMO =   S_ON;
    ////mtp->iStat =  AS_READY;
    na = 0;
    pa_mt = 0;
}

TGEM2MT::TGEM2MT(char ps_mode, long n_nodes)
{
    mtp=&mt[0];
    set_def(0);
    na = 0;
    pa_mt = 0;
    mtp->PsMode = ps_mode;
    mtp->nC = n_nodes;

    math_transport_defaults();
}

TGEM2MT::~TGEM2MT()
{
    mem_kill(0);
    if( pa_mt )
        delete pa_mt;
}

//Calculate record
void TGEM2MT::RecCalc()
{
    try {

        // use particles
        if(mtp->PsMode == RMT_MODE_W) {
            na->SetGrid(mtp->sizeLc, mtp->grid);   // set up grid structure
            pa_mt = new TParticleArray(mtp->nPTypes, mtp->nProps,
                                       mtp->NPmean, mtp->ParTD, mtp->nPmin, mtp->nPmax, na.get());
            pa_mt->setUpCounters();
        }

        // put HydP
        if(mtp->PsMode != RMT_MODE_S && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B) {
            putHydP(na->pNodT0());
            putHydP(na->pNodT1());
        }

        if(mtp->PsVTK != S_OFF)  {
            u_create_directory(pathVTK+nameVTK+"/");
        }

        if(mtp->iStat != AS_RUN) {  // ???? Could we stop and restart calculations?
            mtp->gStat = GS_GOING;
            mt_reset();
            mtp->gStat = GS_DONE;
        }

        // internal calc
        auto iret = internalCalc();
        if(!iret) {
            mtp->iStat = AS_DONE;
        }
        //else we have a stop point
    }
    catch(TError& xcpt)  {
        mtp->gStat = GS_ERR;
        mtp->iStat = AS_INDEF;
        throw xcpt;
    }
}

// read TGEM2MT structure from file
int TGEM2MT::ReadTask(const std::string& gem2mt_file, const std::string& vtk_dir)
{
    try {
        std::string gem2mt_in = gem2mt_file;
        std::fstream ff(gem2mt_in, std::ios::in );
        ErrorIf(!ff.good(), gem2mt_in, "Fileopen error");

        if(gem2mt_in.rfind(".json") != std::string::npos) {
#ifdef USE_NLOHMANNJSON
            io_formats::NlohmannJsonRead in_format(ff, "", "gem2mt");
            from_text_file( in_format );
        }
#else
            io_formats::SimdJsonRead in_format(ff, "", "gem2mt");
            from_text_file( in_format );
        }
#endif
        else {
            io_formats::KeyValueRead in_format( ff );
            from_text_file( in_format );
        }

        pathVTK = vtk_dir;
        if(!pathVTK.empty()) {
            pathVTK += "/";
        }
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
    }
    return 1;
}

int TGEM2MT::ReadTaskString(const std::string json_string)
{
    if(json_string.empty()) {
        return 1;
    }

    try  {
        std::stringstream ss;
        ss.str(json_string);
#ifdef USE_NLOHMANNJSON
        io_formats::NlohmannJsonRead in_format(ss, "", "gem2mt");
        from_text_file( in_format );
#else
        io_formats::SimdJsonRead in_format(ss, "", "gem2mt");
        from_text_file(in_format);
#endif
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
    }
    return 1;
}

// Write TGEM2MT structure to file
int TGEM2MT::WriteTask(const std::string& gem2mt_file)
{
    try  {
        std::string gem2mt_out = gem2mt_file;
        std::fstream ff(gem2mt_out, std::ios::out);
        ErrorIf(!ff.good(), gem2mt_out, "Fileopen error");

        if(gem2mt_out.rfind(".json") != std::string::npos) {
#ifdef USE_NLOHMANNJSON
            io_formats::NlohmannJsonWrite out_format( ff, "");
            to_text_file(out_format, true, false );
        }
#else
            io_formats::SimdJsonWrite out_format(ff, "", true);
            to_text_file(out_format, true, false);
        }
#endif
        else  {
            io_formats::KeyValueWrite out_format(ff);
            to_text_file(out_format, true, false);
        }
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() <<std:: endl;
    }
    return 1;
}

void TGEM2MT::default_VTK(const std::string& work_path)
{
    std::string folder, fname, ext;
    if(!work_path.empty()) {
        u_splitpath(work_path, folder, fname, ext);
    }
    auto pos = fname.rfind("-");
    if(pos != std::string::npos) {
        fname = fname.substr(0, pos);
    }
    if(pathVTK.empty()) {
        pathVTK = folder+"VTK/";
    }
    if(fname.empty()) {
        fname = "vtk";
    }
    if(nameVTK.empty()) {
        nameVTK = fname;
    }
    prefixVTK = nameVTK;
}


// Set up math transport default values7sizes in constructor
void TGEM2MT::math_transport_defaults()
{
    // from set_def(int q), ask gem2mt users for better defaults
    memset( &mtp->Msysb1, 0, sizeof(double)*20 );
    memset( mtp->size[0], 0, sizeof(float)*8 );

    mtp->nVTKfld = 0;
    mtp->Tau[START_] = 0.;
    mtp->Tau[STOP_] = 1000.;
    mtp->Tau[STEP_] = 1.;
    mtp->ntM =1000;
    mtp->cdv = 1e-9;
    mtp->cez = 1e-12;
    mtp->tf = 1.;

    mtp->Msysb1 = 0.;
    mtp->Vsysb1 = 0.;
    mtp->Mwatb1 = 1.;
    mtp->Maqb1 = 1.;
    mtp->Vaqb1 = 1.;

    // Alloc important arrays
    mtp->DiCp = new long int[ mtp->nC][2];
    defaults_DiCp();

    // Alloc the NodeArray
    na = TNodeArray::create(mtp->nC);
    TNodeArray::na = na.get();

    if(mtp->PsMode == RMT_MODE_W) {
        mtp->nPTypes = 10;
    }

    mtp->PsMO =   S_ON;
    mtp->iStat =  AS_READY;
}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
int TGEM2MT::gem3k_files_read(const std::string& ipm_lst_file, const std::string& dbr_lst_file)
{
    // define name of vtk file
    default_VTK(ipm_lst_file);

    // Prepare the array for initial conditions allocation
    std::vector<long int> nodeType;
    for(int ii=0; ii<mtp->nC; ++ii) {
        nodeType.push_back(mtp->DiCp[ii][0]);
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    // if mtp->iStat == AS_RUN we resume calculation ????
    if(na->GEM_init(ipm_lst_file.c_str(), dbr_lst_file.c_str(), nodeType.data(), mtp->iStat == AS_RUN)) {
        // ???? add logger
        return 1;  // error reading files
    }
    return 0;
}

// Set up NodeArray and ParticleArray classes after reading gems3k files
int TGEM2MT::restore_data_from_gems3k(const std::vector<std::string>& dbr_names)
{
    int ii;
    CalcIPM(NEED_GEM_AIA, 0, mtp->nC); //recalc all nodes ?

    // realloc gem2mt memory  (if not read exported gem2mt)
    if(true) {

        // Restore sizes from gems3k export
        mtp->Lsf = na->pCSD()->nDCs;
        mtp->Nf = na->pCSD()->nIC;
        mtp->FIf = na->pCSD()->nPH;
        mtp->nTai = na->pCSD()->nTp;
        mtp->nPai = na->pCSD()->nPp;

        // allocate memory and setup default values
        mem_new(0);
        init_arrays(true);
        // initialization from scripts in GUI
        if(mtp->HydP) {
            defaults_HydP();
        }
        if(mtp->PvGrid != S_OFF && mtp->grid) {
            defaults_Grid();
        }

        // read names
        std::string name;
        for(ii=0; ii<mtp->nIV; ++ii) {
            if(ii<dbr_names.size()) {
                name = std::to_string(ii)+dbr_names[ii];
                strncpy(mtp->nam_i[ii], name.c_str(), MAXIDNAME);
            }
        }
    }

    for(ii=0; ii<mtp->nTai; ++ii) {
        mtp->Tval[ii] = na->pCSD()->TKval[ii]-C_to_K;
    }
    for(ii=0; ii<mtp->nPai; ++ii) {
        mtp->Pval[ii] = na->pCSD()->Pval[ii]/bar_to_Pa;
    }

    return 0;
}

// Here we read the MULTI structure, DATACH and DATABR strings prepared from GEMS
int TGEM2MT::gems3k_strings(const std::string& dch_json, const std::string& ipm_json,
                            const std::vector<std::string>& dbr_json)
{
    // define name of vtk file
    default_VTK("");

    // Prepare the array for initial conditions allocation
    std::vector<long int> nodeType;
    for(int ii=0; ii<mtp->nC; ++ii) {
        nodeType.push_back(mtp->DiCp[ii][0]);
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(na->GEM_init(dch_json, ipm_json, dbr_json, nodeType.data())) {
        // ???? add logger
        return 1;  // error reading files
    }
    return 0;
}


int TGEM2MT::MassTransInit(const std::string &ipm_lst_file, const std::string &dbr_lst_file)
{
    if(gem3k_files_read(ipm_lst_file, dbr_lst_file)) {
        return 1;
    }

    // get dbr file names
    GEMS3KGenerator generator(ipm_lst_file);
    mtp->nIV = generator.load_dbr_lst_file(dbr_lst_file);

    restore_data_from_gems3k(generator.dbr_names());
    return 0;
}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up NodeArray and ParticleArray classes
int TGEM2MT::MassTransStringInit(const std::string& dch_json, const std::string& ipm_json,
                                 const std::vector<std::string>& dbr_json)
{
    if(gems3k_strings(dch_json, ipm_json, dbr_json)) {
        return 1;
    }

    mtp->nIV = dbr_json.size();
    restore_data_from_gems3k({});
    return 0;
}

void TGEM2MT::setVTKfields(const std::vector<std::pair<int, int>> &vtk_fields)
{
    if(mtp->xVTKfld) {
        delete[]  mtp->xVTKfld;
        mtp->xVTKfld = nullptr;
    }
    mtp->nVTKfld = vtk_fields.size();
    if(mtp->nVTKfld>0) {
        mtp->xVTKfld = new long int[mtp->nVTKfld][2];
        for(int ii=0; ii<mtp->nVTKfld; ++ii) {
            mtp->xVTKfld[ii][0] =vtk_fields[ii].first;
            mtp->xVTKfld[ii][1] =vtk_fields[ii].second;
        }
        mtp->PvnVTK = S_ON;
        mtp->PsVTK = mtp->PvnVTK;
    }
    else {
        mtp->PvnVTK = S_OFF;
        mtp->PsVTK = mtp->PvnVTK;
    }
}

void TGEM2MT::setParticle(long int pndx, long int pmean, long int pmin, long int pmax, const std::array<long int, 6> &pparam)
{
    if(mtp->PsMode == RMT_MODE_W && pndx<mtp->nPTypes) {
        mtp->NPmean[pndx] = pmean;
        mtp->nPmin[pndx] = pmin;
        mtp->nPmax[pndx] = pmax;
        std::copy(pparam.begin(), pparam.end(), mtp->ParTD[pndx]);
    }
}

void TGEM2MT::setHydraulicParameters(long pndx, double Vt, double vp, double eps, double Km, double al, double Dif, double nto)
{
    if(mtp->HydP && pndx<mtp->nC) {
        mtp->HydP[pndx][0] = Vt;
        mtp->HydP[pndx][1] = vp;
        mtp->HydP[pndx][2] = eps;
        mtp->HydP[pndx][3] = Km;
        mtp->HydP[pndx][4] = al;
        mtp->HydP[pndx][5] = Dif;
        mtp->HydP[pndx][6] = nto;
    }
}

void TGEM2MT::setPhaseGroupsID(long gndx, const std::string &ids)
{
    if(mtp->MGPid && gndx<mtp->nPG) {
        strncpy( mtp->MGPid[gndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setUnitsPhaseQuantities(long pndx, char units)
{
    if(mtp->UMGP && pndx<mtp->FIf) {
        mtp->UMGP[pndx] = units;
    }
}

void TGEM2MT::setPhaseGroupsQuantities(long gndx, long pndx, double quantity)
{
    if(mtp->MGPid && gndx<mtp->nPG && pndx<mtp->FIf) {
        mtp->PGT[gndx*mtp->FIf+pndx] = quantity;
    }
}

void TGEM2MT::setICsourceQuantities(long gndx, long indx, double quantity)
{
    if(mtp->MGPid && gndx<mtp->nSFD && indx<mtp->Nf) {
        mtp->BSF[gndx*mtp->Nf+indx] = quantity;
    }
}

void TGEM2MT::setFluxSourceReceive(long pndx, double order, double rate, double quantity, double val)
{
    if(mtp->FDLf && pndx<mtp->nFD) {
        mtp->FDLf[pndx][0] = order;
        mtp->FDLf[pndx][1] = rate;
        mtp->FDLf[pndx][2] = quantity;
        mtp->FDLf[pndx][3] = val;
    }
}

void TGEM2MT::setFluxMGPid(long pndx, const std::string &ids)
{
    if(mtp->FDLmp && pndx<mtp->nFD) {
        strncpy( mtp->FDLmp[pndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setFluxIDs(long pndx, const std::string &ids)
{
    if(mtp->FDLid && pndx<mtp->nFD) {
        strncpy( mtp->FDLid[pndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setGridPoint(long pndx, double x, double y, double z)
{
    if(mtp->grid && pndx<mtp->nC) {
        mtp->grid[pndx][0] = x;
        mtp->grid[pndx][1] = y;
        mtp->grid[pndx][2] = z;
    }
}

void TGEM2MT::setFluxSourceReceive(long pndx, long source, long receive)
{
    if(mtp->FDLi && pndx<mtp->nFD) {
        mtp->FDLi[pndx][0] = source;
        mtp->FDLi[pndx][1] = receive;
    }
}

//==========================================================================================

// free dynamic memory in objects and values
void TGEM2MT::mem_kill(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
            "E05GTrem: Attempt to access corrupted dynamic memory.");

    //- if( mtp->lNam) delete[] mtp->lNam;
    //- if( mtp->lNamE) delete[] mtp->lNamE;
    //- if( mtp->tExpr) delete[] mtp->tExpr;
    //- if( mtp->gExpr) delete[] mtp->gExpr;
    if(mtp->sdref) {
        delete[] mtp->sdref;
        mtp->sdref = nullptr;
    }
    if(mtp->sdval) {
        delete[] mtp->sdval;
        mtp->sdval = nullptr;
    }
    if(mtp->DiCp) {
        delete[] mtp->DiCp;
        mtp->DiCp = nullptr;
    }
    if(mtp->FDLi) {
        delete[] mtp->FDLi;
        mtp->FDLi = nullptr;
    }
    //- if( mtp->PTVm) delete[] mtp->PTVm;
    //- if( mtp->StaP) delete[] mtp->StaP;
    if(mtp->xVTKfld) {
        delete[] mtp->xVTKfld;
        mtp->xVTKfld = nullptr;
    }
    //- if( mtp->xEt) delete[] mtp->xEt;
    //- if( mtp->yEt) delete[] mtp->yEt;
    //- if( mtp->Bn) delete[] mtp->Bn;
    if(mtp->HydP) {
        delete[] mtp->HydP;
        mtp->HydP = nullptr;
    }
    //- if( mtp->qpi) delete[] mtp->qpi;
    //- if( mtp->qpc) delete[] mtp->qpc;
    //- if( mtp->xt) delete[] mtp->xt;
    //- if( mtp->yt) delete[] mtp->yt;
    //- if( mtp->CIb) delete[] mtp->CIb;
    //- if( mtp->CAb) delete[] mtp->CAb;
    if(mtp->FDLf) {
        delete[] mtp->FDLf;
        mtp->FDLf = nullptr;
    }
    if(mtp->PGT) {
        delete[] mtp->PGT;
        mtp->PGT = nullptr;
    }
    if(mtp->Tval) {
        delete[] mtp->Tval;
        mtp->Tval = nullptr;
    }
    if(mtp->Pval) {
        delete[] mtp->Pval;
        mtp->Pval = nullptr;
    }
    if(mtp->nam_i) {
        delete[] mtp->nam_i;
        mtp->nam_i = nullptr;
    }
    //- if( mtp->for_i) delete[] mtp->for_i;
    //- if( mtp->stld) delete[] mtp->stld;
    //- if( mtp->CIclb) delete[] mtp->CIclb;
    //- if( mtp->AUcln) delete[] mtp->AUcln;
    if(mtp->FDLid) {
        delete[] mtp->FDLid;
        mtp->FDLid = nullptr;
    }
    if(mtp->FDLop) {
        delete[] mtp->FDLop;
        mtp->FDLop = nullptr;
    }
    if(mtp->FDLmp) {
        delete[] mtp->FDLmp;
        mtp->FDLmp = nullptr;
    }
    if(mtp->MGPid) {
        delete[] mtp->MGPid;
        mtp->MGPid = nullptr;
    }
    if(mtp->UMGP) {
        delete[] mtp->UMGP;
        mtp->UMGP = nullptr;
    }
    //- if( mtp->SBM) delete[] mtp->SBM;
    if(mtp->BSF) {
        delete[] mtp->BSF;
        mtp->BSF = nullptr;
    }
    if(mtp->MB) {
        delete[] mtp->MB;
        mtp->MB = nullptr;
    }
    if(mtp->dMB) {
        delete[] mtp->dMB;
        mtp->dMB = nullptr;
    }
    if(mtp->DDc) {
        delete[] mtp->DDc;
        mtp->DDc = nullptr;
    }
    if(mtp->DIc) {
        delete[] mtp->DIc;
        mtp->DIc = nullptr;
    }
    if(mtp->DEl) {
        delete[] mtp->DEl;
        mtp->DEl = nullptr;
    }
    if(mtp->for_e) {
        delete[] mtp->for_e;
        mtp->for_e = nullptr;
    }
    //- if( mtp->xIC) delete[] mtp->xIC;
    //- if( mtp->xDC) delete[] mtp->xDC;
    //- if( mtp->xPH) delete[] mtp->xPH;
    if(mtp->grid) {
        delete[] mtp->grid;
        mtp->grid = nullptr;
    }
    if(mtp->NPmean) {
        delete[] mtp->NPmean;
        mtp->NPmean = nullptr;
    }
    if(mtp->nPmin) {
        delete[] mtp->nPmin;
        mtp->nPmin = nullptr;
    }
    if(mtp->nPmax) {
        delete[] mtp->nPmax;
        mtp->sdref = nullptr;
    }
    if(mtp->ParTD) {
        delete[] mtp->ParTD;
        mtp->ParTD = nullptr;
    }
    if(mtp->BM) {
        delete[] mtp->BM;
        mtp->BM = nullptr;
    }
    if(mtp->BdM) {
        delete[] mtp->BdM;
        mtp->BdM = nullptr;
    }
    if(mtp->FmgpJ) {
        delete[] mtp->FmgpJ;
        mtp->FmgpJ = nullptr;
    }
    if(mtp->BmgpM) {
        delete[] mtp->BmgpM;
        mtp->BmgpM = nullptr;
    }
    // work
    //- if( mtp->An) delete[] mtp->An;
    //- if( mtp->Ae) delete[] mtp->Ae;
    if(mtp->gfc) {
        delete[] mtp->gfc;
        mtp->gfc = nullptr;
    }
    if(mtp->yfb) {
        delete[] mtp->yfb;
        mtp->yfb = nullptr;
    }
    if(mtp->tt) {
        delete[] mtp->tt;
        mtp->tt = nullptr;
    }
    //- if( mtp->etext) delete[] mtp->etext;
    //- if( mtp->tprn) delete[] mtp->tprn;
    //- FreeNa();
    //- freeNodeWork();
}

// realloc dynamic memory
void TGEM2MT::mem_new(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
            "E04GTrem: Attempt to access corrupted dynamic memory.");

    //- mtp->xIC = new long int[mtp->nICb];
    //- mtp->xDC = new long int[mtp->nDCb];
    //- mtp->xPH = new long int[mtp->nPHb];

    if( mtp->PvGrid == S_OFF )
    { if(mtp->grid) delete[] mtp->grid;
        mtp->grid = 0;
    }
    else
        mtp->grid = new double[ mtp->nC][3];

    if(mtp->PsMode == RMT_MODE_W) {
        mtp->NPmean = new long int[ mtp->nPTypes];
        mtp->nPmin = new long int[ mtp->nPTypes];
        mtp->nPmax = new long int[ mtp->nPTypes];
        mtp->ParTD = new long int[mtp->nPTypes][6];
    }
    else {
        if(mtp->NPmean) delete[] mtp->NPmean;
        if(mtp->nPmin) delete[] mtp->nPmin;
        if(mtp->nPmax) delete[] mtp->nPmax;
        if(mtp->ParTD) delete[] mtp->ParTD;
        mtp->NPmean = nullptr;
        mtp->nPmin = nullptr;
        mtp->nPmax = nullptr;
        mtp->ParTD = nullptr;
    }

    mtp->nam_i= new char[mtp->nIV][ MAXIDNAME ];
    //- mtp->PTVm = new double[ mtp->nIV][5];
    if(!mtp->DiCp) { // could be allocated in constructor
        mtp->DiCp = new long int[ mtp->nC][2];
    }
    //- mtp->StaP = new double[ mtp->nC ][4];

    if(mtp->PvnVTK == S_OFF) {
        if(mtp->xVTKfld) {
            delete[] mtp->xVTKfld;
            mtp->xVTKfld = nullptr;
        }
    }
    else {
        if(!mtp->xVTKfld) { // could be allocated in function setVTKfields
            mtp->xVTKfld = new long int[mtp->nVTKfld][2];
        }
    }

    //- mtp->stld = new char[ mtp->nIV ][EQ_RKLEN];
    mtp->Tval  = new double[ mtp->nTai ];
    mtp->Pval  = new double[ mtp->nPai ];
    //- mtp->Bn = new double[ mtp->nIV][ mtp->Nb ];
    //- mtp->SBM = new char [ mtp->Nb][MAXICNAME+MAXSYMB];

    if(mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B) {
        mtp->HydP = new double[mtp->nC][SIZE_HYDP];
    }
    else {
        if(mtp->HydP) delete[] mtp->HydP;
        mtp->HydP = nullptr;
    }

    //-if( mtp->PvICi == S_OFF )
    //-   {
    //-    if(mtp->CIb) delete[] mtp->CIb;
    //-    if(mtp->CIclb) delete[] mtp->CIclb;
    //-    mtp->CIb = 0;
    //-    mtp->CIclb = 0;
    //-   }
    //-   else
    //-   {
    //-    mtp->CIb = new double[ mtp->nIV][mtp->Nb];
    //-    mtp->CIclb =  new char[ mtp->Nb ];
    //-   }

    //- if( mtp->PvAUi == S_OFF )
    //-    {
    //-     if(mtp->CAb) delete[] mtp->CAb;
    //-      if(mtp->for_i) delete[] mtp->for_i;
    //-      if(mtp->AUcln) delete[] mtp->AUcln;
    //-      if(mtp->An) delete[] mtp->An;
    //-      mtp->CAb = 0;
    //-      mtp->for_i = 0;
    //-      mtp->AUcln = 0;
    //-      mtp->An = 0;
    //-      mtp->Lbi = 0;
    //-    }
    //-    else
    //-    {
    //-      mtp->CAb = new double[ mtp->nIV][ mtp->Lbi ];
    //-      mtp->for_i = new char[ mtp->Lbi][ MAXFORMUNITDT ];
    //-      mtp->AUcln = new char[ mtp->Lbi ];
    //-      mtp->An = new double[ mtp->Lbi][ mtp->Nb ];
    //-   }

    if(mtp->PvFDL == S_OFF) {
        if(mtp->FDLi) delete[] mtp->FDLi;
        if(mtp->FDLf) delete[] mtp->FDLf;
        if(mtp->FDLid) delete[] mtp->FDLid;
        if(mtp->FDLop) delete[] mtp->FDLop;
        if(mtp->FDLmp) delete[] mtp->FDLmp;
        mtp->FDLi = nullptr;
        mtp->FDLf = nullptr;
        mtp->FDLid = nullptr;
        mtp->FDLop = nullptr;
        mtp->FDLmp = nullptr;
        mtp->nFD = 0;
    }
    else  {
        mtp->FDLi = new long int[ mtp->nFD][2];
        mtp->FDLf = new double[ mtp->nFD][4];
        mtp->FDLid= new char[ mtp->nFD][MAXSYMB];
        mtp->FDLop= new char[ mtp->nFD][MAXSYMB];
        mtp->FDLmp = new char[ mtp->nFD][MAXSYMB];
        for(long int ii=0; ii<mtp->nFD; ++ii) {
            fillValue(mtp->FDLid[ii], '\0', MAXSYMB);
            fillValue(mtp->FDLop[ii], '\0', MAXSYMB);
            fillValue(mtp->FDLmp[ii], '\0', MAXSYMB);
        }
    }

    if(mtp->PvPGD == S_OFF) {
        if(mtp->PGT) delete[] mtp->PGT;
        if(mtp->MGPid) delete[] mtp->MGPid;
        if(mtp->UMGP) delete[] mtp->UMGP;
        mtp->PGT = nullptr;
        mtp->MGPid = nullptr;
        mtp->UMGP = nullptr;
        mtp->nPG = 0;
    }
    else  {
        mtp->PGT  =  new double[ mtp->FIf*mtp->nPG ];
        mtp->MGPid = new char[ mtp->nPG ][MAXSYMB];
        mtp->UMGP = new char[ mtp->FIf ];
        for(long int ii=0; ii<mtp->nPG; ++ii) {
            fillValue(mtp->MGPid[ii], '\0', MAXSYMB);
        }
        for(long int ii=0; ii<mtp->FIf; ++ii) {
            mtp->UMGP[ii] = ' ';
        }
    }

    if( mtp->PvSFL == S_OFF )
    { if(mtp->BSF) delete[] mtp->BSF;
        mtp->BSF = 0;
    }
    else
        mtp->BSF = new double[ mtp->nSFD*mtp->Nf ];

    if(mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF) {
        mtp->MB =  new double[mtp->nC*mtp->Nf];
        mtp->dMB = new double[mtp->nC*mtp->Nf];
    }
    else {
        if(mtp->MB) delete[] mtp->MB;
        if(mtp->dMB) delete[] mtp->dMB;
        mtp->MB = nullptr;
        mtp->dMB = nullptr;
    }
    if( mtp->PvDDc == S_OFF )
    {
        if(mtp->DDc) delete[] mtp->DDc;
        mtp->DDc = 0;
    }
    else
        mtp->DDc = new double[mtp->Lsf];

    if( mtp->PvDIc == S_OFF )
    {
        if(mtp->DIc) delete[] mtp->DIc;
        mtp->DIc = 0;
    }
    else
        mtp->DIc = new double[ mtp->Nf ];

    if( mtp->nEl <= 0  )
    {
        if(mtp->DEl) delete[] mtp->DEl;
        if(mtp->for_e) delete[] mtp->for_e;
        //-   if(mtp->Ae) delete[] mtp->Ae;
        mtp->DEl = 0;
        mtp->for_e = 0;
        //-   mtp->Ae = 0;
        mtp->nEl = 0;
    }
    else
    {
        mtp->DEl = new double[ mtp->nEl ];
        mtp->for_e = new char[mtp->nEl][MAXFORMUNITDT];
        //-     mtp->Ae = new double[ mtp->nEl*mtp->Nb ];
    }

    //----------------------------------------------------------------
    //- if( mtp->Nqpt > 0  )
    //-  mtp->qpi   = new double[mtp->Nqpt];
    //- else
    //- { if(mtp->qpi) delete[] mtp->qpi; mtp->qpi = 0;}

    //- if( mtp->Nqpg > 0  )
    //-  mtp->qpc   = new double[mtp->Nqpg];
    //- else
    //- { if(mtp->qpc) delete[] mtp->qpc; mtp->qpc = 0;}

    //- if( mtp->PvMSt == S_OFF )
    //- { if(mtp->tExpr) delete[] mtp->tExpr; mtp->tExpr = 0;}
    //- else
    //-    mtp->tExpr = new char[4096];

    //-if( mtp->PvMSg == S_OFF )
    //-   {
    //-    if(mtp->lNam) delete[] mtp->lNam;
    //-    if(mtp->gExpr) delete[] mtp->gExpr;
    //-    if(mtp->xt) delete[] mtp->xt;
    //-    if(mtp->yt) delete[] mtp->yt;
    //-    mtp->lNam = 0;
    //-    mtp->gExpr = 0;
    //-    mtp->xt = 0;
    //-    mtp->yt = 0;
    //-   }
    //-   else
    //-   {
    //-        mtp->lNam = new char[ mtp->nYS][ MAXGRNAME];
    //-        mtp->gExpr = new char[2048];
    //-        mtp->xt   = new double[ mtp->nC];
    //-        mtp->yt   = new double[ mtp->nC*mtp->nYS];
    //-   }

    //-if( mtp->PvEF == S_OFF )
    //-   {
    //-    if(mtp->lNamE) delete[] mtp->lNamE;
    //-    if(mtp->xEt) delete[] mtp->xEt;
    //-    if(mtp->yEt) delete[] mtp->yEt;
    //-    mtp->lNamE = 0;
    //-    mtp->xEt = 0;
    //-    mtp->yEt = 0;
    //-   }
    //-   else
    //-   {
    //-     mtp->lNamE = new char[ mtp->nYE][ MAXGRNAME];
    //-     mtp->xEt   = new double[ mtp->nE ];
    //-     mtp->yEt   = new double[ mtp->nE*mtp->nYE ];
    //-   }

    if( mtp->Nsd > 0 )
    {
        mtp->sdref = new char[ mtp->Nsd ][ V_SD_RKLEN ];
        mtp->sdval = new char[ mtp->Nsd ][ V_SD_VALEN ];
    }
    else
    {
        if(mtp->sdref) delete[] mtp->sdref;
        if(mtp->sdval) delete[] mtp->sdval;
        mtp->sdref = 0;
        mtp->sdval = 0;
    }
    //- mtp->etext = new char[4096];
    //- mtp->tprn = new char[2048];
    //mtp->gfc = (double *)aObj[ o_mtgfc].Free();
    //mtp->yfb = (double *)aObj[ o_mtyfb].Free();
    //mtp->tt = (double *)aObj[ o_mttt].Free();
}

//=============================================================

// Conversion of concentration units to moles
double TGEM2MT::Reduce_Conc( char UNITP, double Xe, double DCmw, double Vm,
                            double R1, double Msys, double Mwat, double Vaq, double Maq, double Vsys )
{
    double Xincr = 0.;
    switch( UNITP )
    {  // Quantities
    case QUAN_MKMOL: /*'Y'*/
        Xincr = Xe / 1e6;
        goto FINISH;
    case QUAN_MMOL:  /*'h'*/
        Xincr = Xe / 1e3;
        goto FINISH;
    case QUAN_MOL:   /*'M'*/
        Xincr = Xe;
        goto FINISH;
    }
    if( DCmw > 1e-12 )
        switch( UNITP )
        {
        case QUAN_MGRAM: /*'y'*/
            Xincr = Xe / DCmw / 1e3;
            goto FINISH;
        case QUAN_GRAM:  /*'g'*/
            Xincr = Xe / DCmw;
            goto FINISH;
        case QUAN_KILO:  /*'G'*/
            Xincr = Xe * 1e3 / DCmw;
            goto FINISH;
        }
    /* Concentrations */
    if( fabs( R1 ) > 1e-12 )
        switch( UNITP )
        { // mole fractions relative to total moles in the system
        case CON_MOLFR:  /*'n'*/
            Xincr = Xe * R1;
            goto FINISH;
        case CON_MOLPROC:/*'N'*/
            Xincr = Xe / 100. * R1;
            goto FINISH;
        case CON_pMOLFR: /*'f'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow(10., -Xe )* R1;
            goto FINISH;
        }
    if( fabs( Vsys ) > 1e-12 && Vm > 1e-12 )
        switch( UNITP )   /* Volumes */
        {
        case CON_VOLFR:  /*'v'*/
            Xincr = Xe * Vsys * 1e3 / Vm;
            goto FINISH;
        case CON_VOLPROC:/*'V'*/
            Xincr = Xe * Vsys * 10. / Vm;
            goto FINISH;
        case CON_pVOLFR: /*'u'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vsys / Vm;
            goto FINISH;
        }
    if( fabs( Msys ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP ) // Mass fractions relative to mass of the system
        {
        case CON_WTFR:   /*'w'*/
            Xincr = Xe * Msys * 1e3 / DCmw;
            goto FINISH;
        case CON_WTPROC: /*'%'*/
            Xincr = Xe * Msys *10. / DCmw;
            goto FINISH;
        case CON_PPM:    /*'P'*/
            Xincr = Xe * Msys / 1e3 / DCmw;
            goto FINISH;
        }
    if( fabs( Mwat ) > 1e-12 )
        switch( UNITP ) /* Molalities */
        {
        case CON_MOLAL:  /*'m'*/
            Xincr = Xe * Mwat;
            goto FINISH;
        case CON_MMOLAL: /*'i'*/
            Xincr = Xe / 1e3 * Mwat;
            goto FINISH;
        case CON_pMOLAL: /*'p'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Mwat;
            goto FINISH;
        }
    if( fabs( Vaq ) > 1e-12 )
        switch( UNITP )  /* Molarities */
        {
        case CON_MOLAR:  /*'L'*/
            Xincr = Xe * Vaq;
            goto FINISH;
        case CON_MMOLAR: /*'j'*/
            Xincr = Xe * Vaq / 1e3;
            goto FINISH;
        case CON_pMOLAR: /*'q'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vaq;
            goto FINISH;
            /* g/l, mg/l, mkg/l */
        case CON_AQGPL: /*'d'*/
            Xincr = Xe * Vaq / DCmw;
            goto FINISH;
        case CON_AQMGPL: /*'e'*/
            Xincr = Xe * Vaq / DCmw / 1e3;
            goto FINISH;
        case CON_AQMKGPL: /*'b'*/
            Xincr = Xe * Vaq / DCmw / 1e6;
            goto FINISH;
        }
    if( fabs( Maq ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP )     /* Weight concentrations */
        {
        case CON_AQWFR:  /*'C'*/
            Xincr = Xe * Maq * 1e3 / DCmw;
            goto FINISH;
        case CON_AQWPROC:/*'c'*/
            Xincr = Xe * 10. * Maq / DCmw;
            goto FINISH;
        case CON_AQPPM:  /*'a'*/
            Xincr = Xe * Maq / 1e3 / DCmw;
            goto FINISH;
        }
    /* Error */
FINISH:
    return Xincr;
}

//---------------------------------------------------------------------------

