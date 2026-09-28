###  1. Task definition

DiCp Distributing must be setup before reading gems3k files and allocation.

```
/// Constructor defines type of mas transport and number of nodes;
explicit TGEM2MT(char ps_mode, long int  n_nodes);
/// nC:  Input number of local equilibrium cells (nodes)
long int nNodes() const;

/// Get the full name of this GEM2MT task
std::string name() const;
/// Set the full name of this GEM2MT task
void setName(const std::string& task_name);
/// Get the comment of this GEM2MT task
std::string comment() const;
/// Set the comment of this GEM2MT task
void setComment(const std::string& task_notes);
   
/// DiCp:  Change array of indexes of initial system variants for distributing to nodes [nC]
void setDistributing(long int node_ndx, long int sys_ndx);

/// DiCp: The second column DiCp[1] contains the node type for each node:
/// 0:   normal node;
/// Boundary condition nodes:
/// 1:   Dirichlet source (constant composition source);
/// -1:  Dirichlet sink;
/// 2:   Neumann source (constant gradient source);
/// -2:  Neumann sink;
/// 3:   Cauchy source (constant flux source);
/// -3:  Cauchy sink;
/// 4:   Input time-depended function (TBD).
void setNodeType(long int node_ndx, long int type);
   
if(!mtp->DiCp) { // could be allocated in constructor
     mtp->DiCp = new long int[ mtp->nC][2];
}
void defaults_DiCp();

```
   
### 2. Controls on operation

Type flux phase is used for default initialization MGPid, PGT, UMGP, FDLmp, FDLid and must be defined before allocation. 

```
/// PvMSt,    // ? Use math script for start setup (+ -)?    callback to update internal
/// PvMSg,    // ? Use math script for graphic presentation (+ -)?  callback to collect graphic data
/// PvMSc,    // ?  Use math script for control on time steps (+ -)?  CalcControlScript
   
/// PsSIA: Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
void setSIA(char flag);
/// PsMO: Use non stop debug output for nodes (+ -) (default +)
void setOutput(bool enable);
/// PsVTK: Use non stop debug output nodes to VTK format(+ -) (default -)
void setOutVTK(bool enable);
/// PsMPh: Type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
void setTypeFluxPhase(char flag);

```

### 3. Initial scalars and iterators 

Some of scalars are used for default initialization and better define them before allocation.

```
/// Tau:   Physical time iterator (start,end,step)
void setTau(double start, double end, double step);
/// sizeLc:  Spatial dimensions of the medium defines topology of nodes ( x y z )
void setSpatialDimensions(double x, double y, double z);

/// tf:  Advection/diffusion mass transport: time step reduction factor (usually 1)
void setTimeStepReductionFactor(double val);
/// Vt:  Initial total node volume (m^3)
void setInitialTotalNodeVolume(double val);
/// vp:  Fluid advection velocity (m/sec)
void setFluidAdvectionVelocity(double val);
/// eps:  Initial node effective porosity (0 < eps < 1), usually 1
void setInitialNodeEffectivePorosity(double val);
/// Km:  Initial effective permeability, m2, usually 1
void setInitialEffectivePermeability(double val);
/// al:  Initial value of specific longitudinal dispersivity (m), usually 1e-3
void setInitialDispersiviSome of them used for default array initialization, better define before allocationty(double val);
/// Dif:  Initial general aqueous medium diffusivity (m2/sec), usually 1e-9
void setInitialDiffusivity(double val);
/// nto:  Initial tortuosity factor, usually 1
void setInitialTortuosityFactor(double val);
/// cdv:   Cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
void setCutofffICamount(double val);
/// cez:   Cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
void setCutoffMinimalAmountsIC(double val);

```

### 4. gems3k sizes

Some sizes and arrays are defined after reading gem3k files. 

```
mtp->Lsf = na->pCSD()->nDCs;
mtp->Nf = na->pCSD()->nIC;
mtp->FIf = na->pCSD()->nPH;
mtp->nTai = na->pCSD()->nTp;
mtp->nPai = na->pCSD()->nPp;
mtp->Tval;
mtp->Pval;

/// FIf:  Number of phases in (DATABR) for setting box-fluxes
long int nPhases() const;
/// Nf:  Number of ICs in (DATABR) for setting box-fluxes
long int nElements() const;
/// Lsf: of DCs in phases-solutions in Multi (DATACH) for setting box-fluxes
long int nComponents() const;

```


### 5. Allocation and setup flags, sizes

Sizes must be defined before allocation, flags sets when size more then 0;

1. Use export to VTK format 

```
    PvnVTK != S_OFF   // Use selected fields to VTK format (+ -)
    PsVTK             // Use non stop output from nodes to VTK format files (+ -)
    long int nVTKfld; // Number of selected fields to VTK format  (by default 0)
    /// xFlds: Set list of selected fields and indexes to VTK format
    void setVTKfields(const std::vector<std::pair<int, int>>& vtk_fields);
    /// PsVTK: Use non stop debug output nodes to VTK format(+ -) (default -)
    void setOutVTK(bool enable);
```
   

2. If PsMode not in { `S`, `F`, `B` } initial hydraulic parameters are used.
   
```
   mtp->HydP = new double[mtp->nC][SIZE_HYDP];
   /// HydP:  Initial hydraulic parameters in nodes: Vt, vp, eps, Km, al, Dif,  nto
   /// @param pndx: index in array
   /// @param Vt: initial total volume of the node, m3 (for porosity)
   /// @param vp: initial advection velocity, m/s
   /// @param eps: initial effective porosity
   /// @param Km: initial effective permeability
   /// @param al: initial specific longitudinal dispersivity
   /// @param Dif: initial general diffusivity
   /// @param nto: initial tortuosity factor
   void setHydraulicParameters(long int pndx, double Vt, double vp, double eps, double Km, double al, double Dif, double nto);
   
   // Default initialization HydP
   void defaults_HydP();
```

3. PvFDL != S_OFF  // Use flux definition list (+ -)

```
   long int mtp->nFD; // total number of MGP flux definitions, incl elemental source fluxes (by defaults 0)
   mtp->FDLi = new long int[ mtp->nFD][2];
   mtp->FDLf = new double[ mtp->nFD][4];
   mtp->FDLid= new char[ mtp->nFD][MAXSYMB];
   mtp->FDLop= new char[ mtp->nFD][MAXSYMB];  // only allocated in code
   mtp->FDLmp = new char[ mtp->nFD][MAXSYMB];
   /// nFD: Number of MGP fluxes defined in the megasystem, nFD >= 0
   long int nMGPfluxes() const;
   /// nFD: Number of MGP fluxes defined in the megasystem, nFD >= 0
   void setNumberMGPfluxes(long int num);
   /// FDLi: Set Source/Receive box index in the flux definition
   void setFluxSourceReceive(long int  pndx, long int  source, long int  receive);
   /// FDLf: Set the flux defnition: flux order, flux rate, MGP quantity
   void setFluxSourceReceive(long int  pndx, double order, double rate, double quantity, double val);
   /// FDLmp: [nFD] ID of MGP to move in this flux
   void setFluxMGPid(long int  pndx, const std::string& ids);
   /// FDLid: Set IDs of fluxes
   void setFluxIDs(long int  pndx, const std::string& ids);
   
   // Default initialization FDLf, FDLi
   void defaults_FDLi_FDLf();
   // Default initialization MGPid, PGT, UMGP, FDLmp, FDLid
   void defaults_MGPid_PGT_FDLmp_FDLid(bool mode);
```

4. PvPGD != S_OFF  // Use phase groups definitions (+ -)

```
   long int mtp->nPG // number of mobile phase groups MGP (by defaults 0)
      mtp->PGT  =  new double[ mtp->FIf*mtp->nPG ];
      mtp->MGPid = new char[ mtp->nPG ][MAXSYMB];
      mtp->UMGP = new char[ mtp->FIf ];
   /// nMGP: Number of mobile groups of phases, nMGP >= 0
   long int nPhaseGroups() const;
   /// Define the number of mobile groups of phases, nMGP >= 0
   void setNumberPhaseGroups(long int num);
   /// MGPid: ID list of mobile phase groups
   void setPhaseGroupsID(long int  pndx, const std::string& ids);
   /// UMGP: [nFi] units for setting phase quantities in MGP (see PGT )
   void setUnitsPhaseQuantities(long int  gndx, char units);
   /// PGT: Quantities of phases in MGP [Fi][nPG]
   /// @param gndx: phase groups index
   /// @param pndx: phase index
   void setPhaseGroupsQuantities(long int gndx, long int pndx, double quantity);
  
   // Default initialization MGPid, PGT, UMGP, FDLmp, FDLid
   void defaults_MGPid_PGT_FDLmp_FDLid(bool mode);
```

5. PvSFL != S_OFF  // Use source fluxes and elemental stoichiometries for them (+ -)

We have not example of using this mode.  (To be done)

```
   long int nSFD;  // number of elemental source flux definitions (by defaults 0)
     mtp->BSF = new double[ mtp->nSFD*mtp->Nf ];
   /// nSFD:  Number of IC source flux compositions defined in megasystem, nSFD >= 0
   long int nICsourceFluxes() const;
   /// Define the number of elemental source flux definitions, nSFD >= 0
   void setNumberICsourceFluxes(long int num)
   /// BSF: table of bulk compositions of elemental fluxes [nSFD][Nf]
   /// @param gndx: groups index
   /// @param indx: element index
   void setICsourceQuantities(long int gndx, long int indx, double quantity);
   
   // Default initialization BSF (to be done; temporally all 0)
   void defaults_BSF();
```

6.  mtp->PsMode == `W`

```
  long int mtp->nPTypes; // Number of allocated particle types (by default 10)
    mtp->NPmean = new long int[ mtp->nPTypes];
    mtp->nPmin = new long int[ mtp->nPTypes];
    mtp->nPmax = new long int[ mtp->nPTypes];
    mtp->ParTD = new long int[mtp->nPTypes][6];
  /// nPTypes:  Number of allocated particle types < 20
  long int nParticleTypes();
  /// Define number of allocated particle types (<20)
  void setNumberParticleTypes(long int num);
  /// Set of particle statistic property [nPTypes]
  /// @param pndx: index in array
  /// @param NPmean: Array of initial mean particle type numbers per node
  /// @param nPmin: Minimum average total number of particles of each type per one node
  /// @param nPmax: Maximum average total number of particles of each type per one node
  /// @param ParTD: Array of particle type definitions at t0 or after interruption
  void setParticle(long int pndx, long int pmean, long int pmin, long int pmax, const std::array<long int, 6>& pparam);
  
  // Default  particle array setup
  void defaults_particle_setup();
```

7. PvGrid != S_OFF // Use array of grid point locations

We have not example using grid, in examples used default grid, generated from sizeLc.

```
   mtp->grid = new double[ mtp->nC][3];
   /// sizeLc:  Spatial dimensions of the medium defines topology of nodes ( x y z )
   void setSpatialDimensions(double x, double y, double z)
   /// PvGrid: Use array of grid point locations (+ -) (default -)
   void useArrayGridPoints(bool enable);
   
   /// Set grid point location, size is nC [grid]
   /// @param pndx: index in array
   /// @param x: Array of initial mean particle type numbers per node
   /// @param y: Minimum average total number of particles of each type per one node
   /// @param z: Maximum average total number of particles of each type per one node
   void setGrid(long int pndx, double x,  double y,  double z);
   
   // Set default grid coordinate array use predefined sizeLc
   void defaults_Grid();
```
  
8. mtp->PvDDc != S_OFF

Used in `Trans1D` for different transport function, the array is allocated but not used in mas transport calculations

```
   mtp->DDc = new double[mtp->Lsf]; // only allocated in code
   /// PvDDc: Use diffusion coefficients for DC - DDc vector (+ -) (default -)
   void usePvDDc(bool enable);
```

9. mtp->PvDIc != S_OFF 

Used in `Trans1D` for different transport function, the array is allocated but not used in mas transport calculations

```
   mtp->DIc = new double[ mtp->Nf ]; // only allocated in code
   /// PvDIc: Use diffusion coefficients for IC - DIc vector (+ -) (default -)
   void usePvDIc(bool enable);
```

10. Arrays are not used in mas transport calculations

```
   mtp->nEl > 0
   mtp->nEl  // only for allocation DEl and for_e
   mtp->DEl = new double[ mtp->nEl ];  // only allocated in code
   mtp->for_e = new char[mtp->nEl][MAXFORMUNITDT]; // only allocated in code
```





 
  
