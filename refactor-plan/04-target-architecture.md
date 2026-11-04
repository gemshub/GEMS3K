GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 04 — Target architecture

Written after `09-depth-check.md`. The target keeps every contract in `02-api-contract.md`. It is reached in
small steps (05 for Stage A, 06 for Stage B). Proposed comments and documentation text use ASD-STE100 style
(short sentences, one instruction per sentence, approved verbs).

## 1. Design goals (in priority order)

1. **Same answers.** Stage A changes no floating-point operation order (brief §4). Every A step is checked by
   STD and RAW freezes row for row, including `# dec`.
2. **Explicit state.** Every value that one call passes to another goes through a parameter, a return value or a
   named object. No thread-local or static mutable state in the solver (today 00 §4 lists H1–H14).
3. **One definition per fact.** One declaration of each `BASE_PARAM` field; one copy of each dense kernel; one
   status-mapping function; one snapshot helper; one table per character code.
4. **Modules with one job each**, visible in folders and CMake targets, so that a new activity model, a new solver
   leg or a GPU batch path plugs in at one place.
5. **Diagnostics unchanged** (RUN/SET/EFF/KEY/CERT/DECIDE, probes, spdlog channels), but written by one component.

## 2. Modules, folders and CMake targets

```
GEMS3K/                         (repository root: unchanged top level)
└─ GEMS3K/                      (installed include dir name stays "GEMS3K")
   ├─ core/        gems3k_core      errors (TError, Error/ErrorIf), logging set-up, constants (m_const_base),
   │                                version, string helpers (v_service, v_detail), cancel flag
   ├─ data/        gems3k_data      DATACH, DATABR, MULTI (struct kept), BASE_PARAM + base_param_fields.def,
   │                                settings resolvers (optima_dimreduce_passes …), status codes, status mapping
   ├─ io/          gems3k_io        io_template/keyvalue/simdjson/nlohmann policies, DCH/DBR/IPM formats,
   │                                GemDataStream, .lst generator (gems3k_impex), GemsSettings (jsonconfig)
   ├─ thermo/      gems3k_thermo    ThermoSource interface: ThermoFunSource (target), GridSource (legacy until
   │                                B14), CompareSource (10 §2)
   ├─ models/      gems3k_models    TSolMod base + factory; fluid/ (CG, CORK, PR78, PRSV, SRK, STP), solution/
   │                                (Ideal, Margules, …, Berman, CEF, MBW), aqueous/ (DH, Davies, Helgeson, Karpov,
   │                                Shvarov, ELVIS, Pitzer, SIT, EUNIQUAC), other/; kinetics/ (TKinMet); sorption/
   │                                → becomes the ThermoSol library (phase models, from TActivity) and the
   │                                ThermoKin library (kinetics, from TKinetics): thermosol/, thermokin/ (11)
   ├─ numerics/    gems3k_numerics  dense kernels used by the solver: Householder QR (IPM), dense LU with
   │                                refinement, Jacobi symmetric eigen floor, two-phase LP simplex, Karpov simplex,
   │                                interpolation; vendored tnt/jama (internal)
   ├─ kernel/      gems3k_kernel    GEMS3K IPM solver: chemistry evaluation (concentrations, chemical potentials,
   │                                GX, TotalGibbsEnergy), MBR, InteriorPointsMethod, PSSC, StabilityIndexes,
   │                                fill-out, certification, determinacy, leveled start, kinetics step, scaling
   ├─ optima/      gems3k_optima    (USE_OPTIMA_SOLVER) problem build, objective/Hessian, pre-solve, seeds,
   │                                retry/repair tiers as named functions, HOP, native finish/rescue, targets
   ├─ run/         (in gems3k_node) GEM_run layer pipeline: rescue, cap, dry redo, tries, single; RunContext
   ├─ node/        gems3k_node      TNode public API, pack/unpack, name maps, TActivity/TKinetics front ends (11)
   ├─ nodearray/   gems3k_nodearray TNodeArray
   ├─ system/      gems3k_system    SystemDefinition, SystemBuilder, system.yaml reader/writer, derivations
   │                                (stoichiometry, codes, model layout, bulk) — systems without GEMSGUI (13)
   └─ diag/        gems3k_diag      trace writers (RUN/SET/EFF/KEY/CERT/DECIDE), probes, gems3k_timing
```

CMake: each folder is an `OBJECT` library with explicit source lists (no `GLOB`), and with
`target_link_libraries` that state the allowed dependencies (§4). The installed artefacts stay the same:
`GEMS3K::GEMS3K` (shared) and `GEMS3K::GEMS3K-static`, built from the object libraries, with headers installed
flat under `include/GEMS3K/` as today (the source folders are an internal layout; `install(FILES …)` lists the
public headers of 02 §5). Compile flags, definitions and the order of source files inside one translation unit
are kept in Stage A (A-rule 2, 05 §0). The qmake files are generated from the same lists, or dropped (08 Q4).

## 3. Key interfaces (new code in `namespace gems3k`)

### 3.1 `RunContext` — replaces the thread-local flags (H1, H2)

```cpp
// RunContext holds the state of one TNode::GEM_run() call.
// Each layer gets the context as a parameter. No layer keeps a copy.
struct RunContext {
    long requestedMode;                 // NodeStatusCH at entry
    IterationBudget optimaBudget;       // was g_optimaIterBudget, g_optimaCapAllowance
    ReentryGuards guards;               // was g_optimaCapCall, g_optimaDryRedoCall
    AttemptSwitches attempt;            // was g_optimaFixAccept, g_optimaLSEscapeOff, g_optimaLineSearchRetry,
                                        //     g_optimaEscRetryNext, g_optimaCgSeedArmed, g_ipmLevelOff, g_optimaDryOff
    AttemptMemory memory;               // was g_optimaPreSolveEmpty, g_lvlFullCalls, g_optimaTargetWarned,
                                        //     g_optimaDryFailed, g_ipmNoProgArmed/Fired, g_aqueousAbsent
    QuickTryWatch quickTry;             // was g_optimaQuickTry
    DiagnosticsSink* diag;              // trace writer of this call (thread owns it)
};
```

`BeginCalculation()` becomes `RunContext::beginAttempt()` and keeps today's scope (per attempt) in Stage A.
Helper outputs (`g_cgConverged`, `g_cgPseudo`, `g_cgGap`, `g_cgDual`, `g_coldLvlUsed`, `g_lvlShortcut`,
`g_lvlCgStart`, `g_ps0GaveUp`) become return structs (C++ Core Guidelines F.21).

### 3.2 Layer pipeline

```cpp
// A layer is one retry policy. It calls next() to run the inner layers.
long runRescue (TNode&, RunContext&, Next next);   // pa_FailRescue
long runCapped (TNode&, RunContext&, Next next);   // pa_OptimaCallCap
long runDryRedo(TNode&, RunContext&, Next next);   // pa_OptimaDryRetry
long runTries  (TNode&, RunContext&, Next next);   // quick try, retries, side entries
long runSingle (TNode&, RunContext&);              // unpack, dispatch, pack
```

The order is fixed in one place (`run/pipeline.cpp`) and printed in the EFF line. Each re-entry that today goes
through `GEM_run()` again is listed in 01 §2 and kept with the same path (brief §3.7).

### 3.3 Solver entry points

```cpp
struct SolveOutcome { long itFIA, itIPM; double seconds; bool soft; /* MK/PZ */ };
SolveOutcome solveNative(TMultiBase&, RunContext&, bool runKinetics);          // CalculateEquilibriumState
SolveOutcome solveOptima(TMultiBase&, RunContext&, OptimaLeg leg);             // AOP/SOP/ROP, HOP leg, finish
SolveOutcome solveHybrid(TMultiBase&, RunContext&, bool warmNative);           // HOP/SHP
```

`OptimaLeg { Plain, Reference, HopLeg, NativeFinish, NativeRescue }` replaces the member flag
`optima_hop_leg` and the `referenceMode`/`runKinetics` booleans (H11).

### 3.4 Optima adapter split

`CalculateEquilibriumStateOptima` (4 400 lines) becomes a driver of ~15 named functions, one per tier, in the same
order, each taking `OptimaRun&` (problem, state, options, watch, flags) and returning a `TierResult`. The objective
lambda is built once by `makeObjective(view, cfg)` and shared by the pre-solve and the full solve (today two copies,
NOTES ipm_optima part 2). Extraction keeps the operation order: same statements, same order, same compiler
inlining scope (A-rule 3).

### 3.5 Views on MULTI (ownership stays, access becomes typed)

`MULTI` stays (GEMSGUI reads it). New code reads it through small views:
`SystemView` (const: N, L, L1, A, codes, names), `AmountsView` (Y, X, XF, XFA, B), `PotentialsView` (G0, G, F,
F0, lnGam, U), `BoundsView` (DLL, DUL, PLL, PUL). Views are non-owning (`std::span`), built per call. They make
the read/write sets of 09 §4 visible in function signatures (I.4).

### 3.6 Workspace (H3, allocation)

`SolverWorkspace` owned by `TMultiBase`: the scratch now held in function-static thread-locals
(`SolveIpmAugmented` vectors, Eigen matrices), in `Build_compressed_xAN` (rebuilt every call) and in per-call
`new[]` (simplex NEW/WORK, TNT Array2D). In Stage A it only moves the storage (same values, same order); in
Stage B it allows reuse across calls.

### 3.7 Diagnostics

`DiagnosticsSink` (one per thread, created by `RunContext`) owns the trace `FILE*` handles and the warn-once sets
(H3, H4). The record writers move to `diag/` unchanged; `native_trace_decide()` stays as a forwarding function so
the ~150 call sites do not change in Stage A. Environment variables that change results
(`GEMS3K_OPTIMA_NET_RESUME`, `GEMS3K_FILLOUT_BUDGET`) are read once into `RunConfig` and printed in EFF (they
already are for the first); names and meaning stay.

### 3.8 Settings

`BASE_PARAM` gets one declaration file `base_param_fields.def` (X-macro, 08 Q3):
`GEMS3K_PA_FIELD(type, name, default, tag, comment, flags)`. It generates the struct members (same order), the
`MULTI_DYNAMIC_FIELDS` enum, the descriptor rows, `writeField`, the read `case`, both counts
(`std::size(table)`), the SET line and the defaults of `pa_p_` for the members after the positional prefix. The
first 38 members stay hand-written (GEMSGUI binary layout) with `static_assert(offsetof(...))`. Access:
`const BASE_PARAM& settings()` for readers; `LocalSettings` stays the only way to write during a call.

### 3.9 Phase models (preparation for 07)

`PhaseModel` = interface around one `TSolMod` with an explicit input/output block instead of aliased MULTI
pointers: inputs x, T, P, parameters; outputs lnγ (and optional dlnγ/dn). In Stage A it is a thin wrapper that still
writes through today's aliases (no FP change). The factory `SolModCreate` becomes a registration table
(one entry per model: code, constructor, name), so adding a model touches one file (today ≥ 8 sites).
The full design — ThermoSol and ThermoKin as libraries built on ThermoFun v2 principles, with `TActivity` and
`TKinetics` as the starting point — is in `11-thermosol-thermokin.md`; standard states from ThermoFun only in
`10-thermofun-tp.md`; the performance techniques behind 06 in `12-performance-techniques.md`.

## 4. Dependency rules (checked by CMake target links)

```
core ← data ← io
core ← numerics
ThermoFun ← thermo (ThermoSource)
ThermoFun, numerics ← ThermoSol (phase models), ThermoKin (kinetics)     no GEMS3K headers inside them
data, numerics, thermo, ThermoSol, ThermoKin ← kernel
kernel ← optima          (optional; Optima::Optima)
kernel, optima, diag ← node(run) ← nodearray
data, io, thermo, ThermoSol, ThermoKin, ChemicalFun ← system ← node (GEM_init(const SystemDefinition&))
core ← diag ← (kernel, optima, node)     diag only writes; it never changes results
```

No module may include a header of a module to its right. `kernel` does not include `node.h` (today
`TMultiBase::node1` is used for names, ThermoFun and T/P: these become `SystemView` fields and a `ThermoSource`
interface).

## 5. Ownership rules

- `TNode` owns `DATACH`, `DATABR`, `TMultiBase` (unique ownership; `std::unique_ptr` in new code, R.20).
- `TMultiBase` owns the MULTI arrays. Allocation stays in the virtual `alloc_*`/`multi_realloc` (GEMSGUI allocates
  through them); the `new[]`/`delete[]` pairs stay in Stage A.
- Settings: one `BASE_PARAM` per node set, shared read-only; a call that changes settings uses `LocalSettings`.
  Sharing is visible (`std::shared_ptr<const BASE_PARAM>` in B).
- Phase models are owned by `TMultiBase` and **do not own** their inputs; in Stage B they stop aliasing MULTI.
- Optima `Problem`/`State` are owned by one `OptimaRun`; the objective captures `OptimaRun&` only (not
  `&fixedGrad` from a dead scope).
- Warm-start state (U, Y, prev_scfact, floor seeds, prevBulkMol_, warmInputChanged_) becomes one `WarmState`
  object on `TMultiBase`; `TNodeArray` work nodes copy it explicitly (A11).

## 6. Naming

- Keep public names (02). Keep the fixed terms "GEMS3K IPM solver", "warm start", "cold start".
- New types PascalCase, new functions camelCase, files snake_case, namespace `gems3k` for new code.
- DECIDE record names, env variable names and logger names never change.
- A file is named after what it holds (today `ms_multi_diff.cpp` holds thermodynamic loading and
  `ipm_simplex.cpp` holds the native solve entry).

## 7. Guidance: adopted and rejected

Rule titles of the C++ Core Guidelines checked against the official source
(github.com/isocpp/CppCoreGuidelines, `CppCoreGuidelines.md`; isocpp.github.io was blocked by this session's
network policy).

| Source | Item | Decision | Where used |
|---|---|---|---|
| [C++ Core Guidelines](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md) | I.2 Avoid non-const global variables; I.3 Avoid singletons; CP.2 Avoid data races; CP.3 Minimize explicit sharing of writable data; CP.1 Assume multi-threaded | **Adopted** | RunContext, DiagnosticsSink, model statics (A6–A9) |
| same | F.21 Return a struct for multiple "out" values; F.20 prefer return values | **Adopted** | helper outputs (g_cg*) |
| same | I.4 Make interfaces precisely and strongly typed; Enum.3 Prefer class enums | **Adopted** for new code (OptimaLeg, views); char-coded enums kept for files | 3.3, 3.5 |
| same | P.11 Encapsulate messy constructs | **Adopted** | snapshot helper, status mapping, X-macro |
| same | R.1 RAII; R.3 raw pointer is non-owning; R.11 avoid new/delete; R.20 unique_ptr/shared_ptr for ownership | **Adopted for new code; deferred** for MULTI arrays (GEMSGUI allocates through `alloc_*`) | 5 |
| same | C.41 constructor creates a fully initialized object; NR.5 no two-phase initialization | **Rejected for Stage A** (MULTI is allocated in two phases by design: read dims, then `multi_realloc`); revisit with GEMSGUI | — |
| same | C.131 Avoid trivial getters and setters | **Rejected** for TNode (public API, 02) | — |
| same | F.8 Prefer pure functions; Con.2 const member functions | **Adopted** for getters (A5: pure `Get_GibbsEnergy` variant) | A5 |
| same | ES.45 Avoid magic constants | **Adopted** gradually: named constants for thresholds that appear twice or more (1e-12 index cutoff, 1e6 bound sentinel) | B |
| same | Per.1/Per.2/Per.6 no optimization without reason or measurement; Per.7 design to enable optimization; Per.16 compact data; Per.19 access memory predictably | **Adopted**: Stage B steps cite 03 measurements | 06 |
| same | I.27 Pimpl for stable ABI | **Rejected for now**: GEMSGUI subclasses TMultiBase and reads MULTI; Pimpl would break it. Revisit for TNode after 5.0 | — |
| [UK AQuA Book](https://www.gov.uk/guidance/the-aqua-book) and [Duck Book (QA of code for analysis and research)](https://best-practice-and-impact.github.io/qa-of-code-guidance/intro.html), [checklists](https://best-practice-and-impact.github.io/qa-of-code-guidance/checklists.html) | proportionality; verification and validation; reproducible, peer-reviewed changes; version control; automated tests | **Adopted**: each step has a stated verification; freeze pairs are the V&V record; proportionate effort (a pure move needs freeze + instruction count, a numerical change needs the scoring of 06) | 05 §0, 06 §0 |
| [xSDK community policies](https://xsdk.info/policies) | CMake options, 64-bit build, debug build, sustainable compatibility changes, installation | **Adopted in part**: CMake targets per module, no hard-coded paths, debug build option; **rejected**: Spack packaging requirement (owner's conda workflow) | A1 |
| [BSSw PSIP](https://bssw-psip.github.io), [PSIP tools](https://betterscientificsoftware.github.io/PSIP-Tools) | small iterative improvement goals with tracking cards | **Adopted** as format: each 05/06 step is a card (what, why, risk, check, size, deps) | 05, 06 |
| Wilson et al., [Best Practices for Scientific Computing (2014)](https://doi.org/10.1371/journal.pbio.1001745); [Good enough practices (2017)](https://doi.org/10.1371/journal.pcbi.1005510) | make incremental changes; do not repeat yourself; plan for mistakes (assertions, tests); optimize only after it works | **Adopted** | A2, A3, A12, A14 |
| de-RSE [Foundational competencies](https://arxiv.org/abs/2311.11457) and [GI/de-RSE guideline](https://www.se-rwth.de/publications/GI-and-DE-RSE-template-guideline-for-the-efficient-development-of-research-software.pdf); US-RSE (no guideline found; US-RSE24 cites the German paper) | quality management proportional to the software's application class; testing and documentation | **Adopted** as framing for 08 (the owner decides the class: GEMS3K is a long-lived library used by others → highest class) | 08 |
| [Agner Fog, optimization manuals](https://www.agner.org/optimize/) | keep data used together close; sequential access; vectorization needs layout | **Adopted for Stage B** | B4–B6 |
| Data-oriented design ([Acton, CppCon 2014](https://isocpp.org/blog/2015/01/cppcon-2014-data-oriented-design-and-c-mike-acton)) | structure of arrays; design around the data transforms | **Adopted in part**: MULTI is already SoA; batching (07) uses SoA over nodes; **rejected**: ECS-style rewrite | 07 |
| [Kokkos Kernels batched dense](https://kokkos.org/kokkos-kernels/docs/API/batched/dense-index.html), SYCL ([batched iterative on Intel GPUs](https://ar5iv.labs.arxiv.org/html/2308.08417)) | performance-portable batched small solves | **Deferred to 07** (candidate backend for node batches); not a Stage A/B dependency | 07 |
| std::execution ([P2300, in C++26 WD since 2024-07](https://herbsutter.com/2024/07/02/trip-report-summer-iso-c-standards-meeting-st-louis-mo-usa/)), std::simd ([P1928](https://isocpp.org/files/papers/P1928R15.pdf), C++26 target, adoption not confirmed by search) | standard parallel and SIMD | **Rejected for now**: consumers build with C++20 on MSVC/Apple clang; revisit when all three compilers ship them | 07 |
| Reaktoro ([EquilibriumSensitivity](https://reaktoro.org/api/classReaktoro_1_1EquilibriumSensitivity.html); local `../reaktoro`) | one folder per layer, tests and bindings beside each file; ActivityModel as a function object; autodiff | **Adopted**: folder per layer, tests beside code, PhaseModel as function-like interface; **deferred**: autodiff types (FP order change → Stage B/07) | 2, 3.9, 07 |
| Cantera ([vcs_MultiPhaseEquil](https://cantera.org/dev/cxx/da/df6/classCantera_1_1vcs__MultiPhaseEquil.html)) | phase objects shared by the mixture and the solver; one `equilibrate()` front end choosing solvers | **Adopted** (front end = run pipeline; solver legs behind it) | 3.2–3.3 |
| MAGEMin ([Riel et al. 2022](https://doi.org/10.1029/2022GC010427)) | C library parallel over single-point calls; LP + PGE + local minimisation | **Adopted** as model for node-level parallelism (no shared state); method ideas → 07 | 07 |
| PHREEQC/IPhreeqc ([Charlton & Parkhurst 2011](https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqc.pdf)) | move all globals and statics into one instance class; small string API | **Adopted** (exactly the H1–H6 removal) | A6–A9 |
| OSQP code generation ([Banjac, Stellato et al. 2017](https://web.stanford.edu/~boyd/papers/osqp_embedded.html)), FORCES | parametric problem → generated solver; factorisation caching; warm start | **Adopted as idea** for 07 (generated model code with derivatives); **rejected** for the kernel (GEMS3K problems change structure per call) | 07 |
| CEA, IPOPT (TNLP callback class), Clarabel | — | not checked online (no search hit); IPOPT's callback-object pattern matches the Optima adapter's objective; no decision needed | — |
| ThermoFun v2 (local `../thermofun`, `docs/explanation/legacy-and-v2.md`, `autodiff.md`) | named methods and parameters with units; legacy positional format kept; staged deprecation; autodiff two-pass | **Adopted** for the BASE_PARAM/IPM-file migration (02 §6) and as model for derivatives in 07 | 02, 07 |
| ORCHESTRA (local `../orchestra_cpp`: Expander, Parser→ExpressionNode, ExpressionGraph, Calculator) | models as data: text objects expanded, parsed to expression trees, evaluated in 5 stages; per-node memory; FD Jacobian | **Adopted: the idea** "model = data + explicit per-node memory block" (PhaseModel input/output block, 3.9). **Rejected: run-time interpretation** (virtual call per node; activity evaluation is already the largest cost, 03 §3) and the text Expander (GEMS3K has JSON input). Fits closed-form excess models (ideal, Margules, regular, Redlich–Kister, Van Laar, Guggenheim, DH/Davies, SIT, most of Pitzer); does not fit iterative EoS (CG, cubic EoS), ELVIS (integration), EUNIQUAC tables, CEF multisite without extension. Useful form: declarative description → generated C++ with analytic/automatic derivatives (07 §3) | 3.9, 07 |

## 8. What stays as it is (on purpose)

- `MULTI` as a C struct of raw arrays (GEMSGUI dependency; SoA layout is good for the hot loops).
- The char-coded enums in files and in `MULTI` (file format); typed enums only at new interfaces.
- The Optima fork and its feature macros (owner's decision space; 08 Q8).
- JAMA/TNT (used on cold paths; replacing them changes FP order → Stage B with scoring).
