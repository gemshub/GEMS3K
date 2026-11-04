GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 13 — GEMS3K without the GUI: define a chemical system by code or in YAML

**Goal (owner):** a user can make a chemical system without GEMSGUI — through a code API, as in Reaktoro (which has
no GUI), or in a YAML file — and choose the elements, species, phases and phase models there.

**Owner's choices:** the **C++ builder API is the core**. YAML is its serialised form: load and save give the same
object. Python gets the builder through xgems. **Species and standard states come from ThermoFun v2 databases**
(10); phase models come from ThermoSol records (11).

## 1. Today

Every GEMS3K system is made in GEMSGUI and exported. GEMS3K reads only the exported documents:
- `TNode::GEM_init(lst)` reads the DCH, IPM and DBR files named in a `.lst` file (node_copy.cpp l.107);
- `TNode::GEM_init(dch_json, ipm_json, dbr_json, thermo)` (node_copy.cpp l.234) and xgems
  `ChemicalEngine::initializeFromJsonStrings` (ChemicalEngine.hpp l.242) read the same documents from strings.

No code builds DATACH, MULTI and DATABR from a list of species. GEMSGUI does that work in these modules:

| GEMSGUI module | What it decides | Lines |
|---|---|---|
| `Modules/m_icomp.cpp` (IComp) | elements and the charge component: symbol, molar mass, valence, class | 240 |
| `m_dcomp*.cpp`, `m_reacdc.cpp`, `m_reac2.cpp` (DComp, ReacDC) | species: formula, class code, standard-state data and T-P method | ≈ 6 200 |
| `m_phase.cpp`, `m_phase2.cpp` (Phase) | phases: kind code, species list, mixing-model codes (`sMod`), parameters, sorption sites | 3 697 |
| `m_compos.cpp` (Compos) | recipes: amounts of formulas or compounds that give the bulk composition | 1 344 |
| `Submods/ms_rmults.cpp`, `ms_system.cpp`, `ms_sysbc.cpp` (RMULTS, SysEq) | which elements, species and phases are in the system; bulk; bounds DUL/DLL | 3 830 |
| `Submods/ms_mtparm.cpp` (MTPARM) | T-P point or grid | 1 056 |
| `m_param.cpp` | settings (`BASE_PARAM`) | — |
| `v_thermofun.cpp` + the DCH/IPM/DBR writers | export | 963 + |

## 2. How other codes do it

| Code | Database | Choosing elements, species, phases | Attaching phase models | Bulk and conditions | API or input file |
|---|---|---|---|---|---|
| **Reaktoro** (local `../reaktoro`) | `Database` objects; YAML databases (`embedded/databases/reaktoro/*.yaml`: `Species: {Name, Formula, Elements, AggregateState, StandardThermoModel}`); also PHREEQC, SUPCRT, ThermoFun, NASA databases | code: `Phases(db)` with `AqueousPhase(speciate("H O C Na Cl"))`, `GaseousPhase(...)`, `MineralPhase(...)`, `MineralPhases()`, `exclude(tags)` (Core/Phases.hpp l.45–480) | `.set(ActivityModelPitzer())`, `chain(...)` of model generators | `ChemicalState`, `EquilibriumConditions`, `EquilibriumSpecs` | C++ and Python API; no system file |
| **Cantera** | YAML input file (format known; no local repository) | `phases:` entries name their `elements`, `species` and `thermo` model | the `thermo:` key of the phase names the model; species carry their thermo data | `state: {T, P, X}` in the phase entry, or API | YAML file + API |
| **PHREEQC** (local `../solvers/phreeqc3`) | keyword database (`SOLUTION_MASTER_SPECIES`, `SOLUTION_SPECIES`, `PHASES`, `EXCHANGE_MASTER_SPECIES`, `SURFACE_MASTER_SPECIES`, `RATES`) | the problem keywords pick what is present: `SOLUTION`, `EQUILIBRIUM_PHASES`, `SOLID_SOLUTIONS`, `SURFACE`, `EXCHANGE`, `GAS_PHASE`, `KINETICS` | the database chooses the aqueous model (Debye–Hückel, Pitzer, SIT); solid solutions by keyword | in the keywords | input file; IPhreeqc runs the same strings |
| **Thermochimica / Equilipy** (local `../solvers`) | ChemSage `.dat` (CALPHAD) | all phases of the database, optionally a subset | in the database | Python `read_dat`, `equilib_single`, `equilib_batch` with composition, T, P | Python API |
| **MAGEMin** (local `../solvers/MAGEMin`) | built-in databases chosen by name (`*_database` folders) | the database fixes the phases | built in | bulk + T, P through the C/Julia API | API |
| **CEA** (local `../solvers/cea`) | `thermo.lib` | `reactants` and `only` datasets | ideal | `problem` dataset (T, P, ratios) | input file |
| **ThermoFun v2** (local `../thermofun`) | YAML/JSON records: `thermohub: 2.0.0-draft.1`, `thermodataset`, `substances: {formula, class, aggregate_state, charge, standard_state_properties, methods}` (examples/v2/minimal.yaml) | — (a database, not a system) | — | — | file format to follow |

**Lessons taken:**
- From Reaktoro: make phases from a database by an element list (`speciate`), then attach a model to each phase.
  The system is a code object; files are optional.
- From Cantera: one readable YAML file where each phase names its model and lists its species.
- From PHREEQC: the database is separate from the problem.
- From MAGEMin and Equilipy: a named database plus a bulk composition is enough for most users.
- From ThermoFun v2: SI inside, units on read, named records, version key at the top.

## 3. Design

### 3.1 The definition object

`gems3k::SystemDefinition` is a plain value object (copyable, comparable, no solver state):

| Part | Content | Today's GEMSGUI equivalent |
|---|---|---|
| database | one or more ThermoFun v2 datasets (file or string), ThermoSol records, ThermoKin records | project database |
| elements | symbols; the charge component `Zz` added when there are charged species | IComp |
| species | ThermoFun substance ids; the DC class code comes from the ThermoFun `class` and `aggregate_state`; optional bounds DUL/DLL | DComp/ReacDC selection |
| phases | name, kind (aqueous, gas, fluid, solid solution, single solid, sorption, …), species list **or** a speciate rule (element list + filters), model (ThermoSol record id + parameters with units), sorption sites (densities, capacities) | Phase |
| bulk | element amounts, or a recipe of formulas/substances with amounts and units (GEMSGUI Compos) | Compos, SysEq |
| conditions | T and P, one point or a list (no grid: standard states come from ThermoFun, 10) | MTPARM |
| settings | `pa_` fields by name (05 A3); defaults from `pa_p_` | Project settings |

### 3.2 The builder (C++ API, Reaktoro-like)

```cpp
// Build a system from a ThermoFun database. Each call adds to the definition.
// The builder checks the definition before it makes a node.
gems3k::SystemBuilder b( ThermoFun::Database("aq17-gem.json") );
b.elements("H O C Ca Na Cl")
 .aqueous("aq_gen").speciate().model("thermosol:Davies")       // species of the aqueous phase from the elements
 .gas("gas_gen").species({"CO2", "H2O", "CH4"}).model("thermosol:Ideal")
 .solids({"Calcite", "Halite", "Portlandite"})
 .bulk().add("H2O@", 1000.0, "g").add("CaCO3", 0.01, "mol").add("NaCl", 0.1, "mol")
 .conditions(298.15, 1e5)
 .setting("pa_PC", 2);
gems3k::SystemDefinition def = b.definition();   // the value object
TNode node( def );                               // or node.GEM_init( def )
```

Checks of the builder (each gives a clear message, never a silent default):
- every species exists in the database;
- every element of every species is in the element list;
- the phase kind fits the species classes;
- the model exists in ThermoSol and has its parameters;
- the bulk is charge-balanced or the charge component is present;
- the units are known.

### 3.3 What the builder must calculate (the work GEMSGUI does now)

| Result | How | Existing code to reuse |
|---|---|---|
| stoichiometry matrix A, element order | parse each formula | ChemicalFun `FormulaToken`, `DBElements`, `formulasProperties` (`ChemicalFun/FormulaParser/ChemicalData.h`) |
| molar masses `ICmm`, `DCmm` | element masses × A | ChemicalFun `DBElements`; GEMS3K restores `DCmm` from A today (ms_multi_format.cpp) |
| DC class codes (`ccDC`), phase codes (`ccPH`), `nDCinPH`, `LO`, `PG` | from ThermoFun class/aggregate state and the phase kind | the code tables in m_const_base.h |
| mixing-model layout (`sMod`, `LsMod`, `LsMdc`, `IPx`, `PMc`, `DMc`, `MoiSN`, `DQFc`) | from the ThermoSol record; multisite moieties from formulas | ChemicalFun `MoietyParser.h`; tsolmod4rkt `SolModFactory`; replaced by ThermoSol records after B15 (S6) |
| kinetics arrays | from ThermoKin records | 11 §3 |
| bounds DUL/DLL, `RLC`/`RSC` codes | defaults as GEMSGUI, or user values | GEMSGUI `ms_system.cpp` rules |
| `bIC` | recipe × stoichiometry; units converted | GEMSGUI `m_compos.cpp` rules |
| settings | names → `BASE_PARAM` fields | A3 `base_param_fields.def` |

### 3.4 Two output paths

- **Path 1 (first; no new reader):** the builder writes the existing DCH, IPM and DBR JSON documents and the
  ThermoFun JSON, then calls `GEM_init(dch_json, ipm_json, dbr_json, thermo)`. Every check that GEMS3K does on
  read still runs, and the documents can be saved and compared with a GEMSGUI export.
- **Path 2 (later):** `TNode::GEM_init(const SystemDefinition&)` fills DATACH, MULTI and DATABR directly
  (`dbr_dch_api::datach_realloc`, `databr_realloc`, `TMultiBase::multi_realloc`), with no text step.

### 3.5 YAML form (`system.yaml`)

The YAML file is the serialised `SystemDefinition`. `load(save(def)) == def` is a test. Units follow ThermoFun v2
(SI inside, a unit given where the value is not SI).

```yaml
gems3k: 1.0.0-draft.1
system:
  name: calcite-nacl
  database:
    thermofun: [aq17-gem.yaml]          # ThermoFun v2 datasets (files or dataset names)
    thermosol: [aqueous-models.yaml]    # ThermoSol records (phase models)
  elements: [H, O, C, Ca, Na, Cl]       # Zz is added when charged species are present
  phases:
    - name: aq_gen
      kind: aqueous
      speciate: true                    # all aqueous species of the database made of the elements
      model: {id: Davies}
    - name: gas_gen
      kind: gas
      species: [CO2, H2O, CH4]
      model: {id: Ideal}
    - name: Calcite
      kind: single_solid
      species: [Calcite]
  bulk:
    recipe:
      - {formula: H2O, amount: {value: 1000, unit: g}}
      - {formula: CaCO3, amount: {value: 0.01, unit: mol}}
      - {formula: NaCl, amount: {value: 0.1, unit: mol}}
  conditions: {T: {value: 25, unit: Cel}, P: {value: 1, unit: bar}}
  settings: {pa_PC: 2}
```

### 3.6 Python (xgems)

xgems binds `SystemDefinition` and `SystemBuilder` with pybind11, and `ChemicalEngine` takes a definition beside
its file and JSON-string forms (02 §8 migration rules; the old forms stay).

```python
import xgems
b = xgems.SystemBuilder(xgems.ThermoFunDatabase("aq17-gem.json"))
b.elements("H O C Ca Na Cl").aqueous("aq_gen").speciate().model("Davies")
b.solids(["Calcite", "Halite"])
b.bulk().add("H2O@", 1000.0, "g").add("CaCO3", 0.01, "mol")
engine = xgems.ChemicalEngine(b.definition())
```

## 4. Validation

1. **Every corpus project again, without GEMSGUI.** A converter (S3) reads the legacy export of a project
   (DCH, IPM, DBR, ThermoFun JSON) and writes `system.yaml`; no GEMSGUI is needed for this. The builder then makes DCH/IPM/DBR, which are compared with
   GEMSGUI's export field by field: names, A, codes, model arrays, bulk, bounds, settings. Every difference is either
   fixed or recorded with its reason.
2. **Freeze both** (GEMSGUI export and headless definition) in all six modes; score them as 10 §3 b. The
   expected result is identical rows; a difference means a derivation of §3.3 is not yet the same as GEMSGUI's.
3. **Examples and tests in GEMS3K:** small systems written in C++, Python and YAML (CalciteIC, Kaolinite, one solid
   solution, one sorption phase, one gas mixture); a YAML round-trip test; error-message tests for each builder
   check.

## 5. Steps (track S — additive: existing paths do not change, so no Stage A freeze risk)

| Step | What | Depends on | Check |
|---|---|---|---|
| S1 | `SystemDefinition` value object + YAML load/save (round trip) + schema document | A3 (settings by name) | round-trip tests |
| S2 | `SystemBuilder` with the checks and the derivations of §3.3; output path 1 (JSON documents → `GEM_init`) | S1, A15 (ThermoSource) | examples of §4.3 solve; documents compared with GEMSGUI exports for the same systems |
| S3 | converter legacy export (DCH/IPM/DBR + ThermoFun JSON) → `system.yaml` (tool `gems3k-convert`); whole-corpus round trip legacy → YAML → builder → legacy (§4.1–§4.2) | S2 | field comparison report; freezes |
| S4 | Python binding in xgems; `ChemicalEngine(definition)` | S2 | xgems tests, CI on Linux, Windows and macOS |
| S5 | output path 2: direct `GEM_init(const SystemDefinition&)` | S3 | the same freezes as S3, row for row |
| S6 | ThermoSol and ThermoKin records instead of the legacy model and kinetics arrays | B15 (11) | S3 comparison again |
| S7 | named results file (`results.yaml` / `.json`, §7.3) written by `TNode`, xgems and the tools | S1 | schema tests; values equal to the DBR values of the same run |
| S8 | GEMSGUI exports `system.yaml` + a ThermoFun v2 dataset of only the species used (§7.2); legacy export kept as an option for one release cycle | S3, GEMSGUI PR | S3 round trip on GEMSGUI's own test projects |

New module: `system/` (CMake target `gems3k_system`): builder, YAML reader and writer, derivations. It depends on
`data`, `io`, `thermo`, ThermoSol/ThermoKin and ChemicalFun (04 §2, §4).

## 6. Open questions (08 Q30–Q34)

- Q30: sorption and surface complexation (sites, densities, electrostatic models) in S2, or later? Default: in S2,
  with GEMSGUI's arrays.
- Q31: should `system.yaml` also hold processes and sweeps (titrations, T series)? Default: no; a separate
  `run.yaml` later.
- Q32: where is the schema documented? Default: `Docs/system-yaml.md` plus a JSON Schema file used by the reader.
- Q33: should GEMSGUI write `system.yaml` too? **Answered by the owner: yes** — it becomes GEMSGUI's export (§7).
- Q34: results file format — YAML or JSON, one file per node or one per run? Default: JSON for programs and YAML
  for people from the same writer; one file per run with one entry per node and step.

## 7. User-friendly files (owner decision)

### 7.1 What is wrong with the export today

GEMSGUI writes the solver's internal memory, not a user format. Example from the corpus,
`Resources/gems3k/j_10TH_G_seawater_1_0_1_25_0` (12 elements, 122 species, 87 phases):

| Problem | Evidence |
|---|---|
| six files for one system | `-dat.lst`, `-dbr.lst`, `-dch.json`, `-ipm.json`, `-dbr-0-0000.json`, `-fun.json`: 3.9 MB in total |
| phase models as character codes and positional counts | `"sMod": ["ZNNSNNNN", "INNINNNN", "GNNSNNNN"]`, `"LsMod": [242, 4, 6, 0, 0, 0, 1, 0, 3]` |
| interaction parameters keyed by index lists | `"IPxPH": [2, 8, -1, -10, 2, 8, -1, -11, …]` |
| the stoichiometry and the index maps as flat number arrays | `"A"`, `"xic"`, `"xdc"`, `"xph"`, `"ccIC"`, `"ccDC"`, `"ccPH"`, `"nDCinPH"` |
| the system file is mostly a T-P table | `-dch.json` 866 KB: `nTp` = 81 temperatures × 122 species × 7 properties (removed by 10) |
| the ThermoFun file holds the whole dataset | `-fun.json` 2.9 MB for 122 species used |
| settings mixed with allocation flags | `pa_*` values next to `PV`, `PSOL`, `PAalp`, `PSigm`, `Lads`, `FIa`, `FIat` |
| results by index | DBR `bIC`, `xDC`, `xPH`, `gam` are arrays in DATACH/DATABR index order; the names are in another file |
| a setting placed in the wrong block is ignored without a message | trap 2 (`<END_DIM>`) |

### 7.2 The new export

- **GEMSGUI exports `system.yaml`** (§3.5) **plus a ThermoFun v2 dataset that holds only the species the system
  uses** (with their methods and sources). These two files are what users read, edit, share and keep in version
  control.
- **DCH/IPM/DBR become the internal legacy format.** The builder makes them in memory (path 1, then path 2).
  GEMS3K still reads them (all existing projects and corpora keep working), and can write them on request
  (`--legacy`) for old tools and for debugging. Users do not edit them any more.
- **Migration** (ThermoFun `legacy-and-v2.md` pattern):
  - release N: GEMSGUI writes both forms, `system.yaml` first; the converter `gems3k-convert` (S3) turns existing
    exports into `system.yaml`;
  - release N + 1: GEMSGUI writes only `system.yaml` by default; reading legacy files logs one `info` line that
    names the converter;
  - legacy reading stays at least until 5.0 (02 §8).
- **Node arrays** (TNodeArray, transport): `system.yaml` defines the chemistry once; the node states are a list of
  bulk compositions and conditions (`nodes.yaml` or a table), not one DBR file per node.

### 7.3 The results file

`results.yaml` (people) and `results.json` (programs) come from one writer. Everything is named and has units:

```yaml
gems3k-results: 1.0.0-draft.1
system: calcite-nacl              # the system.yaml it belongs to
run:
  - node: 0
    status: {code: OK_GEM_AIA, mode: AIA, iterations: 27, stale: false}   # stale = true on ERR rows (trap 6)
    conditions: {T: {value: 298.15, unit: K}, P: {value: 1.0e5, unit: Pa}}
    system_properties: {G: {value: -1.234e6, unit: J}, V: {value: 1.0e-3, unit: m3}, mass: {value: 1.01, unit: kg}}
    aqueous: {pH: 8.27, pe: 3.1, Eh: {value: 0.18, unit: V}, ionic_strength: {value: 0.105, unit: molal}}
    phases:
      - name: aq_gen
        amount: {value: 55.6, unit: mol}
        mass: {value: 1.006, unit: kg}
        saturation_index: 0.0
        species:
          - {name: Ca+2, amount: {value: 1.2e-4, unit: mol}, molality: 1.2e-4, activity_coefficient: 0.38, log10_activity: -4.34}
      - name: Calcite
        amount: {value: 9.88e-3, unit: mol}
        saturation_index: 0.0
    elements:
      - {name: Ca, total: {value: 0.01, unit: mol}, dual_potential: {value: -1.0e5, unit: J/mol}}
```

Rules:
- one entry per node and step;
- phases that are absent are listed with amount 0 and their saturation index, so a user sees why;
- ERR rows say `stale: true` instead of silently repeating old values (trap 6, 02 §3);
- the writer reads the same getters as xgems (02 §2.1), so the file and the API agree;
- precision is configurable; the default writes 17 significant digits in JSON, 6 in YAML.
