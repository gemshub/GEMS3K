GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 02 — API contract

This file lists everything outside GEMS3K that depends on GEMS3K, as found in the consumers' code, and the rules
the restructuring keeps. "Kept" means: same name, same signature, same behaviour, same file format, same order.
A change needs a migration (§6) with a compatibility layer and a removal date.

## 1. Consumers and how they build

| Consumer | Branch / head | Language, standard | Builds against | Uses |
|---|---|---|---|---|
| xgems | develop_optima / c1ea07f | C++20 + pybind11; CI ubuntu, windows, macos | installed GEMS3K, `#include <GEMS3K/node.h>` | TNode API, DATACH/DATABR fields, status codes, `gems3k_update_loggers` |
| GEMSGUI | develop_gems3kv5 / 1fa83316 | C++20, Qt6; CI linux, macos (2), windows | `GEMS3K::GEMS3K(-static)` | subclasses `TMultiBase`, `TNode`, `TNodeArray`; defines `BASE_PARAM::read/write`; binds BASE_PARAM fields by pointer |
| gems-benchmark (private, not read) | — | — | — | per brief: freeze tools read the `RUN/SET/EFF/KEY/CERT/DECIDE` trace and `# dec` lines; CTest 44 tests |
| node-gem, gemcalc, tools (in repo) | — | C++ (qmake and CMake) | source tree | TNode API |
| tsolmod4rkt (in repo) + Python module | — | C++ | source tree | `TSolMod`, `SolutionData`, `datach_api`, `io_*`, `jsonconfig`, `num_methods`, `v_service`, `verror` |

Compile-time switches that change class or struct layout (consumers must use the same value):
`USE_OPTIMA_SOLVER` (adds members to `TMultiBase` (ms_multi.h l.1394–1405) and methods to `TNode` (node.h l.258,
l.573)); `NO_NODEARRAYLEVEL` (removes transport fields from `DATABR`, databr.h l.74); `USE_NLOHMANNJSON` (I/O
classes). CMake exports `USE_OPTIMA_SOLVER` as PUBLIC; the qmake files do not — a qmake consumer can silently
get a different `TMultiBase` layout (ODR violation). Rule: keep the switches, and add a run-time check
(`gems3k_build_features()` returning the switches, compared by `GEM_init`) in A1.

## 2. Symbols used, per consumer

### 2.1 `TNode` (node.h)

| Group | Members used (count of call sites) | Who |
|---|---|---|
| Life cycle | `GEM_init` (lst; json strings), `GEM_run`, `GEM_read_dbr`, `GEM_write_dbr`, `GEM_print_ipm` | xgems, gemcalc, node-gem, tools |
| Data access | `pCNode` (xgems 30), `pCSD` (21), `cTK`, `cP`, `cVs`, `cMs`, `Get_TK`, `Get_P` | all |
| Index/name maps | `check_Phase_xCH` (24), `check_DC_xCH` (12), `check_IC_xCH`, `DC_xDB_to_xCH` (10), `DC_name_to_xCH`, `DC_name_to_xCH_map`, `IC_name_to_xCH`, `Ph_name_to_xCH`, `IC_name_to_xDB`, `Ph_name_to_xDB`, `DC/IC/Ph_xCH_to_xDB` | xgems, gemcalc, node-gem |
| Phase properties | `Ph_Mass` (14), `Ph_Mole` (10), `Ph_Volume` (8), `Ph_HeatCapacityCp`, `Ph_Entropy`, `Ph_Enthalpy`, `Ph_Present`, `Ph_SatInd`, `Ph_Moles`, `Ph_GibbsEnergy`, `phEnthalpies` | xgems, gemcalc, node-gem |
| Species | `Get_nDC`, `Get_cDC`, `Get_aDC`, `Get_gDC`, `Get_muDC`, `DC_n/a/c/g/mu`, `DC_G0`, `DC_V0`, `DC_H0`, `DC_S0`, `DC_Cp0`, `DC_A0`, `DC_U0`, `DCmm`, `DCaJI`, `Set_DC_G0` | xgems, gemcalc |
| Input setters | `Set_bIC`, `Set_TK`, `Set_P`, `Set_dul`, `Set_dll`, `check_TP` | xgems |
| Targets (Optima) | `Set_pH_target`, `Set_Eh_target`, `Get_ControlCondition_titrant`, `Clear_ControlConditions`, `Get_pH`, `Get_Eh` | xgems, tools |
| Diagnostics | `GEM_trace_regimes` (+ `TraceRegime`), `GEM_rescue_step`, `code_error_IPM`, `description_error_IPM`, `GEM_CalcTime`, `GEM_Iterations` | xgems, tools |
| Transport | `GEM_from_MT`, `GEM_from_MT_time`, `GEM_to_MT`, `GEM_step_MT`, `GEM_restore_MT` | node-gem |
| Virtual hooks (GEMSGUI `TNodeGUI`) | `packDataBr()`, `write_ThermoFun_format_stream`, `init_into_gems3k()` | GEMSGUI |

### 2.2 `DATACH` / `DATABR` fields (datach.h, databr.h)

xgems reads `DATABR.xDC, bIC, NodeStatusCH, dul, dll, Ms, pe, pH, Vs, IterDone, IC, Hs, Gs, Eh` and
`DATACH.ccPH, nDCinPH, nPS, PHNL, ICNL, DCNL, nPH, nIC, nDC, ICmm, DCmm`. GEMSGUI and node-gem use many more. The
structs, field names, types and order are kept (binary DBR/DCH files and TNodeArray copies depend on them).

### 2.3 Status codes (databr.h l.139–148)

`NEED_GEM_{AIA=1, SIA=5, AOP=10, SOP=14, ROP=18, HOP=22, SHP=26}`, `OK = NEED + 1`, `BAD = NEED + 2`,
`ERR = NEED + 3`, `T_ERROR_GEM = 9`. Used by xgems (all seven modes) and GEMSGUI (six modes). Kept; a new mode
gets a new block of four, and must be listed in `packDataBr`'s mode reset (node.cpp l.1738–1758) and in
`TNodeArray::CalcIPM_Node`'s success test (today it accepts only AIA/SIA, 09 §3 fact 2).

### 2.4 `TMultiBase` virtual interface (GEMSGUI `TMulti`, Modules/Submods/ms_multi_new.h l.95)

Overridden by GEMSGUI and therefore API (signatures kept exactly): `base_param()`, `set_def`, `multi_realloc`,
`testMulti`, `DC_LoadThermodynamicData`, `get_PAalp_PSigm`, `STEP_POINT`, all 23 `alloc_*`, `loadData`,
`GEM_IPM_Init_gui1/2`, `initalizeGEM_IPM_Data_GUI`, `calculateActivityCoefficients_scripts`, `multiConstInit_PN`,
`GetName`; plus GEMSGUI calls `GetPM()` (raw `MULTI*`) and reads `MULTI` fields directly. `TNodeArray` virtuals
`CalcIPM_List`, `CalcIPM_Node`, `pVisor_Message` (GEMSGUI `TNodeArrayGUI`).

**Consequence:** in Stage A `MULTI` keeps its field names and types, and `TMultiBase` keeps every virtual
listed here. Internal code may stop using `MULTI` directly (through views, 04 §3), but the struct stays.

### 2.5 `BASE_PARAM` (ms_multi.h l.63–539) — the positional contract

Four separate dependencies on member order:

1. GEMSGUI `BASE_PARAM::read/write` (GEMSGUI Modules/m_param.cpp l.236–247): `readArray(&PC, 10)` (10 `long`
   from `PC`), `readArray(&DG, 28)` (28 `double` from `DG`), 4 bytes. The declarations are in GEMS3K
   (ms_multi.h l.537–538). → the **first 38 members** keep type, order and contiguity.
2. GEMSGUI positional aggregate initializers `SPP_SETTING pa_` (m_param.cpp l.49) and `BASE_PARAM dfBase[5]`
   (l.93), and GEMS3K `pa_p_` (ms_multi_diff.cpp l.6–23) → `BASE_PARAM` stays an **aggregate** and the order of the
   initialized prefix is kept.
3. GEMSGUI `ResetTail` (m_param.cpp l.195–235) and 42 data objects `o_pax*` (Modules/v_mod.h, positional enum;
   `Resources/data/vis_od.ini.json`) bound with `SetPtr(&pa.p.Field)` → field **names** are kept.
4. The IPM file: tag names `<pa_…>` / JSON keys, the 164-field dynamic descriptor order (writing order), and the
   rule that new fields go after `tprn` and after every other member.

### 2.6 Other kept interfaces

| Interface | Kept as | Notes |
|---|---|---|
| `TSolMod(SolutionData*)`, `SolutionData` layout, `TSolMod` virtuals (`PTparam`, `MixMod`, `ExcessProp`, `IdealProp`, …) | as is in A | tsolmod4rkt + its Python module build TSolMod outside TMultiBase |
| `dbr_dch_api::*` free functions | as is | tsolmod4rkt, tests |
| `GemsSettings` / `gemsSettings().gems3k_update_loggers(...)` | as is | xgems |
| `SetCalculationCancel`, `CalculationCancelled` (verror.h l.83–84) | as is | callers' cancel button |
| spdlog logger names `gems3k`, `ipm`, `tnode`, `solmod`, `kinmet` | as is | users set levels by name |
| Environment variables (probes and traces, 09 §1 trap 14; `REAKTORO_TRACE_*` fallbacks) | as is | gems-benchmark scripts |
| Trace record formats `RUN`, `SET`, `EFF`, `KEY`, `CERT`, `DECIDE <name> k=v…` | **byte for byte** | freeze `# set/# eff/# key/# dec`; the column order is part of the contract |
| File formats: DCH/IPM/DBR key-value text, JSON (simdjson/nlohmann), binary (`GemDataStream`), `.lst` | byte for byte on write, tolerant on read | GEMSGUI writes them; corpora are stored in them |
| CMake package `GEMS3K` with targets `GEMS3K::GEMS3K`, `GEMS3K::GEMS3K-static`; headers installed as `<GEMS3K/*.h>` | kept; new targets added beside | xgems, GEMSGUI |

## 3. Behaviour that is part of the contract (keep in Stage A, even where it looks wrong)

- `GEM_run` returns the old status for a non-NEED code (trap 5).
- On `ERR_*` the DATABR values are stale (trap 6) and `TNodeArray` copies them back.
- `DC_G0()` returns the grid value, not `pm.G0` (trap 4); `Get_solver_G0()` returns `pm.G0`.
- `DC_name_to_xCH()` returns the first match (trap 8).
- `Get_GibbsEnergy()` recomputes and rewrites X/XF/XFA; `packDataBr()` evaluates absent-phase indexes (activity
  calls) on every call.
- `TNodeArray` reports Optima-mode OK codes as errors (09 §3 fact 2) — a bug, but changing it changes
  GEMSGUI/node-gem behaviour; it is fixed in Stage B step B1 with an owner decision (08 Q7).
- `copyMULTIData` copies 39 long ints and 53 doubles by layout (09 §3 fact 1) — kept exactly in A (as a named
  list that reproduces the same fields), fixed in B.

## 4. What may change freely (internal)

Everything not listed in §2–§3: file names and folders under `GEMS3K/`, non-virtual private and protected
`TMultiBase` methods, the free functions in `ipm_*.cpp`, thread-locals, the Optima adapter's internals, the
order of functions in a file, the qmake files (if the owner agrees to drop them, 08 Q4). Moving a function to
another translation unit is free **only if** the freeze and the instruction counts show no change (05 §0, A-rule 3).

## 5. Headers: public versus internal

Today every `*.h` is installed (`install(DIRECTORY … FILES_MATCHING PATTERN "*.h")`), including `tnt*`, `jama*`,
`gems3k_timing.h`, `lazy_file_sink.h`, `ipm_optima.h`. Proposed public set (installed, documented, stable):
`node.h`, `nodearray.h`, `datach.h`, `databr.h`, `datach_api.h`, `ms_multi.h` (because GEMSGUI subclasses
`TMultiBase` and reads `MULTI`), `s_solmod.h` (tsolmod4rkt), `s_kinmet.h`, `verror.h`, `v_detail.h`, `jsonconfig.h`,
`gems3k_impex.h`, `gdatastream.h`, `io_*.h` (tsolmod4rkt), `m_const_base.h`, `gems3k_version.h`. Internal (not
installed after the removal date): `tnt*.h`, `jama_*.h`, `gems3k_timing.h`, `lazy_file_sink.h`, `ipm_optima.h`,
`activities.h`, `kinetics.h`, `num_methods.h`, `v_service.h`, `s_sorpmod.h`. Before the removal date the internal
headers stay installed with a `#pragma message` deprecation note (opt-out macro `GEMS3K_NO_DEPRECATION_NOTE`).

## 6. Migration rules (for any change to §2)

1. **Add before remove.** The new form comes first; the old form stays as a thin forwarder with
   `[[deprecated("use X; removed in GEMS3K 5.0")]]` (or a `#pragma message` for macros and headers).
2. **Removal date** = the first release after **both** xgems and GEMSGUI have shipped a version that no longer
   uses the old form, and not earlier than **12 months** after the deprecation release. Proposed: deprecate in
   the release that ends Stage A (call it 4.5), remove in 5.0, not before 2027-12 (08 Q5).
3. **Files are versioned by content, not by date.** Readers stay tolerant (unknown tags are skipped, as today,
   but now with one `info` line naming the tag); writers keep the current format until the owner switches the
   default. A named-field format for `BASE_PARAM` follows the ThermoFun v2 pattern (thermofun
   `docs/explanation/legacy-and-v2.md`): v2 format with names and units; legacy positional format kept for read
   and write; a converter tool; legacy read logs a message naming the converter from release 1.x on.
4. **One pull request per consumer repo per change**, prepared by the change author, merged after the GEMS3K
   change is released (owner rule: no unasked pushes; every cross-repo change is listed in 05 per step).
5. **Contract tests in GEMS3K**: compile-time checks for §2.5.1 (`static_assert` on `offsetof` and types of the
   first 38 members), a test that writes and re-reads an IPM file and compares it byte for byte, a test that the
   trace record formats are unchanged on a fixed input (golden files), and a test that every `NEED_*` code is
   handled in `packDataBr` and `TNodeArray`.

## 7. Cross-repository touch list (per Stage A step; details in 05)

| Step | GEMS3K | xgems | GEMSGUI | gems-benchmark |
|---|---|---|---|---|
| A1 build targets | yes | CMake `find_package` unchanged | unchanged | freeze build script: new build folder name only |
| A2 layout asserts | yes | — | — | — |
| A3 BASE_PARAM single declaration | yes | — | **later PR**: generate `ResetTail`, `o_pax*` and `vis_od` entries from the same list (optional, 08 Q6) | — |
| A5 pure getters | yes (old names forward) | optional switch to the pure names | optional | — |
| A10 file moves | yes | include paths unchanged (installed headers keep names) | unchanged | tools that grep source paths, if any |
| B1 TNodeArray OK test | yes | — | review `TNodeArrayGUI` (overrides) | new freeze rows for node arrays |
