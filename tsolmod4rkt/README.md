## tsolmod4rkt

 The **tsolmod4rkt** library is the isolation of TSolMod (with all its subclasses and constructors) into a separate library.  Where TSolMod instances are initialized after reading the ipm_dat (and minor parts of dch_dat, dbr_dat) JSON documents sent from GEMSGUI into a truncated variant of  MULTI data structure, that only contain arrays used in TSolMod.

### Installing the library

The library is installed with GEMS3K when option BUILD_SOLMOD was set.
To build GEMS3K and install it in your home directory or in the system directory (as in the example below), a typical sequence of commands can be executed in the terminal:

```sh
cd ~/GEMS3K
mkdir build
cd build
cmake .. -DBUILD_SOLMOD=ON
make
sudo make install
```

## C++ API

```cpp

SolModFactory task(input_system_file_list_name);

auto& phase2m = task.SolPhase("Plagioclase");
std::map<std::string, double> x2m = {
    {"Albite", 0.187},
    {"Anorthite", 3.5e-09},
    {"Sanidine", 0.813}};

phase2m.SetMoleFractions(x2m);
phase2m.SolModActivityCoeffs();
auto ln_gamma = phase2.GetlnActivityCoeffs();

```


#### Python API

Added  Python API for using tsolmod4rkt (TSolMod) activity models from Reaktoro (via wrappers on these API methods) or other chemical solver codes, or simply in Python files or Jupyter notebooks.

```python

from solmod4rkt import *

task = SolModFactory("Thermo-time-in/series1-dat.lst")

phase = task.SolPhase("Plagioclase");
wx = {'Albite': 0.186993363098213, 'Anorthite': 3.45294711467247e-09, 'Sanidine': 0.81300663344884}
phase.SetMoleFractions(wx)
print(phase)

phase.SolModActivityCoeffs()
ln_gamma = phase.GetlnActivityCoeffs()
print("\nln_gamma: ", ln_gamma)

map_ideal = phase.SolModIdealProps()
print("\nSolModIdealProps: ", map_ideal)

map_excess = phase.SolModExcessProps()
print("\nSolModExcessProps: ", map_excess)

```

### Changes into main GEMS3K

1. Implemented DataCH and DataBR API (datach_api.h, datach_api.cpp and datach_formats.cpp). Functions for allocation and reading/writing DATACH and DATABR structures are separated from the TNode class.

2. Added phase and component names to the initial structure of TSolMod class, and implemented printing of TSolMod structure to JSON and key-value formats for comparing results.

3. Added subdirectory *tsolmod4rkt* with the TSolModMulti class to initialize and manage Phase models and the SolModCalc class as C++ API for phase models.





