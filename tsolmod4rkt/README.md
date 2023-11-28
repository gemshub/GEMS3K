## tsolmod4rkt

**tsolmod4rkt** is the variant of TSolMod library of mixing and activity models used in GEMS codes (with all its subclasses and constructors) into a separate library with C++ and Python APIs. This allows using TSolMod models not only in GEMS3K and xGEMS code, but in any other C++ code such as Reaktoro (https://reaktoro.org) or independently from Python or Jupyter notebooks for testing or parameterization of solution models. 

Compared with the TSolMod used in GEMS3K, the **tsolmod4rkt** code uses a truncated variant of the MULTI work data structure that only contain arrays used in various mixing and activity models. TSolMod instances for all solution phases are initialized upon reading the ipm_dat (and minor parts of dch_dat, dbr_dat) JSON documents exported from GEMSGUI (GEM-Selektor). Other API methods connect TSolMod with the current data in chemical solver (e.g. concentrations, activity coefficients) or expose these data and bulk phase properties for collecting into tables or plotting. This makes it easy, especially using Python API, to test and fine tune parameters of various models, or develop ones newly added to TSolMod.

### Building and installation

The tsolmod4rkt library is built and installed along with GEMS3K when an option BUILD_SOLMOD=ON is set for CMake. as shown below. To do that in linux and to install in user's home directory or in the system directory as shown below, a typical sequence of commands can be executed in the terminal or in bash script (assuming that the source code tree of GEMS3K resides in /home/username/GEMS3K folder):

```sh
cd ~/GEMS3K
sudo ./install-dependencies.sh
mkdir build
cd build
cmake .. -DBUILD_SOLMOD=ON
make
sudo make install
```
This will build and install GEMS3K and tsolmod4rkt libraries into /usr/local/lib and /usr/local/include system folders. If you wish instead to install into your home directory then run the above script after repacing the 'cmake' command with the following one:

```sh
cmake .. -DBUILD_SOLMOD=ON -DCMAKE_INSTALL_PREFIX=/home/username/local
```

This will install both libraries into the indicated folder; sudo is not required for the 'make install' command then. 

In order to test that the built tsolmod4rkt library works, run the C++ example as follows:

```sh
cd build/bin$ 
./tsolmod4rkt_test
```
The last four lines of the console output should look like this:

```
...
Endmember activity coefficients:
   Albite 1.57948
   Anorthite 5.13304
   Sanidine 0.0537313
$
```


## C++ API use examples

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

## Python API use examples

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

### Changes made in the main branch of GEMS3K to retain consistency

1. Implemented DataCH and DataBR API (datach_api.h, datach_api.cpp and datach_formats.cpp). There, the functions for allocation and reading/writing DATACH and DATABR structures are separated from the TNode class.

2. Added lists of names of solution phases and their emdmembers (components) to the initial structure of TSolMod class, and implemented printing of TSolMod structure to JSON and key-value formats for comparing results.

3. Added a subdirectory *tsolmod4rkt* with source code of the SolModFactory class to initialize and manage solution phase models from the imported GEMS3K JSON documents, and the array of SolModEngine classes with C++ API methods for phase mixing and activity models (derived from TSolMod class instances).





