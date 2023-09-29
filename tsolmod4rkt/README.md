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

        std::shared_ptr<TSolModMulti> multi(new TSolModMulti());
        if( multi->GEM_init(input_system_file_list_name) )  {
            std::cout << "error occured during reading the files" << std::endl;
            return 1;
        }

        auto& phase = multi->get_phase("Plagioclase");
        std::map<std::string, double> wx = {
            {"Albite", 0.186993363098213},
            {"Anorthite", 3.45294711467247e-09},
            {"Sanidine", 0.81300663344884}};

        phase.SetMoleFractionsWx(wx);
        phase.SolModActCoeff();
        auto ln_gamma = phase.GetlnGamma();
```


#### Python API

Added  Python API for using tsolmod4rkt (TSolMod) activity models from Reaktoro (via wrappers on these API methods) or other chemical solver codes, or simply in Python files or Jupyter notebooks.

```python

from solmod4rkt import *

task = TSolModMulti()
task.GEM_init("Thermo-time-in/series1-dat.lst")

phase = task.get_phase("Plagioclase");
wx = {'Albite': 0.186993363098213, 'Anorthite': 3.45294711467247e-09, 'Sanidine': 0.81300663344884}
phase.SetMoleFractionsWx(wx, 0.)
phase.SolModActCoeff()
ln_gamma = phase.GetlnGamma()
print("\nln_gamma: ", ln_gamma)

map_ideal = phase.SolModIdealProp()
print("\nIdealProp: ", map_ideal)

map_excess = phase.SolModExcessProp()
print("\nExcessPropp: ", map_excess)

```

### Changes into main GEMS3K

1. Implemented DataCH and DataBR API (datach_api.h, datach_api.cpp and datach_formats.cpp). Functions for allocation and reading/writing DATACH and DATABR structures are separated from the TNode class.

2. Added phase and components names to initals structure of TSolMod class, and implemented printing of TSolMod structure to JSON and key-value formats for comparing results.

3. Added subdirectory *tsolmod4rkt* with the TSolModMulti class to initialize and manage Phase models and the SolModCalc class as C++ API for phase models.


### To Do

1. Research addition parameters for specific models (AddSolutionData struct), which could be added to the main SolutionData and printing function
2. Add to SolModCalc class functions for using data as input or result if needed. For example some models use not only thermodynamic and mole fractions as input.

```
addsd.arZ = pm.EZ+jb;
addsd.arM = pm.Y_m+jb;
addsd.ardenW = pm.denW;
addsd.arepsW = pm.epsW;
addsd.arG0 = pm.G0+jb;
addsd.arFWGT = pm.FWGT+k;
addsd.arX = pm.X+jb;
```

3. After the first testing and discussion, it would be need to increase the number of `sets` and `gets` functions.

4. Change TSOLMOD_MULTI (TSolModMulti) to a truncated MULTI data structure variant. Now full. It is more easily delete than add. 


