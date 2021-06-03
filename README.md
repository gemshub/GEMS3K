# GEMS3K

Numerical kernel solver of the GEM-Selektor v.3 package for geochemical modelling. 
Just extended with an alternative Activity C++ API as the basis for a new Python API (based on Pybind11) on top of it, see xGEMS repository.

### Briefly about GEMS3K

The code Implements the improved GEM IPM-3 algorithm with excellent mass balance precision and fast convergence to Gibbs energy minimum even in very complex non-ideal chemical systems with two-sided metastability constraints (learn more on GEMS3K web page).

The code is written in C/C++. Using compiler directives, the GEMS3K code can be compiled as a standalone program e.g. 'gemcalc'; as a static or dynamic library for coupling with the mass transport simulator or another code; or as part of the GEM-Selektor v3 code together with GUI and databases.

Input: The standalone GEMS3K code needs to be initialized by reading a set of GEMS3K I/O files (since v.3.8.1, also a set of JSON documents) that can be exported from GEM-Selektor code for any chemical system (SysEq record) or from any GEM2MT task.

* Version: currently 3.8.1.

### License

GEMS3K is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GEMS3K is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with GEMS3K code. If not, see http://www.gnu.org/licenses/. 

### How to clone (download) GEMS3K source code

In your home directory, make a folder named e.g. ~/gitGEMS3 with a subfolder gitGEMS3/standalone.

Change into ~/gitGEMS3/standalone and clone this repository from https://bitbucket.org/gems4/gems3k.git using git, a preinstalled free git client SourceTree or SmartGit (the best way on Windows). 

Alternatively on Mac OS X or linux or Windows10, open a terminal and type in the command line (do not forget a period) to download the actual "trunk" branch:
~~~
git clone https://bitbucket.org/gems4/gems3k.git . 
~~~

To switch to another branch (e.g. devEJDB), use a git client or open a terminal, cd to ~/gitGEMS3/standalone, and type in the command line
~~~
git checkout -b branches/devEJDB --track origin/branches/devEJDB
git pull origin branches/devEJDB
~~~

To switch back to trunk, type
~~~
git checkout trunk
~~~

### How to build GEMS3K library and examples ###

The most common way to build GEMS3K on Linux or MacOS is using cmake ([https://cmake.org/](https://cmake.org/)). Please, make sure that you have cmake installed in your system. 

To build GEMS3K and install it in your home directory or in the system directory (as in the example below), a typical sequence of commands can be executed in the terminal:
~~~
cd ~/gitGEMS3/standalone
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
~~~

The same will be done by executing the install.sh script instead (check that this file has executable status, if not, run a command "chmod +x ./install.sh"): 
~~~
cd ~/gitGEMS3/standalone
sudo ./install.sh
~~~
(will ask for typing sudo password). Edit this script if you need to install GEMS3K libraries to different place than /usr/local/lib, or you would like to use different options.

For debugging purposes and for building examples, call cmake as follows: 
~~~
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/usr/local
~~~

Building demos and examples using cmake

* TBD (in progress), use QtCreator so far.

* For using qmake or QtCreator for building the GEMS3K chemical solver library and examples, please consult this web page: http://gems.web.psi.ch/GEMS3K/techinfo.html

#### Attention: 

Since version 3.7.0, the GEMS3K "ipm.dat" files exported by earlier versions of GEM-Selektor may need a small modification: the "<ID_key> " needs to be entered at the beginning of an exported "ipm.dat file". This can be done using any plain-txt editor (TextEdit, nano, SublimeText, VScode).

Introducing the <ID_key> key was necessary to improve the GEMS3K I/O consistent use of the legacy key-value and JSON data formats. However, to ensure the readability of old GEMS3K file sets (exported from GEM-Selektor before version 3.7), we have made a workaround patch that skips the string after (missing) <ID_key>. This patch is available since the git commit "ffb5283" made on January 15, 2021. 
 
Examples of the beginning of "ipm.dat" file valid since v.3.7.0: 

* KeyValue format with comments:
~~~
#  GEMS3K v.3.8.1 c.0aa600e 
# Comments can be marked with # $ ; as the first character in the line
# IPM text input file for the internal GEM IPM-3 kernel data
# (should be read after the DCH file and before DBR files)

# ID key of the initial chemical system definition
<ID_key>  "Kaolinite G  pHtitrKaS   0    0       1       25      0   "

## (1) Flags that affect memory allocation
# PE: Flag for using electroneutrality condition in GEM IPM calculations (1 or 0)
<pa_PE>  1
# PV: Flag for the volume balance constraint (on Vol IC) for indifferent equilibria at P_Sat (0 or 1)
..........
~~~

* KeyValue format without comments:
~~~
<ID_key>  "Kaolinite G  pHtitrKaS   0    0       1       25      0   "
<pa_PE>  1
..........
~~~

* JSON format (available since version 3.8.1):
~~~
[
   {
      "set": "Kaolinite",
      "ipm": 
	  {
         "ID_key": "Kaolinite G  pHtitrKaS   0    0       1       25      0",
         "pa_PE": 1,
...........          
~~~

More on GEMS3K I/O API, latest KeyValue and JSON formats will be added here soon. 

####TBD

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines
