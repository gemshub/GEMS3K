## GEM2MT Module Operation

The **GEM2MT** module is a tool for automation of one-dimensional (1-D) reactive mass transport simulations coupled with GEM calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps. **GEM2MT** can be used to simulate the open system evolution and fluid-rock reactions occurring along a fluid flow path (e.g. vein, fracture, porous media, etc.). This type of simulation solves simple reactive mass transport problems, which requires the definition of fluxes between boxes of certain fluid/rock mass or volume ratios (or alternatively porosity and permeability) where fluid is moved sequentially or simultaneously from one box to the other.

The **GEM2MT** module can simulate reactive transport in three main modes:

* Sequential reactor chain and waves ('S');
* Flow-through and box-flux sequence ('F', 'B');
* One-dimensional reactive transport with advection/dispersion/diffusion ('A', 'C', 'W')

### Building and installation

The *GEM2MT* library is built and installed along with GEMS3K when an option BUILD_GEM2MT=ON is set for CMake. as shown below. To do that in linux and to install in user's home directory or in the system directory as shown below, a typical sequence of commands can be executed in the terminal or in bash script (assuming that the source code tree of GEMS3K resides in /home/username/GEMS3K folder):

```sh
cd ~/GEMS3K
sudo ./install-dependencies.sh
mkdir build
cd build
cmake .. -DBUILD_GEM2MT=ON
make -j 4
sudo make install
```
This will build and install GEMS3K and *GEM2MT* libraries into /usr/local/lib and /usr/local/include system folders (this needs sudo before make install).

### Tutorial in C++

#### Sequential reactor chains and waves ('S' mode)

Such simulations are useful for predicting the development of rock alteration zones upon progressive infiltration of a fluid supplied from a constant source (e.g. a large fracture). This model is a good approximation for steady-state irreversible processes such as weathering, metasomatism, evaporation of seawater etc.

In this mode, we can simulate a "wave" of a fluid that will pass through each of the rock nodes once per step (i.e., a wave consists of ca. 100 substeps). Many "waves" are possible. 

```
 // Here would be example with parameters (new API)
```

#### Flow-through box-flux models (F, B modes)

In the 'F' mode, the 1-D reactive transport column is represented by many equal-size boxes that may have different rock compositions, temperatures and pressures, if desired. The first box must be coded as a constant-flux (Cauchy) source (3), and the last box - as a constant-flux sink (-3); all the boxes in-between are normal local-equilibrium systems (1). 

At each unit time step, a given mass of reacted aqueous fluid is moved into the next box replacing there the equivalent mass of fluid that moves simultaneously to the next box and so on (fresh fluid from the constant-flux source moves from box 0 into box 1, and the same mass of reacted fluid moves away from the last box). In the 'F' mode, the transport is performed by direct iteration; in the 'B' mode, this simultaneous transport is conducted using the built-in integrator of the ODE system (the Bulirsch-Stoer algorithm) of mass balance equations for all elements in all boxes.
After the transport step, the equilibrium speciation is computed in all nodes by the GEM algorithm, and the next iteration of transport occurs followed by the next equilibrium in all boxes, and so on, until the maximum number of steps or the total time of simulation is reached. After each equilibration, the masses, volumes or concentrations of interest (both in aqueous and solid parts) are **sampled and plotted**  (would change to call back function).


In this mode, we can simulate the simultaneous transport of (reacted) fluid in one direction through all nodes from one node to the next. 

```
 // Here would be example with parameters (new API)
```

#### One-dimensional reactive transport ('A', 'C', 'W' modes)

The GEM2MT module allows the user to set up and run simplified  simulations of 1-D reactive transport in porous media using the finite-difference advection-dispersion algorithm (classical 'A' mode and implicit Crank-Nicolson 'C' mode) for advection-dominated cases. For the diffusion-dominated transport, the random-walk algorithm ('W' mode, also accounts for advection) can be used; the pure diffusion transport method ('D' mode TBD) is under construction in this version of GEM-Selektor. Theoretical background for all these transport modeling methods can be found in numerous textbooks.

The non-iterative coupling of transport and chemistry is implemented in all cases (so-called "operator splitting" approach): the mass transport model integration over all nodes alternates with GEM calculation of equilibrium states in all nodes after each time step, until the maximum time or number of time steps is reached. The profiles of amounts of phases or aqueous species of interest **can be plotted after each time step (or exported to VTK files)**   (would change to call back function)..


```
 // Here would be example with parameters (new API)
```
