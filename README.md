# MolecCav : Molecule trapped in optical Cavity
*Developed by*:  
Mathéo Segaud $^{1,\dagger}$  
*Under the supervision of:*  
Research Professor : David Lauvergnat $^{1,2,\dagger}$

${^1}$ Université Paris-Saclay, 3 rue Joliot Curie, Bâtiment Breguet, 91190 Gif-sur-Yvette  


## Contents
- [Quick introduction and Objectives](#Quick-introduction-and-Objectives) 
- [Setting up environment](#Setting-up-environment)
- [MolecCav Installation](#MolecCav-Installation)
- [General view and Structure](#General-view-and-Structure)
    - [```APP```](#APP)
    - [```Ext_Lib```](#Ext_Lib)
    - [```OBJ```](#OBJ)
    - [```OUT```](#OUT)
    - [```SRC```](#SRC)
        - [```MC_cavity_mode_m.f90```](#MC_cavity_mode_m.f90)
        - [```MC_operator_1D_m.f90```](#MC_operator_1D_m.f90)
        - [```MC_total_hamiltonian_m.f90```](#MC_total_hamiltonian_m.f90)
    - [```TESTS```](#TESTS)
        - [```test_cavity_mode.f90```](#test_cavity_mode.f90)
        - [```test_construct_op.f90```](#test_construct_op.f90)
    - [```The data files```](#The-data-files)
- [Compiling the library](#Compiling-the-library)
- [Running the tests](#Running-the-tests)
- [Running the application](#Running-the-application)
- [Authors](#authors)
- [Reference](#reference)


# Quick Introduction and Objectives
Theoretical studies in photochemical process and spectroscopy have great interest in describing quantically the dynamics of a light-interacting molecular system using the most efficient algorithms possible in order to be able to support computations with a high potential number of degrees of freedom (DOF).  

The aim is to provide numerical and methodological methods for a quantum description of an optical cavity-caged molecule, including the molecular modes with non-adiabatic couplings, and a quantum description of the electric field of the cavity modes and its couplings with the molecule. Numerical developments of a analytically-based computation of the cavity-matter couplings are carried out. This project aims to provide a Modern Fortran library to compute the effects of the cavity and the couplings between a single-mode cavity and the molecule, neglecting the magnetic contribution of the field, and based on analytical properties of the cavity model. The molecular part of the system is assumed to be computed by the code calling the library since it is supposed to be coupled to quantum dynamics programs.  

This library is designed, in a first time, to compute the effects of the interactions between **one** molecule and a **single** mode of the cavity it is confined in. The system is considered non-relativistic and in the dipole approximation. As a first approach, the cavity will be considered having infinite dimensions over the y and z-axis and small over the x-axis. No quantisation of the cavity field will thus be observed on the two former axis, and the considered single stationary cavity mode will resonate only along the latter one. As a consequence, the photon displacement coordinate $$\overrightarrow{x} \in R^3$$ can be reduced to a one-dimensional coordinate along this very axis $$x \in R$$. Within the framework of the Extended Jaynes-Cummings Model (EJCM), this cavity mode will be modeled by a quantum HO of mass m = 1 a.u., eigenpulsation $$\omega_c$$ and force constant $$k_c = m\omega_c$$, each eigenstate $$\ket{n}, \forall n\in N$$ of which i.e. each excitation quanta of which corresponding to the number of photons n resonating in the mode. As for the total system Hamiltonian, it is expressed as written below **for a given electronic state** as well as the expression of the photon field Hamiltonian which comes from the HO character of the electric field model (using atomic units and m = 1 a.u. property), and the expression of the matter and the cavity-matter coupling Hamiltonians which are written directly following the EJCM formalism, adapted to a single-mode cavity.  

Yet, one is actually not concerned with the expression of $$\hat{H}_{mat}$$ since it is not computed in the library and left to the quantum chemistry code calling MolecCav, nor is the action of the matter total dipole moment $$\hat{\mu}_{M}$$.  

```math
\hat{H}_{tot}(x,\overrightarrow{R}) = \hat{H}_{mat}(\overrightarrow{R}) + \hat{H}_{cav}(x) + \hat{H}_{CM}(x,\overrightarrow{R})
\hat{H}_{mat}(\overrightarrow{R}) = \hat{T}_N(\overrightarrow{R}) + E_{elec}(\overrightarrow{R})
\hat{H}_{cav}(x) = \hbar \omega_c\left(\hat{a}^{\dag}\hat{a}+\frac{1}{2}\right) = \frac{\hat{p}^2_x}{2m} + \frac{1}{2}k_c\hat{x}^2 \stackrel{\bra{x}}{=} -\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + \frac{1}{2}k_cx^2
\qquad \qquad = \omega_c\left(\hat{a}^{\dag}\hat{a}+\frac{1}{2}\right) = \frac{1}{2}\left(\hat{p}_x^2 + \omega_c^2\hat{x}^2\right) \stackrel{\bra{x}}{=} \frac{1}{2}\left(-\frac{\partial^2}{\partial x^2} + \omega_c^2x^2\right)
\hat{H}_{CM}(x,\overrightarrow{R}) = \lambda\omega_c\hat{\mu}_{M}(\overrightarrow{R})\hat{x}
```

More detailed information about the model can be found in the Manual "MolecCav_Manual.pdf"  


# Setting up environment
This library is designed to work within a GNU/Linux operating system and the compilation and execution files (makefile and run.sh) were written to use the gfortran compiler for Fortran90.  


# MolecCav installation
One can either install the library by downloading its files (if they are not interested into having a Git control of their version) or by Cloning it (otherwise) from its GitHub repository.


## Downloading the files
    1. Choose directory on your computer where the MolecCav directory is to be created.
    2. Open the GitHub repository of MolecCav and click on "<> Code" on top right of the page > Download zip
    3. Place it into the chosed directory and extract the files


## Cloning the files
    1. Choose directory on your computer where the MolecCav directory is to be created.
    2. Open the GitHub repository of MolecCav and click on "<> Code" at the top right of the project page
    3. In the menue that appeared, select HTTPS to make the URL appear and copy it.
    4. Open the terminal of the computer and ```bash cd``` to the chosed directory
    5. execute the command : ```bash git clone <the copied link>``` 


# General view and Structure
The MolecCav library is structured as follows :  

```.
MolecCav
├── APP
│   └── App_MolecCav.f90
│
├── Ext_Lib
│   ├── QDUtilLib
│   ├── QDUtilLib_loc
│   ├── QDUtilLib-main
│   ├── .gitignore
│   ├── cleanlib
│   └── get_Lib.sh
│
├── OBJ
│   └── obj
│       ├── *.o
│       └── *.mod
│    
├── OUT
│   └── *.log
│
├── SRC
│   ├── MC_cavity_mode_m.f90
│   ├── MC_operator_1D_m.f90
│   └── MC_total_hamiltonian_m.f90
│
├── TESTS
│   ├── test_action_op.f90
│   ├── test_cavity_mode.f90
│   └── test_construct_op.f90
│
├── .gitignore
├── data_app.nml
├── data_tests.nml
├── data.nml
├── makefile
├── MolecCav_Manual.pdf
├── README.md
└── run.sh
```  

More details can be found about each of these files and directories in the manual. Here will be only a general presentation.


## APP
This directory contains the source ```.f90``` file of the application program. It is written to provide an illustration of what the library can realise and how to set-up the simulations, as well as an additional test. As it has to work without the library being coupled to a program that manage the matter description, this simulation takes place in a somewhat particuliar model. It describes a diatomic molecule with harmonic electronic potential trapped in the cavity (so that it can be described using the procedures supposed to describe the cavity and the matter system has only **one** degree of freedom to take into account : $$\overrightarrow{R} = R$$), the dipole moment of the molecule $$\hat{\mu}_M(R)$$ is simplified thanks to a Taylor development and reduced by some approximation (cf. the Manual), and the total system wavefunction ($$matter system wavefunction \otimes photonic system wavefunction$$) is arbitrarily initiallized.  

```math
\mu(R) &= \mu_0 + \frac{\partial\mu}{\partial R}(R-R_0) 
\mu(R) &= \frac{\partial\mu}{\partial R}.R
\mu(R) &= Cte.R
```  

## Ext_Lib
This directory contains the external libraries used by MolecCav. Up to now, the QDUtilLib [[1]](#reference) library designed by David Lauvergnat.  


## OBJ
This directory contains all the object ```.o``` and ```.mod``` files generated as the library is compiled.  


## OUT
This directory contains all the output ```.log``` files generated as the tests or the application are executed.  


## SRC
This directory contains all the source ```.f90``` files of the modules that composes the library.  


### MC_cavity_mode_m.f90
This module contains the procedures to build and describe the modes of the cavity.  


### MC_operator_1D_m.f90
This module contains the procedures to build the quantum mechanics operators relative to a cavity mode, and to compute their action on a photon system vector.  


### MC_total_hamiltonian_m.f90
This module contains the procedures to compute the action of the total hamiltonian $$\hat{H}_{tot}(x,R)$$ over a matrix representing the total system's wavefunction.  


## TESTS
This directory contains all the source ```.f90``` files of the programs designed to test automatically that the library is working as expected. For now, it contains the two following tests.  


## The data files
Three are provided. Please do not modify the two dedicated respectively to the execution of the tests and the application. The third one is the one supposed to be modified and used by the user for specific uses.


### test_cavity_mode.f90
This program test that the corresponding module initializes well the cavity modes from the data file it is supposed to read.  


### test_construct_op.f90
This program test that the corresponding module initializes well all the operators for a cavity mode and that their matrices are consistent with the theoretical ones.  


# Compiling the library
The library can be built using either the makefile or the run.sh shell script. Again, more details are in the manual. Both support the following commands to manage the library :  

    - "getlib" : install the external library if not already downloaded.  
    - "lib" : Compiles the external library's modules and build the .a static library file if needed, then does the same for the MolecCav library.  
    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and link them with the static library file into executables.
    - "clean" : Deletes all the object, executable and output files from MolecCav (not those of the external libraries).
    - "cleanall" : Performs the same cleaning as "clean" in addition also to delete the .mod and static library files, and run the cleanlib shell script from the external QDUtilLib library that takes the same actions upon its files.

They can be executed by commad-line :  
```bash
make <name of the command>
./run.sh <name of the command>
```


# Running the tests
The tests can be compiled and executed using either the makefile or the run.sh shell script. The shell script allows to choose to compile and execute only one of the tests and to choose the data file to read if one is not interested in the default one. It uses these commands on both the makefile and the script :

    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and link them with the static library file into executables.
    - "ut" : Performs the same actions as "all" in addition also to execute the tests, using the indicated data file and direct the output into the corresponding .log output file. Then it grabs the final sentence of these files for each test and display them on the screen. In both makefile and run.sh, it is the default command if no arguments are passed when they are called.  


# Running the application
The application can be compiled and executed using either the makefile or the run.sh shell script. The shell script allows to choose data file to read. It uses these commands on both the makefile and the script :  

    - "all" : Performs the same actions as "lib"" in addition also to compile the source files of the tests and link them with the static library file into executables.
    - "app" : Performs the same actions as "all" in addition also to compile the application source file, link it with the library into an executable and execute it using the indicated data file and direct the output into the corresponding .log output file.

## Authors
*Developed by*: \
Mathéo Segaud $^{1,\dagger}$ \
*Under the supervision of:* \
Research Professor : David Lauvergnat $^{1,2,\dagger}$

${^1}$ Université Paris-Saclay, 3 rue Joliot Curie, Bâtiment Breguet, 91190 Gif-sur-Yvette \

## Reference
[1] [https://github.com/lauvergn/QDUtilLib)  
cf. manual  
