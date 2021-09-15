
# This program

"PyBindE" is a python program for Molecular Mechanics Poisson-Boltzmann Surface Area (MMPBSA) calculations in protein-protein and protein-ligand systems.


## Basic Usage

  The main program is named pybinde_2021-09.py and receives as an input argument a .gro file (a simple the coordinate file for a GROMACS simulation) and the staring and ending atom numbers for the 2 sections of our system that we wish to use to calculate the binding free energy (i.e. starting and ending atom numbers for the protein and ligand).
  
  Example: pybinde_2021-09.py file.gro 1 1000 1001 1100
  
  By running pybinde_2021-09.py with the --help or -h flag, you can see additional optional arguments.
  
  In order to correctly work, the program requires 5 databases.
  - 2 databases named DataBaseT.crg and DataBase.siz, which contained partial atomic charges and atomic radii.
  - 1 database with information regarding atom types (the same information present in forcefield.rtp files).
  - 1 database containing c6 and c12 Lennard-Jones nonbonded parameters (present also in the forcefield_nb.itp file).
  - 1 file named 'classifier.config' with atomic radii.
  
  The program utilises a working DelPhi4Py implementation (to make use of its Poisson-Boltzmann solver), which is provided in the folder 'delphi4py'.

  Several settings may be easily changed or adjusted on the main program's file.
  The current settings are specific for the test system provided, (named 'frame00100.gro'), which is a MD simulation frame for a Beta-2-Microglobin (B2M) dimer.

  Different systems may have amino acid residues not found on the provided databases, missing atomic radii and charges in the SASA classifier file and in the DelPhi4Py databases.
  These can be manually added on the appropriate files, taking care to respect each file's specific formatting.

## Installation & Dependencies

This program does not require installation. However, given its dependencies, it should be used in a Linux-based system.
The dependencies for this program are:

* python2.6>= & python3.8>=
* numpy
* scipy1.7>=
* pandas1.3>=
* freesasa

In newer Ubuntu versions (21.xx), it may be required the installation of a package containing "libgfortran.so.4", for example gfortran-7 using the command:
$ sudo apt install gfortran-7

## License

  This program is distributed under a [LGPL-3.0](./LICENSE), however delphi4py depends on
  DelPhi which is proprietary. To use DelPhi the user is required to
  download the DelPhi license
  [here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

## Contacts

  Please submit a github issue to report bugs and to request new features.
  Alternatively you may find the developer [here](mailto:jnvitorino@fc.ul.pt).
