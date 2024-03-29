# PyBindE

**PyBindE** is a python program for Molecular Mechanics Poisson-Boltzmann Surface Area (MMPBSA) calculations in protein-protein and protein-ligand systems, using GROMOS54a7 parameters.

## Basic Usage

The main program is named _pybinde_2021-09.py_.

In order to work correctly, the program requires 5 databases.

- 2 databases named DataBaseT.crg and DataBase.siz, which contain partial atomic charges and atomic radii, respectively.
- 1 database with information regarding atom types (the same information present in forcefield.rtp files).
- 1 database containing c6 and c12 Lennard-Jones nonbonded parameters (present also in the forcefield_nb.itp file).
- 1 file named 'classifier.config' with atomic radii.

The program utilises a working _DelPhi4Py_ implementation (to make use of its Poisson-Boltzmann solver), which is provided in the folder _'delphi4py'_.

Several settings may be easily changed using the program's flags. To see all the flags, run the program with the flag '-h' or '--help'.

There is an _'example'_ folder which contains a test case for basic usage. A _'frame00100.gro'_ file is provided, (which is a MD simulation frame for a Beta-2-Microglobin dimer), as well as a bash script _'example_run.sh'_ with the a basic usage example.

Different systems may have amino acid residues not found on the provided databases, missing charges or atomic radii in the classifier file and in the DelPhi4Py databases (DataBaseT.crg and DatabaseT.siz). These can be manually added on the appropriate files, taking care to respect each file's specific formatting (be extra careful with extra spaces).

**IMPORTANT NOTE**: it's **crucial** that your system (present in the .gro files you intend to use to run MMPBSA calculations) is **centered in the simulation box** in every frame, with no breaks in the Periodic Boudary Condition, otherwise the results cannot be trusted.

## Installation & Dependencies

This program does not require installation. However, given its dependencies, it should be used in a Linux-based system.
The dependencies for this program are:
|**Dependency**|**Suggested command for installation**|
|------------|-------------------------------------------|
|* python2.6>= & python3.8>=|sudo apt install python3.9|
|* numpy | pip3 install --upgrade numpy|
|* scipy1.7>= | pip3 install --upgrade scipy|
|* pandas1.3>= | pip3 install --upgrade pandas|
|* freesasa | pip3 install --upgrade freesasa|
|* delphi4py | pip3 install --upgrade delphi4py|

In newer Ubuntu versions (20.xx or above), it may be required the installation of a package containing "libgfortran.so.4", for example gfortran-7 (suggested command for installation: sudo apt install gfortran-7).

## License

This program is distributed under a [LGPL-3.0](./LICENSE), however delphi4py depends on
DelPhi which is proprietary. To use DelPhi the user is required to
download the DelPhi license
[here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

## Contacts

Please submit a github issue to report bugs and to request new features.
Alternatively you may find the developer [here](mailto:jnvitorino@fc.ul.pt).
