
# This program

A python program for Molecular Mechanics Poisson-Boltzmann Surface Area (MMPBSA) calculations in protein-protein and protein-ligand systems.


## Basic Usage

  The main program is named mmpbsa_protocol.py and receives as an input argument a .gro file (a simple the coordinate file for a GROMACS simulation).

  In order to correctly work, the program requires:
  - 3 databases in order to calculate Van der Waals and Electrostactic energies. An example of these is provided on the 'databases' folder.
  - 1 file named 'classifier.config' with atomic radii, to calculate the Solvent Acessible Surface Area (SASA).
  - A working DelPhi4Py program, which is provided in the folder 'delphi4github'.

  Several settings may be easily changed or adjusted on the main program's file.
  The current settings are specific for the test system provided, (named 'frame00100.gro'), which is a MD simulation frame for a Beta-2-Microglobin (B2M) dimer.
  
  Different systems may have amino acid residues not found on the provided databases, missing atomic radii and charges in the SASA classifier file and in the DelPhi4Py databases.
  These can be manually added on the appropriate files, taking care to respect each file's specific formatting.

## Installation & Dependencies

This program does not require installation. However, given its dependencies, it should be used in a Linux-based system.
The dependencies for this program are:

* python2.6>= & python3.8>=
* numpy
* scipy
* pandas
* freesasa


## License

  This program is distributed under a [LGPL-3.0](./LICENSE), however delphi4py depends on
  DelPhi which is proprietary. To use DelPhi the user is required to
  download the DelPhi license
  [here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

## Contacts

  Please submit a github issue to report bugs and to request new features.
  Alternatively you may find the developer [here](mailto:jnvitorino@fc.ul.pt).