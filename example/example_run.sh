#!/bin/bash -e

# Welcome!! This an easy-to-understand example of how to run PyBindE.

# Our example system is a B2m dimer where MD simulations were run. 
# Trajectory (.xtc) files were extracted to frames (.gro) where the calculations will be run.

# The program requires as a first argument a .gro file path. 
# Secondly, the atom number for the start and end of each system component.
# You can explore aditional optional inputs using the flag --help or -h.  Example: pybinde_2021-09.py --help


# This next line should work to calculate energies, for a single frame, with a dielectric constant of 4, and saving them to a "energy_file.txt" file under the "saved_files" folder.
python3 ../pybinde_2021-09.py frame00100.gro    1 1149 1150 2298        -ep 4    -sav saved_files/   -ensav saved_files/energy_file.txt
#         call to program    | gro file      | 2 ranges for objects   | epsin |   pdb saving path  |   energy file path
