#!/usr/bin/python3
# version 2021.09.25

import time
from pybinde_classes import pybinde
import argparse
import os

initial_time = time.perf_counter()

parser = argparse.ArgumentParser(description='====< Calculates MM-PBSA for a given .gro file >====')
parser.add_argument('gro',type=str,metavar='gro', help='input (.gro) file path')
parser.add_argument('atom_ranges',type=int, metavar='atom_ranges', help='Range of atom nrs that defines start and end of obj1 and obj2. Example: 1 1901 1902 1973', nargs=4)
parser.add_argument('-dbs','--databases',type=str, metavar='',required=False, help='(opt.) path for all databases')


parser.add_argument('-db_crg','--db_charges',type=str, default=os.path.dirname(os.path.realpath(__file__)) + "/databases/DataBaseT.crg", metavar='',required=False, help='(opt.) charges database path')
parser.add_argument('-db_siz','--db_radii',type=str, default=os.path.dirname(os.path.realpath(__file__)) + "/databases/DataBaseT.siz", metavar='',required=False, help='(opt.) radii database path')
parser.add_argument('-db_rtypes','--db_restypes',type=str, default=os.path.dirname(os.path.realpath(__file__)) + "/databases/residues_rtp.txt", metavar='',required=False, help='(opt.) residues_rtp database path')
parser.add_argument('-db_nb','--db_nonbonded',type=str, default=os.path.dirname(os.path.realpath(__file__)) + "/databases/nonbonded_itp.txt", metavar='',required=False, help='(opt.) nonbonded database path')
parser.add_argument('-db_cls','--db_classifier',type=str, default=os.path.dirname(os.path.realpath(__file__)) + "/databases/classifier.config", metavar='',required=False, help='(opt.) classifier.config databases path')

parser.add_argument('-sav','--pdb_saving_dir',type=str, metavar='',required=False, default='', help='(opt.) saving directory for generated PDB files')
parser.add_argument('-ensav','--energies_saving_path',type=str, metavar='./',required=False, default='energies.txt', help='(opt.) saving path/filename for generated energy files')
parser.add_argument('-ep','--epsin'    ,type=int, metavar='',required=False, default=4, help='(opt.) dielectric constant to be used')
parser.add_argument('-cto','--cutoff'    ,type=float, metavar='',required=False, default=None, help='(opt.) cutoff for the MM calculations')

args = parser.parse_args()


if __name__ == '__main__':
  if args.databases:
    if not args.databases.endswith("/"):
      args.databases += "/"
    args.db_charges = args.databases + "DataBaseT.crg"
    args.db_radii = args.databases + "DataBaseT.siz"
    args.db_restypes = args.databases + "residues_rtp.txt"
    args.db_nonbonded = args.databases + "nonbonded_itp.txt"
    args.db_classifier = args.databases + "classifier.config"


  pybinde_object = pybinde(args.gro,args.atom_ranges,args.epsin,args.db_charges,args.db_radii,args. db_restypes,args.db_nonbonded,args.db_classifier)

  pybinde_object.read_gro_df()
  pybinde_object.construct_objects()

  pybinde_object.create_pdb(args.pdb_saving_dir)

  pybinde_object.calculate_pair_distances(cutoff_n=None)

  pybinde_object.calculate_MM_energy()
  pybinde_object.calculate_SA_energy()
  pybinde_object.calculate_PB_energy()

  pybinde_object.calculate_binding_energy()

  final_time = time.perf_counter()
  elapsed_time=final_time-initial_time

  print("Running Time (s):                  ", round(elapsed_time,2))
  pybinde_object.save_results(args.energies_saving_path,elapsed_time)