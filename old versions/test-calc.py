
# from Bio import PDB
# from Bio.PDB.SASA import ShrakeRupley
# # import numpy as np
# # import os

# parser = PDB.PDBParser()
# sasa_a=[]
# a=parser.get_structure("Cyclodextrine_ligand_names_good.pdb","Cyclodextrine_ligand_names_good.pdb")

# for struct in a:
#   sr = ShrakeRupley()
#   sr.compute(struct, level="S")
#   sasa_a.append(round(struct.sasa, 2))

# print(sasa_a)

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
p = PDBParser(QUIET=1)
# This assumes you have a local copy of 1LCD.pdb in a directory called "PDB"
struct = p.get_structure("Cyclodextrine_ligand_names_good.pdb","Cyclodextrine_ligand_names_good.pdb")
sr = ShrakeRupley()
sr.compute(struct, level="S")
print(round(struct.sasa, 2))