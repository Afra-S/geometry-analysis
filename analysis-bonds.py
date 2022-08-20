from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math
from Bio.PDB import PDBParser


type = input("Enter type of trajectory [pdb/xtc] ")
if(type == 'pdb'):
    with open("list_structures.temp") as file:
        lines = file.readlines()
        for i in lines:
            #parser = PDBParser(PERMISSIVE=1)
            traj = i.split()[1]
    #structure = parser.get_structure("lig", traj)

    #parser = PDBParser(PERMISSIVE=1)

    #structure = parser.get_structure("lig", traj)
        import file0
        file0.file0(traj)
elif(type == 'xtc'):
    name = input("Enter name trajctory: ")
    topname = input("Enter name topology with extension: ")
    traj = md.load(name+'.xtc', top=topname)
    import bonds_afra
    bonds_afra.bonds(traj)
    topology = traj.topology
    nres = traj.n_residues
    print(traj)
    print(nres)
    print(topology)
