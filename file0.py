from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import subprocess
import math
import mdtraj as md
import numpy as np
import math
import pickle
from Bio.PDB import PDBParser
from Bio import PDB
import pandas as pd
import itertools

alldis1 = []
alldis2 = []
alldis3 = []
alldis4 = []
alldis5 = []


# fileout0 = open('dist_C2_ave_fluc.dat', 'w')
nres = []
resname = []


def file0(traj):

    with open("list_structures.temp") as file:
        lines = file.readlines()
    for i in lines:
        parser = PDBParser(PERMISSIVE=1)
        traj = i.split()[1]
        structure = parser.get_structure("lig", traj)

        residues = [r for r in structure.get_residues()]

        for i in range(len(residues)):
            i = i+1
        nres.append(i)
        for i in residues:
            name = i.get_resname()
            resname.append(name)

        for each in itertools.combinations(residues, 2):
            seq1 = each[0].get_id()[1]
            seq2 = each[1].get_id()[1]
            one_c = each[0]["C2"].get_coord()
            two_c = each[1]["C2"].get_coord()
            one_c1 = each[0]["C1'"].get_coord()
            #one_P = each[0]["P"].get_coord()
            two_P = each[1]["P"].get_coord()
            one_O = each[0]["O2'"].get_coord()
            two_N = each[1]["N3"].get_coord()
            two_O2 = each[1]["O2'"].get_coord()

            if (seq2-1) == seq1:
                distance_C2 = float('{}'.format(np.linalg.norm(one_c-two_c)))
                distance_C1 = float('{}'.format(np.linalg.norm(one_c1-two_c)))
                #distance_p = float('{}'.format(np.linalg.norm(one_P-two_P)))
                distance_O = float('{}'.format(np.linalg.norm(one_O-two_N)))
                distance_O2 = float('{}'.format(np.linalg.norm(one_O-two_O2)))

                # print(distance)
                alldis1.append(distance_C2)
                alldis2.append(distance_C1)
                # alldis3.append(distance_p)
                alldis4.append(distance_O)
                alldis5.append(distance_O2)

    sep1 = np.array_split(alldis1, len(lines))
    vert1 = np.vstack(sep1)
    data1 = vert1.T

    sep2 = np.array_split(alldis2, len(lines))
    vert2 = np.vstack(sep2)
    data2 = vert2.T

    sep3 = np.array_split(alldis4, len(lines))
    vert3 = np.vstack(sep3)
    data3 = vert3.T

    sep4 = np.array_split(alldis5, len(lines))
    vert4 = np.vstack(sep4)
    data4 = vert4.T

    # print(data)

    np.savetxt('dist_C2test.dat', data1, fmt='%1.2f')
    np.savetxt('dist_C1test.dat', data2, fmt='%1.2f')
    np.savetxt('dist_Otest.dat', data3, fmt='%1.2f')
    np.savetxt('dist_O2test.dat', data4, fmt='%1.2f')


#df = pd.read_csv('dist_C2test.dat', skiprows=0, delimiter="\s+", header=None)
# print(df)
# np.savetxt('dist_C1.dat', bond1[0:ll+1, 0:k1+1]*10, fmt='%4.2f')
# np.savetxt('dist_P.dat', bond2[0:ll+1, 0:k2+1]*10, fmt='%4.2f')
# np.savetxt('dist_OP_O2p_bis.dat', bond3[0:ll+1, 0:k3+1]*10, fmt='%4.2f')
# np.savetxt('dist_O2p_base.dat', bond4[0:ll+1, 0:k4+1]*10, fmt='%4.2f')
# np.savetxt('dist_OP_O2p_min.dat', bond5[0:ll+1, 0:k5+1]*10, fmt='%4.2f')

# fileout0.close()
