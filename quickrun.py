from Bio.PDB import *
#import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from ring_atom_name_getter import *

#-alternate_3_letter_codes pdb_sugar
#-out:level 100

#-beta
#-auto_detect_glycan_connections
#-alternate_3_letter_codes pdb_sugar



import os
import numpy as np
import pandas as pd
import copy

import gzip
import shutil

BOND_CUT = 1.85
INTERACT = 4.5

in_dir = "./diffdock_xtal/fix/"
#in_dir = "./af3_diffdock/align/"
in_dir = "./rfaa_out/fix/"
#in_dir = "./chai_out/fix/"
xtal_dir = './diff_combo/'
ls = os.listdir(in_dir)
ls.sort()

s = 'PDB,dice,fnat,RiRMS,lRMS,prot-carb_clash,carb-carb_clash\n'

a1 = 1.5
a2 = 8.5

ope = []
for ii in ls:
    if 'DS_' in ii:
        continue;
    if '.pdb' not in ii:
        continue;
    #if '_' in ii:
    #    continue;


    #ii = '7TOH_rfaa.pdb'
    #print(ii)
    if '_rfaa' in ii:
        continue;

    sl = True
    if '_' in ii:
        sl = False
    sl = False

    #try:
    if True:
        d, f, rirms, lrms, ab_clash,aa_clash = calc_metrics(in_dir + ii,xtal_dir + 'good/' + ii[:4] + '.pdb',same_ligand=sl,is_align=True,is_same_num=False)

        s += str(ii) + ',' + str(round(d,3)) + ',' + str(round(f,3)) + ',' + str(round(rirms,3)) + ',' + str(round(lrms,3)) + ',' + str(ab_clash) + ',' + str(aa_clash) + '\n'


    #except:
    #    ope.append(ii)

    #break
print(ope)
print(s)
