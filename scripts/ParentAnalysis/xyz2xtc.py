#xyz2xtc
from __future__ import print_function
import mdtraj as md
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
#from matplotlib.colors import colorConverter
#from sklearn.decomposition import PCA
#from itertools import combinations
#from msmbuilder.cluster import MiniBatchKMedoids
import numpy as np
traj_path = '../TrajectoriesPreparation/Subsection/TDD_r2p_1_Subsection/'
top_path = ''
path = './'
traj_origin = md.load(traj_path + '2.84_2.73_bondlength.xyz', top = top_path +'./dimerization-low-hp-C2-opt.pdb') 
traj_origin[0:10000:1].save_xtc(path + '2.84_2.73_bondlength.xtc')
