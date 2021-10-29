#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
traj_path = '../TrajectoriesPreparation/Subsection/TDD_r2pA_Subsection/'
top_path = ''
outpath = './'
traj_origin = md.load(traj_path + '2.768_2.636_bondlength.xyz', top = top_path +'./dimerization-low-hp-C2-opt.pdb') 
traj_origin[0:10000:1].save_xtc(outpath + '2.768_2.636_bondlength.xtc')


# In[ ]:




