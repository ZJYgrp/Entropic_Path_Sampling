#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import math

#①combination
dele_path = '../Truncation/TDD_r2pA_dele/'
inpath = '../Truncation/TDD_r2pA_dele/'

rand_path = './TDD_r2pA_Subsection/'
outpath = './TDD_r2pA_Subsection/'

file_list = os.listdir(inpath)
combi_file= outpath + 'r_pA_combination.xyz'
combi_list = []
for file in file_list:
    filename = inpath + file
    with open(filename,'r') as fr:
        lines_list = fr.readlines()
        file_len = len(lines_list)
        atom_num = eval(lines_list[0].strip())
        for i in range(int(file_len/(atom_num+2))):            
            combi_list.append(lines_list[i*(atom_num+2):(i+1)*(atom_num+2)])
with open(combi_file ,'w+') as fw:
    for chunk in combi_list:
        fw.writelines(chunk)    

#②subsection 
#bond(4,14) or bond(10,12): selective bond that forms in product A or B
atom1_num = 4
atom2_num = 14
#range of bond length
bond_len_max = 2.900
bond_len_min = 2.768



lines_list = []
every_point_list = []
sub_path = './TDD_r2pA_Subsection/'

sub_file = sub_path  + str(bond_len_max) + "_" + str(bond_len_min) + "_bondlength.xyz"

with open(combi_file,'r') as fr:
    for line in fr:
        lines_list.append(line)
    file_len = len(lines_list)
    atom_num = eval(lines_list[0].strip())
    for i in range(int(file_len/(atom_num+2))):  
        #bond(4,14): selective bond that forms in product A or B
        coord1_list = lines_list[i*(atom_num+2)+atom1_num+1].strip().split()
        coord1 = [eval(coord1_list[j]) for j in range(1,4)]
        coord2_list = lines_list[i*(atom_num+2)+atom2_num+1].strip().split()
        coord2 = [eval(coord2_list[j]) for j in range(1,4)]
        length1 = math.sqrt(sum([(coord1[j]-coord2[j])**2 for j in range(3)]))

        if length1 >= bond_len_min and length1 < bond_len_max :
            every_point_list.append(lines_list[i*(atom_num+2):i*(atom_num+2)+atom_num+2])

with open(sub_file ,'w') as fw:
    for chunk in every_point_list:
        fw.writelines(chunk)  
        
#③points calculation
flag = 0
with open(sub_file,'r') as fr:
    lines_list = fr.readlines()
    file_len = len(lines_list)
    atom_num = eval(lines_list[0].strip())
    every_point_list = []
    for i in range(int(file_len/(atom_num+2))):  
        flag+=1
print(flag) 


# In[ ]:




