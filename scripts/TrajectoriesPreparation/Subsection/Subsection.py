#①combination

import os
import math
rand_path = './TDD_r2p_1_Subsection/'
dele_path = '../Truncation/TDD_r2p_1_dele/'
file_list = os.listdir(dele_path)
dele_file= rand_path + 'r_p_total_dele_P1.xyz'
dele_list = []
for file in file_list:
    traj_dele = dele_path + file
    with open(traj_dele,'r') as fr:
        lines_list = fr.readlines()
        file_len = len(lines_list)
        atom_num = 22
        for i in range(int(file_len/(atom_num+2))):            
            dele_list.append(lines_list[i*(atom_num+2):(i+1)*(atom_num+2)])
with open(dele_file ,'w+') as fw:
    for chunk in dele_list:
        fw.writelines(chunk)    

#②subsection        
atom1_num = 4
atom2_num = 14

bond_len_min = 2.84
bond_len_max = 2.90


lines_list = []
every_point_list = []
orign_path ='./TDD_r2p_1_Subsection/'
dele_path = './TDD_r2p_1_Subsection/'

file_list = os.listdir(orign_path)
traj_file = orign_path + 'r_p_total_dele_P1.xyz'
traj_dele = dele_path  + str(bond_len_max) + "_" + str(bond_len_min) + "_bondlength.xyz"
    #print(traj_file)
with open(traj_file,'r') as fr:
    for line in fr:
        lines_list.append(line)
    file_len = len(lines_list)
    atom_num = eval(lines_list[0].strip())
    #print(atom_num)
    for i in range(int(file_len/(atom_num+2))):  
        #print(i)
        
        #bond(6,13)
        coord1_list = lines_list[i*(atom_num+2)+atom1_num+1].strip().split()
        #coord1 = [eval(coord1_list[1]),eval(coord1_list[2]),eval(coord1_list[3])]
        coord1 = [eval(coord1_list[j]) for j in range(1,4)]
        coord2_list = lines_list[i*(atom_num+2)+atom2_num+1].strip().split()
        #coord2 = [eval(coord2_list[1]),eval(coord2_list[2]),eval(coord2_list[3])]
        coord2 = [eval(coord2_list[j]) for j in range(1,4)]
        length1 = math.sqrt(sum([(coord1[j]-coord2[j])**2 for j in range(3)]))
        #print(length1)
        #P1 or P2
        #coord3_list = lines_list[i*(atom_num+2)+atom3_num+1].strip().split()
        #coord3 = [eval(coord3_list[1]),eval(coord3_list[2]),eval(coord3_list[3])]
        #coord4_list = lines_list[i*(atom_num+2)+atom4_num+1].strip().split()
        #coord4 = [eval(coord4_list[1]),eval(coord4_list[2]),eval(coord4_list[3])]
        #length2 = math.sqrt((coord3[0]-coord4[0])**2 + (coord3[1]-coord4[1])**2 + (coord3[2]-coord4[2])**2)
        if length1 >= bond_len_min and length1 < bond_len_max :
            every_point_list.append(lines_list[i*(atom_num+2):i*(atom_num+2)+atom_num+2])

with open(traj_dele ,'w') as fw:
    for chunk in every_point_list:
        fw.writelines(chunk)  
 
#③points calculation
orign_path = './TDD_r2p_1_Subsection/'
file_list = os.listdir(orign_path)
flag = 0
traj_file= orign_path +  str(bond_len_max) + "_" + str(bond_len_min) + "_bondlength.xyz"
with open(traj_file,'r') as fr:
    lines_list = fr.readlines()
    file_len = len(lines_list)
    atom_num = eval(lines_list[0].strip())
    every_point_list = []
    for i in range(int(file_len/(atom_num+2))):  
        flag+=1
print(flag) 
