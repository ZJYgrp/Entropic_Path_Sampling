import os
import math
  
atom1_num = 6
atom2_num = 13
#bond(4,14) or bond(10,12)
atom3_num = 4
atom4_num = 14

orign_path_start = '../Processing/TDD_r2p'
dele_path = './TDD_r2p_1_dele/'
for product in ['A/', 'B/']:
    orign_path = orign_path_start + product
    file_list = os.listdir(orign_path)
    for file in file_list:
        traj_file = orign_path + file
        traj_dele = dele_path  + file
        #print(traj_file)
        with open(traj_file,'r') as fr:
            lines_list = fr.readlines()
            file_len = len(lines_list)
            atom_num = eval(lines_list[0].strip())
            every_point_list = []
            for i in range(int(file_len/(atom_num+2))):   
                #bond(6,13)
                coord1_list = lines_list[i*(atom_num+2)+atom1_num+1].strip().split()
                coord1 = [eval(coord1_list[j]) for j in range(1,4)]
                coord2_list = lines_list[i*(atom_num+2)+atom2_num+1].strip().split()
                coord2 = [eval(coord2_list[j]) for j in range(1,4)]
                length1 = math.sqrt(sum([(coord1[j]-coord2[j])**2 for j in range(3)]))
                #P1 or P2
                coord3_list = lines_list[i*(atom_num+2)+atom3_num+1].strip().split()
                coord3 = [eval(coord3_list[j]) for j in range(1,4)]
                coord4_list = lines_list[i*(atom_num+2)+atom4_num+1].strip().split()
                coord4 = [eval(coord4_list[j]) for j in range(1,4)]
                length2 = math.sqrt(sum([(coord3[j]-coord4[j])**2 for j in range(3)]))
                if length1 <= 5.0:
                    every_point_list.append(lines_list[i*(atom_num+2):(i+1)*(atom_num+2)])
                #print(i)
                #print(every_point_list)
                if length2 <= 1.58:
                     break
        with open(traj_dele ,'w') as fw:
            for chunk in every_point_list:
                fw.writelines(chunk)  
