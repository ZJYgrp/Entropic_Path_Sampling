import os
import pandas as pd
from glob import glob
def Total_Entropy(filename):
    total = []
    with open (filename,'r') as fr:
        line_list = fr.readlines()
        for line in line_list:
            words = line.split()
            #print(words)
            if words[0] == 'TOTAL' and words[1] == 'CONFIGURATIONAL':
                words[4]= float(words[4])
                file = filename.split('/')
                #total_entropy = '%s\t%s\t%s\t%s : %.2f'%(file[1],words[0],words[1],words[2],words[4])
                total_entropy = '%s\t : %.2f'%(file[1],words[4])
                #print(total_entropy)
    MI_total_entropy = outpath + 'MI_total_entropy.txt'
#make new file called total_outputs,and move PARENT-master_MIST.txt in every output file into total_ouputs file.
inpath = './PARENT-master/total_outputs/'
outpath = './PARENT-master/output_analysis/'

# Reset outpath for new analysis files
os.system('rm ' + outpath + '*')

for filename in glob(inpath+'*.txt'):
    Total_Entropy(filename) 

D1_D2_conf_file = outpath + 'D1_D2_conf.txt'
file_list = os.listdir(inpath)
for file in file_list:
    D2_MI_total = 0
    D1_entropy_total = 0
    conf_entropy_total = 0
    filename = inpath +file
    with open (filename,'r') as fr:
        line_list = fr.readlines()
        for line in line_list:
            words = line.split()
            #print(words)
            if words[1] == '1D' :
                D1_entropy= float(words[5])
                D1_entropy_total = D1_entropy + D1_entropy_total
            if words[1] == '2D' and words[3] =='MUTUAL':
                D2_MI= float(words[6])
                D2_MI_total = D2_MI + D2_MI_total
                file = filename.split('/')
            if words[1] == 'CONFIGURATIONAL':
                conf_entropy= float(words[4])
                conf_entropy_total = conf_entropy + conf_entropy_total
                #total_entropy = '%s\t%s\t%s\t%s : %.2f'%(file[1],words[0],words[1],words[2],words[4])
        D1_D2_conf_total  = "%s D1_entopy_total %.2f D2_MI_total %.2f conf_entropy_total %.2f" %(file[1],float(D1_entropy_total),float(D2_MI_total),float(conf_entropy_total))
    with open (D1_D2_conf_file,'a') as fw:
        fw.writelines(D1_D2_conf_total)
        fw.writelines('\n') 
title = []
D1_entropy = []
D2_MI = []
conf_entropy = []
with open (D1_D2_conf_file,'r') as fr:
    line_list = fr.readlines()
    for line in line_list:
        words = line.split()
        title.append(words[0])
        D1_entropy.append(words[2])
        D2_MI.append(words[4])
        conf_entropy.append(words[6])
data = {"total_outputs":title,"D1_entropy": D1_entropy,"D2_MI": D2_MI,"conf_entropy":conf_entropy}
f1 = pd.DataFrame(data,columns = ["total_outputs","D1_entropy","D2_MI","conf_entropy"])     
f1.to_excel(outpath + 'total_outputs.xls',sheet_name='data',index=None)         
