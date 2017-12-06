#!/usr/bin/python
import numpy as np
filename='abcbc_23.nb'

def mv_atom_center(cell,data,tot_num,i_n):  # i_n for the atom need to move center
    tmp = cell[0][0]*0.5
    dis_x = data[i_n][1]-tmp
    tmp = cell[1][1]*0.5
    dis_y = data[i_n][2]-tmp
    tmp = cell[2][2]*0.5
    dis_z = data[i_n][3]-tmp
    for i in range(tot_num):
        data[i][1]-=dis_x
        data[i][2]-=dis_y
        data[i][3]-=dis_z
    return data

def get_xsf_cart_in(cell,data,tot_num,file_out) :
    fout=open(file_out,'w')
    fout.write('CRYSTAL\nPRIMVEC\n')
    for i in range(3):
        for j in range(3):
            fout.write( '%12.8f  ' % cell[i][j] )
        fout.write('\n')
    fout.write(str('PRIMCOORD\n'+str( tot_num)+' 1\n'))

    for i in range(tot_num):
        fout.write('%s  %12.8f  %12.8f  %12.8f\n'%(data[i][0],data[i][1],data[i][2],data[i][3] ))  
    fout.close()
    return
#cell=[[50.0000,0.00000000,0.00000000],
#  [0.00000000,50.00000000,0.00000000],  
#  [0.00000000,0.00000000,50.00000000]]  

cell=[[14.46210000,0.00000000,0.00000000],
  [0.00000000,15.02946000,0.00000000],  
  [0.00000000,0.00000000,23.61660000]]  
ln_num=[]
with open(filename,'r') as fin:
    pos=fin.tell()
    for num,ln in enumerate(fin):
        if 'neighbors' in ln:
            ln_num.append(num)
    fin.seek(pos)
    lines=fin.readlines()

for i in range(len(ln_num)):
    tot_num=int(lines[ln_num[i]-1].split()[2])+1
    data=np.zeros((tot_num,4))
    for j in range(tot_num):
        for k in range(1,4):
            data[j][k]=float(lines[ln_num[i]+5+j].split()[k])
    mv_atom_center(cell,data,tot_num,0)
    file_out='abcbc_23-'+str(i)+'.xsf'
    get_xsf_cart_in(cell,data,tot_num,file_out) 
