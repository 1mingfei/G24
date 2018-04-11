#!/usr/bin/python
import numpy as np

# get xsf file format from a cartesian input
def get_xsf_cart_in(flag,cell,data,tot_num,file_out) :
    fout=open(file_out,'w')
    fout.write('CRYSTAL\nPRIMVEC\n')
    for i in range(3):
        for j in range(3):
            fout.write( '%12.8f  ' % cell[i][j] )
        fout.write('\n')
    fout.write(str('PRIMCOORD\n'+str( tot_num)+' 1\n'))
    if flag=='F':
        tmp=np.dot(data[:,1:4],cell)
        for i in range(tot_num):
            fout.write('%s  %12.8f  %12.8f  %12.8f\n'%(data[i][0],tmp[i][0],tmp[i][1],tmp[i][2]))  
    elif flag=='C':
        for i in range(tot_num):
            fout.write('%s  %12.8f  %12.8f  %12.8f\n'%(data[i][0],data[i][1],data[i][2],data[i][3]))  

    fout.close()
    return


