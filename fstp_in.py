#!/usr/bin/python
import numpy as np
import linecache

def get_fstp(filename):
    g=[]   # to store g2 g4 function parameters
    with open(filename,'r') as fin:
        for num,line in enumerate(fin,1):
            if 'SYMMFUNC' in line:
                myline = linecache.getline(filename, (num+1))
                tot_num=int(myline.split()[0])
                g_tmp=[]
                for i in range(tot_num):
                    myline = linecache.getline(filename, (num+2+i))
                    tmp=myline.split()
                    for j in range(len(tmp)):
                        g_tmp.append(tmp[j].split('=')[-1])
                    g.append(g_tmp)
                    g_tmp=[]
    return g

