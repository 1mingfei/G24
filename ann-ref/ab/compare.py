#!/usr/bin/python
import numpy as np

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def test_structures(filename):
    with open(filename+'.debug','r') as fin1:
        lines1=fin1.readlines()
        data1=np.zeros((len(lines1),1))
        for i in range(len(lines1)):
            data1[i]=float(lines1[i].split()[0])
    with open('debug.'+str(filename)+'.own','r') as fin2:
        lines2=fin2.readlines()
        data2=np.zeros((len(lines2),1))
        for i in range(len(lines2)):
            data2[i]=float(lines2[i].split()[0])
    
    a=rmse(data1,data2) 
    flag=0
    if (a > 0.05) :
        print "%s failed the test with RMSE: %.04g" %(filename,a)
        flag=0
    else:
        print "%s passed the test with RMSE: %.04g" %(filename,a)
        flag=1
    return flag

filename='ab'
test_structures(filename)






