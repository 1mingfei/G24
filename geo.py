#!/usr/bin/python
import formatxfer as fx
import numpy as np
import outputs as op

#cell_expd_cart(cell,data,tot_num,x,y,z)  expand unit cell by x,y,z
# mv_atom_center(cell,data,tot_num,in)   in for the atom need to move center
# calc_dist_cart(lst,n1,n2,H) Cart calculated distance of n1 and n2
# calc_dist_direct(lst,n1,n2,H) Direct calculated distance of n1 and n2
# nb_lst(cell,data,tot_num,i_n,Rc)
# expd_prmt(cell,Rc=6.5)
#get cosine of theta(angle of j-i-k)
# get_cos_cart(lst,i,j,k,H)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# get_cos_direct(lst,i,j,k,H)  !!!!!!!!!!!!still have some problem

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if vector==[0,0,0]:
        return 1
    else:
        return vector / np.linalg.norm(vector)

#Cartesian inputs gives calculated distance of n1 and n2
def calc_dist_cart(lst,n1,n2,H):
    if  lst[n2][1]-lst[n1][1]>0.5*H[0][0]:
        a=(lst[n1][1]-lst[n2][1]+H[0][0])**2
    elif  lst[n2][1]-lst[n1][1]<(-0.5*H[0][0]):
        a=(lst[n1][1]-lst[n2][1]-H[0][0])**2
    else:
        a=(lst[n1][1]-lst[n2][1])**2

    if  lst[n2][2]-lst[n1][2]>0.5*H[1][1]:
        b=(lst[n1][2]-lst[n2][2]+H[1][1])**2
    elif  lst[n2][2]-lst[n1][2]<(-0.5*H[1][1]):
        b=(lst[n1][2]-lst[n2][2]-H[1][1])**2
    else:
        b=(lst[n1][2]-lst[n2][2])**2

    if  lst[n2][3]-lst[n1][3]>0.5*H[2][2]:
        c=(lst[n1][3]-lst[n2][3]+H[2][2])**2
    elif  lst[n2][3]-lst[n1][3]<(-0.5*H[2][2]):
        c=(lst[n1][3]-lst[n2][3]-H[2][2])**2
    else:
        c=(lst[n1][3]-lst[n2][3])**2

    out=np.sqrt(a+b+c)
    return out

#fractional inputs gives calculated distance of n1 and n2
def calc_dist_direct(lst,n1,n2,H):
    if  lst[n2][1]-lst[n1][1]>=0.5:
        a=(lst[n1][1]-lst[n2][1]+1)
    elif  lst[n2][1]-lst[n1][1]<(-0.5):
        a=(lst[n1][1]-lst[n2][1]-1)
    else:
        a=((lst[n1][1]-lst[n2][1]))

    if  lst[n2][2]-lst[n1][2]>=0.5:
        b=((lst[n1][2]-lst[n2][2]+1))
    elif  lst[n2][2]-lst[n1][2]<(-0.5):
        b=((lst[n1][2]-lst[n2][2]-1))
    else:
        b=((lst[n1][2]-lst[n2][2]))

    if  lst[n2][3]-lst[n1][3]>=0.5:
        c=((lst[n1][3]-lst[n2][3]+1))
    elif  lst[n2][3]-lst[n1][3]<(-0.5):
        c=((lst[n1][3]-lst[n2][3]-1))
    else:
        c=((lst[n1][3]-lst[n2][3]))
    tmp=[a,b,c]
    tmp_cart=np.dot(tmp,H)
    aa,bb,cc=tmp_cart[0],tmp_cart[1],tmp_cart[2]

    return np.sqrt(aa**2+bb**2+cc**2)

#get cell expand parameters 
def expd_prmt(cell,Rc=6.5):
    #if cell[0][0]> (2*Rc+0.5) : x=0
    if cell[0][0]> (2*Rc+0.5) : x=1
    else:
        x= int((2*Rc+0.5)/cell[0][0])+1
    #if cell[1][1]> (2*Rc+0.5) : y=0
    if cell[1][1]> (2*Rc+0.5) : y=1
    else:
        y= int((2*Rc+0.5)/cell[1][1])+1
    #if cell[2][2]> (2*Rc+0.5) : z=0
    if cell[2][2]> (2*Rc+0.5) : z=1
    else:
        z= int((2*Rc+0.5)/cell[2][2])+1
    return [x,y,z]

#get neighbor list for i_n within dist Rc
def nb_lst(cell,data,tot_num,i_n,Rc):
    lst=[]
    for i in range(tot_num):
        a = calc_dist_cart(data,i_n,i,cell) 
        if a <= Rc and (a > 0.5) and (i != i_n):
            lst.append([i,a])
    return lst


# fractional inputs cos_theta for atom i with j and k
def get_cos_direct(lst,i,j,k,H):
    if  lst[j][1]-lst[i][1]>=0.5:
        a=(lst[i][1]-lst[j][1]+1)
    elif  lst[j][1]-lst[i][1]<(-0.5):
        a=(lst[i][1]-lst[j][1]-1)
    else:
        a=((lst[i][1]-lst[j][1]))
    if  lst[j][2]-lst[i][2]>=0.5:
        b=((lst[i][2]-lst[j][2]+1))
    elif  lst[j][2]-lst[i][2]<(-0.5):
        b=((lst[i][2]-lst[j][2]-1))
    else:
        b=((lst[i][2]-lst[j][2]))
    if  lst[j][3]-lst[i][3]>=0.5:
        c=((lst[i][3]-lst[j][3]+1))
    elif  lst[j][3]-lst[i][3]<(-0.5):
        c=((lst[i][3]-lst[j][3]-1))
    else:
        c=((lst[i][3]-lst[j][3]))
    tmp=[a,b,c]
    print(tmp)
    print(H)
    aa=np.dot(tmp,H)
    print(aa)
    if  lst[k][1]-lst[i][1]>=0.5:
        a=(lst[i][1]-lst[k][1]+1)
    elif  lst[k][1]-lst[i][1]<(-0.5):
        a=(lst[i][1]-lst[k][1]-1)
    else:
        a=((lst[i][1]-lst[k][1]))
    if  lst[k][2]-lst[i][2]>=0.5:
        b=((lst[i][2]-lst[k][2]+1))
    elif  lst[k][2]-lst[i][2]<(-0.5):
        b=((lst[i][2]-lst[k][2]-1))
    else:
        b=((lst[i][2]-lst[k][2]))
    if  lst[k][3]-lst[i][3]>=0.5:
        c=((lst[i][3]-lst[k][3]+1))
    elif  lst[k][3]-lst[i][3]<(-0.5):
        c=((lst[i][3]-lst[k][3]-1))
    else:
        c=((lst[i][3]-lst[k][3]))
    tmp1=[a,b,c]
    bb=np.dot(tmp1,H)
    print(bb)

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#    return np.dot(aa,bb)/np.linalg.norm(aa)/np.linalg.norm(bb)

def get_cos_cart(lst,i,j,k,H):
    if  lst[j][1]-lst[i][1]>0.5*H[0][0]:
        a=(lst[i][1]-lst[j][1]+H[0][0])
    elif  lst[j][1]-lst[i][1]<(-0.5*H[0][0]):
        a=(lst[i][1]-lst[j][1]-H[0][0])
    else:
        a=(lst[i][1]-lst[j][1])
    if  lst[j][2]-lst[i][2]>0.5*H[1][1]:
        b=(lst[i][2]-lst[j][2]+H[1][1])
    elif  lst[j][2]-lst[i][2]<(-0.5*H[1][1]):
        b=(lst[i][2]-lst[j][2]-H[1][1])
    else:
        b=(lst[i][2]-lst[j][2])
    if  lst[j][3]-lst[i][3]>0.5*H[2][2]:
        c=(lst[i][3]-lst[j][3]+H[2][2])
    elif  lst[j][3]-lst[i][3]<(-0.5*H[2][2]):
        c=(lst[i][3]-lst[j][3]-H[2][2])
    else:
        c=(lst[i][3]-lst[j][3])
    aa=[a,b,c]
    #print(aa)
    if  lst[k][1]-lst[i][1]>0.5*H[0][0]:
        a=(lst[i][1]-lst[k][1]+H[0][0])
    elif  lst[k][1]-lst[i][1]<(-0.5*H[0][0]):
        a=(lst[i][1]-lst[k][1]-H[0][0])
    else:
        a=(lst[i][1]-lst[k][1])

    if  lst[k][2]-lst[i][2]>0.5*H[1][1]:
        b=(lst[i][2]-lst[k][2]+H[1][1])
    elif  lst[k][2]-lst[i][2]<(-0.5*H[1][1]):
        b=(lst[i][2]-lst[k][2]-H[1][1])
    else:
        b=(lst[i][2]-lst[k][2])

    if  lst[k][3]-lst[i][3]>0.5*H[2][2]:
        c=(lst[i][3]-lst[k][3]+H[2][2])
    elif  lst[k][3]-lst[i][3]<(-0.5*H[2][2]):
        c=(lst[i][3]-lst[k][3]-H[2][2])
    else:
        c=(lst[i][3]-lst[k][3])
    bb=[a,b,c]
    #print(bb)
    v1_u = unit_vector(aa)
    v2_u = unit_vector(bb)
    return (np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    #return np.dot(aa,bb)/float(np.linalg.norm(aa))/float(np.linalg.norm(bb))

def geo(filename,ftype,n_spe):
    ab=fx.info(filename,ftype,n_spe)
    tot_num=ab.tot_num
    cell= ab.cell
    data= ab.data
    return cell,data,tot_num

def cell_expd_frac_2sides(cell,data,tot_num,xyz):
    #expand x cells in both +/- x-direction
    #expand y cells in both +/- y-direction
    #expand z cells in both +/- z-direction
    x,y,z=int(xyz[0]),int(xyz[1]),int(xyz[2])
    #data_new=np.zeros(int(x*y*z*tot_num) , dtype='S2,float,float,float')
    data_new=np.zeros((int((2*x+1)*(2*y+1)*(2*z+1)*tot_num),4))
    #print(data)
    for i in range(tot_num):
        data_new[i][0]=data[i][0]
        data_new[i][1]=data[i][1]/(2*x+1)+x/(2*x+1)
        data_new[i][2]=data[i][2]/(2*y+1)+y/(2*y+1)
        data_new[i][3]=data[i][3]/(2*z+1)+z/(2*z+1)
    for ix in range(tot_num,tot_num*(x+1)):
        data_new[ix][0]=data_new[ix%tot_num][0]
        tmp=ix//(tot_num)
        data_new[ix][1]=data_new[ix%tot_num][1]+tmp/(2*x+1)
        data_new[ix][2]=data_new[ix%tot_num][2]
        data_new[ix][3]=data_new[ix%tot_num][3]
    for ix in range(tot_num*(x+1),tot_num*(2*x+1)):
        data_new[ix][0]=data_new[ix%tot_num][0]
        tmp=ix//(tot_num)-x
        data_new[ix][1]=data_new[ix%tot_num][1]-(tmp)/(2*x+1)
        data_new[ix][2]=data_new[ix%tot_num][2]
        data_new[ix][3]=data_new[ix%tot_num][3]
    for iy in range(tot_num*(2*x+1),tot_num*(2*x+1)*(y+1)):
        data_new[iy][0]=data_new[iy%(tot_num*(2*x+1))][0]
        data_new[iy][1]=data_new[iy%(tot_num*(2*x+1))][1]
        tmp=iy//(tot_num*(2*x+1))
        data_new[iy][2]=data_new[iy%(tot_num*(2*x+1))][2]+(tmp)/(2*y+1)
        data_new[iy][3]=data_new[iy%(tot_num*(2*x+1))][3]
    for iy in range(tot_num*(2*x+1)*(y+1),tot_num*(2*x+1)*(2*y+1)):
        data_new[iy][0]=data_new[iy%(tot_num*(2*x+1))][0]
        data_new[iy][1]=data_new[iy%(tot_num*(2*x+1))][1]
        tmp=iy//(tot_num*(2*x+1))-y
        data_new[iy][2]=data_new[iy%(tot_num*(2*x+1))][2]-(tmp)/(2*y+1)
        data_new[iy][3]=data_new[iy%(tot_num*(2*x+1))][3]
    for iz in range(tot_num*(2*x+1)*(2*y+1),tot_num*(2*x+1)*(2*y+1)*(z+1)):
        data_new[iz][0]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][0]
        data_new[iz][1]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][1]
        data_new[iz][2]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][2]
        tmp=iz//(tot_num*(2*x+1)*(2*y+1))
        data_new[iz][3]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][3]+(tmp)/(2*z+1)
    for iz in range(tot_num*(2*x+1)*(2*y+1)*(z+1),tot_num*(2*x+1)*(2*y+1)*(2*z+1)):
        data_new[iz][0]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][0]
        data_new[iz][1]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][1]
        data_new[iz][2]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][2]
        tmp=iz//(tot_num*(2*x+1)*(2*y+1))-z
        data_new[iz][3]=data_new[iz%(tot_num*(2*x+1)*(2*y+1))][3]-(tmp)/(2*z+1)
    cell[0][0]*=(2*x+1)
    cell[0][1]*=(2*x+1)
    cell[0][2]*=(2*x+1)
    cell[1][0]*=(2*y+1)
    cell[1][1]*=(2*y+1)
    cell[1][2]*=(2*y+1)
    cell[2][0]*=(2*z+1)
    cell[2][1]*=(2*z+1)
    cell[2][2]*=(2*z+1)
    tot_num*=((2*x+1)*(2*y+1)*(2*z+1))
    return cell,data_new,tot_num


def cell_expd_frac(cell,data,tot_num,xyz):
    #expand x cells in x-direction
    #expand y cells in y-direction
    #expand z cells in z-direction
    x,y,z=int(xyz[0]),int(xyz[1]),int(xyz[2])
    #data_new=np.zeros(int(x*y*z*tot_num) , dtype='S2,float,float,float')
    data_new=np.zeros((int(x*y*z*tot_num),4))
    print(data)
    for i in range(tot_num):
        data_new[i][0]=data[i][0]
        data_new[i][1]=data[i][1]/x
        data_new[i][2]=data[i][2]/y
        data_new[i][3]=data[i][3]/z
    for ix in range(tot_num,tot_num*x):
        data_new[ix][0]=data_new[ix%tot_num][0]
        data_new[ix][1]=data_new[ix%tot_num][1]+(ix//tot_num)/x
        data_new[ix][2]=data_new[ix%tot_num][2]
        data_new[ix][3]=data_new[ix%tot_num][3]
    for iy in range(tot_num*x,tot_num*x*y):
        data_new[iy][0]=data_new[iy%(tot_num*x)][0]
        data_new[iy][1]=data_new[iy%(tot_num*x)][1]
        data_new[iy][2]=data_new[iy%(tot_num*x)][2]+(iy//(tot_num*x))/y
        data_new[iy][3]=data_new[iy%(tot_num*x)][3]
    for iz in range(tot_num*x*y,tot_num*x*y*z):
        data_new[iz][0]=data_new[iz%(tot_num*x*y)][0]
        data_new[iz][1]=data_new[iz%(tot_num*x*y)][1]
        data_new[iz][2]=data_new[iz%(tot_num*x*y)][2]
        data_new[iz][3]=data_new[iz%(tot_num*x*y)][3]+(iz//(tot_num*x*y))/z
    print(data_new)
    cell[0][0]*=x
    cell[0][1]*=x
    cell[0][2]*=x
    cell[1][0]*=y
    cell[1][1]*=y
    cell[1][2]*=y
    cell[2][0]*=z
    cell[2][1]*=z
    cell[2][2]*=z
    tot_num*=(x*y*z)
    return cell,data_new,tot_num

def cell_expd_cart(cell,data,tot_num,xyz):
    #expand total x, y and z in each direction
    #expand x cells in x-direction
    #expand y cells in y-direction
    #expand z cells in z-direction
    x,y,z=int(xyz[0]),int(xyz[1]),int(xyz[2])
    #data_new=np.zeros(int(x*y*z*tot_num) , dtype='S2,float,float,float')
    data_new=np.zeros((int(x*y*z*tot_num),4))
    print(data_new.shape)
    for i in range(tot_num):
        data_new[i][0]=data[i][0]
        data_new[i][1]=data[i][1]
        data_new[i][2]=data[i][2]
        data_new[i][3]=data[i][3]
    for ix in range(tot_num,tot_num*x):
        data_new[ix][0]=data_new[ix%tot_num][0]
        data_new[ix][1]=data_new[ix%tot_num][1]+cell[0][0]*(ix/(tot_num))
        data_new[ix][2]=data_new[ix%tot_num][2]
        data_new[ix][3]=data_new[ix%tot_num][3]
    for iy in range(tot_num*x,tot_num*x*y):
        data_new[iy][0]=data_new[iy%(tot_num*x)][0]
        data_new[iy][1]=data_new[iy%(tot_num*x)][1]
        data_new[iy][2]=data_new[iy%(tot_num*x)][2]+cell[1][1]*(iy/(tot_num*x))
        data_new[iy][3]=data_new[iy%(tot_num*x)][3]
    for iz in range(tot_num*x*y,tot_num*x*y*z):
        data_new[iz][0]=data_new[iz%(tot_num*x*y)][0]
        data_new[iz][1]=data_new[iz%(tot_num*x*y)][1]
        data_new[iz][2]=data_new[iz%(tot_num*x*y)][2]
        data_new[iz][3]=data_new[iz%(tot_num*x*y)][3]+cell[2][2]*(iz/(tot_num*x*y))

    cell[0][0]+=(cell[0][0]*(x-1))
    cell[0][1]+=(cell[0][1]*(x-1))
    cell[0][2]+=(cell[0][2]*(x-1))
    cell[1][0]+=(cell[1][0]*(y-1))
    cell[1][1]+=(cell[1][1]*(y-1))
    cell[1][2]+=(cell[1][2]*(y-1))
    cell[2][0]+=(cell[2][0]*(z-1))
    cell[2][1]+=(cell[2][1]*(z-1))
    cell[2][2]+=(cell[2][2]*(z-1))
    tot_num*=(x*y*z)
    return cell,data_new,tot_num

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

def mv_atom_center_frac(cell,data,tot_num,i_n):  # i_n for the atom need to move center
    dis_x = data[i_n][2]-0.5
    dis_y = data[i_n][2]-0.5
    dis_z = data[i_n][3]-0.5
    for i in range(tot_num):
        data[i][1]-=dis_x
        data[i][2]-=dis_y
        data[i][3]-=dis_z
    return data


'''
#check if atom coordinates is out of boundary (not finished yet)
def check_boundary(cell,data,tot_num): 
    for i in range(tot_num):
        while 
        tmp = int(data[i][1]/ cell[0][0])
'''






#testing use only comment all afterwards
a=geo('CONTCAR','vasp',1)
#pm=expd_prmt(a[0],6.5)
#b=cell_expd_frac(a[0],a[1],a[2],pm)
b=cell_expd_frac_2sides(a[0],a[1],a[2],[1,1,3])
#c=mv_atom_center_frac(b[0],b[1],b[2],0)
#lst1=nb_lst(b[0],c,b[2],1,6.5)
#print(lst1)
#op.get_xsf_cart_in('F',a[0],a[1],a[2],'tmp.xsf')
op.get_xsf_cart_in('F',b[0],b[1],b[2],'tmp.xsf')

#print(get_cos_direct(a[1],1,0,2,a[0]))  wrong!!!
#print(get_cos_cart(a[1],1,0,2,a[0]))

#print(get_cos_direct(c,133,170,174,b[0]))
#print(np.arccos(get_cos_cart(c,0,0,1,b[0])) /float(np.pi) *  180.0)



