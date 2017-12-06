#!/usr/bin/python
import geo
import Gfx
import fstp_in
import numpy as np
import outputs as op

# original : a ; after expand : b ; after move select atom to the center : c ;
name_in='ab'
name=name_in+'.xsf'
a=geo.geo(name)

pm=geo.expd_prmt(a[0],6.5)
b=geo.cell_expd_cart(a[0],a[1],a[2],pm)

fpt1=fstp_in.get_fstp('Al.fingerprint.stp')



score = 0.0

#!!! examples to remind myself
#!!! G=2 type2=Al  eta=1.428426  Rs=0.0000  Rc=6.5000
#!!! G=4 type2=Al type3=Al   eta=0.000357 lambda=-1.0  zeta=1.0  Rc=6.5000
#!!! examples end

for i_i in range(a[2]):
    c=geo.mv_atom_center(b[0],b[1],b[2],i_i)
    lst1=geo.nb_lst(b[0],c,b[2],i_i,6.5)
    d=[] # to store nb within Rc
    print i_i,len(lst1)

    for i in range(len(fpt1)):
        score = 0.0


        if fpt1[i][0]=='2' :
            type2 = str(fpt1[i][1])
            eta   = float(fpt1[i][2])
            Rs    = float(fpt1[i][3])
            Rc    = float(fpt1[i][4])

            #tmp test use
            #Rij=3.0
            #score += Gfx.G2_func(Rij,Rc,Rs,eta)
            #print score
            #score = 0.0
    
            #going seriously
            for i_j in range(len(lst1)):
                Rij = lst1[i_j][1]
                score += Gfx.G2_func(Rij,Rc,Rs,eta)
                d.append(c[lst1[i_j][0]])

        elif fpt1[i][0]=='4' :
            type2 = str(fpt1[i][1])
            type3 = str(fpt1[i][2])
            eta   = float(fpt1[i][3])
            lamda = float(fpt1[i][4])
            zeta  = float(fpt1[i][5])
            Rc    = float(fpt1[i][6])
            #print type2,type3,eta,lamda,zeta,Rc
    
            score = 0.0
            #tmp test use
            #Rij=2.0
            #Rik=3.0
            #Rjk=1.9
            #cos_theta=0.3
    
            #going seriously

            for i_j in range(len(lst1)):
                for i_k in range(i_j+1,len(lst1)):
                    if i_j <> i_k:
                        Rij = lst1[i_j][1]
                        Rik = lst1[i_k][1]
                        Rjk = geo.calc_dist_cart(c,lst1[i_j][0],lst1[i_k][0],b[0])
                        cos_theta=geo.get_cos_cart(c,i_i,lst1[i_j][0],lst1[i_k][0],b[0])
                        score += Gfx.G4_func(Rij,Rik,Rjk,Rc,cos_theta,eta,lamda,zeta)      
        print score
        score = 0.0

    #name_tmp=str(name_in)+'-'+str(i_i)+'.xsf'
    #op.get_xsf_cart_in(b[0],d,len(lst1),name_tmp)

