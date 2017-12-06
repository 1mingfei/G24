#!/usr/bin/python
import numpy as np



def cutoff_func(Rij,Rc):
    if Rij<=Rc:
        a=0.5*(np.cos(np.pi*Rij/Rc)+1)
    else:
        a=0
    return a

def G2_func(Rij,Rc,Rs,eta):
    return np.exp(-eta*(Rij-Rs)**2)*cutoff_func(Rij,Rc)


def G4_func(Rij,Rik,Rjk,Rc,cos_theta,eta,lamda,zeta):
    p1=(1+lamda*cos_theta)**zeta
    e=np.exp(-eta*(Rij**2+Rik**2+Rjk**2))
    f=cutoff_func(Rij,Rc)*cutoff_func(Rik,Rc)*cutoff_func(Rjk,Rc)
    fnl= p1*e*f*2**(1-zeta)

    return fnl

