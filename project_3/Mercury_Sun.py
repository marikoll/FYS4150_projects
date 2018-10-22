#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 12:23:49 2018

@author: maritkollstuen
"""

import numpy as np
# import matplotlib.pyplot as plt
import time



def velocity_verlet(pos0, vel0, N, dt):
    G = 4*np.pi**2
    c2 = 63239.7263**2
    dt2 = dt**2
    
    start = time.time()
    
    pos = pos0
    vel = vel0
    acc = -pos[0]*G/np.linalg.norm(pos)**3*(1+(3*np.linalg.norm(np.cross(pos,vel))**2)/(np.linalg.norm(pos)**2*c2))

    
    pos_p = []
    
    
    posnew = pos + dt*vel + 0.5*dt2*acc
        
    l = np.linalg.norm(np.cross(posnew,vel))
    r = np.linalg.norm(posnew)
    
    accnew = -posnew*G/r**3*(1+(3*l**2)/(r**2*3999262982.498912))
    velnew = vel + 0.5*(accnew+accnew)*dt
    
    posprev = pos
    velprev = vel
    accvel  = acc
        
    pos = posnew
    vel = velnew
    acc = accnew

    for i in range(1,N-1):
        posnew = pos + dt*vel + 0.5*dt2*acc
        
        l = np.linalg.norm(np.cross(posnew,vel))
        r = np.linalg.norm(posnew)
        
        accnew = -posnew*G/r**3*(1+(3*l**2)/(r**2*3999262982.498912))
        velnew = vel + 0.5*(accnew+accnew)*dt
        
        if np.linalg.norm(pos)< r and np.linalg.norm(pos)< np.linalg.norm(posprev):
            pos_p.append(pos)
        
        posprev = pos
        velprev = vel
        accvel  = acc
        
        pos = posnew
        vel = velnew
        acc = accnew
        
        
    stop = time.time()
    print("Verlet: {:.3f} sec".format(stop-start))
    return pos_p

if __name__ == "__main__":    
    pass
    
    
    
    
    