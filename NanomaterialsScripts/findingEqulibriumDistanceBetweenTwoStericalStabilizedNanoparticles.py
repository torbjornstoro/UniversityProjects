# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:21:18 2019

@author: to-bo
"""
import numpy as np
import matplotlib.pyplot as plt
#Constants

r=3.4*10**-9 # radii of particles
A=1.9*1.602*10**-19 #Hamaker constant for material, in Joule from eV
d_s = 18*10**-10 # Thickness of monolayer, 18 Å
s_t = 4.3*10**-10 # Diameter of area covered by head group (thiol)
k_b = 1.38*10**-23#Boltzman constant
T=298 # Temperature in kelvin, assuming room temperature

def energyVanDerWaal(c):
    #total van der Waal potential between nanoparticles
    return -A/12*(4*r**2/(c**2-4*r**2)+4*r**2/c**2 + 2*np.log((c**2-4*r**2)/c**2))

def energySteric(c):
    #steric repulsion potential between nanoparticles between carbon chains covering the nanoparticles
    return(100*r*d_s**2)/(np.pi*(c-2*r)*s_t**2)*k_b*T*np.exp(-(np.pi*(c-2*r))/(d_s))
    
lowVal = energyVanDerWaal(2*r+10**-11)+energySteric(2*r+10**-11) # Set a start value
eqDistance=0 # variable for  holding equlibrium distance

#lists for plotting:

xList=[]
vdwList=[]
stericList=[]
sumList=[]

#Iterate through c, to find the c which makes energyVanDerWaal and energySteric the most similar. starting at 6.9 nm, as anything smaller would not make physical sense
for i in range(698,2000):
    if(energyVanDerWaal(i*10**-11)+energySteric(i*10**-11))<lowVal:
        lowVal=energyVanDerWaal(i*10**-11)+energySteric(i*10**-11)
        eqDistance=i*10**-11
    xList.append(i)
    vdwList.append(energyVanDerWaal(i*10**-11))
    stericList.append(energySteric(i*10**-11))
    sumList.append(energySteric(i*10**-11)+energyVanDerWaal(i*10**-11))
print(eqDistance)
print(lowVal)        
    
plt.figure()
#plt.plot(xList,vdwList)
#plt.plot(xList,stericList)
plt.plot(xList,sumList)
plt.axis([698,2000,-1*10**-21,1*10**-21])