#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 23:25:55 2020

@author: reeseboucher
#"""

import matplotlib.pyplot as plt
import numpy as np

realizations = 0
energies     = []
phi          = []
N            = []
high_list    = []
low_list     = []
allEnergies  = []
allEnergies  = np.array(allEnergies)
phi          = np.array(phi)
high_list    = np.array(high_list)
low_list     = np.array(low_list)


#Input Variables
stateNum        = 1000                    #length of initial lattice
w               = 1                       #Max and abs(min) energy 
realizationsNum = 100                     #Number of separate calulcations
alpha           = 0.1                    #power of r in potential

while realizations < realizationsNum:
    highVal = 0 
    lowVal  = 0
    highInd = 0
    lowInd  = 0
    
    phi      = np.random.uniform(-w,w,stateNum)
    energies = np.zeros(stateNum)
    N        = np.zeros(stateNum)
    
    N[np.where(phi > 0)] = 0
    N[np.where(phi < 0)] = 1


    for i in range(stateNum):    #Modifies and returns Energies array with updated energies
        sigmaE_noPhi = 0
        for j in range(stateNum):
            r_ij  = abs(j - i)
            if r_ij != 0:
                sigmaE_noPhi = -(0.5-N[j])/(r_ij**alpha)+sigmaE_noPhi        
        energies[i] = phi[i] + sigmaE_noPhi
                            
        if energies[i] > 0 and N[i] == 1:
            challenger = energies[i]
            if challenger > highVal:            
                highVal = challenger
                highInd = i

        elif energies[i] < 0 and N[i] == 0:
            challenger = energies[i]
            if challenger < lowVal:
                lowVal = challenger
                lowInd = i

    stop = False    
    while stop == False:
        dN_high = 0
        dN_low  = 0
                      
        if highVal > 0:
            N[highInd] = 0
            dN_high    = -1    #dN_high = N[highInd] - n_high  #always equal to -1 if highVal > 0

        if lowVal < 0:
            N[lowInd] = 1
            dN_low    = 1      #dN_low = N[lowInd] - n_low  #always equal to 1 if lowVal < 0

        
        rDiff = np.sqrt(abs(highInd-lowInd)**2)
        
        for x_idx in range(0,stateNum):    #Tracks position of y and x of q_2
            
            rHigh = abs(highInd - x_idx)          
            rLow  = abs(lowInd - x_idx)
            
            if highInd != -1:
                if rHigh != 0:
                    energies[x_idx] = energies[x_idx] + dN_high/(rHigh**alpha)  
                if rLow == 0 and rDiff != 0:
                    energies[x_idx] = energies[x_idx] + dN_high/(rDiff**alpha)
                
            if lowInd != -1:    
                if rLow != 0:
                    energies[x_idx] = energies[x_idx] + dN_low/(rLow**alpha)                    
                if rHigh == 0 and rDiff != 0:
                    energies[x_idx] = energies[x_idx] + dN_low/(rDiff**alpha)
                                             
        lowVal  = 0
        highVal = 0 
        highInd = -1
        lowInd  = -1
              
        for i in range(stateNum):
            if energies[i] > 0 and N[i] == 1:
                challenger = energies[i]
                if challenger > highVal:             #HighVal is energy of occupied state that gets modified
                    highVal = challenger
                    highInd = i
        
            elif energies[i] < 0 and N[i] == 0:
                challenger = energies[i]
                if challenger < lowVal:
                    lowVal = challenger
                    lowInd = i     

        if highVal == 0 and lowVal == 0:
            stop = True
            allEnergies = np.append(allEnergies,energies)
            
    realizations += 1
    print("===========")
    print(realizations)
    

 
np.save("1-D_1000_100real_r^0.1", allEnergies)

#Density of States Plot
def DOSCalc(spectrum):
        size                =   spectrum.size;
        hist, bin_edges     =   np.histogram(spectrum, bins=1000, density=False);
        bin_width           =   bin_edges[2] - bin_edges[1]
        norm_hist           =   hist/(size*bin_width);
        return norm_hist, bin_edges;
 

DOS, DOS_bin_edges = DOSCalc(allEnergies);

plt.figure(1)
plt.plot(DOS_bin_edges[:-1], DOS, color='green')
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.title("1-D: DOS vs Energy Plot for Interacting System")
plt.grid()
plt.savefig("1-D_1000_100real_r^0.1_DOS.pdf")
plt.show()

