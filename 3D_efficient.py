#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 00:19:42 2020

@author: reeseboucher
"""

import matplotlib.pyplot as plt
import numpy as np

realizations = 0
energies     = []
phi          = []
N            = []
allEnergies  = []
allEnergies  = np.array(allEnergies)
lowInd       = [-1,-1,-1]
highInd      = [-1,-1,-1]

#input variables
stateNum        = 3375                    #length of initial lattice
lattice_x       = 15                       #stateNum = lattice_x*lattice_y
lattice_y       = 15
lattice_z       = 15 
w               = 1                       #Max and abs(min) energy 
realizationsNum = 100                     #Number of separate calculations
alpha           = 1                       #Power of r in potential

while realizations < realizationsNum:
    phi      = np.random.uniform(-w,w,stateNum)
    energies = np.zeros(stateNum)
    N        = np.zeros(stateNum)
    
    phi      = np.reshape(phi,(lattice_z,lattice_y,lattice_x))
    energies = np.reshape(energies,(lattice_z,lattice_y,lattice_x))
    N        = np.reshape(N,(lattice_z,lattice_y,lattice_x))

    N[np.where(phi > 0)] = 0
    N[np.where(phi < 0)] = 1

    highVal = 0
    lowVal  = 0
    for h in range(0,lattice_z):
        print(h)
        for i in range(0,lattice_y):        #i and j track y and x of q_1
            for j in range(0,lattice_x):
                sigmaE_noPhi = 0
                for z_idx in range(0,lattice_z):
                    for y_idx in range(0,lattice_y):
                        for x_idx in range(0,lattice_x):    #Tracks position of y and x of q_2

                            r = np.sqrt((h - z_idx)**2 + (i - y_idx)**2 + (j - x_idx)**2) #magnitude of r (rz^2+ry^2+rx^2)^(1/2)

                            if r != 0:
                                sigmaE_noPhi = -(0.5 - N[z_idx][y_idx][x_idx])/(r**alpha)+sigmaE_noPhi
        
                    energies[h][i][j] = phi[h][i][j] + sigmaE_noPhi
       
                    if energies[h][i][j] > 0 and N[h][i][j] == 1:     #HighVal and LowVal are states that will be modified
                        challenger = energies[h][i][j]
                        if challenger > highVal:             
                            highVal = challenger
                            highInd[0] = h
                            highInd[1] = i
                            highInd[2] = j
        
                    elif energies[h][i][j] < 0 and N[h][i][j] == 0:
                        challenger = energies[h][i][j]
                        if challenger < lowVal:
                            lowVal = challenger
                            lowInd[0] = h
                            lowInd[1] = i
                            lowInd[2] = j

    stop = False    
    while stop == False:
        dN_high = 0
        dN_low  = 0
        
        if highVal > 0:
            N[highInd[0]][highInd[1]][highInd[2]] = 0
            dN_high = -1


        if lowVal < 0:
            N[lowInd[0]][lowInd[1]][lowInd[2]] = 1
            dN_low = 1

        
        rDiff = np.sqrt(abs(highInd[0]-lowInd[0])**2+abs(highInd[1]-lowInd[1])**2+abs(highInd[2]-lowInd[2])**2) #questionable
        
        for z_idx in range(0,lattice_z):
            for y_idx in range(0,lattice_y):
                for x_idx in range(0,lattice_x):    #Tracks position of y and x of q_2
                    
                    rHigh = np.sqrt((highInd[2]-x_idx)**2 + (highInd[1]-y_idx)**2 + (highInd[0]-z_idx)**2)            
                    rLow  = np.sqrt((lowInd[2]-x_idx)**2 + (lowInd[1]-y_idx)**2 + (lowInd[0]-z_idx)**2)
                    if highInd[0]!=-1:
                        if rHigh != 0:
                            energies[z_idx][y_idx][x_idx] = energies[z_idx][y_idx][x_idx] + dN_high/(rHigh**alpha) 
                        if rLow == 0 and rDiff != 0:
                            energies[z_idx][y_idx][x_idx] = energies[z_idx][y_idx][x_idx] + dN_high/(rDiff**alpha)                        
                        
                    if lowInd[0]!=-1:
                        if rLow != 0:
                            energies[z_idx][y_idx][x_idx] = energies[z_idx][y_idx][x_idx] + dN_low/(rLow**alpha)
                            
                        if rHigh == 0 and rDiff != 0:
                            energies[z_idx][y_idx][x_idx] = energies[z_idx][y_idx][x_idx] + dN_low/(rDiff**alpha)        
                                                      
                    
        lowVal  = 0
        highVal = 0
        highInd[0] = -1
        lowInd[0]  = -1
    
        for h in range(lattice_z):              
            for i in range(lattice_x):
                for j in range(lattice_y):
                    if energies[h][i][j] > 0 and N[h][i][j] == 1:
                        challenger = energies[h][i][j]
                        if challenger > highVal:             #HighVal is energy of state that gets modified
                            highVal = challenger
                            highInd[0] = h
                            highInd[1] = i
                            highInd[2] = j
                            
                
                    elif energies[h][i][j] < 0 and N[h][i][j] == 0:
                        challenger = energies[h][i][j]
                        if challenger < lowVal:
                            lowVal = challenger
                            lowInd[0] = h
                            lowInd[1] = i
                            lowInd[2] = j
        print(highVal)
        print(lowVal)
        if highVal == 0 and lowVal == 0:
            stop = True
            allEnergies = np.append(allEnergies,energies)
            
    realizations += 1
    print("===========")
    print(realizations)
    

#Probability Density    
def DOSCalc(spectrum):
        size                =   spectrum.size;
        hist, bin_edges     =   np.histogram(spectrum, bins = 1000, density=False);
        bin_width           =   bin_edges[2] - bin_edges[1]
        norm_hist           =   hist/(size*bin_width);
        return norm_hist, bin_edges;
 
np.save("3d_15x15x15_r^1_100.npy",allEnergies)
DOS, DOS_bin_edges     =       DOSCalc(allEnergies);

plt.figure(1)
plt.plot(DOS_bin_edges[:-1], DOS, color='green')
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.title("3-D DOS vs Energy Plot: Interacting System")
plt.grid()
plt.savefig("3-D_15x15x15_100_r^1.pdf")
plt.show()

