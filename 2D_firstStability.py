#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:32:22 2020

@author: reeseboucher
"""

import matplotlib.pyplot as plt
import numpy as np

realizations = 0
firstStability = 0
energies     = []
phi          = []
N            = []
allEnergies  = []
allEnergies  = np.array(allEnergies)
phi          = np.array(phi)
lowInd       = [-1,-1]
highInd      = [-1,-1]


#input variables
stateNum        = 625                   #sites in lattice
lattice_x       = 25                      #stateNum = lattice_x*lattice_y
lattice_y       = 25
w               = 1                        #Max and abs(min) energy 
realizationsNum = 100                      #Number of separate calculations
alpha           = 0.5                        #Power of r in potential


while realizations < realizationsNum:
    phi      = np.random.uniform(-w,w,stateNum)
    energies = np.zeros(stateNum)
    N        = np.zeros(stateNum)
    
    phi      = np.reshape(phi,(lattice_y,lattice_x))
    energies = np.reshape(energies,(lattice_y,lattice_x))
    N        = np.reshape(N,(lattice_y,lattice_x))

    N[np.where(phi > 0)] = 0
    N[np.where(phi < 0)] = 1

    highVal = 0
    lowVal  = 0
    
    for i in range(0,lattice_y):        #i and j track y and x of q_1
        print(i)
        for j in range(0,lattice_x):
            sigmaE_noPhi = 0
            for y_idx in range(0,lattice_y):
                for x_idx in range(0,lattice_x):    #Tracks position of y and x of q_2
                    r = np.sqrt(((i - y_idx)**2)+((j - x_idx)**2))
                    if r != 0:
                        sigmaE_noPhi = -(0.5 - N[y_idx][x_idx])/(r**alpha)+sigmaE_noPhi

            energies[i][j] = phi[i][j] + sigmaE_noPhi

            if energies[i][j] > 0 and N[i][j] == 1:     #HighVal and LowVal are energies of states that will be modified
                challenger = energies[i][j]
                if challenger > highVal:             
                    highVal = challenger
                    highInd[0] = i
                    highInd[1] = j

            elif energies[i][j] < 0 and N[i][j] == 0:
                challenger = energies[i][j]
                if challenger < lowVal:
                    lowVal = challenger
                    lowInd[0] = i
                    lowInd[1] = j

    stop = False    
    while stop == False:
        dN_high = 0
        dN_low  = 0
                      
        if highVal > 0:
            N[highInd[0]][highInd[1]] = 0
            dN_high = -1
        if lowVal < 0:
            N[lowInd[0]][lowInd[1]] = 1
            dN_low = 1
        
        rDiff = np.sqrt((highInd[0]-lowInd[0])**2+(highInd[1]-lowInd[1])**2) 
        
        for y_idx in range(0,lattice_y):
            for x_idx in range(0,lattice_x):    #Tracks position of y and x of q_2

                rHigh = np.sqrt((highInd[1] - x_idx)**2 + (highInd[0] - y_idx)**2)         
                rLow  = np.sqrt((lowInd[1] - x_idx)**2 + (lowInd[0] - y_idx)**2)

                if highInd[0] != -1:
                    if rHigh != 0:
                        energies[y_idx][x_idx] = energies[y_idx][x_idx] + dN_high/(rHigh**alpha)
                    if rLow == 0 and rDiff!=0:
                        energies[y_idx][x_idx] = energies[y_idx][x_idx] + dN_high/(rDiff**alpha)
 
                
                if lowInd[0] != -1:
                    if rLow != 0:
                        energies[y_idx][x_idx] = energies[y_idx][x_idx] + dN_low/(rLow**alpha)
                    if rHigh == 0 and rDiff != 0:
                        energies[y_idx][x_idx] = energies[y_idx][x_idx] + dN_low/(rDiff**alpha)
                    
        lowVal     =  0
        highVal    =  0
        highInd[0] = -1                         #No energy contribution from high or low electron if there
        lowInd[0]  = -1                         #is none that is modified and highInd[0]=lowInd[0]=-1
        for i in range(lattice_y):
            for j in range(lattice_x):
                if energies[i][j] > 0 and N[i][j] == 1:
                    challenger = energies[i][j]
                    if challenger > highVal:             #HighVal is energy of state that gets modified
                        highVal = challenger
                        highInd[0] = i
                        highInd[1] = j
            
                elif energies[i][j] < 0 and N[i][j] == 0:
                    challenger = energies[i][j]
                    if challenger < lowVal:
                        lowVal = challenger
                        lowInd[0] = i
                        lowInd[1] = j
        print(highVal)
        print(lowVal)
            
            
        if highVal == 0 and lowVal == 0:         #First Order Stability Criterion
            i = 0
            j = 0
            complete = True
            while i < lattice_y and complete == True:        
                while j < lattice_x and complete == True:                   
                    k = 0
                    l = 0
                    while k < lattice_y and complete == True:
                        while l < lattice_x and complete == True:                    
                            if energies[k][l] > energies[i][j]  and N[i][j] == 1 and N[k][l] == 0 and energies[k][l]-energies[i][j] < (1/np.sqrt(((k-i)**2 + (l-j)**2)**alpha)):
                                highVal =  1
                                lowVal  = -1
                                highInd[0] = i
                                highInd[1] = j
                                lowInd[0]  = k
                                lowInd[1]  = l
                                complete   = False
                                firstStability += 1
                            l+=1
                        k+=1
                    j+=1
                i += 1

            print("First Stability Corrections = "+str(firstStability))
            complete=True                        
            if complete == True:
                stop = True
                allEnergies = np.append(allEnergies,energies)
            
    realizations += 1
    print("===========")
    print(realizations)
    


#Probability Density
def DOSCalc(spectrum):
        size                =   spectrum.size;
        hist, bin_edges     =   np.histogram(spectrum, bins=100, density=False);
        bin_width           =   bin_edges[2] - bin_edges[1]
        norm_hist           =   hist/(size*bin_width);
        return norm_hist, bin_edges;
 
#np.save("2-D_150x150_100_r^0.5_first.npy",allEnergies)
DOS, DOS_bin_edges     =       DOSCalc(allEnergies);

#x=np.linspace(-.65,.65,100)
#y=(2/np.pi)*abs(x)
plt.figure(1)
plt.plot(DOS_bin_edges[:-1], DOS, color='green')
#plt.plot(x,y,r)
plt.xlabel('Energy')
plt.ylabel('DOS')
#plt.title("2-D: DOS vs Energy Plot: Interacting System")
plt.grid()
#plt.savefig("2-D-150x150_100_r^0.5_first.pdf")
plt.show()

