#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 14:28:14 2020

@author: reeseboucher
"""

import numpy as np
import scipy as sci
from scipy import stats
from scipy import special
import matplotlib.pyplot as plt

energies  = np.load("1-D_1000_100real_r^0.1_DOS_first.npy") #first
energies2 = np.load("1-D_1000_100real_r^0.1.npy")        #zeroth

def DOSCalc(spectrum):   
#        spectrum = [x for x in spectrum if -0.5 < x < 0.5 ]
#        print(np.log((np.abs(spectrum)).min()))            
#        print(spectrum)
        size            =   len(spectrum)
        hist, bin_edges =   np.histogram(spectrum, bins = 400, density=False);
        bin_width       =   bin_edges[2] - bin_edges[1]
        norm_hist       =   hist/(size*bin_width);
        return norm_hist, bin_edges;


DOS, DOS_bin_edges = DOSCalc(energies);
DOS2, DOS2_bin_edges=DOSCalc(energies2);

#Theory Equation given in Skinner paper
alpha = 0.1
d     = 1
Vo = 1
k  = ((d**2)*(sci.special.gamma(d/2)))/(2*alpha*(np.pi**(d/2))*(Vo**(d/alpha)))
x  = np.linspace(-0.78,0.78,100)
y  = (abs(x)**((d/(alpha))-1))*k


##one dimension test
#x=np.linspace(-1,1,100)
#y=(1/np.pi)*np.log(abs(x))+1.6

##two dimension
#x = np.linspace(-.75,.75,100)
#y = (2/np.pi)*abs(x)

#three dimension
#x = np.linspace(-.65,.65,100)
#y = (3/np.pi)*(x**2)


##Slope and intercept of log-log relation
#slope     = stats.linregress(np.log(abs((DOS2_bin_edges[:-1]))),np.log(DOS))[0]
#intercept = stats.linregress(np.log(abs((DOS2_bin_edges[:-1]))),np.log(DOS))[1]
#
#print(slope)
#print(intercept)

plt.figure(1)
plt.plot(x, y, 'blue', label="Theory")                              #Theory
plt.plot(DOS_bin_edges[:-1], DOS, color='green', label = "First")   #First Stability energies
plt.plot(DOS2_bin_edges[:-1], DOS2, color='red', label = "Zeroth")  #zeroth Stability energies2

##Log plot
#plt.plot(np.log(abs(x)),np.log(abs(y)),'blue', label="Theory")  
#plt.plot(np.log(abs((DOS2_bin_edges[:-1]))), np.log(DOS),"ro", label="Log") 
#plt.xlim(-8,0)


plt.plot()
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.title("1D, 1000, 100 realizations, alpha=0.1")
plt.legend()
plt.grid()
plt.savefig("1-D_1000_100_r^0.1_combo.pdf")
plt.show()


