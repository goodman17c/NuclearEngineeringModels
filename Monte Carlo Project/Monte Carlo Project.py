#Charles Goodman
#5 December 2018
#NE 301
#Monte Carlo

#The code reads in number of energy groups, cross sections for each energy group, and number of particles to simulate.
#The code assumes isotropic scattering in a homogenous medium. All neutrons start in group 1.
# The code outputs a file with statistics from the simulation.

import numpy as np
import random

filename = input("Input file: ")
f = open(filename,"r")
num = int(f.readline()) #number of groups

#import absorption crosssections for each group to a list
abs = np.zeros((num,1))
for i in range(num):
    abs[i] = float(f.readline())

#import scattering cross sections for each each to nested list
scat = np.zeros((num,num))
for i in range(num):
    for j in range(num):
        scat[i,j] = float(f.readline())
#import starting number of simulated particles
n = int(f.readline())
f.close()

#calculate total cross sections for each group
tot = abs + np.sum(scat, axis=1, keepdims=1)
neginvtot = -1/tot
probs = np.cumsum(np.concatenate((abs, scat), axis=1)/tot , axis=1)

#Create statistics to be tracked
rsquaredbar = 0
tracklength = 0
var = 0
scatterbar = 0;
varscatter = 0;

for i in range(n): #for each simulated particle (this loop can be parallelized)
    #starting values for a neutron
    start = np.array([0,0,0]) #location
    group = 1 #energy group
    scatters = 0
    while(group > 0):
        phi = 2*np.pi*random.random()
        theta = np.pi*random.random()
        r = neginvtot[group-1]*np.log(1-random.random())
        react = random.random()
        for j in range(num+1):
            if react < probs[group-1,j]:
                group = j
                break
        if (group > 0):
            scatters += 1
        tracklength += r
        st=np.sin(theta)
        ct=np.cos(theta)
        sp=np.sin(phi)
        cp=np.cos(phi)
        start = start + r*np.array([st*cp,st*sp,ct])
    rsquared = np.sum(np.square(start))
    rsquaredbar = (rsquaredbar*i+rsquared)/(i+1)
    var = ((rsquared-rsquaredbar)**2+i*var)/(i+1)
    scatterbar = (scatters+scatterbar*i)/(i+1)
    varscatter = ((scatters-scatterbar)**2+i*varscatter)/(i+1)

s = np.sqrt(n/(n-1)*var)
error = s/np.sqrt(n)


f = open("soln_"+ filename, "w+")
f.write("Number of particles: %d\n" %n)
f.write("Average R Squared: %E\n" %rsquaredbar)
f.write("Standard Deviation: %E\n" %s)
f.write("Estimated Error: %E\n" %error)
f.write("Average scatter interactions: %f\n" %scatterbar)
f.write("Variance in scatter interactions: %f\n" %varscatter)
f.close();
