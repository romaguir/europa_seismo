import os
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.europa_seismo import rayleigh,minos

#run model
modesfile='modes.out'
eigfile='eig.out'
eps=1e-10
wgrav=1.315
jcom=3
lmin=0
lmax=1000
wmin=0.1
wmax=100.0 #10 s period
nmin=0.0
nmax=0.0
#modelfile='../data/models/icehot_20km_simple.txt'
modelfile='../data/models/aicehot_10km_prem.txt'
minos.main(modelfile,modesfile,eigfile,eps,wgrav,jcom,lmin,lmax,wmin,wmax,nmin,nmax)

f = np.genfromtxt('modes.out')
freq = f[:,4]
phase_vel = f[:,3]
group_vel= f[:,6]

#add Gaussian noise
#sigma = 0.1 #standard deviation of error = 0.1 km / s
#error = np.random.normal(loc = 0.0, scale = sigma, size = len(group_vel))
#group_vel += error

#plot dispersion curves
rayleigh.plot_deck_model(modelfile)
plt.plot(freq,group_vel,label='group velocity',color='blue')
plt.plot(freq,phase_vel,label='phase_velocity',color='purple')

plt.xlabel('frequency (mHz)')
plt.ylabel('velocity (km/s)')
plt.legend()
plt.show()

#save dispersion curve
np.savetxt('disp_icehot20_test.dat',np.c_[freq,group_vel],fmt='%5f')
