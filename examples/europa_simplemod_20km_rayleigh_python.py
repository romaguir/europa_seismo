import os
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.europa_seismo import rayleigh,rayleigh_python

cwd = os.getcwd()

#run model
eps=1e-15
dt=5.0
#npow=10
npow=10
fnyquist=0
nbran=0 #number of branches
cmin=0.135
#cmin=0.1
cmax=500.0
maxlyr=1.0
#modelfile='../data/models/icehot_20km_simple.txt'
modelfile='ice20.deck'

#run model
(modearray,nmodes) = rayleigh_python.rayleigh(eps,npow,dt,fnyquist,nbran,
                                                  cmin,cmax,maxlyr,modelfile)

#plot model
#rayleigh.plot_deck_model(modelfile)

#plot dispersion curves
parray = 1./modearray[2,:nmodes] #list of periods
phase_vel = modearray[3,:nmodes] #phase velocity (km/s)
group_vel = modearray[4,:nmodes] #group velocity (km/s)
freq_mhz = (1./parray)*1000.0

plt.plot(freq_mhz,group_vel,label='group velocity',color='blue',alpha=0.5)
plt.plot(freq_mhz,phase_vel,label='phase_velocity',color='purple',alpha=0.5)
plt.xlabel('frequency (mHz)')
plt.ylabel('velocity (km/s)')
plt.legend()
plt.show()
