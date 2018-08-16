import os
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.europa_seismo import rayleigh,rayleigh_python

cwd = os.getcwd()

h2o_thickness = 1565000.0 - 1438875
regolith_thickness = 2000.0
ice_thicknesses = np.linspace(5000.0,60000.0,100)

for i,ice_thickness in enumerate(ice_thicknesses):
    print i,ice_thickness

    ocean_thickness = h2o_thickness - regolith_thickness - ice_thickness
    layers = rayleigh.make_layers_dict(regolith_thickness,ice_thickness,ocean_thickness,
                                       vp_regolith=2000.,vs_regolith=1000.)

    #write model
    rayleigh.write_deck_model(layers,output_model='testmod.deck',base_model=cwd+'/../data/'+'icehot_20km_simple.txt')

    #run model
    eps=1e-10
    dt=10.0
    npow=10
    fnyquist=0
    nbran=0
    cmin=0.15
    cmax=500.0
    maxlyr=1
    modelfile='testmod.deck'

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

    if i == 0:
        plt.plot(freq_mhz,group_vel,label='group velocity',color='blue')
        plt.plot(freq_mhz,phase_vel,label='phase_velocity',color='purple')
    else:
        plt.plot(freq_mhz,group_vel,color='blue')
        plt.plot(freq_mhz,phase_vel,color='purple')

plt.xlabel('frequency (mHz)')
plt.ylabel('velocity (km/s)')
plt.legend()
plt.show()
