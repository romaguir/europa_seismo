import os
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.europa_seismo import rayleigh,minos

cwd = os.getcwd()

h2o_thickness = 1565000.0 - 1438875
regolith_thickness = 2000.0
ice_thicknesses = np.linspace(5000.0,60000.0,10)

for i,ice_thickness in enumerate(ice_thicknesses):
    print i,ice_thickness

    ocean_thickness = h2o_thickness - regolith_thickness - ice_thickness
    layers = rayleigh.make_layers_dict(regolith_thickness,ice_thickness,ocean_thickness,
                                       vp_regolith=2000.,vs_regolith=1000.)

    #write model
    rayleigh.write_deck_model(layers,output_model='testmod.deck',base_model=cwd+'/../data/'+'icehot_20km_simple.txt')

    #run model
    modesfile='modes.out'
    eigfile='eig.out'
    eps=1e-9
    wgrav=1.5
    jcom=3
    lmin=0
    lmax=1000
    wmin=0.0
    wmax=100.0 #10 s period
    nmin=0.0
    nmax=0.0

    modelfile='testmod.deck'
    minos.main(modelfile,modesfile,eigfile,eps,wgrav,jcom,lmin,lmax,wmin,wmax,nmin,nmax)

    #plot model
    #rayleigh.plot_deck_model(modelfile)

    #plot dispersion curves
    f = np.genfromtxt('modes.out')
    freq = f[:,4]
    phase_vel = f[:,3]
    group_vel= f[:,6]

    if i == 0:
        plt.plot(freq,group_vel,label='group velocity',color='blue')
        plt.plot(freq,phase_vel,label='phase_velocity',color='purple')
    else:
        plt.plot(freq,group_vel,color='blue')
        plt.plot(freq,phase_vel,color='purple')

plt.xlabel('frequency (mHz)')
plt.ylabel('velocity (km/s)')
plt.legend()
plt.show()
