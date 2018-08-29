import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from europa_seismo.europa_seismo import rayleigh,minos

cwd = os.getcwd()

h2o_thickness = 1565000.0 - 1438875
regolith_thickness = 2000.0
ice_thicknesses = np.linspace(3000.0,98000.0,20)

#setup colors
viridis = plt.get_cmap('viridis')
inferno = plt.get_cmap('inferno')
values = range(len(ice_thicknesses))
cNorm = colors.Normalize(vmin=0,vmax=values[-1])
scalarMap1 = cm.ScalarMappable(norm=cNorm,cmap=viridis)
scalarMap2 = cm.ScalarMappable(norm=cNorm,cmap=inferno)

#setup figure
fig,axes = plt.subplots(1,2,figsize=(10,6))

for i,ice_thickness in enumerate(ice_thicknesses):
    print i,ice_thickness
    colorVal1 = scalarMap1.to_rgba(values[i])
    colorVal2 = scalarMap2.to_rgba(values[i])

    ocean_thickness = h2o_thickness - regolith_thickness - ice_thickness
    ocean_depth_km = (ice_thickness + regolith_thickness) / 1000.
    layers = rayleigh.make_layers_dict(regolith_thickness,ice_thickness,ocean_thickness,
                                       vp_regolith=4000.,vs_regolith=2000.)

    #write model
    rayleigh.write_deck_model(layers,output_model='testmod.deck',base_model=cwd+'/../data/models/'+'icehot_20km_simple.txt')

    #run model
    modesfile='modes.out'
    eigfile='eig.out'
    eps=1e-12
    wgrav=1.315
    jcom=3
    lmin=0
    lmax=1000
    wmin=0.1
    wmax=100.0 #10 s period
    nmin=0.0
    nmax=0.0

    modelfile='testmod.deck'
    minos.main(modelfile,modesfile,eigfile,eps,wgrav,jcom,lmin,lmax,wmin,wmax,nmin,nmax)

    #plot model
    #rayleigh.plot_deck_model(modelfile)

    #plot dispersion curves
    f = np.genfromtxt('modes.out')
    freq_mhz = f[:,4]
    phase_vel = f[:,3]
    group_vel= f[:,6]

    axes[0].plot(freq_mhz,group_vel,color=colorVal1,label='{:3.0f}'.format(ocean_depth_km))         
    axes[1].plot(freq_mhz,phase_vel,color=colorVal2,label='{:3.0f}'.format(ocean_depth_km))

axes[0].set_xlabel('frequency (mHz)')
axes[0].set_ylabel('velocity (km/s)')
axes[1].set_xlabel('frequency (mHz)')
axes[1].set_ylabel('velocity (km/s)')
axes[0].set_title('group velocity')
axes[1].set_title('phase velocity')
axes[0].grid()
axes[1].grid()

axes[0].legend()
axes[1].legend()
plt.show()
