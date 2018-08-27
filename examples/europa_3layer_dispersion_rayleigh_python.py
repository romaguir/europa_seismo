import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from europa_seismo.europa_seismo import rayleigh,rayleigh_python

cwd = os.getcwd()

h2o_thickness = 1565000.0 - 1438875
regolith_thickness = 2000.0
ice_thicknesses = np.linspace(3000.0,48000.0,10)
#ice_thicknesses = np.linspace(50000.0,100000.0,100)

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
    eps=1e-10     #1e-10
    dt=5.0        #10.0
    npow=10
    fnyquist=0
    nbran=0 #number of branches
    cmin=0.135     #0.15
    #cmin=0.1
    cmax=500.0
    maxlyr=1
    modelfile='testmod.deck'
    #modelfile='../data/models/icehot_20km_simple.txt'

    #run model
    (modearray,nmodes) = rayleigh_python.rayleigh(eps,npow,dt,fnyquist,nbran,
                                                  cmin,cmax,maxlyr,modelfile)

    #plot model
    print 'OCEAN DEPTH', ocean_depth_km
    #rayleigh.plot_deck_model(modelfile)

    #plot dispersion curves
    parray = 1./modearray[2,:nmodes] #list of periods
    phase_vel = modearray[3,:nmodes] #phase velocity (km/s)
    group_vel = modearray[4,:nmodes] #group velocity (km/s)
    freq_mhz = (1./parray)*1000.0

    axes[0].plot(freq_mhz,group_vel,color=colorVal1,label='{:3.0f}'.format(ocean_depth_km))
    axes[1].plot(freq_mhz,phase_vel,color=colorVal2,label='{:3.0f}'.format(ocean_depth_km))

#plot simple model
modelfile='../data/models/icehot_20km_simple.txt'
(modearray,nmodes) = rayleigh_python.rayleigh(eps,npow,dt,fnyquist,nbran,
                                              cmin,cmax,maxlyr,modelfile)
parray = 1./modearray[2,:nmodes] #list of periods
phase_vel = modearray[3,:nmodes] #phase velocity (km/s)
group_vel = modearray[4,:nmodes] #group velocity (km/s)
freq_mhz = (1./parray)*1000.0

axes[0].plot(freq_mhz,group_vel,color='black',linewidth=3,linestyle='--',label='icehot_20km_simple')
axes[1].plot(freq_mhz,phase_vel,color='black',linewidth=3,linestyle='--',label='icehot_20km_simple')

axes[0].set_xlabel('frequency (mHz)')
axes[0].set_ylabel('velocity (km/s)')
axes[1].set_xlabel('frequency (mHz)')
axes[1].set_ylabel('velocity (km/s)')
axes[0].set_title('group velocity')
axes[1].set_title('phase velocity')
axes[0].legend()
axes[1].legend()
plt.show()
