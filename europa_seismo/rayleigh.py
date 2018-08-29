import os
import numpy as np
import matplotlib.pyplot as plt

def read_deck_model(modelfile):
    f = open(modelfile,'r')
    lines = f.readlines()
    modelname = lines[0].strip()[0]
    nlayers = int(lines[2].strip().split()[0])
    nic = int(lines[2].strip().split()[1])
    noc = int(lines[2].strip().split()[2])
    ncr = int(lines[2].strip().split()[3])
    indices = (nlayers,nic,noc,ncr)
    f.close()
    vals_arr = np.loadtxt(modelfile,skiprows=3)

    return vals_arr,indices

def plot_deck_model(modelfile,planet_radius_km=1565.,ax=None,show=True):
    vals_arr, indices = read_deck_model(modelfile)
    rad = vals_arr[:,0]/1000.0
    rho = vals_arr[:,1]/1000.0
    vp = vals_arr[:,2]/1000.0
    vs = vals_arr[:,3]/1000.0
    depth = planet_radius_km - rad

    if ax == None:
        fig,ax = plt.subplots(1,figsize=(6,6))
        ax.plot(rho,depth,color='red',label='rho')
        ax.plot(vp,depth,color='blue',label='vp')
        ax.plot(vs,depth,color='green',label='vs')
        ax.set_ylim([200.0,0.0])
        ax.set_ylabel('depth')
        ax.set_xlabel('velocity (km/s), density (g/cm3)')
        ax.legend()

    else:
        ax.plot(rho,depth,color='red')
        ax.plot(vp,depth,color='blue')
        ax.plot(vs,depth,color='green')

    if show:
        plt.show()
    else:
        return ax

def make_layers_dict(regolith_thickness,ice_thickness,ocean_thickness,**kwargs):
    vp_ice = kwargs.get('vp_ice',4000.0)
    vs_ice = kwargs.get('vs_ice',2000.0)
    rho_ice = kwargs.get('rho_ice',1000.0)
    vp_ocean = kwargs.get('vp_ocean',1600.0)
    vs_ocean = kwargs.get('vs_ocean',0.0)
    rho_ocean = kwargs.get('rho_ocean',1000.0)
    vp_regolith = kwargs.get('vp_regolith',4000.0)
    vs_regolith = kwargs.get('vs_regolith',2000.0)
    rho_regolith = kwargs.get('rho_regolith',1000.0)
    Qkappa = kwargs.get('Qkappa',40000)
    Qmu = kwargs.get('Qmu',200)
    eta = kwargs.get('eta',1)

    layers = {'regolith':{'thickness':regolith_thickness,'vp':vp_regolith,'vs':vs_regolith,'rho':rho_regolith,'Qkappa':Qkappa,'Qmu':Qmu,'eta':eta},
              'ice': {'thickness':ice_thickness,'vp':vp_ice,'vs':vs_ice,'rho':rho_ice,'Qkappa':Qkappa,'Qmu':Qmu,'eta':eta},
              'ocean': {'thickness':ocean_thickness,'vp':vp_ocean,'vs':vs_ocean,'rho':rho_ocean,'Qkappa':Qkappa,'Qmu':Qmu,'eta':eta}}

    return layers

def write_deck_model(layers_dict,output_model,base_model,**kwargs):

    #read base model
    pts_in_layer = 30
    vals_arr,indices = read_deck_model(base_model)
    nic = indices[1]
    noc = nic + 2 + (pts_in_layer) - 1
    ncr = nic + 4 + (pts_in_layer * 2) - 2
    nlayers = nic + 6 + (pts_in_layer * 3) - 3

    #write header
    f = open(output_model,'w')
    f.write(output_model+'\n')
    f.write('  1    1.00000  1\n')
    f.write('{} {} {} {}\n'.format(nlayers,nic,noc,ncr))
    #f.close()

    #write core and mantle model (i.e. to start of ocean)
    #fmt='%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f'
    #fmt='%10.2f '*9
    fmt = '%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f' 
    np.savetxt(f, vals_arr[0:nic], fmt=fmt)

    #append ocean, ice, and regolith layer on top of the mantle
    fmt = '{:8.0f}{:9.2f}{:9.2f}{:9.2f}{:9.1f}{:9.1f}{:9.2f}{:9.2f}{:9.5f}' 
    #fmt='{:>10.2f} '*9
    ric = vals_arr[:,0][nic]
    roc = ric + layers_dict['ocean']['thickness']
    rcr = roc + layers_dict['ice']['thickness']
    rplanet = rcr + layers_dict['regolith']['thickness']

    for i in range(0,pts_in_layer):
        f.write(fmt.format(ric + (i*((roc-ric)/pts_in_layer)),
                           layers_dict['ocean']['rho'],
                           layers_dict['ocean']['vp'],
                           layers_dict['ocean']['vs'],
                           layers_dict['ocean']['Qkappa'],
                           layers_dict['ocean']['Qmu'],
                           layers_dict['ocean']['vp'],
                           layers_dict['ocean']['vs'],
                           layers_dict['ocean']['eta'])+'\n')
    f.write(fmt.format(roc,
                       layers_dict['ocean']['rho'],
                       layers_dict['ocean']['vp'],
                       layers_dict['ocean']['vs'],
                       layers_dict['ocean']['Qkappa'],
                       layers_dict['ocean']['Qmu'],
                       layers_dict['ocean']['vp'],
                       layers_dict['ocean']['vs'],
                       layers_dict['ocean']['eta'])+'\n')

    for i in range(0,pts_in_layer):
        f.write(fmt.format(roc + (i*((rcr-roc)/pts_in_layer)),
                           layers_dict['ice']['rho'],
                           layers_dict['ice']['vp'],
                           layers_dict['ice']['vs'],
                           layers_dict['ice']['Qkappa'],
                           layers_dict['ice']['Qmu'],
                           layers_dict['ice']['vp'],
                           layers_dict['ice']['vs'],
                           layers_dict['ice']['eta'])+'\n')
    f.write(fmt.format(rcr,
                       layers_dict['ice']['rho'],
                       layers_dict['ice']['vp'],
                       layers_dict['ice']['vs'],
                       layers_dict['ice']['Qkappa'],
                       layers_dict['ice']['Qmu'],
                       layers_dict['ice']['vp'],
                       layers_dict['ice']['vs'],
                       layers_dict['ice']['eta'])+'\n')
    for i in range(0,pts_in_layer):
        f.write(fmt.format(rcr + (i*((rplanet-rcr)/pts_in_layer)),
                           layers_dict['regolith']['rho'],
                           layers_dict['regolith']['vp'],
                           layers_dict['regolith']['vs'],
                           layers_dict['regolith']['Qkappa'],
                           layers_dict['regolith']['Qmu'],
                           layers_dict['regolith']['vp'],
                           layers_dict['regolith']['vs'],
                           layers_dict['regolith']['eta'])+'\n')
    f.write(fmt.format(rplanet,
                       layers_dict['regolith']['rho'],
                       layers_dict['regolith']['vp'],
                       layers_dict['regolith']['vs'],
                       layers_dict['regolith']['Qkappa'],
                       layers_dict['regolith']['Qmu'],
                       layers_dict['regolith']['vp'],
                       layers_dict['regolith']['vs'],
                       layers_dict['regolith']['eta'])+'\n')
    f.close()

