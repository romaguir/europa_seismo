import os
import shutil
import pickle
import subprocess
import matplotlib
import numpy as np
import pymc as pm
import matplotlib.pyplot as plt
from europa_seismo.europa_seismo import rayleigh, minos
from scipy.interpolate import interp1d
from pymc.Matplot import plot

cwd = os.getcwd()

#--------------------------------------------------------------------
# set up run parameters
#--------------------------------------------------------------------
runname = 'testrun'
input_file = '../data/dispersion_curves/group_vel_ice5_Mw3_d20_noisemodelD.dat'
Tmin = 10.0
Tmax = 200.0
nT = 100
periods = np.linspace(Tmin,Tmax,nT)
n_iter = 100
n_burn = 50
gcarc_true = 20.00
km_per_deg = 27.31
dist_km = gcarc_true * km_per_deg
vs_i = 2.0

#--------------------------------------------------------------------
# read 'observation'
#--------------------------------------------------------------------
f_input = np.loadtxt(input_file)
freq_obs = f_input[:,0]    #frequency in mHz
grpvel_obs = f_input[:,1]  #group velocity in km/s
grptt_obs = grpvel_obs * dist_km # convert to group time

#--------------------------------------------------------------------
# setup model parameterization
#--------------------------------------------------------------------
h_r = pm.Uniform('h_r',lower=0.0,upper=2.0) #regolith thickness
h_i = pm.Uniform('h_i',lower=h_r,upper=50.0) #ice thickness
vs_r = pm.Uniform('vs_r',lower=0.5,upper=3.5) #regolith shear velocity
gcarc = pm.Uniform('gcarc',lower=10.0,upper=50.0) # distance from event

#--------------------------------------------------------------------
# define data variance
#--------------------------------------------------------------------
sigma = 50.0
tau = np.power(sigma,-2)

#--------------------------------------------------------------------
# define forward model
#--------------------------------------------------------------------
@pm.deterministic
def mu(h_r=h_r, vs_r=vs_r, h_i=h_i, vs_i=vs_i, gcarc=gcarc):
    
    #write current model
    h2o_thickness = 1565000.0 - 1438875
    ocean_thickness = h2o_thickness - (h_i*1000.0)
    layers = rayleigh.make_layers_dict(regolith_thickness = h_r*1000.0,
                                       ice_thickness = h_i*1000.0,
                                       ocean_thickness = ocean_thickness,
                                       vp_regolith=4000.,
                                       vs_regolith=vs_r*1000.0,
                                       vs_ice=vs_i*1000.0)
    
    rayleigh.write_deck_model(layers,output_model='europa.deck',
                              base_model=cwd+'/../data/models/'+'icehot_20km_simple.txt')

    
    #run current model
    modesfile='modes.out'
    eigfile='eig.out'
    eps=1e-8
    wgrav=1.315
    jcom=3
    lmin=0
    lmax=1000
    wmin=0.1
    wmax=100.0 #10 s period
    nmin=0.0                                                                    
    nmax=0.0

    modelfile='europa.deck'
    minos.main(modelfile,modesfile,eigfile,eps,wgrav,jcom,lmin,lmax,wmin,wmax,nmin,nmax)
    f = np.genfromtxt('modes.out')
    freq_mhz = f[:,4]
    group_vel = f[:,6] #group velocity
    
    #interpolate to observed frequencies
    freq_interp = interp1d(freq_mhz,group_vel,bounds_error=None,fill_value='extrapolate')
    vel_modeled = freq_interp(freq_obs)
    
    # use delta to convert distpersion curves from velocity to time
    dist_km = gcarc * km_per_deg
    t_p = dist_km / vel_modeled
    t_0 = (1./len(t_p))*np.sum(t_p-grptt_obs)
    group_time_optimal = t_p - t_0
    
    return group_time_optimal

#--------------------------------------------------------------------
# define likelyhood
#--------------------------------------------------------------------
y = pm.Normal('y', mu=mu, tau=tau, value=grptt_obs, observed=True)

#--------------------------------------------------------------------
# inference
#--------------------------------------------------------------------
m = pm.Model([h_r,h_i,vs_r,gcarc,tau,y])
mc = pm.MCMC(m, db='pickle', dbname=runname+'.pkl')
mc.sample(iter=n_iter,burn=n_burn)

#--------------------------------------------------------------------
# save models
#--------------------------------------------------------------------
mc.write_csv(runname+'.csv',variables=['h_r','h_i','vs_r','gcarc'])

#--------------------------------------------------------------------
# clean up and organize
#--------------------------------------------------------------------
os.mkdir(runname)
subprocess.call('mv ./*.png',shell=True)
subprocess.call('mv ./*.pkl {}'.format(runname),shell=True)
subprocess.call('mv ./*.csv {}'.format(runname),shell=True)
