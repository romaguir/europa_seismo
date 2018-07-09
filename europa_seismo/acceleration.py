import obspy
import instaseis
import numpy as np
import matplotlib.pyplot as plt

def get_corner_freq(Mo,stress_drop=3.4e6,vs=3.9e3):
    '''
    Returns corner frequency for a given seismic moment, assuming a 
    circular fault patch (e.g., Madariaga [1976]).

    positional arguments--------------------------------------------------
    Mo: scalar_moment [Nm]
    stress_drop: stress drop [Pa]
    vs: shear velocity [m/s]

    returns---------------------------------------------------------------
    corner_freq: corner frequency in Hz
    '''
    corner_freq = 0.42*vs*((stress_drop/Mo)**(1./3.))
    return corner_freq

def get_trise_tdur(corner_freq,scale_factor=8.):
    '''
    Returns the rise time and duration time of a Haskell fault model,
    given the corner frequency. By default, the duration is assumed 
    to be 8 times longer than the rise time (e.g., Allmann and Shearer, 2009)

    positional arguments--------------------------------------------------
    corner_freq: corner frequency [Hz]
    scale_factor: t_dur = scale_factor*t_rise

    returns---------------------------------------------------------------
    t_rise, t_dur: rise time, duration is s
    '''
    t_dur = 1.0/(np.pi*corner_freq)
    t_rise = t_dur/scale_factor
    return t_rise,t_dur
