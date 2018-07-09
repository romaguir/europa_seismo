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
