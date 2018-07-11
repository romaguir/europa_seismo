#!/usr/bin/env python

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

def get_acceleration_spectrum(omega_min,omega_max,seismic_response,
                              corner_freq,t_rise,t_dur,Mo,Q,tt):
    '''
    Returns the acceleration spectrum recorded at a receiver, assuming
    a Haskell fault model (Haskell, 1969)

    positional arguments--------------------------------------------------
    omega_min: minimum frequency of spectrum (Hz)
    omega_max: maximum frequency of spectrum (Hz)
    seismic_response: normalized spectrum of the Green's function
    corner_freq: corner frequency (Hz)
    t_rise: rise time (s)
    t_dur: duration (s)
    Mo: scalar moment (Nm)
    Q: bulk attenuation
    tt: travel time between source and receiver 
   
    returns---------------------------------------------------------------
    A_w = acceleration spectrum 
    '''
    omega_pts = 1000
    omega = np.linspace(omega_min,omega_max,omega_pts)
    A_w = -omega**2*seismic_response*np.exp((-np.pi*corner_freq*tt)/(Q))*\
          Mo*np.abs(np.sinc(omega*2*np.pi*t_rise)*np.sinct_dur(omega*2*np.pi*t_dur))

    return A_w

def get_response_spectrum(data,dt,kind='acceleration'):
