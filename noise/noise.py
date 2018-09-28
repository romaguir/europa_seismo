import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

def make_some_noise(f,npts,dt):
    '''
    Make noise for a given amplitude spectrum

    params:
        f (type:np.array), amplitude spectrum of noise
        npts (type:int), number of points in noise vector
        dt (type:float), time step of noie vector
    '''
    f = np.array(f, dtype='complex')
    Np = (len(f) - 1) // 2
    phases = np.random.rand(Np) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    f[1:Np+1] *= phases
    f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    return np.fft.ifft(f * np.sqrt(npts / 2 / dt)).real

def noise_from_external(input_file,npts=1024,dt=1.0,units='amplitude'):
    '''
    Returns a vector of noise created from the amplitude spectrum given in an external
    file.

    params:
        input_file (type:str), name of external file
        npts (type:int), number of samples in noise vector
        dt (type:float), time step of noise vector
        units (type:str), units used in external file (either 'amplitude', 'power', or 'dB')
        
    returns:
        noise (type:np.array), a vector of noise
    '''
    noisefile = np.loadtxt(input_file)
    #per = noisefile[:,0]
    f_in = noisefile[:,0]
    spec = noisefile[:,1]
    
    if units=='dB':
        spec = 10**(spec/10.0)
    elif units=='power':
        spec = np.sqrt(spec)

    #f_in = 1./per
    f_interp = interp1d(f_in,spec,bounds_error=None,fill_value='extrapolate')
    freqs = np.abs(np.fft.fftfreq(npts, dt))
    f = f_interp(freqs)

    return make_some_noise(f, npts, dt)

def SEASHELS_noise(min_freq, max_freq, samples=1024, sampling_rate=1.0):
    freqs = np.abs(np.fft.fftfreq(samples, 1/sampling_rate))
    f = np.zeros(samples)
    for idx,freq in enumerate(freqs):
        if freq > 0.01 and freq < 20.0:
            f[idx] = 1e-8
        elif freq >= 20.0:
            f[idx] = 5e-8
    return make_some_noise(f)

def nlnm_noise(npts=1024,dt=1.0):
    p, power = get_nlnm() # returns period and power of acceleration PSD
    f_in = 1. / p
    power_in = 10**(power/10.0)
    f_interp = interp1d(f_in,power_in,bounds_error=None,fill_value='extrapolate')

    freqs = np.abs(np.fft.fftfreq(npts, dt))
    f = f_interp(freqs)
    return make_some_noise(f, npts, dt)

def nhnm_noise(npts=1024,dt=1.0):
    p, power = get_nhnm() # returns period and power of acceleration PSD
    f_in = 1. / p
    power_in = 10**(power/10.0)
    f_interp = interp1d(f_in,power_in,bounds_error=None,fill_value='extrapolate')

    freqs = np.abs(np.fft.fftfreq(npts, dt))
    f = f_interp(freqs)
    return make_some_noise(f, npts, dt)

def band_limited_noise(min_freq, max_freq, samples=1024, sampling_rate=1.0):
    freqs = np.abs(np.fft.fftfreq(samples, 1/sampling_rate))
    f = np.zeros(samples)
    idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    f[idx] = 1
    return make_some_noise(f, npts, dt)

