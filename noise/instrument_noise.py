import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

def fftnoise(f,samples=1024,sampling_rate=1.0):
    f = np.array(f, dtype='complex')
    Np = (len(f) - 1) // 2
    phases = np.random.rand(Np) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    f[1:Np+1] *= phases
    f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    #return np.fft.ifft(f).real
    return np.fft.ifft(f * np.sqrt(samples / 2 / (1./sampling_rate))).real #as in Stahlers seismic_noise package

def band_limited_noise(min_freq, max_freq, samples=1024, sampling_rate=1.0):
    freqs = np.abs(np.fft.fftfreq(samples, 1/sampling_rate))
    f = np.zeros(samples)
    idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    f[idx] = 1
    return fftnoise(f)

def SEASHELS_noise(min_freq, max_freq, samples=1024, sampling_rate=1.0):
    #TODO: Figure out if I need to square the values of the PSD.
    freqs = np.abs(np.fft.fftfreq(samples, 1/sampling_rate))
    f = np.zeros(samples)
    for idx,freq in enumerate(freqs):
        if freq > 0.01 and freq < 20.0:
            f[idx] = 1e-8
        elif freq >= 20.0:
            f[idx] = 5e-8
    return fftnoise(f)

#def from_custom_file(min_freq,max_freq,filename,samples=1024,sampling_rate=1.0):
#    noisefile = np.loadtxt(filename)

def plot_noise(noise,sampling_rate):
    time = np.linspace(0,len(noise)*(1./sampling_rate),len(noise))
    
    plt.plot(np.abs(np.fft.fft(noise)))
    plt.show()
