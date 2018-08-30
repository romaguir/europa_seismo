import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from numpy.fft import fft,ifft,fftfreq

def make_noise(npts,dt=1.0):
    time = np.linspace(0,npts*dt,npts)
    noise = np.random.randn(npts)
    noise_w = fft(noise)
    noise /= np.abs(noise)
    w = fftfreq(npts,dt)
    return time,noise,w,noise_w

def plot_noise(time,noise,w,noise_w):
    fig,axes = plt.subplots(2,figsize=(5,10))

    axes[0].plot(time,noise)
    axes[0].set_xlabel('time (s)')

    axes[1].loglog(w,noise_w)
    axes[1].set_xlabel('frequency (hz)')
    plt.show()

time,noise,w,noise_w = make_noise(1000,dt=1.0)
plot_noise(time,noise,w,noise_w)
