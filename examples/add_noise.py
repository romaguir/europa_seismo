import os
import instaseis
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.noise import environment_noise
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.fftpack import fft,ifft,fftfreq
from scipy.interpolate import interp1d
import matplotlib.mlab as mlab

n_samples = 10000
dt = 0.1
freqs = np.abs(fftfreq(n_samples, dt))

#--------------------------------------------------------------------
# Petterson 1993 high noise model example
#--------------------------------------------------------------------
per,power_dB = get_nhnm() #returns period (in s) and power of nlnm (in dB)
f_in = 1/per
power = 10**(power_dB/10.) #convert from decibels to (m/s^2)^2/Hz
power_i = np.interp(freqs,f_in,power)
#plt.loglog(freqs,np.sqrt(power_i),c='r')
#plt.loglog(f_in,10**(power_dB/20.))

#method 1: use fftnoise function
#noise1 = fftnoise(np.sqrt(power_i),n_samples,1./dt)

#method 2: add noise to a normalized white noise spectrum
whitenoise = np.random.random(n_samples)
whitenoise_w = fft(whitenoise)
whitenoise_w /= np.abs(whitenoise_w)
whitenoise_w *=  np.sqrt(power_i)

time = np.linspace(0,n_samples*dt,n_samples)
noise2 = ifft(whitenoise_w * np.sqrt(n_samples / 2 / dt) )
noise2 -= np.mean(noise2)
power1,f1 = mlab.psd(x=noise1, Fs=1./dt, NFFT = 2**14)
power2,f2 = mlab.psd(x=noise2, Fs=1./dt, NFFT = 2**14)

fig,axes = plt.subplots(2,figsize=(8,8))
axes[0].semilogx(1./f2[f2 > 0], 10*np.log10(power2[f2 > 0]),label='noise spectrum')
axes[0].semilogx(per,power_dB,label='peterson 1993 nhnm')
axes[0].set_xlim([0.1,500.0])
axes[0].set_xlabel('period (s)')
axes[0].set_ylabel('power (dB)')
axes[0].legend()
#axes[1].plot(time,noise1)
#axes[1].set_xlabel('time (s)')
#axes[1].set_ylabel('acceleration (m/s^2)')
plt.show()

#--------------------------------------------------------------------
# PPSD noise model of europa
#--------------------------------------------------------------------
HOME=os.getenv("HOME")
ppsd_file = HOME+'/Tools/europa_seismo/data/noise_catalogs/icehot20km_preferred/ppsd_mean.dat'
g = np.loadtxt(ppsd_file)
per = g[:,0]
power_dB = g[:,1] + 50.0
f_in = 1./per

power = 10**(power_dB/10.) #convert from decibels to (m/s^2)^2/Hz
min_power = np.min(power_dB)
get_power = interp1d(f_in,power_dB,bounds_error=False,fill_value=min_power)
power_i = get_power(freqs)
power_i = 10**(power_i/10.) #convert from decibels to (m/s^2)^2/Hz

#generate noise
noise = fftnoise(np.sqrt(power_i),n_samples,1./dt)
power,f = mlab.psd(x=noise, Fs=1./dt, NFFT = 2**14)
fig,axes = plt.subplots(2,figsize=(8,8))

#plot
axes[0].semilogx(1./f[f > 0], 10*np.log10(power[f > 0]),label='noise spectrum')
axes[0].semilogx(per,power_dB,label='europa_icehot20km PPSD')
axes[0].set_xlim([0.1,500.0])
axes[0].set_xlabel('period (s)')
axes[0].set_ylabel('power (dB)')
axes[0].legend()
axes[1].plot(time,noise)
axes[1].set_xlabel('time (s)')
axes[1].set_ylabel('acceleration (m/s^2)')
plt.show()
