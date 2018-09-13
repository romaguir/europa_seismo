from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from europa_seismo.noise.instrument_noise import fftnoise
from scipy.interpolate import interp1d
import numpy as np

def nlnm_noise(min_freq,max_freq,samples=1024,sampling_rate=1.0):
    p, power = get_nlnm() # returns period and power of acceleration PSD
    f_in = 1. / p
    power_in = 10**(power/10.0)
    f_interp = interp1d(f_in,power_in,bounds_error=None,fill_value='extrapolate')

    freqs = np.abs(np.fft.fftfreq(samples, 1./sampling_rate))
    f = f_interp(freqs)
    return fftnoise(f)

def nhnm_noise(min_freq,max_freq,samples=1024,sampling_rate=1.0):
    p, power = get_nhnm() # returns period and power of acceleration PSD
    f_in = 1. / p
    power_in = 10**(power/10.0)
    f_interp = interp1d(f_in,power_in,bounds_error=None,fill_value='extrapolate')

    freqs = np.abs(np.fft.fftfreq(samples, 1./sampling_rate))
    f = f_interp(freqs)
    return fftnoise(f)
