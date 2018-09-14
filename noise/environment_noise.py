from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from europa_seismo.noise.instrument_noise import fftnoise
from scipy.interpolate import interp1d
import numpy as np

def ppsd_noise(min_freq,max_freq,input_file,samples=1024,sampling_rate=1.0):
    noisefile = np.loadtxt(input_file)
    p = noisefile[:,0]
    power = noisefile[:,1] #PSD in db
    f_in = 1./p
    power_in = 10**(power/10.0)
    f_interp = interp1d(f_in,power_in,bounds_error=None,fill_value='extrapolate')

    freqs = np.abs(np.fft.fftfreq(samples, 1./sampling_rate))
    f = f_interp(freqs)
    return fftnoise(f)
