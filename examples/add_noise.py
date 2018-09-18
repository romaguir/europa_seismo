import instaseis
import numpy as np
import matplotlib.pyplot as plt
from europa_seismo.noise.environment_noise import ppsd_noise

ppsd_file = '/Users/rossmaguire/Tools/europa_seismo/data/noise_catalogs/icehot20km_preferred/ppsd_mean.dat'
sampling_rate=20.0
noise = ppsd_noise(min_freq=1/200.0,max_freq=2.0,input_file=ppsd_file,samples=2000,sampling_rate=sampling_rate)
n_w = np.fft.fft(noise)
n_w_dB = 10*np.log10(np.abs(n_w))
w = np.fft.fftfreq(len(noise),1./sampling_rate)

#plt.plot(noise)
plt.plot(w,n_w_dB)
plt.show()
