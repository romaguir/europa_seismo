import matplotlib.pyplot as plt
import numpy as np

#mineos
f = np.genfromtxt('modes.out')
period = f[:,5]
grp_vel = f[:,6]
#plt.plot(1./period, grp_vel, label = 'mineos')
plt.plot(period, grp_vel, label = 'mineos')

#bessel
g = np.loadtxt('aicehot_10km_10deg_bessel_10-300s.dat')
freq = g[:,0]
period = 1./ freq
grp_vel = g[:,1]
#plt.plot(1./period, grp_vel, label = 'bessel')
plt.plot(period, grp_vel, label = 'bessel')

#butter
g = np.loadtxt('aicehot_10km_butter_10-150s.dat')
freq = g[:,0]
period = 1./ freq
grp_vel = g[:,1]
#plt.plot(1./period, grp_vel, label = 'butter')
plt.plot(period, grp_vel, label = 'butter')

plt.legend()
#plt.xlabel('period (s)')
plt.xlabel('frequency (Hz)')
plt.ylabel('group velocity (km/s)')

plt.xlim([10,300])
#plt.xlim([1./150,1./10])

plt.show()
