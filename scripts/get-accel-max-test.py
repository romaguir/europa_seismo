import instaseis
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from europa_seismo import acceleration,utils
from obspy.taup import TauPyModel

mars = TauPyModel('EH45Tcold')
db = instaseis.open_db('/Volumes/G-RAID with Thunderbolt/axisem_databases/Mars_EH45Tcold_5s_noatten_database')

def get_max_accel(Mw,gcarc,Q=600,units='dB',plot=False):
    '''
    returns maximum acceleration for a Haskell source

    params:
    Mw: moment magnitude of earthquake
    gcarc: distance in degrees
    Q: quality factor (defaults to 600)
    units: only dB so far... 0 dB is defined as 1e-4 (i.e., 10 mgal)
    '''

    #find scalar moment
    Mo = utils.Mw_to_Mo(Mw)

    #instaseis source (dipslip at equator)
    source = instaseis.Source(latitude=0.0,longitude=0.0,depth_in_m=300.0,m_rr=3.98e13,m_pp=-3.98e13)

    #receiver at 60 degrees distance
    receiver = instaseis.Receiver(latitude=0.0,longitude=gcarc)

    #get seismograms
    st = db.get_seismograms(source,receiver,components='Z',kind='displacement')
    freqs,d_w = signal.periodogram(st[0].data,st[0].stats.sampling_rate)

    #find corner frequency
    corner_freq = acceleration.get_corner_freq(Mo=Mo)
    #print 'corner_freq =', corner_freq

    #Haskell parameters
    t_rise,t_dur = acceleration.get_trise_tdur(corner_freq=corner_freq,scale_factor=8.)
    #print 'trise, t_dur = ', t_rise,t_dur

    #Haskell source 
    sinc1 = np.sinc((freqs*t_rise)/2.)
    sinc2 = np.sinc((freqs*t_dur)/2.)
    haskell = sinc1*sinc2

    #traveltime (use first arriving phase)
    arrs = mars.get_travel_times(source.depth_in_m/1000.0,gcarc)
    #print arrs
    tt = arrs[0].time

    #attenuation factor
    atten = np.exp((-np.pi*corner_freq*tt)/(Q))

    #acceleration spectrum
    A_w = (freqs**2)*d_w*atten*Mo*haskell
    A_w_max = np.max(A_w)
    #print 'MAXIMUM ACCELERATION', A_w_max

    #plot
    if plot:
       plt.loglog(freqs,A_w)
       plt.show()

    if units=='dB':
        return 10*np.log10((A_w_max)/(1e-4))
    else:
        raise ValueError('units not implemented')

Mws = np.arange(0,10,0.25)
gcarcs = np.arange(0,181,1)
x = []
y = []
amax = []

for gcarc in gcarcs:
    print gcarc
    for Mw in Mws:
        dB = get_max_accel(Mw=Mw,gcarc=gcarc)
        x.append(gcarc)
        y.append(Mw)

        if dB >= 0.0:
            amax.append(dB)
        else:
            amax.append(0.0)

plt.scatter(x,y,c=amax,edgecolor=None)
plt.colorbar()
plt.show()

np.savetxt('AMAX_pts.dat',np.c_[x,y,amax],fmt='%4f')
