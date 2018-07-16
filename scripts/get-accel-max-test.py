import instaseis
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from europa_seismo import acceleration,utils
from obspy.taup import TauPyModel

mars = TauPyModel('EH45Tcold')
db = instaseis.open_db('/home/romaguir/Documents/axisem_databases/Mars_EH45Tcold_5s_noatten_database')
print db

vs_surface = 3.0
mars_radius = 3389.5 
mars_circ = 2*np.pi*mars_radius
km_per_deg = mars_circ / 360.0

def get_max_accel(Mw,gcarc,Q=57822,units='dB',plot=False):
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
    print Mo

    #instaseis source (dipslip at equator)
    source = instaseis.Source(latitude=0.0,longitude=0.0,depth_in_m=300.0,m_rr=3.98e13,m_pp=-3.98e13)

    #receiver at 60 degrees distance
    receiver = instaseis.Receiver(latitude=0.0,longitude=gcarc)

    #get seismograms
    st = db.get_seismograms(source,receiver,components='Z',kind='displacement')
    st_accel = db.get_seismograms(source,receiver,components='Z',kind='acceleration')

    st.resample(2.0)
    st_accel.resample(2.0)

    #freqs,d_w = signal.periodogram(st[0].data,st[0].stats.sampling_rate)
    freqs,d_w = signal.welch(st[0].data,st[0].stats.sampling_rate)
    freqs,a_w = signal.welch(st_accel[0].data,st[0].stats.sampling_rate)

    #find corner frequency
    corner_freq = acceleration.get_corner_freq(Mo=Mo)
    #print 'corner_freq =', corner_freq

    #Aki and Richards 10.36
    s_w = Mo/((1.0 + ((freqs/corner_freq)**2))**(1.5))
    plt.loglog(freqs,s_w)
    plt.show()

    #Haskell parameters
    t_rise,t_dur = acceleration.get_trise_tdur(corner_freq=corner_freq,scale_factor=8.)
    print 'trise, t_dur = ', t_rise,t_dur

    #Haskell source 
    sinc1 = np.sinc((2*np.pi*freqs*t_rise)/2.)
    sinc2 = np.sinc((2*np.pi*freqs*t_dur)/2.)
    haskell = np.abs(sinc1*sinc2)
    
    s1_arg = ((2*np.pi*freqs*t_rise)/(2.))
    s2_arg = ((2*np.pi*freqs*t_dur)/(2.))
    haskell_v2 = np.abs( ((np.sin(s1_arg)/s1_arg)) * ((np.sin(s2_arg)/s2_arg)) )
    haskell_v2[0] = 1.0

    #traveltime (use first arriving phase)
    #arrs = mars.get_travel_times(source.depth_in_m/1000.0,gcarc)
    #print arrs
    #tt = arrs[0].time
    d_km = gcarc*km_per_deg
    tt = d_km * (1./vs_surface)

    #attenuation factor
    atten = np.exp((-np.pi*corner_freq*tt)/(Q))
    print Q,tt,atten

    #acceleration spectrum
    #A_w = (freqs**2)*d_w*atten*Mo*haskell
    A_w = ((2*np.pi*freqs)**2)*d_w*atten*Mo
    A_w_max = np.max(A_w)
    print 'MAXIMUM ACCELERATION', A_w_max

    a_w = a_w*Mo*atten*haskell_v2
    a_w_max=np.max(a_w)

    #plot
    if plot:
       #plt.loglog(freqs,A_w)
       fig,axes = plt.subplots(nrows=2,ncols=2,figsize=(14,8))
       axes[0,0].loglog(freqs,d_w)
       axes[0,0].set_xlabel('Frequency (Hz)')
       axes[0,0].set_ylabel('displacement PSD m^2/Hz)')
       axes[0,1].loglog(freqs,A_w,label='calculated')
       axes[0,1].loglog(freqs,a_w,label='instaseis')
       axes[0,1].set_xlabel('Frequency (Hz)')
       axes[0,1].set_ylabel('acceleration PSD (m/s^2)^2/Hz)')
       axes[0,1].legend()
       axes[1,1].semilogx(freqs,haskell)
       axes[1,1].semilogx(freqs,haskell_v2)
       plt.show()

    if units=='dB':
        #return 10*np.log10((A_w_max)/(1e-4))
        return 10*np.log10((a_w_max)/(1e-4))
    else:
        raise ValueError('units not implemented')

Mws = np.arange(0,10,0.5)
gcarcs = np.arange(0,182,2)
x = []
y = []
amax = []

'''
for gcarc in gcarcs:
    print gcarc
    for Mw in Mws:
        dB = get_max_accel(Mw=Mw,gcarc=gcarc)
        x.append(gcarc)
        y.append(Mw)

        if dB >= 0.0:
            amax.append(dB)
        else:
            amax.append(-10.0)

plt.scatter(x,y,c=amax,edgecolor=None)
plt.colorbar()
plt.show()
np.savetxt('AMAX_pts.dat',np.c_[x,y,amax],fmt='%4f')
'''

dB = get_max_accel(Mw=7.0,gcarc=60.0,plot=True)
print dB
