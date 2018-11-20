import instaseis
import numpy as np
import matplotlib.pyplot as plt
from obspy import geodetics
from scipy.signal import hilbert
from europa_seismo.europa_seismo.utils import gauss_filter
from obspy.geodetics.base import gps2dist_azimuth

kind='displacement'
components='Z'
evlo = 0.0
evla = 0.0
evdp = 0.0
stlo = 30.0
stla = 0.0
Tmin = 10.0 #minimum period
Tmax = 200.0 #maximum period
nT = 100
noise_level = 0.10
#insta_db_name = '/home/romaguir/Documents/10s_PREM_ANI_FORCES'
insta_db_name = '/media/romaguir/wdseis/instaseis_databases/aicehot_5km_prem/aicehot_5km_prem_database/'
#window_min = 600.0
#window_max = 1600.0
window_min = 400.0
window_max = 3600.0
m_rr = 3.98e13
m_pp = 3.98e13
planet_radius = 1565.0
#planet_radius = 6371.0
distaz = gps2dist_azimuth(lat1=evla, lon1=evlo, lat2=stla, lon2=stlo, a=planet_radius*1000.0, f=0.0) # no flattening
gcarc_m = distaz[0]
dist_km = gcarc_m / 1000.0


#get instaseis seismogram
instaseis_db = instaseis.open_db(insta_db_name)
source = instaseis.Source(latitude=evla,
                          longitude=evlo,
                          depth_in_m=evdp*1000.,
                          m_rr = m_rr,
                          m_pp = m_rr)
receiver = instaseis.Receiver(latitude=stla,
                              longitude=stlo)
stream = instaseis_db.get_seismograms(source=source,
                                      receiver=receiver,
                                      components=components,
                                      kind = kind,
                                      remove_source_shift=True)

stream = stream.slice(stream[0].stats.starttime + window_min,
                      stream[0].stats.starttime + window_max)
#stream.plot()

time = np.linspace(window_min,
                   window_min+(stream[0].stats.npts * stream[0].stats.delta),
                   stream[0].stats.npts)

#make ensemble of dispersion curves w/ different noise realizations
def mft(stream,Tmin,Tmax,nT,noise_level,a=1.2,b=2.0,plot=True):
    stream_copy = stream.copy()
    per_picks = [] #period at which grp vel is picked
    vel_picks = [] #measured group velocity
    periods = np.linspace(Tmin,Tmax,nT)
    sigmax = np.max(np.abs(stream[0].data))
    noise = np.random.random(len(stream[0].data))
    noise = (noise*2.0) -1.
    noise *= (noise_level * sigmax)
    stream_copy[0].data += noise

    if plot:
        fig,ax = plt.subplots(1,sharex=True,figsize=[10,2])
        ax.plot(time,stream[0].data,c='k',label='no noise')
        ax.plot(time,stream_copy[0].data,c='b',label='noise level = {}'.format(noise_level))
        ax.set_xlabel('time (s)')
        plt.legend()
        plt.show()

    #do mft analysis
    for i,period in enumerate(periods):
        #print i,period
        #stream_copy.plot()

        Tstart = period / a
        Tend = period + (period / b)
        freqmin = 1. / Tend
        freqmax = 1. / Tstart

        if freqmin > 0:
            tr = stream_copy[0].filter('bandpass',
                freqmin=freqmin,freqmax=freqmax,corners=4,
                zerophase=True)
        else:
            tr = stream_copy[0].filter('lowpass',
                freq=freqmax)

        dcol = tr.data
        env = np.abs(hilbert(dcol.real))
        #gabor_matrix[:,i] = env

        if np.max(env) > 0.0:
            per_picks.append(period)
            vel_picks.append(dist_km / time[np.argmax(env)])

        #plt.plot(time,tr.data)
        #plt.scatter(period,np.max(np.abs(env)))
        #plt.show()
        stream_copy = stream.copy()
        stream_copy[0].data += noise
        
    #plt.plot(per_picks,vel_picks)
    #plt.show()
    return per_picks,vel_picks

print 'STARTING MFT'
#mft(stream=stream,Tmin=Tmin,Tmax=Tmax,nT=100,noise_level=0.10,a=1.2,b=2.0,plot=True)

nrealizations = 200
data_matrix = np.zeros((nT,nrealizations))
for i in range(0,nrealizations):
    print i
    per_picks,vel_picks = mft(stream=stream,Tmin=Tmin,Tmax=Tmax,nT=nT,noise_level=noise_level,a=1.2,b=2.0,plot=False)
    data_matrix[:,i] = vel_picks

np.save('data_matrix_nl{}'.format(noise_level),data_matrix)
