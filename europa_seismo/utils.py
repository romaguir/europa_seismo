import obspy
import numpy as np
from scipy import signal

def Mw_to_Mo(Mw,units='Nm'):
    '''
    Converts moment magnitude to scalar moment

    positional arguments-------------------------------------
    Mw: The moment magnitude
    units (optional): units of the scalar moment. default is 'Nm'.
                      also can use 'dynecm'
    '''
    if units=='Nm':
        Mo = 10.0**(1.5*Mw + 9.1)
    elif units=='dynecm':
        Mo = 10.0**(1.5*Mw + 16.1)

    return Mo

def Mo_to_Mw(Mo,units='Nm'):
    '''
    Converts scalar moment to moment magnitude

    positional arguments-------------------------------------
    Mo: The scalar moment
    units (optional): units of the scalar moment. default is 'Nm'.
                      also can use 'dynecm'
    '''
    if units=='Nm':
        Mw = (2./3.)*(np.log10(Mo)-9.1)
    elif units=='dynecm':
        Mw = (2./3.)*(np.log10(Mo)-16.1)

    return Mw

def get_MO(MT):
    '''
    Scalar Moment M0 in Nm

    MT: moment tensor [m_rr,m_tt,m_pp,m_rt,m_rp,m_tp]
    '''
    m_rr = MT[0]
    m_tt = MT[1]
    m_pp = MT[2]
    m_rt = MT[3]
    m_rp = MT[4]
    m_tp = MT[5]

    MO = (np.sqrt(m_rr**2 + m_tt**2 + m_pp**2 +
          2.*m_rt**2 + 2.*m_rp**2 + 2.*m_tp**2)) * (1./np.sqrt(2.))

    return MO
    #return (self.m_rr ** 2 + self.m_tt ** 2 + self.m_pp ** 2 +
    #        2 * self.m_rt ** 2 + 2 * self.m_rp ** 2 +
    #        2 * self.m_tp ** 2) ** 0.5 * 0.5 ** 0.5

def phase_window(data,evdp,gcarc,sampling_rate,phase,
                 window_start,window_end,taup_model,origin_time=0.0):
    '''
    Returns a window around a specified phase

    
    positional arguments-------------------------------------
    data: the seismic time series
    evdp: earthquake depth (km)
    gcarc: great circle distance (degrees)
    sampling_rate: sampling rate in Hz
    phase: seismic phase around which to window
    window_start: window start time relative to phase arrival (s)
    window_end: window end relative to phase arrival (s)
    taup_model: obspy TauPyModel object, used to predict phase arrivals
    origin_time (optional): origin time of earthquake (default = start of time series)
    '''
    
    arrs = taup_model.get_travel_times(source_depth_in_km = evdp,
                                       distance_in_degree = gcarc,
                                       phase_list = [phase]) 
    phase_arr = arrs[0]
    phase_tt = phase_arr.time
    phase_idx = int(phase_tt*sampling_rate) + int(origin_time*sampling_rate) 
    window_start_idx = phase_idx - int(window_start*sampling_rate)
    window_end_idx = phase_idx + int(window_end*sampling_rate)

    return data[window_start_idx:window_end_idx]

def gauss_filter(data,sampling_rate,w_0,alpha,return_filter=False):
    F_data = np.fft.fft(data)
    freqs = np.fft.fftfreq(len(data),d = 1/sampling_rate)
    F_filter = np.zeros(len(F_data))
    for i in range(0,len(freqs)):
        if freqs[i] > 0.0:
            fact = -alpha*((freqs[i]-w_0)/(freqs[i]))**2
            F_filter[i] = np.exp(fact)
    nfreqs = len(freqs)

    #remove negative frequencies
    F_data[int(nfreqs/2):] = 0
    F_filter[int(nfreqs/2):] = 0
    
    F_filtered = F_filter * F_data
    data_filtered = np.fft.ifft(F_filtered)
    if return_filter == False:
        return data_filtered
    else:
        return data_filtered,F_filter[0:int(nfreqs/2)]

def bessel_filter(data,freqmin,freqmax,sampling_rate,corners=4,zerophase=False):
    fe = 0.5 * sampling_rate
    low = freqmin / fe
    high = freqmax / fe
    z,p,k = signal.iirfilter(corners,[low,high],btype='bandpass',ftype='bessel',output='zpk')
    sos = signal.zpk2sos(z,p,k)

    if zerophase:
        firstpass = signal.sosfilt(sos, data)
        return signal.sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return signal.sosfilt(sos, data)

#def stream_filter_bessel(stream,freqmin,freqmax,corners=4,zerophase=False):
def stream_filter_bessel(stream,center_freq,gamma=0.5,corners=4,zerophase=False):
    filtered_stream = obspy.Stream() 
    for tr in stream.copy():
        fe = 0.5 * tr.stats.sampling_rate

        freqmin = center_freq - (gamma * center_freq)
        freqmax = center_freq + (gamma * center_freq)

        low = freqmin / fe
        high = freqmax / fe
        z,p,k = signal.iirfilter(corners,[low,high],btype='bandpass',ftype='bessel',output='zpk')
        sos = signal.zpk2sos(z,p,k)

        if zerophase:
            first_pass = signal.sosfilt(sos,tr.data)
            second_pass = signal.sosfilt(sos,first_pass[::-1])[::-1] 
            tr.data = second_pass
            filtered_stream += tr
        else:
            filtered_data = signal.sosfilt(sos,tr.data)
            tr.data = filtered_data
            filtered_stream += tr

    return filtered_stream
