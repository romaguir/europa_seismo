import numpy as np

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
