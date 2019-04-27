import sys
#sys.path.append("../aln/")
import warnings

#import h5py
import scipy
import numpy as np
#import matplotlib.pyplot as plt
#import defaultParameters as dp
#import timeIntegration as ti
#from matplotlib.colors import LogNorm

def analyse_run(measure = 'domfr', result = [], dt = 0.1):
    '''
    Analysis routine for bifurcation diagrams with a 'rect' stimulus that is used
    to detect bistability.
    
    measure:    Pass the measure you want to compute in string form, such as "domfr_power_exc"
    result:     Timeseries of a successful simulation.
    dt:         Integration timestep of simulations in ms

    '''
    t = result['t']
    
    down_window = (2000<t) & (t<3000) # time period in ms where we expect the down-state
    up_window = (5000<t) & (t<6000) # and up state
    
    if measure.endswith('inh'):
        rate = result['rate_inh']
    else:
        rate = result['rate_exc']
    
        
    if measure.startswith('domfr_power'):
        if np.any((rate > 0)):
            spectrum_windowsize =  0.5 # in seconds 
            f, Pxx_spec = scipy.signal.welch(
            rate[down_window],
            1000/dt,
            window='hanning',
            nperseg=int(spectrum_windowsize * 1000 / dt)-1,
            scaling='spectrum')
            f = f[f < 70]
            Pxx_spec = Pxx_spec[0:len(f)]
            #domfr = f[Pxx_spec.argmax()] if max(Pxx_spec) > 1 else 0
            #print("domfr: {} max_spec: {}".format(domfr, max(Pxx_spec)))
            return np.max(Pxx_spec)
        else: 
            return 0.0
    elif measure.startswith('domfr'):
        if np.any((rate > 0)):
            spectrum_windowsize =  0.5 # in seconds 
            f, Pxx_spec = scipy.signal.welch(
            rate[down_window],
            1000/dt,
            window='hanning',
            nperseg=int(spectrum_windowsize * 1000 / dt)-1,
            scaling='spectrum')
            f = f[f < 70]
            Pxx_spec = Pxx_spec[0:len(f)]
            domfr = f[Pxx_spec.argmax()] if max(Pxx_spec) > 1 else 0
            #print("domfr: {} max_spec: {}".format(domfr, max(Pxx_spec)))
            return domfr
        else: 
            return 0.0
        
    elif measure.startswith('max'):
        return np.max(rate[up_window])
    
    elif measure.startswith('min'):
        return np.min(rate[up_window])    
    
    elif measure.startswith('updowndiff'):
        up_state_rate = np.mean(rate[up_window])
        down_state_rate = np.mean(rate[down_window])
        up_down_difference = up_state_rate - down_state_rate
        return up_down_difference

    elif measure.startswith('spectrum'):
        if np.any((rate > 0)):
            spectrum_windowsize = 1.0 
            f, Pxx_spec = scipy.signal.welch(
            rate[t>1000],
            1000/dt,
            window='hanning',
            nperseg=int(spectrum_windowsize * 1000 / dt),
            scaling='spectrum')
            f = f[f < 70]
            Pxx_spec = Pxx_spec[0:len(f)]
            Pxx_spec /= np.max(Pxx_spec)
            return f, Pxx_spec

    elif measure.startswith('amplitudes'):
        a = rate[t>1000]
        
        these_maxima = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
        these_maxima = these_maxima[1:-1] # cut first and last peak because theyre not peaks
        these_minima = np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
        these_minima = these_minima[1:-1] # cut first and last peak because theyre not peaks

        a = a[1:-1]
        #return these_maxima, these_minima
        return a[these_maxima], a[these_minima]





