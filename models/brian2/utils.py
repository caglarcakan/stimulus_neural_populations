# helper functions for signal generation stuff, power spectra, measures, etc.



from __future__ import print_function
import numpy as np
from numpy import zeros
from math import sqrt
import matplotlib.pyplot as plt
import scipy.signal

# try to import numba
# or define dummy decorator
try:
    from numba import njit, jit
except:
    def njit(func):
        return func
    jit = njit

#
#
#
@njit
def ou_x(runtime, dt, tau, mean, sigma_stat, X0, rands):
    '''
    generate OU process. [cf. https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process]
    parameters: tau, mean, sig_stat [sig_stat is the std of the stationary OU process]
    simulating the ou process according to the Langevin Eq.:
    dX = 1/tau*(mu - X)*dt + sigL * dW(t).
    sigL = sigma_stat*sqrt(2/tau) '''

    sigL = sigma_stat*np.sqrt(2./tau)
    steps = int(runtime/dt)
    x = np.zeros(steps+1)
    # optimizations
    sigL_sqrt_dt = sigL * sqrt(dt)
    dt_tau = dt/tau
    x[0] = X0
    for i in xrange(steps):
        x[i+1] = x[i] + dt_tau * (mean - x[i]) + sigL_sqrt_dt * rands[i]
    return x

def generate_OUinput(params):
    simtime = params['runtime']
    dt = params['min_dt']
    tau =params['ou_tau']
    sigma =params['ou_sigma']
    mu =params['ou_mean']
    if params['ou_stationary']:
        X0 = mu
    else:
        X0 = params['ou_X0']
    rands = np.random.randn(int(simtime/dt))
    ou_trace = ou_x(simtime,dt,tau,mu,sigma,X0,rands)
    return ou_trace


# computes the filtered trace x_filtered of an input trace x
# using a gaussian or a biexponential filter function
def x_filter(x, params):
    runtime = params['runtime']
    dt      = params['min_dt']
    tsteps   = int(runtime/dt)
    t = np.linspace(0., runtime, tsteps+1)
    if params['filter_type']=='gauss':
        # sigma of the filter function in [ms]
        gauss_filter_sigma = params['filter_gauss_sigma'] # [ms]
        sigma3 = gauss_filter_sigma*3
        # gauss_filter_sigma translates to the number of
        # points sigmaN which fullfill the equation sigmaN*dt=sigma
        sigmaN = int(gauss_filter_sigma/dt)
        # take as many elements for the discrete gaussian kernel N_tot
        # such that N_tot*dt=sigma*6 = sigma*3 (left side) + sigma*3 (right side)
        N_tot = sigmaN*3*2+1 #symmtery
        filter_t = np.linspace(-sigma3, sigma3, N_tot)
        filter_function = np.exp((-filter_t**2)/ (gauss_filter_sigma**2*2))
        # normalize filter_function
        filter_function /=np.sum(filter_function)
#        x_filtered_wrong = np.convolve(x,filter_function, mode='same') # wrong: boundary effects        
        x_filtered_inner = np.convolve(x,filter_function, mode='valid')
        x_filtered = np.concatenate((x_filtered_inner[::-1][-(N_tot//2):], 
                                     x_filtered_inner, 
                                     x_filtered_inner[::-1][:(N_tot-N_tot//2-1)]))
        assert len(x_filtered) == len(t) 
    # stuff below should implement a biexponential filter
    # elif params['filter_type']=='bi_exp':
    #     tau_r = params['filter_bi_exp_tau_r']
    #     tau_d = params['filter_bi_exp_tau_d']
    #     filter_function = np.exp(-t/tau_d)-np.exp(-t/tau_r)
    #     plt.plot(filter_function)
    #     plt.show()
    #     # normalize filter_function
    #     filter_function /=np.sum(filter_function)
    #     x_filtered = np.convolve(x, filter_function, mode='full')
    #     # make x and filtered x equal in length
    #     x_filtered = x_filtered[:len(x)]
    else:
        raise NotImplementedError('{} filter not implemented!'.format(params['filter_type']))
    return x_filtered

# for getting e.g. a ramp or quadratic increase in an input (sub)interval
def get_changing_input(runtime, start_change, dt,
              start_val, end_val, type_of_input='ramp',
              quad_factor = None, duration_change=None):
    steps=int(runtime/dt)
    t = np.linspace(0.,runtime,steps+1)
    # print(type_of_input)
    if type_of_input=='ramp':
        idx_s = int(start_change/dt)
        idx_d = int(duration_change/dt)
        input = np.ones_like(t)*start_val
        input[idx_s:idx_s+idx_d]=np.linspace(start_val,end_val,idx_d)
        input[idx_s+idx_d:]=end_val
    elif type_of_input=='quad':
        idx_s = int(start_change/dt)
        input = start_val+quad_factor*(t-start_change)**2/1000.**2
        input[:idx_s]=start_val
        idx_up=np.where(input>end_val)
        input[idx_up]=end_val
    else:
        raise NotImplementedError('Input type {} is not implemented'.format(type))
    return input

# efficient interpolating and lookup functions for reduced models
@njit
def interpolate(xi, yi, range_x, range_y):
    #no problems here
    if xi < range_x[0]:
        x_floor_id = 0
        x_dist_id  = 0
    elif xi >= range_x[-1]:
        x_floor_id = -1
        x_dist_id  =  0
    else:
        x_nearest_index = np.argmin(np.abs(range_x-xi))
        if (xi - range_x[x_nearest_index]) > 0:
            x_floor_id = x_nearest_index
            x_dist_id  = (xi - range_x[int(x_floor_id)])/(range_x[int(x_floor_id+1)]-range_x[int(x_floor_id)])
        else:
            x_floor_id = x_nearest_index-1
            x_dist_id  = (xi - range_x[int(x_floor_id)])/(range_x[int(x_floor_id)]-range_x[int(x_floor_id-1)])
    if yi < range_y[0]:
        y_floor_id = 0
        y_dist_id  = 0
    elif yi >= range_y[-1]:
        y_floor_id = -1
        y_dist_id  = 0
    else:
        y_nearest_index = np.argmin(np.abs(range_y-yi))
        if (yi - range_y[y_nearest_index]) > 0:
            y_floor_id = y_nearest_index
            y_dist_id  = (yi -range_y[int(y_floor_id)])/(range_y[int(y_floor_id+1)]-range_y[int(y_floor_id)])
        else:
            y_floor_id = y_nearest_index-1
            y_dist_id  = (yi - range_y[int(y_floor_id)])/(range_y[int(y_floor_id)]-range_y[int(y_floor_id-1)])

    weights = np.zeros(4)
    weights[0] = float(x_floor_id)
    weights[1] = float(x_dist_id)
    weights[2] = float(y_floor_id)
    weights[3] = float(y_dist_id)
    return weights

    #return (float(x_floor_id), float(x_dist_id), float(y_floor_id), float(y_dist_id))
@njit
def look_up_function_1d(Z, weights):
        #shapex, shapey, shapez = Z.shape
        xid1 = weights[0]
        dxid = weights[1]
        yid1 = weights[2]
        dyid = weights[3]
        #section both have definite values namely either first or last one
        if xid1 == -1 and yid1 == -1:
            return Z[int(xid1),int(yid1)]
        elif xid1 != -1 and yid1 == -1:
            return Z[int(xid1), int(yid1)]*(1-dxid)+Z[int(xid1+1), int(yid1)]*dxid
        elif xid1 == -1 and yid1 != -1:
            return Z[int(xid1), int(yid1)]*(1-dyid)+Z[int(xid1), int(yid1+1)]*dyid
        else:
            return Z[int(xid1),int(yid1)]*(1-dxid)*(1-dyid) +\
                   Z[int(xid1+1),int(yid1)]*dxid*(1-dyid)   +\
                   Z[int(xid1),int(yid1+1)]*(1-dxid)*dyid   +\
                   Z[int(xid1+1),int(yid1+1)]*dxid*dyid

# use xi/yi -> mu/sigma to have correct warnings
@njit
def interpolate_xy(xi, yi, rangex, rangey):
    weights = np.zeros(4)
    dimx = rangex.size
    dimy = rangey.size
    # determine weights for x-coordinate
    if xi <= rangex[0]:
        idx = 0
        distx = 0
    elif xi >= rangex[-1]:
        idx = -1
        distx = 0
    else:
        for i in xrange(dimx-1):
            if rangex[i] <= xi < rangex[i+1]:
                idx = i
                distx = (xi-rangex[i])/(rangex[i+1]-rangex[i])
    # determine weights for y-coordinate
    if yi <= rangey[0]:
        idy = 0
        disty = 0
    elif yi >= rangey[-1]:
        idy = -1
        disty = 0
    else:
        for i in xrange(dimy-1):
            if rangey[i] <= yi < rangey[i+1]:
                idy = i
                disty = (yi-rangey[i])/(rangey[i+1]-rangey[i])
    weights[0] = float(idx)
    weights[1] = float(distx)
    weights[2] = float(idy)
    weights[3] = float(disty)
    return weights

@njit
def lookup_xy(table, weights):
    idx = weights[0]
    distx = weights[1]
    idy = weights[2]
    disty = weights[3]
    return table[int(idx),int(idy)]*(1-distx)*(1-disty) +\
           table[int(idx+1),int(idy)]*distx*(1-disty)   +\
           table[int(idx),int(idy+1)]*(1-distx)*disty   +\
           table[int(idx+1),int(idy+1)]*distx*disty



# function for outside grid warnings.
# if mu/sigma get smallel/larger than min/max
# of the precalculated mui-sig rectangle a warning is shown
#@njit
def outside_grid_warning(xi, yi, rangex, rangey, when):
    # mu warnings
    if(xi < rangex[0]):
        print('--- OUTSIDE-GRID-WARNING: mu too low: ', xi, ' at time: ', when, 's')
    elif(xi > rangex[-1]):
        print('--- OUTSIDE-GRID-WARNING: mu too high: ', xi, ' at time: ', when, 's')
    # sigma warnings
    if(yi < rangey[0]):
        print('--- OUTSIDE-GRID-WARNING: sigma too low: ', yi, ' at time: ', when, 's')
    elif(yi > rangey[-1]):
        print('--- OUTSIDE-GRID-WARNING: sigma too high: ', yi, ' at time: ', when, 's')



# functions to compute synaptic mean and std of the synaptic input
@njit
def get_mu_syn(K, J, mu_ext, delay_type, step, r, r_d, n_d, taud, dt):
    # no delay
    if delay_type == 0:
        r_rec = r[step]
    # general constant delay
    elif delay_type == 1:
        r_rec = r[step-n_d]
    # exponential delay distribution
    elif delay_type == 2:
        r_d[step+1] = r_d[step] + dt*(r[step]-r_d[step])/taud
        r_rec = r_d[step]
    # exponential delay distribution + constant delay
    elif delay_type == 3:
        r_d[step+1] = r_d[step] + dt*(r[step-n_d]-r_d[step])/taud
        r_rec = r_d[step]
    else:
        raise NotImplementedError
    return mu_ext[step]+K*J*r_rec

@njit
def get_sigma_syn(K, J, sigma_ext, delay_type, step, r, r_d, n_d, taud, dt):
    #general constant delay
    if delay_type == 0:
        r_rec = r[step]
    elif delay_type == 1:
        r_rec = r[step-n_d]
    #exponential delay distribution
    elif delay_type == 2:
        r_d[step+1] = r_d[step] + dt*(r[step]-r_d[step])/taud
        r_rec = r_d[step]
    #exponential delay distribution + constant delay
    elif delay_type == 3:
        r_d[step+1] = r_d[step] + dt*(r[step-n_d]-r_d[step])/taud
        r_rec = r_d[step]
    else:
        raise NotImplementedError
    return sqrt(sigma_ext[step]**2 + K*J**2*r_rec)


# sample (one population) connection matrix with exactly k presynaptic contacts per neuron
def fixed_connectivity(n, k):
    prelist = np.zeros(k * n, dtype = int)
    postlist = np.zeros_like(prelist)
    for j in xrange(n):
        presynapses = choose_k_from_n(n, k)
        prelist[j * k:(j + 1) * k] = presynapses
        postlist[j * k:(j + 1) * k] = j * np.ones(k, dtype = int)
    
    return prelist, postlist

# chooses excatly k random numbers from 0 to n-1
@jit
def choose_k_from_n(n, k):
    # use vaguely estimated metric of when sorting random numbers is better
    if float(k) / float(n) > 0.125:
        ans = np.argsort(np.random.rand(n))[:k]
        return ans

    nums = range(n)
    swaps = (np.random.rand(k) * xrange(n, n - k, -1)).astype('int') + xrange(k)
    for i in xrange(k):
        # swap with some random element from here to end - these swap positions precalculated
        nums[i], nums[swaps[i]] = nums[swaps[i]], nums[i]

    ans = nums[:k]
    return ans


class SubplotRect(): 
    
    def __init__(self, rows, cols, current=1, create_axes=False):
        self.rows = rows
        self.cols = cols
        self.current = current
        
        # creating the axes
        if create_axes:
            for i in range(self.rows):
                for j in range(self.cols):
                    self.nextcol()
                self.nextrow()

    def next(self):
        self.current = (self.current + 1) % (self.cols*self.rows)
        
    def last(self):
        self.current = (self.current - 1) % (self.cols*self.rows)
            
    def nextcol(self, sharex=None):
        self.current = self.current % self.cols + 1
        self.current_axes(sharex)
    
    def nextrow(self, sharex=None):
        self.current = (self.current+self.cols -1) % (self.cols*self.rows) +1
        self.current_axes(sharex)     

    def first(self, sharex=None):
        self.current = 1
        self.current_axes()
    
    def current_axes(self, sharex=None):
        plt.subplot(self.rows, self.cols, self.current, sharex=sharex)
        

# returns the power spectral density for real-valued data x sampled with dt [ms] 
# using welch's method with half overlapping windows (with winsize # of data points)
# note that we use boxcar instead of the scipy default 'hanning'
def powerspec(x, dt, winsize, window='boxcar'):
    noverlap = min(winsize, len(x)) // 2 # make code compatible with short time series x
    freqs, psd_r = scipy.signal.welch(x, fs=1000.0/dt, nperseg=winsize, noverlap=noverlap, 
                                  return_onesided=True, window=window)                        
    psd_r *= 0.5 # because of onesided spectrum for the real data
    return freqs, psd_r


# computes the average of y with a moving rectangular window yielding y values 
# on the coarser x grid with equidistant spacing dx (center value interpretation)
def kernelavg(x, y, dx):
    
    dx_coarse = dx
    dx_fine = x[1]-x[0]
    assert abs(round(dx_coarse / (2*dx_fine)) - dx_coarse/(2*dx_fine)) < 1e-15 # enforce even multiple
    M = int(round(dx_coarse / (2*dx_fine)))

    # moving window
    x_avg = x[M-1:-M] + dx_fine/2. # begin after the first M intervals => t_mwin[0] = t[M-1 + 1/2]
    y_avg = np.convolve(y, np.ones(2*M)/float(2*M), mode='valid')
    # scipy version (equivalent)
    #y_avg = scipy.signal.convolve(y, scipy.signal.boxcar(2*M)/float(2*M), mode='valid')

    return x_avg, y_avg


# calculates non-overlapping averages of N consecutive elements of y assuming N to 
# be even and x equidistant. the elements of x are interpreted as centers of the 
# small intervals. the resampled center x is shifted by dx/2 to make its elements centers 
# of the respective enlarged intervals
def resample(x, y, N):
    assert N % 2 == 0
    dx = x[1] - x[0]
    
    x_resampled = x[N/2-1 : -N/2 : N] + dx/2.
    y_truncated_rowise = y[:len(x_resampled)*N].reshape(-1, N)
    y_resampled = np.mean(y_truncated_rowise, axis=1)
    
    return x_resampled, y_resampled

# resampling by dx (has to be a multiple of x[1]-x[0])
# uses resample() from above
def resample_dx(x, y, dx):
    dx_fine = x[1] - x[0]
    dx_coarse = dx
    assert abs(round(dx_coarse / dx_fine) - dx_coarse/dx_fine) < 1e-10
    N = int(round(dx_coarse / dx_fine))
    return resample(x, y, N)

# return left/right x,y arrays for step like line plotting from center 
# interpreted data with interval length dx
def steplike(x, y, dx):
    
    x_left = x-dx/2.0
    x_right = x+dx/2.0

    x_steps = np.zeros(2*len(x))
    x_steps[::2] = x_left
    x_steps[1::2] = x_right

    y_steps = np.repeat(y, 2)
    
    return x_steps, y_steps
    
def rectify_trace(t, f, name_f):
    negative_f = f < 0
    if negative_f.any():
        print('WARNING: rectifying negative values {}({}) = {}'.format(name_f, t[negative_f], f[negative_f]))
    f[negative_f] = 0.0
    
# todo: recursive comparison and fine-grained warning/verbose output control
def compare_dicts(d1, d2, ignore):
    for k1 in d1.keys():
        if k1 not in ignore:
            if k1 not in d2.keys():
                print('DEBUG: second dict does not contain: {}'.format(k1))
            else:
                if d1[k1] != d2[k1]:
                    print('DEBUG: the dictionaries differ in key {}: {} vs. {}'.format(k1, d1[k1], d2[k1]))
                
    for k2 in d2.keys():
        if k2 not in ignore and k2 not in d1.keys():
            print('DEBUG: first dict does not contain: {}'.format(k2))

# deprecated but maybe interesting stuff:

# resample the rate by increasing the bin size to dt assuming an equidistant t grid
def resample_rate_deprecated(t, rate, dt):    
    dt_data = t[1]-t[0] # equidistancy: dt_data is independent of t
    t_points = len(t)
    binsize = int(dt//dt_data)
    bins = int(t_points/binsize)  + (1 if t_points%binsize != 0 else 0)
    rate_resampled = zeros(bins)
    t_smoothed = zeros(bins)
    for k in range(bins):
        idx_left = k*binsize
        idx_right = (k+1)*binsize
        diff_idx_right = max(0, idx_right-t_points) # nonzero only if last interval
        if diff_idx_right > 0:
            idx_right -= diff_idx_right 
        binsize_k = idx_right-idx_left # the binsize of the last bin can be shorter
        rate_resampled[k] = rate[idx_left:idx_right].sum()/float(binsize_k)
        t_smoothed[k] = 0.5*(t[idx_left]+t[idx_right-1]) # center point of the bins
    return t_smoothed, rate_resampled
   

# integrate the rate over subintervals of size dt assuming an equidistant t grid 
# this decreases the sample size to dt 
def integrate_rate_subints_deprecated(t, rate, dt):    
    dt_data = t[1]-t[0] # equidistancy: dt_data is independent of t
    t_points = len(t)
    subint_size = int(dt//dt_data)
    subints = int(t_points/subint_size)  + (1 if t_points%subint_size > 1 else 0) # K*subint_size+1 is the equidist subint case, therefore modulo > 1
    rate_integral = zeros(subints)
    t_smoothed = zeros(subints)
    for k in range(subints):
        idx_left = k*subint_size
        idx_right = (k+1)*subint_size
        diff_idx_right = max(0, idx_right-t_points) # nonzero only if last interval
        if diff_idx_right > 0:
            idx_right -= diff_idx_right 
#        subint_size_k = idx_right-idx_left # the subint_size of the last bin can be shorter
        rate_integral[k] = sum(0.5*(rate[idx_left:idx_right-(0 if k<subints-1 else 1)]+rate[idx_left+1:idx_right+(1 if k<subints-1 else 0)]))*dt_data
        t_smoothed[k] = 0.5*(t[idx_left]+t[idx_right-(0 if k<subints-1 else 1)]) # center point of the subints
    return t_smoothed, rate_integral   


# if model_dt != min_dt interpolate input
def interpolate_input(external_input,params,model):
    mu_ext = external_input[0]
    sigma_ext = external_input[1]
    t_ext = params['t_ext']
    if model == 'reduced':
        if params['uni_dt'] > params['min_dt']:
            steps = int(params['runtime']/params['uni_dt'])
            t_interpolated = np.linspace(0., params['runtime'],steps+1)
            mu_ext = np.interp(t_interpolated,t_ext,mu_ext)
            sigma_ext = np.interp(t_interpolated,t_ext,sigma_ext)
    elif model == 'fp':
        if params['fp_dt'] > params['min_dt']:
            steps = int(params['runtime']/params['fp_dt'])
            t_interpolated = np.linspace(0., params['runtime'],steps+1)
            mu_ext = np.interp(t_interpolated, t_ext, mu_ext)
            sigma_ext = np.interp(t_interpolated, t_ext, sigma_ext)
    elif model == 'net':
        if params['net_dt'] > params['min_dt']:
            steps = int(params['runtime']/params['net_dt'])
            t_interpolated = np.linspace(0., params['runtime'], steps+1)
            mu_ext = np.interp(t_interpolated, t_ext, mu_ext)
            sigma_ext = np.interp(t_interpolated, t_ext, sigma_ext)
    return [mu_ext, sigma_ext]


