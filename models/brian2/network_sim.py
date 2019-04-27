from __future__ import print_function
from __future__ import division
from brian2 import *
import time
import os
import gc
import shutil
from utils_net import smooth_trace, fixed_connectivity





cpp_default_dir = 'brian2_compile'


def network_sim(signal, params, rec = False, standalone_dir = cpp_default_dir):

    # todo only one network simulation for a,b = const. and a(t),b(t)
    # combine PIF, EIF, LIF in one simulation code

    if params['brian2_standalone']:
        # build on run = False (changed for brian2_rc3)
        set_device(params['brian2_device'], build_on_run=False)
        device.insert_code('main', 'srand('+str(int(time.time())+os.getpid())+');')
        #standalone_dir = 'standalone_{}_pid{}'.format(time.strftime("%Y-=%m-%dT%H:%M:%S"), os.getpid())
    else:
        prefs.codegen.target = 'numpy'

    #stripe on brian units within the copied dict for the simulation so that brian can work with them
    dt_sim =        params['net_dt']*ms
    C =             params['C']*pF
    gL =            params['gL']*nS
    taum =          params['taum']*ms
    EL =            params['EL']*mV
    Ew =            params['Ew']*mV
    VT =            params['VT']*mV
    negVT = -VT
    deltaT =        params['deltaT']*mV
    Vcut =          params['Vcut']*mV
    tauw =          params['tauw']*ms
    Vr =            params['Vr']*mV
    t_ref =         params['t_ref']*ms
    net_w_init =  params['net_w_init']*pA
    rates_dt = params['net_record_dt'] *ms
    runtime = params['runtime']*ms
    time_r = np.arange(0., runtime/ms, rates_dt/ms)
    time_v = np.arange(0., runtime/ms, dt_sim/ms)
    #ratio used for smoothing the spikehistogram to match resolution
    ratesdt_dtsim=int(rates_dt/dt_sim)

    #N_total = params['N_total']

    N_e = params['N_e']
    N_i = params['N_i']
    N_total = N_e + N_i

    K_ee = params['K_ee']
    K_ii = params['K_ii']
    K_ei = params['K_ei']
    K_ie = params['K_ie']

    J_ee = params['J_ee']*mV/ms
    J_ii = params['J_ii']*mV/ms
    J_ei = params['J_ei']*mV/ms
    J_ie = params['J_ie']*mV/ms

    tau_se = 2.0 * ms
    tau_si = 5.0 * ms

    # set and rescale c's
    cee = params['cee'] * mV/ms  / np.abs(J_ee)
    cie = params['cie'] * mV/ms  / np.abs(J_ie)
    cii = params['cii'] * mV/ms  / np.abs(J_ii)
    cei = params['cei'] * mV/ms  / np.abs(J_ei)

    #print("cee={} cie={} cii={} cei={}".format(cee, cie, cii, cei))

    print("Simulating {} excitatory and {} inhibitory neurons (total {}) ...".format(N_e, N_i, N_total))
    # garbage collect so we can run multiple brian2 runs
    gc.collect()

    # seed our random number generator!  necessary on the unix server
    np.random.seed()

    mu_e_ext_array = signal[0][0] # [mV/ms]
    sigma_e_ext_array = signal[0][1] # [mV/sqrt(ms)]
    mu_i_ext_array = signal[1][0] # [mV/ms]
    sigma_i_ext_array = signal[1][1] # [mV/sqrt(ms)]


    # what is recorded during the simulation
    record_spikes = params['net_record_spikes']
    record_synapses = record_spikes
    record_w = params['net_record_w']
    record_v_example_trace = 0
    record_all_v_at_times = True if params['net_record_all_neurons'] else False
    record_v_stats = params['net_record_v_stats']
    record_w_stats = params['net_record_w_stats']

    # simulation timestep
    simclock = Clock(dt_sim)

    w_refr = ' (unless refractory)' if 'net_w_refr' not in params or params['net_w_refr'] else ''

    a = params['a']
    b = params['b']

    # convert to array if adapt params are scalar values
    if type(a) in [int,float]:
        a = np.ones_like(mu_e_ext_array)*a
    if type(b) in [int,float]:
        b = np.ones_like(mu_e_ext_array)*b

    # decide if there's adaptation
    have_adap = True if  (a.any() > 0.0)  or (b.any() > 0.0) else False
    #print("> adaptation: %s" % have_adap)

    # convert numpy arrays to TimedArrays
    a = TimedArray(a*nS, dt_sim)
    b = TimedArray(b*pA, dt_sim)

    #transform the external input into TimedArray
    mu_e_ext = TimedArray(mu_e_ext_array*(mV/ms), dt_sim)
    sigma_e_ext = TimedArray(sigma_e_ext_array*(mV/sqrt(ms)), dt_sim)

    mu_i_ext = TimedArray(mu_i_ext_array*(mV/ms), dt_sim)
    sigma_i_ext = TimedArray(sigma_i_ext_array*(mV/sqrt(ms)), dt_sim)    

    #print("inputs: exc = {} mV/ms inh = {} mV/ms".format(mu_e_ext_array[0], mu_i_ext_array[0]))

    # -----------------------
    # -   MODEL EQUATIONS   -
    # -----------------------

    #get the model specific term EIF/PIF
    if params['neuron_model'] == 'EIF' :
        model_term = '((EL - v) + deltaT * exp((negVT + v) / deltaT)) / taum'
    elif params['neuron_model'] == 'PIF':
        model_term = ''
    elif params['neuron_model'] == 'LIF':
        model_term = '(EL - v) / taum'
    else:
        mes = 'The model "{}" has not been implemented yet. For options see params dict.'.format(params['neuron_model'])
        raise NotImplementedError(mes)

    model_eqs_e = '''
        dv/dt = %s %s + mu_e_ext(t) + J_ee * see + J_ei * sei + sigma_e_ext(t) * xi_exc : volt (unless refractory)
        %s
        dsee/dt = -see/tau_se : 1 (unless refractory) 
        dsei/dt = -sei/tau_si : 1 (unless refractory) 
        ''' % (model_term,'- (w / C)' if have_adap else '', ('dw/dt = (a(t) * (v - Ew) - w) / tauw : amp %s' % w_refr)
               if have_adap else '')
    model_eqs_i = '''
        dv/dt = %s + mu_i_ext(t) + J_ii * sii + J_ie * sie + sigma_i_ext(t) * xi_inh : volt (unless refractory)
        dsii/dt = -sii/tau_si : 1 (unless refractory) 
        dsie/dt = -sie/tau_se : 1 (unless refractory) 
        ''' % (model_term)     




    # -------------
    # -- NETWORK --
    # -------------
    # initialize Neuron group
    GE = NeuronGroup(N = N_e, model = model_eqs_e,
                    threshold = 'v > Vcut',clock = simclock,
                    reset = 'v = Vr%s' % ('; w += b(t)' if have_adap else ''),
                    refractory = t_ref, method = params['net_integration_method'])
    GI = NeuronGroup(N = N_i, model = model_eqs_i,
                    threshold = 'v > Vcut',clock = simclock,
                    reset = 'v = Vr',   
                    refractory = t_ref, method = params['net_integration_method'])

    # initialize PopulationRateMonitor
    rate_monitor_e = PopulationRateMonitor(GE, name = 'aeif_ratemon_e')
    rate_monitor_i = PopulationRateMonitor(GI, name = 'aeif_ratemon_i')

    # intitialize net
    Net = Network(GE, rate_monitor_e, GI, rate_monitor_i)

    if rec:
        #print('building synapses...')
        start_synapses = time.time()
        J = params['J']*mV
        K = params['K']

        # ----------------
        # -   SYNAPSES   -
        # ----------------

        #synapses object
        #synapses from G --> G! only one population
        #this only specifies the dynamics of the synapses. they get acctually created when the .connect method is called

        Syn_EE = Synapses(GE,GE, on_pre = 'see += (1-see) * cee ')
        Syn_II = Synapses(GI,GI, on_pre = 'sii += (1-sii) * cii ')
        Syn_EI = Synapses(GE,GI, on_pre = 'sie += (1-sie) * cie ')
        Syn_IE = Synapses(GI,GE, on_pre = 'sei += (1-sei) * cei ')
        

        #sparsity_E = float(K)/N_total
        #assert 0 <= sparsity <= 1.0
        #sparsity = float(K)/N_total
        #assert 0 <= sparsity <= 1.0        

        # connectivity type
        if params['connectivity_type'] == 'binomial':
            Syn_EE.connect(True, p = K_ee/N_e)
            Syn_II.connect(True, p = K_ii/N_i)
            Syn_EI.connect(True, p = K_ei/N_e)
            Syn_IE.connect(True, p = K_ie/N_i)

        elif params['connectivity_type'] == 'fixed':
            prelist_ee, postlist_ee = fixed_connectivity(N_e, K_ee)
            prelist_ii, postlist_ii = fixed_connectivity(N_i, K_ii)
            prelist_ei, postlist_ei = fixed_connectivity(N_i, K_ei)
            prelist_ie, postlist_ie = fixed_connectivity(N_e, K_ie)
            Syn_EE.connect(i=prelist_ee, j=postlist_ee)
            Syn_II.connect(i=prelist_ii, j=postlist_ii)
            Syn_EI.connect(i=prelist_ei, j=postlist_ei)
            Syn_IE.connect(i=prelist_ie, j=postlist_ie)

        elif params['connectivity_type'] == 'random':
            print("Creating random network: p_ee={} p_ie={} p_ii={} p_ei={}".format(K_ee/N_e, K_ii/N_e, K_ei/N_i, K_ie/N_i))
            Syn_EE.connect(p=K_ee/N_e)
            Syn_II.connect(p=K_ii/N_i)
            Syn_EI.connect(p=K_ei/N_i)
            Syn_IE.connect(p=K_ie/N_e)
        else:
            raise Exception('Synapses are not connected.')

        # --------------
        # -   DELAYS   -
        # --------------
        # no delay; nothing has to be implemented
        if params['delay_type'] == 0:
            pass
        # constant delay
        elif params['delay_type'] == 1:
            Syn_EE.delay = '{} * ms'.format(params['const_delay_e'])
            Syn_II.delay = '{} * ms'.format(params['const_delay_i'])
            Syn_EI.delay = '{} * ms'.format(params['const_delay_i'])
            Syn_IE.delay = '{} * ms'.format(params['const_delay_e'])
        # exponentially distributed delays
        elif params['delay_type'] == 2:
            #taud is the mean=standard deviation of the distribution
            Syn_EE.delay = '-log(rand()) * {} * ms'.format(params['taud'])
            Syn_II.delay = '-log(rand()) * {} * ms'.format(params['taud'])
            Syn_EI.delay = '-log(rand()) * {} * ms'.format(params['taud'])
            Syn_IE.delay = '-log(rand()) * {} * ms'.format(params['taud'])
        # exp. delay dist. + const. delay
        elif params['delay_type'] == 3:
            Syn_EE.delay = ('({} -log(rand()) * {}) * ms'.format(format(params['const_delay']), params['taud']))
            Syn_II.delay = ('({} -log(rand()) * {}) * ms'.format(format(params['const_delay']), params['taud']))
            Syn_EI.delay = ('({} -log(rand()) * {}) * ms'.format(format(params['const_delay']), params['taud']))
            Syn_IE.delay = ('({} -log(rand()) * {}) * ms'.format(format(params['const_delay']), params['taud']))
        # add synapses to the network
        else:
            raise NotImplementedError
        Net.add(Syn_EE)
        Net.add(Syn_II)
        Net.add(Syn_IE)
        Net.add(Syn_EI)

        #print('build synapses time: {}s'.format(time.time()-start_synapses))

    # ----------------------------------------
    # - INITIAL VALUES & BOUNDARY CONDITIONS -
    # ----------------------------------------

    #initial distribution of the network simulation
    if params['net_v_init'] == 'delta':
        GE.v = np.ones(len(GE)) * params['net_delta_peak']*mV
        GI.v = np.ones(len(GI)) * params['net_delta_peak']*mV
    elif params['net_v_init'] == 'normal':
        GE.v = params['net_normal_sigma'] * np.random.randn((len(GE))) * mV + params['net_normal_mean'] * mV
        GI.v = params['net_normal_sigma'] * np.random.randn((len(GI))) * mV + params['net_normal_mean'] * mV
    elif params['net_v_init'] == 'uniform':
        len_interval = Vcut - Vr
        GE.v = np.random.rand(len(GE)) * len_interval + Vr
        GI.v = np.random.rand(len(GI)) * len_interval + Vr

    # s init test
    GE.see = np.ones(len(GE)) * params['net_s_init']
    GE.sei = np.ones(len(GE)) * params['net_s_init']
    GI.sie = np.ones(len(GI)) * params['net_s_init']
    GI.sii = np.ones(len(GI)) * params['net_s_init']

    # initial distribution of w_mean
    if have_adap:
        # standart deviation of w_mean is set to 0.1
        GE.w = 0.1 * np.random.randn(len(GE)) * pA + net_w_init

    # include a lower bound for the membrane voltage
    if 'net_v_lower_bound' in params and params['net_v_lower_bound'] is not None:
        #new in Brian2.0b4: custom_operation --> run_regularly
        V_lowerbound = GE.run_regularly('v = clip(v, %s * mV, 10000 * mV)'
                                       % float(params['net_v_lower_bound']),
                                when = 'end', order = -1, dt = dt_sim)
        print('Lower bound active at {}'.format(params['net_v_lower_bound']))
        Net.add(V_lowerbound)

    # make sure that synapse activity is between 0 and 1

    # --------------------
    # -     MONITORS     - 
    # --------------------



    if record_v_example_trace:
        print("recoding v test {}".format(record_v_example_trace))
        # record v from params['net_record_v_trace'] neurons
        v_monitor_example_trace = StateMonitor(GE, 'v',
                                 record = range(min(params['net_record_example_v_traces']
                                                    ,N_total)))
        Net.add(v_monitor_example_trace)


    if record_all_v_at_times:
        # define Clock which runs on a very course time grid (memory issue)
        clock_record_all = Clock(params['net_record_all_neurons_dt']*ms)
        v_monitor_record_all = StateMonitor(GE, 'v', record=True, clock = clock_record_all)
        Net.add(v_monitor_record_all)

    if record_v_stats:
        # create statistics neuton 'group' with 1 neuron
        eqs_Gvstats = '''
                        v_mean : volt
                        v_var  : volt**2
                    '''
        Gvstats = NeuronGroup(1, eqs_Gvstats, clock = simclock)
        # connect the two neuron groups with Synapses
        eqs_Sstats = '''
                        v_mean_post = v_pre/N_pre             : volt    (summed)
                        v_var_post  = (v_pre-v_mean)**2/N_pre : volt**2 (summed)
                     '''
        Svstats = Synapses(GE, Gvstats, eqs_Sstats)
        Svstats.connect(True)
        VstatsMon = StateMonitor(Gvstats, ['v_mean', 'v_var'], record = True)
        # maybe add this a bit later ... together with all the other stuff
        Net.add(Svstats)
        Net.add(Gvstats)
        Net.add(VstatsMon)

    if record_w_stats and have_adap:
        eqs_Gwstats = '''
                      w_mean : amp
                      w_var  : amp**2
                      '''
        Gwstats = NeuronGroup(1, eqs_Gwstats, clock = simclock)
        eqs_Stats = '''
                    w_mean_post = w_pre/N_pre               : amp (summed)
                    w_var_post = (w_pre-w_mean)**2/N_pre    : amp**2 (summed)
                    '''
        Swstats = Synapses(GE, Gwstats, eqs_Stats)
        Swstats.connect(True)
        WstatsMon = StateMonitor(Gwstats, ['w_mean', 'w_var'], record=True)
        Net.add(Swstats)
        Net.add(Gwstats)
        Net.add(WstatsMon)


    if record_spikes > 0:
        record_spikes_group_e = Subgroup(GE, 0, min(record_spikes, N_e))
        record_spikes_group_i = Subgroup(GI, 0, min(record_spikes, N_i))
        spike_monitor_e = SpikeMonitor(record_spikes_group_e, name = 'aeif_spikemon_e')
        spike_monitor_i = SpikeMonitor(record_spikes_group_i, name = 'aeif_spikemon_i')

        Net.add(spike_monitor_e, record_spikes_group_e, spike_monitor_i, record_spikes_group_i)


    if record_w > 0 and have_adap:
        record_w_group = Subgroup(GE, 0, min(record_w, N_total))
        w_monitor = StateMonitor(record_w_group, 'w', record = range(params['net_record_w']),
                                 dt = dt_sim, name = 'aeif_wmon')
        Net.add(w_monitor, record_w_group)


    if record_synapses:
        record_synapses_e = Subgroup(GE, 0, min(record_spikes, N_e))
        see_monitor = StateMonitor(record_synapses_e, 'see', record = range(record_spikes), 
                                dt = dt_sim, name = 'see_mon')
        sei_monitor = StateMonitor(record_synapses_e, 'sei', record = range(record_spikes), 
                                dt = dt_sim, name = 'sei_mon')
        record_synapses_i = Subgroup(GI, 0, min(record_spikes, N_i))
        sii_monitor = StateMonitor(record_synapses_i, 'sii', record = range(record_spikes), 
                                dt = dt_sim, name = 'sii_mon')
        sie_monitor = StateMonitor(record_synapses_i, 'sie', record = range(record_spikes), 
                                dt = dt_sim, name = 'sie_mon')
        Net.add(see_monitor, record_synapses_e)
        Net.add(sei_monitor, record_synapses_e)
        Net.add(sii_monitor, record_synapses_i)
        Net.add(sie_monitor, record_synapses_i)
        

    # ------------------------
    # -    RUN SIMULATION    -
    # ------------------------

    #print('Running network...')
    start_time = time.time()
    Net.run(runtime, report = 'text')

    if params['brian2_standalone']:
        project_dir = standalone_dir + '/sim' + str(os.getpid())
        device.build(directory = project_dir, compile=True, run=True)


    # -----------------
    # -    RESULTS    -
    # -----------------

    # extract results

    # function for binnung the population rate trace
    def binning(arr, N):
        if int(N) in [0, 1]:
            return arr
        else:
            len_return = int(len(arr)/N)
            return np.array([np.mean(arr[k:min(k+N,len(arr))]) for k in range(0, len(arr), N)])

    # unbinned quantites
    net_rates_e = rate_monitor_e.rate/Hz
    net_t_e = rate_monitor_e.t/ms

    net_rates_i = rate_monitor_i.rate/Hz
    net_t_i = rate_monitor_i.t/ms




    if record_v_example_trace > 0:
        v_neurons = v_monitor_example_trace.v/mV
        t_v = v_monitor_example_trace.t/ms

    # old way of saving wm
    if record_w > 0 and have_adap:
        net_w = w_monitor.w / pA
        net_wt = w_monitor.t/ms

    if record_spikes > 0:
        # multiply by 1 like this to ensure brian extracts the results before we delete the compile directory
        net_spikes = spike_monitor_e.it
        i, t = net_spikes
        i = i * 1; t = t * 1
        net_spikes = [i, t]
        print("E spikes: {} I spikes: {}".format(spike_monitor_e.num_spikes, spike_monitor_i.num_spikes))

    if record_synapses > 0:
        net_see = see_monitor.see 
        net_sii = sii_monitor.sii 
        net_sei = sei_monitor.sei 
        net_sie = sie_monitor.sie 
        net_t = see_monitor.t / ms

    if record_v_stats:
        v_mean = VstatsMon.v_mean[0]/mV
        v_var = VstatsMon.v_var[0]/mV**2
        v_std = np.sqrt(v_var)

    if record_w_stats and have_adap:
         w_mean = WstatsMon.w_mean[0]/pA
         w_var = WstatsMon.w_var[0]/pA**2
         w_std = np.sqrt(w_var)

    if record_all_v_at_times:
        v_all_neurons = v_monitor_record_all.v/mV
        t_all_neurons = v_monitor_record_all.t/ms



    run_time = time.time() - start_time
    print('runtime: %1.1f' % run_time)




    #for smoothing function net_rates do: helpers.smooth_trace(net_rates, int(rates_dt / dt_sim))
    # smooth out our hyper-resolution rate trace manually cause brian2 can't do it
    results_dict = {'brian_version':2, 't':time_r}
    results_dict['r'] = smooth_trace(net_rates_e, ratesdt_dtsim)
    results_dict['r_e'] = smooth_trace(net_rates_e, ratesdt_dtsim)
    results_dict['r_i'] = smooth_trace(net_rates_i, ratesdt_dtsim),

#     results_dict = {'brian_version':2, 'r':net_rates, 't':net_t}
    # print(len(results_dict['t']))
    # time binnig

    if params['brian2_standalone']:
        shutil.rmtree(project_dir)
        device.reinit()

    if record_synapses > 0:
        results_dict['s_t'] = net_t
        results_dict['see'] = net_see
        results_dict['sei'] = net_sei
        results_dict['sie'] = net_sie
        results_dict['sii'] = net_sii


    if record_v_example_trace > 0:

        results_dict['v'] = v_neurons
        results_dict['t_v'] = t_v


    if record_w > 0 and have_adap:
        results_dict['net_w'] = net_w
        results_dict['net_wt'] = net_wt

    # if record_w > 0 and have_adap:
    #     results_dict['net_w_samples'] = net_w[:min(10, np.size(net_w, 0)), :]
    #     results_dict['wm'] = np.mean(net_w, 0)
    #     # also include
    #     results_dict['w_std'] = np.std(net_w, 0)
    #     results_dict['net_w_dt'] = rates_dt

    if record_spikes > 0:
        results_dict['net_spikes'] = net_spikes
        results_dict['modelname'] = 'net'
    if record_v_stats:
        results_dict['v_mean'] = v_mean
        results_dict['v_var'] = v_var
        results_dict['v_std'] = v_std
    if record_w_stats and have_adap:
        # maybe change this back again ....
        results_dict['wm'] = smooth_trace(w_mean, ratesdt_dtsim)
        results_dict['w_var'] = smooth_trace(w_var, ratesdt_dtsim)
        results_dict['w_std'] = smooth_trace(w_std , ratesdt_dtsim)
    if record_all_v_at_times:
        results_dict['v_all_neurons'] = v_all_neurons
        results_dict['t_all_neurons'] = t_all_neurons
    return results_dict
