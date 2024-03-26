from ..imports import *
import time

#changed nr of tested durations to 4(from 10): 0.01, 0.04, 0.07, 0.1

def full_injection_recovery(self,
                            nfake=10,
                            pool=False,
                            minimum_planet_radius=0.5 * u.R_earth,
                            maximum_planet_radius=3 * u.R_earth,
                            minimum_period=0.5 * u.d,
                            maximum_period=10 * u.d,
                            ld=[0.385, 0.304],
                            clean_kw = {'dust_removal':False, 'bad_weather_removal':False, 'cosmics_removal':True,
                                        'cosmic_boxsize':0.08,'cosmic_nsigma':3},
                            detrend_method = "gp",
                            detrend_bin = 20 * u.minute,
                            detrend_kw = {'do_first_sigma_clip':True, 'do_second_sigma_clip':True,
                                     'running_mean_boxsize':0.08, 'nsigma':3, 'plot':False},
                            bls_kw = {"minimum_period":0.5, "maximum_period":10,
                                      'transit_durations':np.linspace(0.01, 0.1, 4), 'plot':False, 'verbose': False},
                            bls_bin=7.5 * u.minute,
                            recovery_kw = {'condition_on_epoch':1 * u.hour},
                            plot=False,
                            verbose=False,
                            time_this_process=False,
                            svname="injected_planets.csv",
                            ):

    lcs_with_transits, planets = self.inject_lots_of_transits(nfake=nfake, pool=pool,
                                                              minimum_planet_radius=minimum_planet_radius,
                                                              maximum_planet_radius=maximum_planet_radius,
                                                              minimum_period=minimum_period,
                                                              maximum_period=maximum_period,
                                                              ld=ld,
                                                              fname=svname)

    bls_kw['verbose'] = verbose
    clean_lcs, detrend_lcs, bls_lcs = [],[],[]

    for i, lc in enumerate(lcs_with_transits):
        try:
            if time_this_process:
                t0 = time.time()

            print(f"{i+1}/{len(lcs_with_transits)}...")
            planets = pd.read_csv(svname)
            planets.loc[i, 'injected'] = 1.0
            planets.loc[i, 'observed'] = int(lc.was_planet_observed())
            clean_targ, detrend_targ, bls_targ, planets = self.single_injection_recovery(lc=lc, planets=planets, i=i,
                                                                                         clean_kw=clean_kw, detrend_method=detrend_method,
                                                                                         detrend_bin=detrend_bin, detrend_kw=detrend_kw,
                                                                                         bls_kw=bls_kw, bls_bin=bls_bin,
                                                                                         recovery_kw=recovery_kw, plot=plot,
                                                                                         verbose=verbose)
            planets.to_csv(svname, index=False)
            if time_this_process:
                t1 = time.time()
                print(f"Time to inject-recover: {t1-t0}")
            clean_lcs.append(clean_targ)
            detrend_lcs.append(detrend_targ)
            bls_lcs.append(bls_targ)
        except Exception as e:
            print(e)


    # print summary
    print(f"Planets recovered: {100 * len(planets.loc[planets['recovered'] == 1.0]) / len(planets['recovered'])}%")
    print("But not all of those planets transited during the observation window...")
    if len(planets['recovered'][planets['observed'] == 1.0]) > 0:
        print(f"""Observed Planets recovered: {100 * len(planets.loc[(planets['recovered'] == 1.0) &
                (planets['observed'] == 1.0)]) / len(planets['recovered'][planets['observed'] == 1.0])}%""")
    else:
        print("None of the planets injected were observed!")

    return lcs_with_transits, clean_lcs, detrend_lcs, bls_lcs, planets

def single_injection_recovery(self, lc, planets, i, clean_kw, detrend_bin, detrend_kw, bls_kw, bls_bin, recovery_kw,
                              detrend_method, predetrend_bls=True, plot=False, verbose=False,):

    if plot:
        ax = self.plot(color='C0', label='raw data')
        lc.plot(ax=ax, ylims=[0.9, 1.1], color='C1', label='injected transit')
        plt.show()

    # clean
    clean_targ = lc.clean(**clean_kw)

    if detrend_bin is not None:
        # bin before detrending
        bin_targ = clean_targ.bin(dt=detrend_bin)
    else:
        bin_targ = clean_targ

    if predetrend_bls:
        # do an initial search for transits to mask during detrending
        orig_bin_targ = bin_targ._create_copy()
        removed_nans = bin_targ.remove_nans()
        bls_targs = removed_nans.find_transits(**bls_kw)
        if len(bls_targs) > 0:
            in_transit = bls_targs[0].metadata['BLS_transits_ind']
            bls_targs[0].masks['transit'] = np.zeros(bls_targs[0].ntime)
            bls_targs[0].masks['transit'][in_transit] = 1.0
            init_transits_masked = bls_targs[0].clean(zero_flux_removal=False, nan_flux_removal=False,
                                                      bad_weather_removal=False, threshold_removal=False,
                                                      dust_removal=False, cosmics_removal=False)
            bin_targ = init_transits_masked

    # detrend the light curve
    if detrend_method == "gp":
        # gp detrend
        gp_targ = bin_targ.gp_detrend(**detrend_kw)

        # we can also predict the GP for 7.5 min binning and use that for the BLS
        bin_targ = clean_targ.bin(dt=bls_bin)
        gp_model = 1 + \
                   gp_targ.metadata['gp'].predict(y=gp_targ.metadata['data_to_condition_gp'], t=bin_targ.time.value)[0]
        detrended_targ = bin_targ._create_copy()
        detrended_targ.timelike['flux'] = detrended_targ.timelike['flux'] / gp_model

    elif detrend_method == "lsq":
        detrended_targ = bin_targ.lsq_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)

    elif detrend_method == "mcmc":
        detrended_targ = bin_targ.mcmc_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)

    elif detrend_method == "ridge":
        if predetrend_bls:
            detrended_targ = bin_targ.lsq_ridge_detrend_each_night(orig_lc=orig_bin_targ,**detrend_kw)
        else:
            detrended_targ = bin_targ.lsq_ridge_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)
    elif detrend_method==None:
        print("No detrending method selected!")
        detrended_targ = clean_targ._create_copy().bin(dt=bls_bin)
    else:
        print(f"{detrend_method} is not recognised. Please choose one of: 'gp', 'lsq', 'ridge' or 'mcmc'")
        return None, None, None, None


    if plot:
        ax = bin_targ.plot(color='C0', label='clean lc')
        detrended_targ.plot(ax=ax, ylims=[0.9, 1.1], color='C1', label=f'{detrend_method}-detrended')
        plt.show()

    # search for transit
    removed_nans = detrended_targ.remove_nans()
    bls_targs = removed_nans.find_transits(**bls_kw)

    if len(bls_targs) == 0:
        if verbose:
            print("No transit found!\n")
    else:
        recovered = False
        for bls_targ in bls_targs:
            # determine whether the injected planet was recovered
            if bls_targ.metadata['BLS_transits_found'] == True:

                if verbose:
                    print("Transit was found - checking if it matches the injected transit!")

                rec = bls_targ.was_injected_planet_recovered(**recovery_kw)
                bls_targ.metadata['recovery'] = rec

                # loop over all detected transits (in this case most likely 1):
                for r in range(len(bls_targ.metadata['BLS_transits_params']['depth'])):
                    # !!! the following assumes we have only injected 1 planet at a time !!!:
                    if rec[0][r]:
                        recovered = True
                if recovered:
                    # total_recovered += 1
                    planets.loc[i, 'recovered'] = 1.0
                    planets.loc[i, 'log_Prec'] = bls_targ.metadata['BLS_transits_params']['period'][0].to_value('d')
                    planets.loc[i, 'rec_depth'] = bls_targ.metadata['BLS_transits_params']['depth'][0]
                    planets.loc[i, 'rec_duration'] = bls_targ.metadata['BLS_transits_params']['duration'][0].to_value('d')
                    planets.loc[i, 'rec_epoch'] = bls_targ.metadata['BLS_transits_params']['epoch'][0].to_value('d')
                    planets.loc[i, 'snr'] = np.max(bls_targ.metadata['BLS_transits_params']['snr'])#[0]

                    if verbose:
                        print("Planet was successfully recovered")
                else:
                    # planets.loc[i, 'snr'] = np.max(bls_targ.metadata['BLS_transits_params']['snr'])
                    if verbose:
                        print("Planet was not successfully recovered")
        else:
            if verbose:
                print("No transit found!\n")


    # bls_lcs.append(bls_targ)

    return clean_targ, detrended_targ, bls_targs, planets

def was_injected_planet_recovered(self, condition_on_depth=None, condition_on_overlap=None, condition_on_epoch=None,
                                  condition_on_period=None, condition_on_snr=None):
    """
    Returns a list of booleans whether each transit injected into the light curve was recovered by BLS based on user-
    defined conditions.
    :param self:
    :param condition_on_depth: Fraction from 0-1 of injected depth to be recovered (e.g. 0.5 means the recovered depth
    must be at least 0.5*injected depth)
    :param condition_on_overlap: Fraction from 0-1 of injected depth to be recovered (e.g. 0.5 means the recovered
    depth must be at least 50% of injected depth)
    :param condition_on_epoch: Time value. The recovered epoch must be within injected epoch +/- condition_on_epoch to
    be recovered.
    :param condition_on_period: Fraction from 0-1 of injected period that needs to be recovered (e.g. 0.1 means the
    recovered period must match the injected period to within 10%)
    :param condition_on_snr: SNR threshold for a transit to be recovered successfully.
    :return:
    """
    recovered_all_planets = []
    recovered_all_transits = []

    injected_params = self.metadata['injected_planet']
    recovered_params = self.metadata['BLS_transits_params']

    # loop over all planet transits injected
    for planet in range(len(injected_params['depth'])):

        transit_t = injected_params['epoch'][planet].to_value('d')
        all_epochs = []
        while transit_t < self.time.value[-1]:
            all_epochs.append(transit_t)
            transit_t += injected_params['period'][0].to_value('d')

        # loop over all recovered transits
        for transit in range(len(recovered_params['depth'])):
            recovered = True

            if len(all_epochs) == 0:
                # if the planet occurred after the end of the observation window:
                recovered = False
                recovered_all_transits.append(recovered)
                continue

            if condition_on_depth is not None:
                if recovered_params['depth'][transit] < (condition_on_depth * injected_params['depth'][planet]):
                    recovered = False

            if condition_on_overlap is not None:
                recovered = False
                # for a in all_epochs:
                closest_transit = find_nearest(all_epochs, recovered_params['epoch'][transit].to_value('d'))
                inj_transit_start = all_epochs[closest_transit]*u.d - (0.5*injected_params['duration'][planet])
                inj_transit_end = all_epochs[closest_transit]*u.d + (0.5 * injected_params['duration'][planet])
                overlap = min(inj_transit_end, recovered_params['epoch_end'][transit]) - \
                              max(inj_transit_start, recovered_params['epoch_start'][transit])

                if overlap < (condition_on_overlap * injected_params['duration']):
                    recovered = False

            if condition_on_epoch is not None:
                # print(all_epochs, recovered_params['epoch'][transit], transit_t, self.time.value[-1])
                # print(self.metadata)
                closest_transit = find_nearest(all_epochs, recovered_params['epoch'][transit].to_value('d'))
                if abs(recovered_params['epoch'][transit] - all_epochs[closest_transit]*u.d) > \
                        condition_on_epoch:
                    recovered = False

            if condition_on_period is not None:
                period_sway = condition_on_period * injected_params['period'][planet]
                if (recovered_params['period'][transit] < injected_params['period'][planet] - period_sway) or \
                        (recovered_params['period'][transit] > injected_params['period'][planet] + period_sway):
                    recovered = False

            if condition_on_snr is not None:
                if recovered_params['snr'][transit] < condition_on_snr:
                    recovered = False

            recovered_all_transits.append(recovered)
        recovered_all_planets.append(recovered_all_transits)

    return recovered_all_planets

def was_planet_observed(self, fraction_overlap=0.5, planet_i=0):

    # set times of the first transit in the observation span
    duration = self.metadata['injected_planet']['duration'][planet_i]
    period = self.metadata['injected_planet']['period'][planet_i]
    transit_mid = self.metadata['injected_planet']['epoch'][planet_i]
    transit_start = transit_mid - (0.5*duration)
    transit_end = transit_mid + (0.5*duration)

    observed = 0

    days = self.split_lightcurve(split_every=0.5*u.d)[0]
    for day in days:
        while transit_start <= day[-1]:
            if transit_end >= day[0]:
                overlap = min(transit_end, day[-1]) - max(transit_start, day[0])
                if overlap >= (fraction_overlap * duration):
                    observed = 1
            transit_start = transit_start + period
            transit_end = transit_end + period
    return bool(observed)

def split_lightcurve(self, split_every=0.5*u.d):
    t = self.time.value * u.d
    start = t[0]
    nextdays = t[np.absolute(t - start) > split_every]
    split = []

    while nextdays != []:
        start = nextdays[0]
        ind_st = np.where(t == start)[0][0]
        split.append(ind_st)
        time = t[ind_st:]
        nextdays = time[np.absolute(time - start) > split_every]

    times = np.split(t, split)

    return times, split