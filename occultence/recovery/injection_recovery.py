from ..imports import *

def full_injection_recovery(self,
                            nfake=10,
                            pool=False,
                            minimum_planet_radius=0.5 * u.R_earth,
                            maximum_planet_radius=3 * u.R_earth,
                            minimum_period=0.5 * u.d,
                            maximum_period=0.5 * u.d,
                            clean_kw = {'dust_removal':False, 'bad_weather_removal':True, 'cosmics_removal':True,
                                        'cosmic_boxsize':0.08,'cosmic_nsigma':3},
                            gp_bin = 20 * u.minute,
                            gp_kw = {'do_first_sigma_clip':True, 'do_second_sigma_clip':True,
                                     'running_mean_boxsize':0.08, 'nsigma':3, 'plot':False},
                            bls_kw = {"minimum_period":0.5, "maximum_period":10,
                                      'transit_durations':np.linspace(0.01, 0.1, 10), 'plot':True},
                            bls_bin=7.5 * u.minute,
                            recovery_kw = {'condition_on_epoch':1 * u.hour},
                            plot=False,
                            ):

    lcs_with_transits, planets = self.inject_lots_of_transits(nfake=nfake, pool=pool,
                                                              minimum_planet_radius=minimum_planet_radius,
                                                              maximum_planet_radius=maximum_planet_radius,
                                                              minimum_period=minimum_period,
                                                              maximum_period=maximum_period, )

    total_injected = len(lcs_with_transits)
    total_recovered = 0

    clean_lcs, gp_lcs, bls_lcs = [],[],[]

    for i, lc in enumerate(lcs_with_transits):
        if plot:
            ax = self.plot()
            lc.plot(ax=ax, ylims=[0.9, 1.1])
            plt.show()
        # clean
        clean_targ = lc.clean(**clean_kw)
        clean_lcs.append(clean_targ)
        # bin 20 mins for the GP
        bin_targ = clean_targ.bin(dt=gp_bin)
        # gp detrend
        gp_targ = bin_targ.gp_detrend(**gp_kw)
        gp_lcs.append(gp_targ)

        # we can also predict the GP for 7.5 min binning and use that for the BLS
        bin_targ = clean_targ.bin(dt=bls_bin)
        gp_model = 1 + \
                   gp_targ.metadata['gp'].predict(y=gp_targ.metadata['data_to_condition_gp'], t=bin_targ.time.value)[0]
        gp_smaller_binning = bin_targ._create_copy()
        gp_smaller_binning.timelike['flux'] = gp_smaller_binning.timelike['flux'] / gp_model

        if plot:
            ax = bin_targ.plot()
            gp_smaller_binning.plot(ax=ax, ylims=[0.9, 1.1])
            plt.show()

        # search for transit
        bls_targ = gp_smaller_binning.find_transits(**bls_kw)

        # determine whether the injected planet was adequately recovered
        if bls_targ.metadata['BLS_transits_found'] == True:
            recovered = False
            print("Transit found!")
            # print(bls_targ.metadata, "\n")
            rec = bls_targ.was_injected_planet_recovered(**recovery_kw)
            # loop over all recovered transits (in this case most likely 1):
            for r in range(len(bls_targ.metadata['BLS_transits_params']['depth'])):
                if rec[0][r]:
                    recovered = True
            if recovered:
                total_recovered += 1
                planets['recovered'][i] = 1.0
                planets['log_Prec'][i] = bls_targ.metadata['BLS_transits_params']['period'][0].to_value('d')
                planets['rec_depth'][i] = bls_targ.metadata['BLS_transits_params']['depth'][0]
                planets['rec_duration'][i] = bls_targ.metadata['BLS_transits_params']['duration'][0].to_value('d')
                planets['rec_epoch'][i] = bls_targ.metadata['BLS_transits_params']['epoch'][0].value
                planets['snr'][i] = bls_targ.metadata['BLS_transits_params']['snr'][0]
        else:
            print("No transit found!\n")

        bls_lcs.append(bls_targ)

    return lcs_with_transits, clean_lcs, gp_lcs, bls_lcs, planets

def was_injected_planet_recovered(self, condition_on_depth=None, condition_on_overlap=None, condition_on_epoch=None,
                                  condition_on_period=None):
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
    :return:
    """
    recovered_all_planets = []
    recovered_all_transits = []

    injected_params = self.metadata['injected_planet']
    recovered_params = self.metadata['BLS_transits_params']

    for planet in range(len(injected_params['depth'])):
        for transit in range(len(recovered_params['depth'])):
            recovered = True

            if condition_on_depth is not None:
                if recovered_params['depth'][transit] < (condition_on_depth * injected_params['depth'][planet]):
                    recovered=False

            if condition_on_overlap is not None:
                inj_transit_start = injected_params['epoch'][planet] - (0.5*injected_params['duration'][planet])
                inj_transit_end = injected_params['epoch'][planet] + (0.5 * injected_params['duration'][planet])
                overlap = min(inj_transit_end, recovered_params['epoch_end'][transit]) - \
                          max(inj_transit_start, recovered_params['epoch_start'][transit])

                if overlap < (condition_on_overlap * injected_params['duration']):
                    recovered = False

            if condition_on_epoch is not None:
                if abs((recovered_params['epoch'][transit] - injected_params['epoch'][planet]).to_value('jd')) > \
                        condition_on_epoch.to_value('d'):
                    recovered = False

            if condition_on_period is not None:
                period_sway = condition_on_period * injected_params['period'][planet]
                if (recovered_params['period'][transit] < injected_params['period'][planet] - period_sway) or \
                        (recovered_params['period'][transit] > injected_params['period'][planet] + period_sway):
                    recovered=False

            recovered_all_transits.append(recovered)
        recovered_all_planets.append(recovered_all_transits)

    return recovered_all_planets

