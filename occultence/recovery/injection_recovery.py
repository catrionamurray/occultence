from ..imports import *
import time
from multiprocessing import Pool
import warnings
from astropy.utils.exceptions import AstropyWarning
import pickle as pkl

warnings.filterwarnings("ignore", category=AstropyWarning,
                       message=".*Input data contains invalid values.*")


#changed nr of tested durations to 4(from 10): 0.01, 0.04, 0.07, 0.1

def full_injection_recovery(self,
                            nfake=10,
                            pool_bls=False,
                            poolkw={'chunksize': 1},
                            ncores=4,
                            minimum_planet_radius=0.5 * u.R_earth,
                            maximum_planet_radius=3 * u.R_earth,
                            minimum_period=0.5 * u.d,
                            maximum_period=10 * u.d,
                            ld=[0.385, 0.304],
                            normalize_each_night=False,
                            clean_kw={'dust_removal': False, 'bad_weather_removal': False, 'cosmics_removal': True,
                                      'cosmic_boxsize': 0.08, 'cosmic_nsigma': 3},
                            flare_kw={},
                            detrend_method="gp",
                            detrend_bin=20 * u.minute,
                            detrend_kw={'do_first_sigma_clip': True, 'do_second_sigma_clip': True,
                                        'running_mean_boxsize': 0.08, 'nsigma': 3, 'plot': False},
                            bls_kw={"minimum_period": 0.5, "maximum_period": 10,
                                    'transit_durations': np.linspace(0.01, 0.1, 4), 'plot': False, 'verbose': False},
                            bls_bin=7.5 * u.minute,
                            recovery_kw={'condition_on_epoch': 1 * u.hour, 'min_n_transits': 1,
                                         'meet_condition': 'any'},
                            plot=False,
                            plot_dir="/plots/",
                            save_plots=False,
                            plot_kw={'ylims': [0.95, 1.05]},
                            verbose=False,
                            time_this_process=False,
                            save_lcs= False,
                            svname="injected_planets.csv",
                            min_save=True,
                            planets=None,
                            lcs_with_transits=None,
                            ):


    if save_lcs:
        self.save(fname=f"{plot_dir}/{self.name}.pkl")

    if os.path.exists(f"../../Results/{self.name}_BLS.pkl"):
        print("Deleting BLS file...")
        os.remove(f"../../Results/{self.name}_BLS.pkl")

    if planets is None or lcs_with_transits is None:
        lcs_with_transits, planets = self.inject_lots_of_transits(nfake=nfake,
                                                                  minimum_planet_radius=minimum_planet_radius,
                                                                  maximum_planet_radius=maximum_planet_radius,
                                                                  minimum_period=minimum_period,
                                                                  maximum_period=maximum_period,
                                                                  ld=ld,
                                                                  fname=svname)

    bls_kw['verbose'] = verbose
    clean_lcs, bin_lcs, detrend_lcs, bls_lcs = [], [], [], []

    if time_this_process:
        t0 = time.time()

    planets = pd.read_csv(svname)

    bls_kw['save_plots'] = save_plots
    bls_kw['plot_dir'] = plot_dir
    bls_kw['min_save'] = min_save
    detrend_kw['save_plots'] = save_plots
    detrend_kw['plot_dir'] = plot_dir
    detrend_kw['min_save'] = min_save
    flare_kw['save_plots'] = save_plots
    flare_kw['plot_dir'] = plot_dir
    flare_kw['min_save'] = min_save
    # plot_kw['min_save'] = min_save

    if pool_bls:
        import occultence.recovery.single_inj_rec as sir
        print(f"Pooling {nfake} injection-recoveries with {ncores} cores and kw={poolkw}.")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # results = pool.starmap(sir.single_injection_recovery,
            #                        [(self, lc, planets, i, clean_kw, detrend_bin, detrend_kw, bls_kw, bls_bin,
            #                          recovery_kw, detrend_method, time_this_process) for i, lc in enumerate(lcs_with_transits)],
            #                        **poolkw)

            # run = []
            for i, lc in enumerate(lcs_with_transits):
                observed = int(lc.was_planet_observed())
                planets.loc[i, 'injected'] = 1.0
                planets.loc[i, 'observed'] = observed

                if normalize_each_night:
                    lc = lc.normalize_each_night()

                clean_targ, bin_targ = self.single_clean_bin(lc, clean_kw, detrend_bin)
                removed_nans = bin_targ.remove_nans()
                bin_lcs.append(removed_nans)
                clean_lcs.append(clean_targ)

            # pool the pre-detrend BLS
            pool = Pool(ncores)
            bin_lcs_pre = pool.starmap(sir.single_predetrend, [(lc, i, bls_kw,) for lc in bin_lcs])

            # for now I cannot pool the GP-detrending - I believe it's an issue with george or NaNs
            to_run = []
            for i, (bin_lc, clean_lc, orig_bin_lc) in enumerate(zip(bin_lcs_pre, clean_lcs, bin_lcs)):
                if planets.loc[i, 'observed'] == 1:
                    detrend_lc = bin_lc.single_detrend(clean_targ=clean_lc, orig_bin_targ=orig_bin_lc,
                                                       detrend_kw=detrend_kw,
                                                       detrend_method=detrend_method, bls_bin=bls_bin,
                                                       plot_dir=plot_dir,
                                                       save_plots=save_plots,
                                                       i=i,
                                                       plotkw=plot_kw,
                                                       )
                    detrend_lcs.append(detrend_lc)
                    to_run.append(i)
                else:
                    print("Planet was not observed, skipping...")

            # pool the post-detrend BLS
            pool = Pool(ncores)
            results = pool.starmap(sir.single_bls, [(i, lc, bls_kw, recovery_kw, planets,
                                                     time_this_process, verbose, plot_dir,
                                                     save_plots, plot_kw) for i, lc in zip(to_run, detrend_lcs)])

            bls_lcs = [r[0] for r in results]
            planets_list = [r[1] for r in results]

        for i, p in enumerate(planets_list):
            planets.loc[i] = p.iloc[i]

        planets.to_csv(svname, index=False)

    else:
        bls_meta = {}
        for i, lc in enumerate(lcs_with_transits):
            # try:
            print(f"{i + 1}/{len(lcs_with_transits)}...")
            observed = int(lc.was_planet_observed())
            planets.loc[i, 'injected'] = 1.0
            planets.loc[i, 'observed'] = observed

            if observed:
                clean_targ, detrend_targ, bls_targ, planets = self.single_injection_recovery(lc=lc,
                                                                                        planets=planets,
                                                                                        i=i,
                                                                                        normalize_each_night=normalize_each_night,
                                                                                        flare_kw=flare_kw,
                                                                                        clean_kw=clean_kw,
                                                                                        detrend_method=detrend_method,
                                                                                        detrend_bin=detrend_bin,
                                                                                        detrend_kw=detrend_kw,
                                                                                        bls_kw=bls_kw,
                                                                                        bls_bin=bls_bin,
                                                                                        recovery_kw=recovery_kw,
                                                                                        plot=plot,
                                                                                        plot_dir=plot_dir,
                                                                                        save_plots=save_plots,
                                                                                        time_this_process=time_this_process,
                                                                                        verbose=verbose,
                                                                                        plotkw=plot_kw,
                                                                                        save_lcs=save_lcs,
                                                                                        min_save=min_save,)

                planets.to_csv(svname, index=False)
                clean_lcs.append(clean_targ)
                detrend_lcs.append(detrend_targ)
                bls_lcs.append(bls_targ)
                bls_meta[i] = bls_targ[0].metadata['BLS_transits_params']
                bls_meta[i]['frac_inj_dur'] = bls_targ[0].metadata['BLS_transits_params']['dt'] / bls_targ[0].metadata['injected_planet']['duration'][0]
                bls_meta[i]['injected'] = bls_targ[0].metadata['injected_planet']

            else:
                print("Planet was not observed, skipping...")
                planets.to_csv(svname, index=False)
                clean_lcs.append(None)
                detrend_lcs.append(None)
                bls_lcs.append(None)

            # except Exception as e:
            #     print(e)

    if time_this_process:
        t1 = time.time()
        print(f"Time to inject-recover {nfake} planets: {t1 - t0}")

    pkl.dump(bls_meta, open(f"{plot_dir}/{self.name}_BLS.pkl", 'wb'))

    # print summary
    print(f"Planets recovered: {100 * len(planets.loc[planets['recovered'] == 1.0]) / len(planets['recovered'])}%")
    print("But did all of those planets transited during the observation window?...")
    if len(planets['recovered'][planets['observed'] == 1.0]) > 0:
        print(f"""Observed Planets recovered: {(100 * len(planets.loc[(planets['recovered'] == 1.0) &
                                                                      (planets['observed'] == 1.0)]) / len(planets['recovered'][planets['observed'] == 1.0])):.1f}%""")
    else:
        print("None of the planets injected were observed!")

    return lcs_with_transits, clean_lcs, detrend_lcs, bls_lcs, planets


def was_injected_planet_recovered(self, min_n_transits=1, condition_on_depth=None, condition_on_overlap=None,
                                  condition_on_epoch=None, condition_on_period=None, condition_on_snr=None,
                                  condition_on_blspower=None, verbose=False, **kw):
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
    :param condition_on_blspower: BLS power threshold for a transit to be recovered. If objective was 'snr' then this is
    equivalent to the phase-folded SNR.
    :return:
    """
    recovered_all_planets = []
    recovered_all_transits = []

    injected_params = self.metadata['injected_planet']
    recovered_params = self.metadata['BLS_transits_params']

    if verbose:
        print(f"Injected planet parameters: {injected_params}")
        print(f"Recovered planet parameters: {recovered_params}")

    # loop over all planet transits injected
    for planet in range(len(injected_params['depth'])):

        transit_t = injected_params['epoch'][planet].to_value('d')
        all_epochs = []
        while transit_t < self.time.value[-1]:
            all_epochs.append(transit_t)
            transit_t += injected_params['period'][0].to_value('d')

        n_transits = len(recovered_params['depth'])
        if n_transits < min_n_transits:
            recovered = [False]*n_transits
            if verbose:
                print(f"We detect {n_transits} but the minimum required # of transits = {min_n_transits}")

        # loop over all recovered transits
        for transit in range(n_transits):
            recovered = True

            if len(all_epochs) == 0:
                # if the planet occurred after the end of the observation window:
                recovered = False
                recovered_all_transits.append(recovered)
                if verbose:
                    print(f"The planet transited after the end of the observation window")
                continue

            if condition_on_depth is not None:
                if recovered_params['depth'][transit] < (condition_on_depth * injected_params['depth'][planet]):
                    recovered = False
                    if verbose:
                        print(f"""
                        The recovered planet's depth ({recovered_params['depth'][transit]}) < {(condition_on_depth * injected_params['depth'][planet])}
                        """)

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

                    if verbose:
                        print(f"""
                        The recovered planet's overlaps with the true transit < {(condition_on_overlap * injected_params['duration'])}
                        """)

            if condition_on_epoch is not None:
                # print(all_epochs, recovered_params['epoch'][transit], transit_t, self.time.value[-1])
                # print(self.metadata)
                closest_transit = find_nearest(all_epochs, recovered_params['epoch'][transit].to_value('d'))
                if abs(recovered_params['epoch'][transit] - all_epochs[closest_transit]*u.d) > \
                        condition_on_epoch:
                    recovered = False

                    if verbose:
                        print(f"""
                        The recovered planet's epoch > {condition_on_epoch} from the true transit.
                        """)

            if condition_on_period is not None:
                period_sway = condition_on_period * injected_params['period'][planet]
                if (recovered_params['period'][transit] < injected_params['period'][planet] - period_sway) or \
                        (recovered_params['period'][transit] > injected_params['period'][planet] + period_sway):
                    recovered = False

                    if verbose:
                        print(f"""
                        The recovered planet's period > {period_sway} from the true transit.
                        """)

            if condition_on_snr is not None:
                if recovered_params['snr'][transit] < condition_on_snr:
                    recovered = False

                    if verbose:
                        print(f"""
                        The recovered planet's SNR < {condition_on_snr}.
                        """)

            if condition_on_blspower is not None:
                if recovered_params['power'][transit] < condition_on_blspower:
                    recovered = False

                    if verbose:
                        print(f"""
                        The recovered planet's BLS power < {condition_on_blspower}.
                        """)


            recovered_all_transits.append(recovered)
        recovered_all_planets.append(recovered_all_transits)

    return recovered_all_planets

def was_planet_observed(self, fraction_overlap=0.5, n_transits=1, planet_i=0):

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
                    observed += 1
            transit_start = transit_start + period
            transit_end = transit_end + period
    return bool(observed >= n_transits)

# def split_lightcurve(self, split_every=0.5*u.d):
#     t = self.time.value * u.d
#     start = t[0]
#     nextdays = t[np.absolute(t - start) > split_every]
#     split = []
#
#     while nextdays != []:
#         start = nextdays[0]
#         ind_st = np.where(t == start)[0][0]
#         split.append(ind_st)
#         time = t[ind_st:]
#         nextdays = time[np.absolute(time - start) > split_every]
#
#     times = np.split(t, split)
#
#     return times, split

def split_lightcurve(self, split_every=0.5 * u.d):
    """Vectorized version - 100-1000x faster than original."""

    # Get time array with units
    t = self.time.value * u.d

    # Handle unit conversion
    if hasattr(split_every, 'to'):
        split_val = split_every.to(u.d).value
    else:
        split_val = split_every

    # Vectorized approach: compute all time differences at once
    time_vals = t.value
    time_diffs = np.diff(time_vals)

    # Find where consecutive points are more than split_every apart
    gap_indices = np.where(time_diffs > split_val)[0]

    # Split indices are right after each gap
    split = (gap_indices + 1).tolist()

    # Split the time array
    times = np.split(t, split)

    return times, split

