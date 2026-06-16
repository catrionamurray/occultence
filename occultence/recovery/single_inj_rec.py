from ..imports import *
import time
import pickle as pkl

def single_injection_recovery(self, lc, planets, i, normalize_each_night, clean_kw, flare_kw, detrend_bin, detrend_kw,
                              bls_kw, bls_bin, recovery_kw, detrend_method, time_this_process, plot_dir, save_plots,
                              predetrend_bls=True, plot=False, verbose=False, plotkw={'ylims': [0.95, 1.05]},
                              save_lcs=False, min_save=True):

    name = lc.name

    # plot injected LC
    if plot:
        ax = self.plot(color='C0', label='raw data')
        lc.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label='injected transit')
        try:
            lc.plot(quantity="transit_model", color='k', ax=ax, linestyle="-")
        except Exception as e:
            print(e)
        if save_plots:
            plt.savefig(f"{plot_dir}/{name}_injected_transit")
        else:
            plt.show()

    # save injected LC
    if save_lcs:
        self.save(fname=f"{plot_dir}/{name}.pkl")

    # normalize
    if normalize_each_night:
        lc = self.normalize_each_night(lc)
        if save_lcs:
            if not min_save:
                # svname = svname + "_norm"
                lc.save(fname=f"{plot_dir}/{lc.name}.pkl")

    # clean
    if time_this_process:
        t0 = time.time()
    clean_targ = lc.clean(**clean_kw)
    if save_lcs:
        # svname = svname + "_clean"
        clean_targ.save(fname=f"{plot_dir}/{clean_targ.name}.pkl")
    if time_this_process:
        t1 = time.time()
        print(f"Time to clean LC: {t1-t0}")

    # flares
    if time_this_process:
        t0 = time.time()
    clean_targ = clean_targ.detect_flares_sclip(i_lc=i, **flare_kw)
    if save_lcs:
        if not min_save:
            # svname = svname + "_flares"
            clean_targ.save(fname=f"{plot_dir}/{clean_targ.name}.pkl")
    if time_this_process:
        t1 = time.time()
        print(f"Time to remove flares: {t1-t0}")

    if detrend_bin is not None:
        # bin before detrending
        bin_targ = clean_targ.bin(dt=detrend_bin)
        bin_targ.telescope = np.array([bin_targ.telescope[0]] * bin_targ.ntime)
        bin_targ.filter = np.array([bin_targ.filter[0]] * bin_targ.ntime)
    else:
        bin_targ = clean_targ

    if time_this_process:
        t0 = time.time()
    if predetrend_bls:
        # do an initial search for transits to mask during detrending
        orig_bin_targ = bin_targ._create_copy()
        removed_nans = bin_targ.remove_nans()
        bls_targs = removed_nans.find_transits(i_lc=i, **bls_kw)
        if len(bls_targs) > 0:
            in_transit = bls_targs[0].metadata['BLS_transits_ind']
            bls_targs[0].masks['transit'] = np.zeros(bls_targs[0].ntime)
            bls_targs[0].masks['transit'][in_transit] = 1.0
            init_transits_masked = bls_targs[0].clean(zero_flux_removal=False, nan_flux_removal=False,
                                                      bad_weather_removal=False, threshold_removal=False,
                                                      dust_removal=False, cosmics_removal=False)
            bin_targ = init_transits_masked
    if time_this_process:
        t1 = time.time()
        print(f"Time to BLS before detrending: {t1-t0}")

    # detrend the light curve
    if time_this_process:
        t0 = time.time()

    if detrend_method == "gp":
        # gp detrend
        gp_targ = bin_targ.gp_detrend(i_lc=i, **detrend_kw)

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

    if save_lcs:
        # svname = svname + "_detrended"
        detrended_targ.save(fname=f"{plot_dir}/{detrended_targ.name}.pkl")

    if time_this_process:
        t1 = time.time()
        print(f"Time to detrend: {t1-t0}")

    if plot:
        ax = bin_targ.plot(color='C0', label='clean lc')
        detrended_targ.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label=f'{detrend_method}-detrended')
        if save_plots:
            plt.savefig(f"{plot_dir}/{detrended_targ.name}")
        else:
            plt.show()

    # search for transit
    if time_this_process:
        t0 = time.time()
    removed_nans = detrended_targ.remove_nans()

    if save_lcs:
        # svname = svname + "_bls"
        removed_nans.save(fname=f"{plot_dir}/{removed_nans.name}.pkl")


    bls_targs = removed_nans.find_transits(i_lc=i, **bls_kw)

    bls_meta = {}

    if save_lcs:
        # svname = svname + "_bls"
        if not min_save:
            for b in bls_targs:
                b.save(fname=f"{plot_dir}/{b.name}.pkl")

    if time_this_process:
        t1 = time.time()
        print(f"Time to BLS search: {t1-t0}")

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

                bls_meta[i] = bls_targ.metadata['BLS_transits_params']
                rec = bls_targ.was_injected_planet_recovered(**recovery_kw)
                bls_targ.metadata['recovery'] = rec

                recovered_n_transits = len(bls_targ.metadata['BLS_transits_params']['depth'])

                if 'meet_condition' in recovery_kw:
                    if recovery_kw['meet_condition'] == "any":
                        min_n_transits = 1
                    if recovery_kw['meet_condition'] == "all":
                        min_n_transits = recovered_n_transits

                if 'min_n_transits' in recovery_kw:
                    min_n_transits = recovery_kw['min_n_transits']

                n_recovered = 0
                # loop over all detected transits (in this case most likely 1):
                for r in range(recovered_n_transits):
                    # !!! the following assumes we have only injected 1 planet at a time !!!:
                    if rec[0][r]:
                        n_recovered += 1
                        # recovered = True

                planets.loc[i, 'n_recovered'] = n_recovered

                if n_recovered >= min_n_transits:
                    # total_recovered += 1
                    planets.loc[i, 'recovered'] = 1.0
                    planets.loc[i, 'log_Prec'] = bls_targ.metadata['BLS_transits_params']['period'][0].to_value('d')
                    planets.loc[i, 'rec_depth'] = bls_targ.metadata['BLS_transits_params']['depth'][0]
                    planets.loc[i, 'rec_duration'] = bls_targ.metadata['BLS_transits_params']['duration'][0].to_value('d')
                    planets.loc[i, 'rec_epoch'] = bls_targ.metadata['BLS_transits_params']['epoch'][0].to_value('d')
                    planets.loc[i, 'bls_power'] = np.max(bls_targ.metadata['BLS_transits_params']['power'])  # [0]
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

    # pkl.dump(bls_meta, open(f"{plot_dir}/{self.name}_BLS.pkl", 'wb'))
    # bls_lcs.append(bls_targ)

    return clean_targ, detrended_targ, bls_targs, planets


def single_clean_detrend(self, lc, clean_kw, detrend_bin, detrend_kw, detrend_method, bls_kw, bls_bin,
                         time_this_process=False, predetrend_bls=True, plot=False, verbose=False, save_plots=False,
                         plot_dir="", i=0, plotkw={'ylims':[0.95, 1.05]}):
    if plot:
        ax = self.plot(color='C0', label='raw data')
        lc.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label='injected transit')
        if save_plots:
            plt.savefig(f"{plot_dir}/{self.name}_injected_transit")
        else:
            plt.show()

    # *** clean *****************
    if time_this_process:
        t0 = time.time()
    clean_targ = lc.clean(**clean_kw)
    if time_this_process:
        t1 = time.time()
        print(f"Time to clean LC: {t1 - t0}")
    # ****************************

    # ***  bin  *****************
    if detrend_bin is not None:
        # bin before detrending
        bin_targ = clean_targ.bin(dt=detrend_bin)
        bin_targ.telescope = np.array([bin_targ.telescope[0]] * bin_targ.ntime)
        bin_targ.filter = np.array([bin_targ.filter[0]] * bin_targ.ntime)
    else:
        bin_targ = clean_targ
    # ****************************

    # *** pre-detrend BLS ********
    if time_this_process:
        t0 = time.time()
    if predetrend_bls:
        # do an initial search for transits to mask during detrending
        orig_bin_targ = bin_targ._create_copy()
        removed_nans = bin_targ.remove_nans()
        bls_targs = removed_nans.find_transits(i_lc=i, **bls_kw)
        if len(bls_targs) > 0:
            in_transit = bls_targs[0].metadata['BLS_transits_ind']
            bls_targs[0].masks['transit'] = np.zeros(bls_targs[0].ntime)
            bls_targs[0].masks['transit'][in_transit] = 1.0
            init_transits_masked = bls_targs[0].clean(zero_flux_removal=False, nan_flux_removal=False,
                                                      bad_weather_removal=False, threshold_removal=False,
                                                      dust_removal=False, cosmics_removal=False)
            bin_targ = init_transits_masked
    if time_this_process:
        t1 = time.time()
        print(f"Time to BLS before detrending: {t1 - t0}")
    # ****************************

    # ***  detrend  **************
    if time_this_process:
        t0 = time.time()

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
            detrended_targ = bin_targ.lsq_ridge_detrend_each_night(orig_lc=orig_bin_targ, **detrend_kw)
        else:
            detrended_targ = bin_targ.lsq_ridge_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)
    elif detrend_method == None:
        print("No detrending method selected!")
        detrended_targ = clean_targ._create_copy().bin(dt=bls_bin)
    else:
        print(f"{detrend_method} is not recognised. Please choose one of: 'gp', 'lsq', 'ridge' or 'mcmc'")
        return None, None, None, None

    if time_this_process:
        t1 = time.time()
        print(f"Time to detrend: {t1 - t0}")

    if plot:
        ax = bin_targ.plot(color='C0', label='clean lc')
        detrended_targ.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label=f'{detrend_method}-detrended')
        if save_plots:
            plt.savefig(f"{plot_dir}/{detrended_targ.name}")
        else:
            plt.show()

    # ****************************

    return clean_targ, detrended_targ

def single_clean_bin(self, lc, clean_kw, detrend_bin, time_this_process=False, plot=False, verbose=False, save_plots=False,
                         plot_dir="", plotkw={'ylims': [0.95, 1.05]}):
    if plot:
        ax = self.plot(color='C0', label='raw data')
        lc.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label='injected transit')
        if save_plots:
            plt.savefig(f"{plot_dir}/{self.name}_injected_transit")
        else:
            plt.show()

    # *** clean *****************
    if time_this_process:
        t0 = time.time()
    clean_targ = lc.clean(**clean_kw)
    if time_this_process:
        t1 = time.time()
        print(f"Time to clean LC: {t1 - t0}")
    # ****************************

    # ***  bin  *****************
    if detrend_bin is not None:
        # bin before detrending
        bin_targ = clean_targ.bin(dt=detrend_bin)
        bin_targ.telescope = np.array([bin_targ.telescope[0]] * bin_targ.ntime)
        bin_targ.filter = np.array([bin_targ.filter[0]] * bin_targ.ntime)
    else:
        bin_targ = clean_targ
    # ****************************

    return clean_targ, bin_targ

def single_predetrend(self, i, bls_kw, time_this_process=False, predetrend_bls=True, verbose=False):
    # *** pre-detrend BLS ********
    if time_this_process:
        t0 = time.time()
    if predetrend_bls:
        # do an initial search for transits to mask during detrending
        bin_targ = self._create_copy()
        removed_nans = bin_targ.remove_nans()
        bls_targs = removed_nans.find_transits(i_lc=i, **bls_kw)
        if len(bls_targs) > 0:
            in_transit = bls_targs[0].metadata['BLS_transits_ind']
            bls_targs[0].masks['transit'] = np.zeros(bls_targs[0].ntime)
            bls_targs[0].masks['transit'][in_transit] = 1.0
            init_transits_masked = bls_targs[0].clean(zero_flux_removal=False, nan_flux_removal=False,
                                                      bad_weather_removal=False, threshold_removal=False,
                                                      dust_removal=False, cosmics_removal=False)
            bin_targ = init_transits_masked
    if time_this_process:
        t1 = time.time()
        print(f"Time to BLS before detrending: {t1 - t0}")
    # ****************************
    return bin_targ
def single_detrend(self, clean_targ, orig_bin_targ, detrend_kw, detrend_method, bls_bin,
                         time_this_process=False, predetrend_bls=True, plot=False, verbose=False, save_plots=False,
                         plot_dir="", i=0, plotkw={'ylims': [0.95, 1.05]}):

    # ***  detrend  **************
    if time_this_process:
        t0 = time.time()

    if detrend_method == "gp":
        # gp detrend
        gp_targ = self.gp_detrend(**detrend_kw)

        # we can also predict the GP for 7.5 min binning and use that for the BLS
        bin_targ = clean_targ.bin(dt=bls_bin)
        gp_model = 1 + \
                   gp_targ.metadata['gp'].predict(y=gp_targ.metadata['data_to_condition_gp'], t=bin_targ.time.value)[0]
        detrended_targ = bin_targ._create_copy()
        detrended_targ.timelike['flux'] = detrended_targ.timelike['flux'] / gp_model

    elif detrend_method == "lsq":
        detrended_targ = self.lsq_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)

    elif detrend_method == "mcmc":
        detrended_targ = self.mcmc_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)

    elif detrend_method == "ridge":
        if predetrend_bls:
            detrended_targ = self.lsq_ridge_detrend_each_night(orig_lc=orig_bin_targ, **detrend_kw)
        else:
            detrended_targ = self.lsq_ridge_detrend_each_night(**detrend_kw)
        detrended_targ = detrended_targ.bin(dt=bls_bin)
    elif detrend_method == None:
        print("No detrending method selected!")
        detrended_targ = clean_targ._create_copy().bin(dt=bls_bin)
    else:
        print(f"{detrend_method} is not recognised. Please choose one of: 'gp', 'lsq', 'ridge' or 'mcmc'")
        return None, None, None, None

    if time_this_process:
        t1 = time.time()
        print(f"Time to detrend: {t1 - t0}")

    if plot:
        ax = self.plot(color='C0', label='clean lc')
        detrended_targ.plot(ax=ax, ylims=plotkw['ylims'], color='C1', label=f'{detrend_method}-detrended')
        if save_plots:
            plt.savefig(f"{plot_dir}/{detrended_targ.name}")
        else:
            plt.show()

    # ****************************

    return detrended_targ

def single_bls(i, detrended_targ, bls_kw, recovery_kw, planets, time_this_process=False,
               verbose=False, save_plots=False, plot_dir="", plotkw={'ylims': [0.95, 1.05]}):
    # search for transit
    if time_this_process:
        t0 = time.time()
    removed_nans = detrended_targ.remove_nans()
    bls_targs = removed_nans.find_transits(i_lc=i, **bls_kw)
    if time_this_process:
        t1 = time.time()
        print(f"Time to BLS search: {t1-t0}")

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
                    planets.loc[i, 'bls_power'] = np.max(bls_targ.metadata['BLS_transits_params']['power'])
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

    return bls_targs, planets

def find_transits_wrapper(lc, dict_args):
    return lc.find_transits(**dict_args)

# def single_predetrend_wrapper(lc, dict_args):
#     return lc.single_predetrend(**dict_args)