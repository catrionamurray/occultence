from ..imports import *
from .gp import *
from .sigma_clipping import first_sigma_clip, second_sigma_clip

def gp_detrend(self, do_first_sigma_clip=True, do_second_sigma_clip=True, nsigma=3, running_median_boxsize=0.04,
               use_uncertainties=True, rotation_period=None, rotation_amp=None, plot=True, figsize=(12, 4),
               verbose=False, plot_dir="", save_plots=False, plotkw={'ylims': [0.95, 1.05]}, i_lc=0, **kw):

    detrended_lightcurve = self._create_copy()

    if type(nsigma) == int:
        nsigma_lower, nsigma_upper = nsigma, nsigma
    else:
        nsigma_lower, nsigma_upper = nsigma[0], nsigma[1]

    x, y, yerr = self.time, self.flux, self.uncertainty
    if type(x) == astropy.time.core.Time:
        x = x.value

    if plot:
        plt.figure(figsize=figsize)

    if do_first_sigma_clip:
        if plot:
            plt.plot(x, y, '.', alpha=0.5, label='Before First Sigma-Clip')


        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="After First Sigma-Clip")
    if do_second_sigma_clip:
        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="Before Second Sigma-Clip")

        y = second_sigma_clip(x=x, y=y, dy=yerr, nsigma_lower=nsigma_lower, nsigma_upper=nsigma_upper,
                              running_median_boxsize=running_median_boxsize, use_uncertainties=use_uncertainties)
        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="After Second Sigma-Clip")

    if plot:
        plt.legend()

    # remove NaNs - the GP will not work with NaNs!
    cond_nans = ~np.isnan(y)
    # x = x[cond_nans]
    # yerr = yerr[cond_nans]
    # y = y[cond_nans]

    gp_mu, gp_var, gp_jitter, gp_mu_og, gp_var_og, gp_kernel, gp_func = gp(x[cond_nans], y[cond_nans] - 1,
                                                                           yerr[cond_nans], x,
                                                                           rotation_period=rotation_period,
                                                                           rotation_amp=rotation_amp,
                                                                           plot=plot, figsize=figsize, verbose=verbose,
                                                                           **kw)

    detrended_lightcurve.timelike['gp_model'] = (gp_mu_og + 1)
    detrended_lightcurve.timelike['gp_model_err'] = np.sqrt(gp_var_og)
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = (detrended_lightcurve.timelike['flux'] / (gp_mu_og + 1)) * 1

    if plot:
        plt.figure(figsize=figsize)
        plt.errorbar(self.time.value, detrended_lightcurve.timelike['flux'], self.uncertainty, fmt=".k", capsize=0,
                     zorder=0, alpha=0.3)
        plt.title("GP-detrended data")
        plt.ylabel("Flux")
        plt.xlabel("Time [d]")
        if save_plots:
            plt.savefig(f"{plot_dir}/{self.name}_gp_detrended")
        else:
            plt.show()
        plt.close()

        if self.ndays > 1:
            ax = self.plot()
            x_pred = np.linspace(np.min(self.time.value), np.max(self.time.value), 1000)
            pred_mu, pred_var = gp_func.predict(y[cond_nans] - 1, x_pred, return_var=True)
            y_pred = pred_mu + 1
            yerr_pred = np.sqrt(pred_var)

            for i, ax_i in enumerate(ax):
                plt.sca(ax_i)
                plt.fill_between(x_pred, y_pred - yerr_pred, y_pred + yerr_pred,
                                 color="orange", alpha=0.5, zorder=2)
                plt.plot(x_pred, y_pred, "orange", lw=1.5, alpha=0.8, zorder=2)
                plt.xlim(self.split_day(i).time[0].value, self.split_day(i).time[-1].value)
            if save_plots:
                plt.savefig(f"{plot_dir}/{self.name}_gp_detrended2")
            else:
                plt.show()

    # store some metadata about the kernel too:
    detrended_lightcurve.metadata['gp'] = gp_func
    detrended_lightcurve.metadata['data_to_condition_gp'] = y[cond_nans] - 1
    detrended_lightcurve.metadata['kernel'] = {}
    for k in gp_kernel:
        detrended_lightcurve.metadata['kernel'][k] = gp_kernel[k]

    detrended_lightcurve._set_name(detrended_lightcurve.name + "_gpdetrend")
    return detrended_lightcurve


def gp_detrend_each_night(self, **kw):
    """
    Perform GP fit individually for each night in the light curve timeseries.
    :param self: LightCurve object
    :param kw: keywords to pass to self.gp_detrend
    :return: detrended LightCurve object
    """
    lsq_days = []
    for i in range(self.ndays):
        lsq_days.append(self.split_day(i).gp_detrend(**kw))

    reconst = lsq_days[0]
    for md in lsq_days[1:]:
        reconst = reconst.concatenate(md)

    return reconst





