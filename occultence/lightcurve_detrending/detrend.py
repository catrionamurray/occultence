from ..imports import *
from .gp import *
from .sigma_clipping import first_sigma_clip, second_sigma_clip

def gp_detrend(self, do_first_sigma_clip=True, do_second_sigma_clip=True, nsigma=3, running_mean_boxsize=0.04,
                        rotation_period=None, rotation_amp=None, plot=True, figsize=(12, 4), verbose=False, **kw):

    detrended_lightcurve = self._create_copy()

    if type(nsigma)==int:
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

        y = first_sigma_clip(y=y,nsigma_lower=nsigma_lower, nsigma_upper=nsigma_upper)

        if plot:
            plt.plot(x,y, '.', alpha=0.5, label="After First Sigma-Clip")
    if do_second_sigma_clip:
        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="Before Second Sigma-Clip")

        y = second_sigma_clip(x=x, y=y, nsigma_lower=nsigma_lower, nsigma_upper=nsigma_upper,
                               running_mean_boxsize=running_mean_boxsize)
        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="After Second Sigma-Clip")

    if plot:
        plt.legend()

    # remove NaNs - the GP will not work with NaNs!
    cond_nans = ~np.isnan(y)
    # x = x[cond_nans]
    # yerr = yerr[cond_nans]
    # y = y[cond_nans]

    gp_mu, gp_var, gp_jitter, gp_mu_og, gp_var_og, gp_kernel, gp_func = gp(x[cond_nans], y[cond_nans]-1, yerr[cond_nans], x,
                                                                  rotation_period=rotation_period,
                                                                  rotation_amp=rotation_amp,
                                                                  plot=plot,figsize=figsize, verbose=verbose, **kw)

    detrended_lightcurve.timelike['gp_model'] = (gp_mu_og+1)
    detrended_lightcurve.timelike['gp_model_err'] = np.sqrt(gp_var_og)
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = (detrended_lightcurve.timelike['flux'] / (gp_mu_og+1)) * 1

    if plot:
        plt.figure(figsize=figsize)
        plt.errorbar(self.time.value, detrended_lightcurve.timelike['flux'], self.uncertainty, fmt=".k", capsize=0,zorder=0, alpha=0.3)
        plt.title("GP-detrended data")
        plt.ylabel("Flux")
        plt.xlabel("Time [d]")

    # store some metadata about the kernel too:
    detrended_lightcurve.metadata['gp'] = gp_func
    detrended_lightcurve.metadata['data_to_condition_gp'] = y[cond_nans]-1
    detrended_lightcurve.metadata['kernel'] = {}
    for k in gp_kernel:
        detrended_lightcurve.metadata['kernel'][k] = gp_kernel[k]

    detrended_lightcurve._set_name(detrended_lightcurve.name + "_gpdetrend")
    return detrended_lightcurve








