from ..imports import *
from .gp import *

def gp_detrend_lightcurve(self, first_sigma_clip=True, second_sigma_clip=True, nsigma=3, running_mean_boxsize=0.04,
                        rotation_period=None, rotation_amp=None, binkw={}, plot=True,figsize=(6, 4)):

    detrended_lightcurve = self._create_copy()

    if type(nsigma)==int:
        nsigma_lower, nsigma_upper = nsigma, nsigma
    else:
        nsigma_lower, nsigma_upper = nsigma[0], nsigma[1]

    x, y, yerr = self.time, self.flux, self.uncertainty
    if first_sigma_clip:
        y = self.first_sigma_clip(y=y,nsigma_lower=nsigma_lower, nsigma_upper=nsigma_upper)
    if second_sigma_clip:
        y = self.second_sigma_clip(x=x, y=y, nsigma_lower=nsigma_lower, nsigma_upper=nsigma_upper,
                               running_mean_boxsize=running_mean_boxsize)

    # removing the NaNs will mess up the shape of the timelike arrays - but the GP may fail with NaNs - TO TEST!
    # cond_nans = ~np.isnan(y)
    # x = x[cond_nans]
    # yerr = yerr[cond_nans]
    # y = y[cond_nans]

    gp_mu, gp_var, gp_jitter, gp_kernel = gp(x, y, yerr, rotation_period=rotation_period,rotation_amp=rotation_amp,
                                             plot=plot,figsize=figsize)

    detrended_lightcurve.timelike['gp_model'] = gp_mu
    detrended_lightcurve.timelike['gp_model_err'] = np.sqrt(gp_var)
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = (detrended_lightcurve.timelike['flux'] / gp_mu) * 1

    # store some metadata about the kernel too:
    detrended_lightcurve.metadata['kernel'] = {}
    for k in gp_kernel:
        detrended_lightcurve.metadata['kernel'][k] = gp_kernel[k]

    return detrended_lightcurve








