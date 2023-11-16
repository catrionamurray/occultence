from ..imports import *
from .gp import *
from .sigma_clipping import first_sigma_clip, second_sigma_clip
from scipy.optimize import least_squares
def lsq_detrend(self, params, degree, **kw):
    """

    :param self:
    :param params:
    :param degree:
    :param kw:
    :return:
    """
    detrended_lightcurve = self._create_copy()
    def poly(p, x, degree):
        y = 0
        if degree > 0:
            for d in range(degree - 1):
                y += p[d + 1] * x ** (d + 1)
        return y

    def fit_all_params(p, x, y, degree):
        total_model = p[0]
        if len(p) > 1:
            for i, x_i in enumerate(x):
                p_i = p[1:][i * degree:(i + 1) * degree]
                total_model += poly(p=p_i, x=(x_i - np.mean(x_i)) / np.std(x_i), degree=degree)
        return total_model - y

    p0 = [1.0]
    for d in range(degree * len(params)):
        p0.append(0.0)

    params_list = [self.timelike[par] for par in params]

    res = least_squares(fun=fit_all_params, x0=p0, f_scale=np.median(self.uncertainty),
                        args=[params_list, self.flux, degree])

    total_model = res.x[0]
    for i, x_i in enumerate(params_list):
        if degree > 0:
            p_i = res.x[1:][i * degree:(i + 1) * degree]
        else:
            p_i = []
        total_model += poly(p=p_i, x=(x_i - np.mean(x_i)) / np.std(x_i), degree=degree)

    detrended_lightcurve.timelike['lsq_model'] = total_model
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = detrended_lightcurve.timelike['flux'] / total_model
    detrended_lightcurve._set_name(detrended_lightcurve.name + "_lsqdetrend")

    return detrended_lightcurve

def lsq_detrend_each_night(self, params, degree, **kw):
    """
    Perform least-squares opt. of polynomials individually for each night in the light curve timeseries.
    :param self: LightCurve object
    :param params: List of parameters to detrend with
    :param degree: Polynomial degree
    :param kw: keywords to pass to self.lsq_detrend
    :return: detrended LightCurve object
    """
    lsq_days = []
    for i in range(self.ndays):
        lsq_days.append(self.split_day(i).lsq_detrend(params=params, degree=degree, **kw))

    reconst = lsq_days[0]
    for md in lsq_days[1:]:
        reconst = reconst.concatenate(md)

    return reconst



def mcmc_detrend(self, params, degree, draws=1000, tune=1000, chains=4, cores=4, norm=True, plot=True):
    """
    Perform a PyMC NUTS detrend of using timeseries
    :param self: LightCurve object
    :param degree: Number of polynomial degrees
    :param params: List of names of parameters to fit (must be in self.timelike)
    :param draws: Number of samples to draw from the posterior
    :param tune: Number of tuning steps
    :param chains: Number of MCMC chains
    :param cores: Number of computer cores to run in parallel
    :param norm: Boolean whether to normalize the parameters, recommended (default = True)
    :param plot: Boolean whether to plot the traces to check for convergence
    :return: LightCurve object detrended
    """
    import pymc as pm
    import arviz as az

    detrended_lightcurve = self._create_copy()
    p, new_x = {}, {}

    with pm.Model() as mod:
        mu = pm.Normal("A", mu=1.0, sigma=1e-2)
        for par in params:
            p[par] = []
            if norm:
                x = (detrended_lightcurve.timelike[par] - \
                     np.nanmean(detrended_lightcurve.timelike[par])) / np.nanstd(detrended_lightcurve.timelike[par])
                new_x[par] = x
            else:
                x = detrended_lightcurve.timelike[par]
                new_x[par] = x

            for d in range(degree):
                p[par].append(pm.Normal(f"{par}_p{d + 1}", mu=0, sigma=1e-2))
                mu += p[par][d] * (x ** (d + 1))

        # Likelihood (sampling distribution) of observations
        Y_obs = pm.Normal("Y_obs", mu=mu, sigma=detrended_lightcurve.uncertainty, observed=detrended_lightcurve.flux)

    with mod:
        trace = pm.sample(draws=draws, tune=tune, chains=chains, cores=cores)
        summary = az.summary(trace, round_to=20)

    if plot:
        az.plot_trace(trace, combined=True)

    p = {}
    mu = summary['mean']['A']
    for par in params:
        p[par] = summary['mean']['A']
        for d in range(degree):
            mu += summary['mean'][f"{par}_p{d + 1}"] * (new_x[par] ** (d + 1))
            p[par] += summary['mean'][f"{par}_p{d + 1}"] * (new_x[par] ** (d + 1))


    detrended_lightcurve.timelike['mcmc_model'] = mu
    detrended_lightcurve.metadata['mcmc_summary'] = summary
    detrended_lightcurve.metadata['mcmc_trace'] = trace
    detrended_lightcurve.metadata['mcmc_newx'] = new_x
    detrended_lightcurve.metadata['mcmc_individual_parameter_models'] = p
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = detrended_lightcurve.timelike['flux'] / mu
    detrended_lightcurve._set_name(detrended_lightcurve.name + "_mcmcdetrend")

    return detrended_lightcurve


def mcmc_detrend_each_night(self, params, degree, **kw):
    """
    MCMC detrend (using PyMC) each night
    :param self:
    :param params:
    :param degree:
    :param kw:
    :return:
    """
    mcmc_days = []
    for i in range(self.ndays):
        mcmc_days.append(self.split_day(i).mcmc_detrend(params=params, degree=degree, **kw))

    reconst = mcmc_days[0]
    for md in mcmc_days[1:]:
        reconst = reconst.concatenate(md)

    return reconst

def gp_detrend(self, do_first_sigma_clip=True, do_second_sigma_clip=True, nsigma=3, running_mean_boxsize=0.04,
               rotation_period=None, rotation_amp=None, plot=True, figsize=(12, 4), verbose=False, **kw):
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

        3

        if plot:
            plt.plot(x, y, '.', alpha=0.5, label="After First Sigma-Clip")
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

    # store some metadata about the kernel too:
    detrended_lightcurve.metadata['gp'] = gp_func
    detrended_lightcurve.metadata['data_to_condition_gp'] = y[cond_nans] - 1
    detrended_lightcurve.metadata['kernel'] = {}
    for k in gp_kernel:
        detrended_lightcurve.metadata['kernel'][k] = gp_kernel[k]

    detrended_lightcurve._set_name(detrended_lightcurve.name + "_gpdetrend")
    return detrended_lightcurve
