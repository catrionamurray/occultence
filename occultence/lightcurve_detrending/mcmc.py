from ..imports import *
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