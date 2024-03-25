# thanks to Soichiro Hattori for this piece of code!
from ..imports import *
from sklearn.model_selection import ParameterGrid
from tqdm import tqdm

def lsq_ridge_detrend_each_night(self, param_list, highest_order=3, orig_lc=None, orders=None, plot=False, verbose=False):

    # X_f = generate_design_matrix(second_night, 4, 4, 4)
    # weights = calc_ridge_coeff(X_f, second_night.flux, 1e-5)
    # model = X_f @ weights

    ridge_days = []
    for i in range(self.ndays):
        if orig_lc is None:
            orig = None
        else:
            orig = orig_lc.split_day(i)

        if orders is None:
            order = orders
        else:
            order = orders[i]
        ridge_days.append(self.split_day(i).lsq_ridge_detrend(param_list=param_list, highest_order=highest_order,
                                                              orig_lc=orig, orders=order, plot=plot, verbose=verbose))

    reconst = ridge_days[0]

    for rd in ridge_days[1:]:
        reconst = reconst.concatenate(rd)

    return reconst

def lsq_ridge_detrend(self, param_list, highest_order=3, orig_lc=None, orders=None, plot=False, verbose=False):

    if orig_lc is None:
        detrended_lightcurve = self._create_copy()
    else:
        detrended_lightcurve = orig_lc._create_copy()

    scores = []
    models = []
    param_grid = {}

    if orders is None:
        for i in param_list:
            param_grid[i] = np.arange(1, highest_order + 2)
        param_grid['precision'] = 10.0 ** np.arange(-5, -2)
        param_grid = ParameterGrid(param_grid)
    else:
        if type(orders) == dict:
            param_grid = [orders]
        else:
            param_grid = orders

    # print(f"Calculating all polynomial combinations for: {param_list}...")
    for params in param_grid:
        X = generate_design_matrix(self=self, N_params=params)
        weights = calc_ridge_coeff(design_matrix=X, flux=self.flux, precision=params['precision'])
        model = X @ weights
        models.append(model)
        # Calculate the score here and append to scores array
        nd = 1
        for k, v in params.items():
            nd += v - 1
        aic = calculate_AIC(flux = self.flux, ypred=model, sigma=self.uncertainty, k=nd)
        scores.append(aic)

    models = np.array(models)
    scores = np.array(scores)

    i_min = np.argmin(scores) # find best/lowest AIC

    if plot:
        fig, ax = plt.subplots(dpi=300, figsize=(10, 10))
        ax.plot(self.time.jd, self.flux, '.k')
        for i in range(0, len(param_grid), 1):
            if len(param_grid) > 30:
                ax.plot(self.time.jd, models[i])#, c=[scores[i]]*len(self.time.jd))
            else:
                ax.plot(self.time.jd, models[i], label=f'{param_grid[i]}, {scores[i]}')
        ax.plot(self.time.jd, models[i_min], color='k', label=f"best AICc, {param_grid[i_min]}")
        ax.legend(fontsize=6)

    if orig_lc is None:
        mu = models[i_min]
    else:
        # if we've already removed transits before detending we will linearly interpolate our model here and detrend the
        # original LC
        mu = np.interp(detrended_lightcurve.time.jd, self.time.jd, models[i_min])

    detrended_lightcurve.timelike['ridge_model'] = mu
    detrended_lightcurve.metadata['ridge_parameter_degrees'] = param_grid[i_min]
    detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
    detrended_lightcurve.timelike['flux'] = detrended_lightcurve.timelike['flux'] / mu
    detrended_lightcurve._set_name(detrended_lightcurve.name + "_ridgedetrend")

    return detrended_lightcurve

def calculate_AIC(flux, ypred, sigma, k):
    log_lik = ln_likelihood(flux, ypred, sigma)
    correction = (2 * k * (k + 1)) / (len(flux) - k - 1) # to change from AIC to AICc (this is mostly negligible for n>>k)
    aic = (2*k) - 2*log_lik + correction
    return aic

def likelihood(flux, ypred, sigma):
    return (1/(np.sqrt(2*np.pi)*sigma)) * np.exp(-0.5 * ((flux-ypred)/sigma)**2)

def ln_likelihood(flux,ypred,sigma):
    return np.sum(np.log(likelihood(flux, ypred, sigma)))

def generate_powers(param, N):
    x = np.vander(param, N, increasing=True)
    return x[:, 1::]  # remove constant offset

def generate_design_matrix(self, N_params):
    """

    :param df: LightCurve object
    :param N_params: Dict with number of polynomial degrees for each timelike parameter
    :return:
    """
    param_list = [np.repeat(1, self.time.size).reshape(-1, 1)]
    for param, N in N_params.items():
        if param != "precision":
            res = generate_powers(param=self.timelike[param], N=N)
            if res.size != 0:
                param_list.append(res)
    X = np.hstack(param_list)
    return X

def calc_ridge_coeff(design_matrix, flux, precision):
    X = design_matrix
    XTX = X.T @ X
    R = precision * np.identity(XTX.shape[0])
    A = XTX + R
    B = X.T @ flux
    weights = np.linalg.solve(A, B)  # solves for x in expression Ax=B
    return weights