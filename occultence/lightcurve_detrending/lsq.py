from ..imports import *
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

    params_list = []
    for par in params:
        if par == 'time':
            params_list.append(self.time.jd)
        else:
            params_list.append(self.timelike[par])

    res = least_squares(fun=fit_all_params, x0=p0, f_scale=np.median(self.uncertainty),
                        args=[params_list, self.flux, degree], **kw)

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
