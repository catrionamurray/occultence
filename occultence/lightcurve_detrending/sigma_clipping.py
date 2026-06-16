from ..imports import *
from ..utils import running_box

def first_sigma_clip(y, nsigma_upper, nsigma_lower, **kw):
    """
    Perform an overall sigma clip.
    :param y: Data to sigma clip
    :param nsigma_upper: Number of sigma to clip above the data.
    :param nsigma_lower: Number of sigma to clip below the data.
    :return: Sigma-clipped data
    """
    return sigma_clip(y, sigma_upper=nsigma_upper, sigma_lower=nsigma_lower, **kw).filled(np.nan)

def second_sigma_clip(x,y, dy, nsigma_upper, nsigma_lower, running_median_boxsize, use_uncertainties=False,
                      plot=False, operation="median"):
    """
    Perform running sigma clip.
    :param x: time data.
    :param y: data to sigma clip.
    :param nsigma_upper: Number of sigma to clip above the data.
    :param nsigma_lower: Number of sigma to clip below the data.
    :param running_mean_boxsize: Size of box for running mean + std dev.
    :return: Sigma-clipped data
    """
    cm = len(y)
    it = True

    # iteratively sigma clip until data is no longer clipped
    while it == True:
        run_med = running_box(x, y, running_median_boxsize, operation=operation)
        if use_uncertainties:
            avg_std = dy
        else:
            run_std = running_box(x, y, running_median_boxsize, operation='std')
            avg_std = np.nanmedian(run_std)
        cond = np.logical_or(y > run_med + (nsigma_upper * avg_std),
                              y < run_med - (nsigma_lower * avg_std))
        y[cond] = np.nan
        cmasked = np.count_nonzero(~np.isnan(y))
        if cmasked - cm == 0:
            it = False
        else:
            cm = cmasked
    return y

def global_and_local_sigma_clip(self, global_sc_kw={'nsigma_upper':3, 'nsigma_lower':3},
                                local_sc_kw={'nsigma_upper': 5, 'nsigma_lower': 5, 'running_median_boxsize': 0.4}):

    clipped_lc = self._create_copy()

    y_clip_global = first_sigma_clip(self.flux, **global_sc_kw)
    y_clip_local = second_sigma_clip(self.time.value, y_clip_global.copy(), self.uncertainty, **local_sc_kw)

    clipped_lc.timelike['flux'] = y_clip_local
    clipped_lc._set_name(clipped_lc.name + "_sigmaclipped")
    return clipped_lc
