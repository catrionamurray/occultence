from ..imports import *
from ..utils import running_box

def first_sigma_clip(y, nsigma_upper, nsigma_lower):
    """
    Perform an overall sigma clip.
    :param y: Data to sigma clip
    :param nsigma_upper: Number of sigma to clip above the data.
    :param nsigma_lower: Number of sigma to clip below the data.
    :return: Sigma-clipped data
    """
    return sigma_clip(y, sigma_upper=nsigma_upper, sigma_lower=nsigma_lower).filled(np.nan)

def second_sigma_clip(x,y, nsigma_upper, nsigma_lower, running_mean_boxsize):
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
        run_med = running_box(x, y, running_mean_boxsize, 'median')
        run_std = running_box(x, y, running_mean_boxsize, 'std')
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