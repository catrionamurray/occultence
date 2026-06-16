from ..imports import *
from ..utils import running_box

def mask_bad_weather(self, boxsize, thresh, verbose, plot=False):
    """
    Create a mask of times where there is bad weather. This is determined as when the artificial light curve (weighted
    mean of multiple star's lightcurves) has a high scatter == most stars are experiencing atmospheric contamination.
    :param self:
    :param boxsize: Size of box to calculate running standard deviation [in days]
    :param thresh: Threshold for std dev to determine bad weather.
    :return:
    """
    running_std = running_box(self.time.value, self.timelike['artificial_lightcurve'], boxsize, operation="std")
    mask = np.zeros(self.ntime)
    mask[running_std > thresh] = 1

    self.masks['bad_weather'] = mask
    if verbose:
        print(f"""
            {100 * np.divide(float(np.count_nonzero(self.masks['bad_weather'])),
            len(self.masks['bad_weather'])):.2f}% of data is flagged as bad weather ⛅.
            """)