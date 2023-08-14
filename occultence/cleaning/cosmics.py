from ..imports import *
from ..utils import running_box

def mask_cosmics(self, boxsize, nsigma):
    """
    Mask the times when there were cosmic ray hits (determined using a running sigma clip)
    :param self:
    :param boxsize:
    :param nsigma:
    :return:
    """
    if type(self.time[0]) == astropy.time.core.Time:
        time = self.time.value
    else:
        time = self.time
    run_med = running_box(time, self.flux, boxsize, 'median')
    run_std = running_box(time, self.flux, boxsize, 'clipped_std')

    self.masks['cosmics'] = np.array((self.flux > run_med + (nsigma * run_std)), dtype=int)
