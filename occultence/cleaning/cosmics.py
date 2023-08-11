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
    run_med = running_box(self.time, self.flux, boxsize, 'median')
    run_std = running_box(self.time, self.flux, boxsize, 'std')
    # this mask will be the wrong size???
    self.masks['cosmics'] = (self.flux > run_med + (nsigma * run_std))
