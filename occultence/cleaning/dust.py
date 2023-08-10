from ..imports import *

def mask_dust(self, starts_of_dust_periods, ends_of_dust_periods):
    """
    Mask entire periods during which dust affected observations [WARNING: this will remove ALL data in these periods,
    not just individual nights]

    :param self:
    :param starts_of_dust_periods: List of astropy.Time objects marking the starts of each dust period
    :param ends_of_dust_periods: List of astropy.Time objects marking the ends of each dust period
    :return:
    """

    self.masks['dust'] = np.zeros(self.ntime)
    for s, e in zip(starts_of_dust_periods, ends_of_dust_periods):
        self.masks['dust'][(self.time > s) & (self.time < e)] = 1
