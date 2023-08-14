from ..imports import *
from .timelike_masks import mask_timelike_threshold
from .bad_weather import mask_bad_weather
from .cosmics import mask_cosmics
from .dust import mask_dust
from ..lightcurve import *

operator_dict = { ">": operator.gt,
                  "<": operator.lt,
                  ">=": operator.ge,
                  "<=": operator.le,
                  "==": operator.eq,
                  "!=": operator.ne}

def clean(self,
          zero_flux_removal = True,
          nan_flux_removal = True,
          bad_weather_removal = True,
          threshold_removal=True,
          dust_removal=True,
          cosmics_removal=True,
          bad_weather_boxsize: float = 0.02,
          bad_weather_threshvalue: float = 0.06,
          thresholds: dict = {},
          threshold_operators: dict = {},
          dust_crossing_events: dict = {'Ganymede': Time([(2458605, 2458654)], format='jd', scale='tdb'),
                                        'Callisto': Time([(2458384, 2458385)], format='jd', scale='tdb')},
          cosmic_boxsize: float = 0.02,
          cosmic_nsigma: int = 4,
          ):
    """

    :param self:
    :param zero_flux_removal: Boolean whether to remove/mask data points where flux == 0.0. (Default: True)
    :param nan_flux_removal: Boolean whether to remove/mask data points where flux is NaN. (Default: True)
    :param bad_weather_removal: Boolean whether to remove/mask data points affected by bad weather. (Default: True)
    :param threshold_removal: Boolean whether to remove/mask data points above/below/equal to a threshold. (Default: True)
    :param dust_removal: Boolean whether to remove/mask data points affected by dust moving. (Default: True)
    :param cosmics_removal: Boolean whether to remove/mask data points affected by cosmic rays. (Default: True)
    :param bad_weather_boxsize: The boxsize (in days) of the bad weather running standard deviation. (Default: 0.02 [d])
    :param bad_weather_threshvalue: The threshold value for flagging bad weather (Default: 0.06 [i.e. 6%]).
    :param thresholds: A dictionary containing the threshold value to apply to the timelike array. (Default: {})
    :param threshold_operators: A dictionary containing the operation to apply to the timelike array and its threshold.
     The options are >, <, >=, <=, == and !=. (Default: {})
    :param dust_crossing_events: A dictionary
    :return:
    """

    self.masks = {}
    self.metadata['thresholds'] = {}

    # mask where flux is nan or zero
    if nan_flux_removal:
        print("Masking NaN fluxes...")
        self.masks['flux_is_nan'] = np.zeros(self.ntime)
        self.masks['flux_is_nan'][np.isnan(self.flux)] = 1
    if zero_flux_removal:
        print("Masking zero fluxes...")
        self.masks['flux_is_zero'] = np.zeros(self.ntime)
        self.masks['flux_is_zero'][self.flux==0] = 1

    ### Bad Weather ###
    if bad_weather_removal:
        if "bad_weather" in self.timelike:
            print("Masking bad weather using 'bad_weather' in timelike dictionary...")
            self.masks['bad_weather'] = self.timelike['bad_weather']
        else:
            print("Masking bad weather...")
            if "artificial_lightcurve" in self.timelike:
                self.mask_bad_weather(bad_weather_boxsize, bad_weather_threshvalue)
            else:
                # If there is no artifical lightcurve stored we cannot check for bad weather!
                warnings.warn(f""" The LightCurve's timelike dictionary does not appear to have the 'artifical_lightcurve' in it.
                 This array is necessary to calculate bad weather, therefore this will ** not ** be applied!""")

    ### Thresholds ###
    # Apply a threshold to each key passed in thresholds dictionary. Threshold_operators tells the function which
    # operation to use e.g. >, ==, <= etc.
    if threshold_removal:
        for thresh in thresholds:
            if thresh in self.timelike:
                print(f"Masking {thresh}{threshold_operators[thresh]}{thresholds[thresh]}...")
                if threshold_operators[thresh] in operator_dict.keys():
                    self.mask_timelike_threshold(timelike_key=thresh,
                                                 threshold=thresholds[thresh],
                                                 op=threshold_operators[thresh])
                else:
                    # If the thresholds_operators dictionary contains a key which isn't valid (in operator_dict)
                    message = f""" The threshold_operators dictionary passed contains the value
                            {threshold_operators[thresh]}. The options are: {", ".join(operator_dict.keys())}. 
                            """
                    cheerfully_suggest(message)

                self.metadata['thresholds'][thresh] = thresholds[thresh]
            else:
                # If the thresholds dict contains a key which isn't in the LightCurve.timelike dictionary.
                message = f"""
                            The thresholds dictionary contains the key {thresh} that doesn't exist in the timelike 
                            dictionary. The options for thresholds to set are: {", ".join(self.timelike.keys())} 
                                """
                cheerfully_suggest(message)

    ### Dust Crossing Issues ###
    if dust_removal:
        print(f"Masking dust-crossing events...")
        starts = [d[0] for d in dust_crossing_events[self.telescope]]
        ends = [d[1] for d in dust_crossing_events[self.telescope]]
        self.mask_dust(starts_of_dust_periods=starts,
                       ends_of_dust_periods=ends)

    # clean before cosmic ray removal because it will reduce the running std dev. (making cosmics easier to spot)
    # CAN'T DO THIS - IT MESSES UP THE SIZES OF THE ARRAYS...
    # WILL CHECK IF COSMIC REMOVAL WORKS ON MASKED ARRAYS
    # self.get_clean_mask()
    # self.get_clean_timelike()

    ### Cosmic Ray Hits ###
    if cosmics_removal:
        print("Masking cosmic ray hits...")
        self.mask_cosmics(boxsize=cosmic_boxsize, nsigma=cosmic_nsigma)

    self.get_clean_mask()
    self.get_clean_timelike()

    # this might make a weird reference loop?
    cleaned_lightcurve = self._create_copy()
    cleaned_lightcurve.timelike = {}
    for k in self.timelike:
        if type(self.timelike[k]) == astropy.time.core.Time:
            cleaned = self.timelike[k][self.masks['total'] == 0]
        else:
            cleaned = self.timelike[k][self.masks['total'] == 0] * 1
        cleaned_lightcurve.timelike[k] = cleaned
    # self.cleaned = cleaned_lightcurve

    return cleaned_lightcurve

def get_clean_mask(self):
    """
    Generate a 'total' mask by combining all masks.
    :param self:
    :return:
    """
    self.masks['total'] = np.zeros(self.ntime)
    for mask in self.masks:
        self.masks['total'] = (self.masks['total']!=0) | (self.masks[mask]!=0)
    self.masks['total'] = np.array(self.masks['total'], dtype=int)

def get_clean_timelike(self):
    """
    Create new clean_timelike dictionary with only clean values from timelike.
    :param self:
    :return:
    """
    if "total" not in self.masks:
        self.get_clean_mask()

    self.clean_timelike = {}
    for t in self.timelike:
        if type(self.timelike[t]) == astropy.time.core.Time:
            self.clean_timelike[t] = self.timelike[t][self.masks['total'] == 0]
        else:
            self.clean_timelike[t] = self.timelike[t][self.masks['total'] == 0] * 1

def apply_masks(self, **arrs):
    """
    Apply the total mask to a series of arrays.
    :param self:
    :param arrs: Arrays to be cleaned.
    :return:
    """
    new_arrs = []
    for arr in arrs:
        if "total" in self.masks:
            new_arrs.append(np.ma.masked_where(self.masks['total'],arr))
    return (*new_arrs,)

@property
def clean_time(self):
    """
    The 1D array of time (with astropy units of time).
    """
    return self.clean_timelike.get("time", None)

@property
def clean_flux(self):
    """
    The 2D array of fluxes (row = wavelength, col = time).
    """
    return self.clean_timelike.get("flux", None)

@property
def clean_uncertainty(self):
    """
    The 2D array of uncertainties on the fluxes.
    """
    return self.clean_timelike.get("uncertainty", None)