from ..imports import *

operator_dict = { ">": operator.gt,
                  "<": operator.lt,
                  ">=": operator.ge,
                  "<=": operator.le,
                  "==": operator.eq,
                  "!=": operator.ne}

def clean(self,
          bwboxsize,
          bwthreshvalue,
          thresholds: dict,
          threshold_operators: dict):
    """

    :param self:
    :param bwboxsize:
    :param bwthreshvalue:
    :param thresholds: A dictionary containing the threshold value to apply to the timelike array.
    :param threshold_operators: A dictionary containing the operation to apply to the timelike array and its threshold.
     The options are >, <, >=, <=, == and !=.
    :return:
    """

    self.masks = {}
    self.thresholds = {}

    # Apply a threshold to each key passed in thresholds dictionary. Threshold_operators tells the function which
    # operation to use e.g. >, ==, <= etc.
    for thresh in thresholds:
        if thresh in self.timelike:
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

            self.thresholds[thresh] = thresholds[thresh]
        else:
            # If the thresholds dict contains a key which isn't in the LightCurve.timelike dictionary.
            message = f"""
                        The thresholds dictionary contains the key {thresh} that doesn't exist in the timelike dictionary.
                        The options for thresholds to set are: {", ".join(self.timelike.keys())} 
                            """
            cheerfully_suggest(message)





# BAD WEATHER REMOVAL
#     if bad_weather:
#         try:
#             splc.get_bad_weather_mask(bwbox, bwthresh, lowfluxthresh, rms)
#         except:
#             print("BAD WEATHER FLAG FAILED:")
#             splc.bw_err = np.nanstd(splc.f) * np.ones(len(splc.t))

        # splc.get_bad_conditions_mask(airmass_thresh, bkg_thresh, fwhm_thresh, saturation_counts,ramove_thresh)
        #
        # print(str(100 * np.divide(float(np.count_nonzero(splc.bw_mask)),
        #                           len(splc.bw_mask))) + "% of data is flagged as bad weather")
        # print(str(100 * np.divide(float(np.count_nonzero(splc.am_mask)),
        #                           len(splc.am_mask))) + "% of data is flagged as high airmass")
        # print(str(100 * np.divide(float(np.count_nonzero(splc.bkg_mask)),
        #                           len(splc.bkg_mask))) + "% of data is flagged as high sky background")
        # print(str(100 * np.divide(float(np.count_nonzero(splc.fwhm_mask)),
        #                           len(splc.fwhm_mask))) + "% of data is flagged as high FWHM")
        # print(str(100 * np.divide(float(np.count_nonzero(splc.sat_mask)),
        #                           len(splc.sat_mask))) + "% of data is flagged as saturated")

        # if do_plot:
        #     print("Plotting ALC... " + svplt)
        #     splc.plot_alc(bwthresh, ccd_temps, svplt, show)

    print("Calculating Running RMS...")
    splc.get_running_rms(running_mean_box)


    ### REMOVE DATA DURING DUST CROSSING ISSUES
    if dust_removal:
        if splc.tel == "Ganymede":
            splc.get_dust_mask(start=2458605, finish=2458654)
        if splc.tel == "Callisto":
            splc.get_dust_mask(start=2458384, finish=2458385)

    splc.get_clean_mask()
    splc.get_clean()


def get_running_rms(self, boxsize):
    self.e = running_box(self.t, self.f, boxsize, 'std')
    # self.e = calculate_running_rms(self.t, self.f, boxsize)


def get_bad_weather_mask(self, bwbox, bwthresh, lowfluxthresh, rms):
    bad_weather, err = bad_weather_mask(self.t, self.cf, bwbox, bwthresh, lowfluxthresh, rms)
    self.bw_mask = bad_weather
    self.bw_err = err
    self.bw = True


def get_bad_conditions_mask(self, airmass_thresh, bkg_thresh, fwhm_thresh, sat_thresh, ramove_thresh):
    am_mask = bad_conditions_mask(self.t, self.airmass, airmass_thresh)
    self.am_mask = am_mask

    bkg_mask = bad_conditions_mask(self.t, self.bkg, bkg_thresh)
    self.bkg_mask = bkg_mask

    fwhm_mask = bad_conditions_mask(self.t, self.fwhm, fwhm_thresh)
    self.fwhm_mask = fwhm_mask

    shift_mask = bad_conditions_mask_posneg(self.t, self.ramove, ramove_thresh)
    self.shift_mask = shift_mask

    try:
        sat_mask = bad_conditions_mask(self.t, self.peak, sat_thresh)
        self.sat_mask = sat_mask
    except:
        print("SATURATION FLAG FAILED")


def get_dust_mask(self, start, finish):
    # dust issue on Ganymede from 20190501 - 20190619
    self.dust_mask[(self.t > start) & (self.t < finish)] = 1


def apply_masks(self, f, am_mask=False, bkg_mask=False, fwhm_mask=False, bw_mask=False, dust_mask=False, sat_mask=False,
                flare_mask=False, shift_mask=False):
    masks_bool = [am_mask, bkg_mask, fwhm_mask, bw_mask, dust_mask, sat_mask, flare_mask, shift_mask]
    masks = [self.am_mask, self.bkg_mask, self.fwhm_mask, self.bw_mask, self.dust_mask, self.sat_mask, self.flare_mask,
             self.shift_mask]
    for m in range(len(masks_bool)):
        if masks_bool[m] == True:
            f = np.ma.masked_where(masks[m], f)
    return f


def get_clean_mask(self):
    self.clean_mask = (self.flare_mask != 0) | (self.bw_mask != 0) | (self.cosmic_mask != 0) | (self.dust_mask != 0) | \
                      (np.isnan(self.f)) | (self.am_mask != 0) | (self.bkg_mask != 0) | (self.fwhm_mask != 0) | \
                      (self.sat_mask != 0) | (self.lowflux_mask != 0) | (self.shift_mask != 0)


def get_clean(self):
    self.clean_t = self.t[self.clean_mask == False]
    self.clean_f = self.f[self.clean_mask == False]
    self.clean_e = self.e[self.clean_mask == False]
    self.clean_airmass = self.airmass[self.clean_mask == False]