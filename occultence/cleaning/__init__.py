from ..imports import *

# ## THRESHOLDS - make these user-defined
# filt = "I+z"
# date = "20*"
# targ = ""
# bwthresh = 0.06  # 0.08
# lowfluxthresh = 0.6  # 0.2
# bwbox = 0.01 * 2
# nsigma = 4
# running_mean_box = bwbox  # 0.04
# airmass_thresh = 2.5
# bkg_thresh = 3000
# fwhm_thresh = 7.5
# saturation_counts = 61000
# ramove_thresh = 5
# binsize = 5
# transit_durations = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
# rms = True
# bad_weather = True
# dust_removal = True
# night_plot = True
# import_global_lc = True
# scale = False
# correct = True
# clip_gp = True
# night_gp = False

from clean import *
from bad_weather import *
from timelike_masks import *
from dust import *
from cosmics import *

