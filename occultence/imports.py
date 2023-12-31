import numpy as np
import copy
import astropy
from astropy.time import Time
from astropy.timeseries import TimeSeries, aggregate_downsample
import astropy.units as u
from .utils import *
# import fitsio
import math
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from astropy.io import fits, ascii
from astropy.stats import sigma_clip
import pickle
import operator
import george
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from astropy.timeseries import BoxLeastSquares

import warnings, textwrap

# copied from chromatic (credit Zach Berta-Thompson):
def custom_formatwarning(message, *args, **kwargs):
    return f"🌈🤖 {textwrap.dedent(str(message)).strip().strip()}\n\n"

original_warning_format = warnings.formatwarning
def cheerfully_suggest(*args, **kwargs):
    warnings.formatwarning = custom_formatwarning
    warnings.warn(*args, **kwargs)
    warnings.formatwarning = original_warning_format
