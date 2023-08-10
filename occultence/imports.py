import numpy as np
import matplotlib
import astropy
from astropy.time import Time

import warnings, textwrap

# copied from chromatic (credit Zach Berta-Thompson):
def custom_formatwarning(message, *args, **kwargs):
    return f"🌈🤖 {textwrap.dedent(str(message)).strip().strip()}\n\n"

original_warning_format = warnings.formatwarning
def cheerfully_suggest(*args, **kwargs):
    warnings.formatwarning = custom_formatwarning
    warnings.warn(*args, **kwargs)
    warnings.formatwarning = original_warning_format

import operator