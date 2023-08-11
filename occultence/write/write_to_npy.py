from ..imports import *
from ..version import version

def to_rainbow_npy(self, filepath, **kw):
    """
    Write a LightCurve to a file in the .lightcurve.npy format.

    Parameters
    ----------

    self : LightCurve
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    assert ".lightcurve.npy" in filepath

    # populate a dictionary containing the four core dictionaries
    dictionary_to_save = self._get_core_dictionaries()

    # save that to a file
    np.save(filepath, [dictionary_to_save, version()], allow_pickle=True)