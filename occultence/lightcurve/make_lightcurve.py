from ..imports import *

class LightCurve:
    def __init__(self,
                 name: str,
                 time: list,
                 flux: list,
                 telescope: str,
                 filter: str,
                 dates: list,
                 uncertainty: list,
                 timelike: dict = {}):

        self.metadata = {'name': name,
                         'telescope': telescope,
                         'filter': filter,
                         'dates': dates}

        if type(time[0]) != astropy.time.core.Time:
            message = f"""
                        Warning! The time array is not an astropy.Time object, therefore there is no info about the 
                        format or scale.
                        We will assume that it is JD and TDB from here on!
                        """
            cheerfully_suggest(message)
            time = Time(time, format='jd', scale='tdb')

        self.timelike = {'time': time * 1,
                         'flux': flux * 1,
                         'uncertainty': uncertainty * 1}

        for k in timelike:
            if np.shape(timelike[k]) == np.shape(self.time):
                self.timelike[k] = timelike[k] * 1
            else:
                message = f"""
                Something doesn't line up!
                The timelike array {k} has a shape of {timelike[k]}.
                The time array has {self.ntime} times.
                """
                cheerfully_suggest(message)

    def __repr__(self):
        return f"<Lightcurve {self.metadata['name']}({self.ntime}t)>"
    def _sort(self):
        i_time = np.argsort(self.time)
        if "original_time_index" not in self.timelike:
            self.timelike["original_time_index"] = np.arange(self.ntime)

        for k in self.timelike:
            if self.timelike[k] is not None:
                self.timelike[k] = self.timelike[k][i_time]
    @property
    def name(self):
        """
        The name of this `Rainbow` object.
        """
        return self.metadata.get("name", None)

    @property
    def dates(self):
        """
        The name of this `Rainbow` object.
        """
        return self.metadata.get("dates", None)

    @property
    def time(self):
        """
        The 1D array of time (with astropy units of time).
        """
        return self.timelike.get("time", None)

    @property
    def flux(self):
        """
        The 2D array of fluxes (row = wavelength, col = time).
        """
        return self.timelike.get("flux", None)

    @property
    def uncertainty(self):
        """
        The 2D array of uncertainties on the fluxes.
        """
        return self.timelike.get("uncertainty", None)

    @property
    def ntime(self):
        """
        The number of times.
        """
        if self.time is None:
            return 0
        else:
            return len(self.time)

    @property
    def dt(self):
        """
        The typical timestep.
        """
        if self.time is None:
            return None
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return np.nanmedian(np.diff(self.time))


