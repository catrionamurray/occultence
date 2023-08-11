from ..imports import *

class LightCurve:
    def __init__(self,
                 name: str,
                 time: list,
                 flux: list,
                 uncertainty: list,
                 timelike: dict = {},
                 metadata: dict = {}):

        self.metadata = {'name': name}
        for m in metadata:
            self.metadata[m] = metadata[m]

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

    def _validate_core_dictionaries(self):
        """
        Do some simple checks to make sure this LightCurve
        is populated with the minimal data needed to do anything.
        It shouldn't be run before the LightCurve is fully
        initialized; otherwise, it might complain about
        a half-populated object.
        """

        # make sure there are some times defined
        if self.ntime is None:
            cheerfully_suggest(
                f"""
            No times are defined for this LightCurve.
            """
            )

        # does the flux have the right shape?
        if (self.shape != np.shape(self.flux)) or (np.shape(self.time)!=np.shape(self.flux) ) :
            message = f"""
            Something doesn't line up!
            The flux array has a shape of {np.shape(self.flux)}.
            The time array has {self.ntime} times.
            """
            if self.shape == np.shape(self.flux)[::-1]:
                cheerfully_suggest(
                    f"""{message}
                    Any chance your flux array is transposed?
                    """
                )
            else:
                cheerfully_suggest(message)

        for n in ["uncertainty"]:
            x = getattr(self, n)
            if x is not None:
                if x.shape != np.shape(self.flux):
                    message = f"""
                    Watch out! The '{n}' array has
                    a shape of {x.shape}, which doesn't match the
                    flux array's shape of {np.shape(self.flux)}.
                    """
                    cheerfully_suggest(message)

        self._sort()

    def _initialize_from_dictionaries(
        self, timelike={}, metadata={}
    ):
        """
        Populate from dictionaries in the correct format.

        Parameters
        ----------
        timelike : dict
            A dictionary containing 1D arrays with the same
            shape as the time axis. It must at least
            contain the keys 'time', 'flux' and 'uncertainty.
        metadata : dict
            A dictionary containing all other metadata
            associated with the dataset, generally lots of
            individual parameters or comments.
        """

        # update the three core dictionaries of arrays
        for k in timelike:
            self.timelike[k] = timelike[k] * 1
        # multiplying by 1 is a kludge to prevent accidental links

        # update the metadata
        self.metadata.update(**metadata)

        # validate that something reasonable got populated
        self._validate_core_dictionaries()


    def _create_copy(self):
        """
        Create a copy of self, with the core dictionaries copied.
        """
        new = type(self)()
        new._initialize_from_dictionaries(
            **copy.deepcopy(self._get_core_dictionaries())
        )
        return new


