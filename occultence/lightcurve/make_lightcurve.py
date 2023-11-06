from ..imports import *

class LightCurve:
    def __init__(self,
                 name: str = None,
                 time: list = None,
                 flux: list = None,
                 uncertainty: list = None,
                 timelike: dict = None,
                 metadata: dict = None,
                 **kw):

        self._core_dictionaries = ["timelike", "metadata"]

        # self.metadata = {'name': name}
        self._set_name(name)
        self.metadata['target'] = name
        self.timelike = {}
        self.masks = {}
        self.plot_method = LightCurve.plot_split
        # for m in metadata:
        #     self.metadata[m] = metadata[m]

        if (type(timelike) == dict):
            if (time is not None):
                timelike['time'] = time
            if (flux is not None):
                timelike['flux'] = flux
            if (uncertainty is not None):
                timelike['uncertainty'] = uncertainty

            self._initialize_from_dictionaries(
                timelike=timelike,
                metadata=metadata)
        # then try to initialize from arrays
        elif (time is not None) and (flux is not None) and (uncertainty is not None):
            self._initialize_from_arrays(
                time=time,
                flux=flux,
                uncertainty=uncertainty,
                **kw,
            )

        if metadata is not None:
            self.metadata.update(**metadata)

        if self.time is not None:
            if type(self.time[0]) != astropy.time.core.Time:
                message = f"""
                            Warning! The time array is not an astropy.Time object, therefore there is no info about the 
                            format or scale.
                            We will assume that it is JD and TDB from here on!
                            """
                cheerfully_suggest(message)
                self.timelike['time'] = Time(self.timelike['time'] * 1,
                                             format='jd',
                                             scale='tdb')

        # for k in timelike:
        #     if np.shape(timelike[k]) == np.shape(self.time):
        #         self.timelike[k] = timelike[k] * 1
        #     else:
        #         message = f"""
        #         Something doesn't line up!
        #         The timelike array {k} has a shape of {timelike[k]}.
        #         The time array has {self.ntime} times.
        #         """
        #         cheerfully_suggest(message)

    def __repr__(self):
        return f"<🌟 Lightcurve {self.metadata['name']} ({self.ntime}t) 🌟>"
    def _sort(self):
        """
        Sort all timelike quantities in light curve by time
        :return:
        """
        i_time = np.argsort(self.time)
        if "original_time_index" not in self.timelike:
            self.timelike["original_time_index"] = np.arange(self.ntime)

        for k in self.timelike:
            if self.timelike[k] is not None:
                self.timelike[k] = self.timelike[k][i_time]

    def _set_name(self, name):
        if hasattr(self, 'metadata'):
            self.metadata['name'] = name
        else:
            self.metadata = {'name':name}

    @property
    def name(self):
        """
        The name of this `Rainbow` object.
        """
        return self.metadata.get("name", None)

    @property
    def shape(self):
        """
        The shape of the flux array (ntimes).
        """
        return (self.ntime,)

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


    def _initialize_from_arrays(
        self, time=None, flux=None, uncertainty=None, **kw
    ):
        """
        Populate from arrays.

        Parameters
        ----------
        time : Quantity, Time, optional
            A 1D array of times, in any unit.
        flux : array, optional
            A 2D array of flux values.
        uncertainty : array, optional
            A 2D array of uncertainties, associated with the flux.
        **kw : dict, optional
            Additional keywords will be interpreted as arrays
            that should be sorted into the appropriate location
            based on their size.
        """
        # store the time, flux + uncertainty
        self.timelike = {'time': time,
                         'flux': flux * 1,
                         'uncertainty': uncertainty * 1}

        for k, v in kw.items():
            if type(v) == astropy.time.core.Time:
                self.timelike[k] = v
            else:
                self.timelike[k] = v * 1

        # validate that something reasonable got populated
        self._validate_core_dictionaries()

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
            if type(timelike[k]) == astropy.time.core.Time:
                self.timelike[k] = timelike[k]
            else:
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

    def _get_core_dictionaries(self):
        """
        Get the core dictionaries of this LightCurve.

        Returns
        -------
        core : dict
            Dictionary containing the keys
            ['timelike', 'metadata']
        """
        return {k: vars(self)[k] for k in self._core_dictionaries}

    def __getattr__(self, key):
        """
        If an attribute/method isn't explicitly defined,
        try to pull it from one of the core dictionaries.

        Let's say you want to get the 2D uncertainty array
        but don't want to type `self.fluxlike['uncertainty']`.
        You could instead type `self.uncertainty`, and this
        would try to search through the four standard
        dictionaries to pull out the first `uncertainty`
        it finds.

        Parameters
        ----------
        key : str
            The attribute we're trying to get.
        """
        if key == "_core_dictionaries":
            try:
                return self.__dict__[key]
            except KeyError:
                pass
            except RecursionError:
                return
        else:
            if key not in self._core_dictionaries:
                for dictionary_name in self._core_dictionaries:
                    try:
                        return self.__dict__[dictionary_name][key]
                    except KeyError:
                        pass
        message = f".{key} does not exist for this LightCurve"
        raise AttributeError(message)

    def remove_nans(self,):
        new_lc = self._create_copy()

        if np.count_nonzero(~np.isfinite(new_lc.flux)) > 0:
            ind_finite = np.where(np.isfinite(new_lc.flux) == True)
            for k in new_lc.timelike.keys():
                new_lc.timelike[k] = new_lc.timelike[k][ind_finite]

        if np.count_nonzero(~np.isfinite(new_lc.uncertainty)) > 0:
            ind_finite = np.where(np.isfinite(new_lc.uncertainty) == True)
            for k in new_lc.timelike.keys():
                new_lc.timelike[k] = new_lc.timelike[k][ind_finite]

        return new_lc

    def plot(self, **kw):
        return self.plot_method(self,**kw)

    def plot_all(self, ax=None, figsize=(12,4), ylims=[0.98,1.02], color='C0', **kw):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        ax.plot(self.time.value, self.flux, '.', color=color **kw)
        ax.errorbar(self.time.value, self.flux, self.uncertainty, fmt='.', color=color, alpha=0.1, **kw)
        ax.set_ylim(ylims[0], ylims[1])
        ax.set_ylabel("Flux")
        ax.set_xlabel("Time [d]")
        ax.legend()
        return ax

    def plot_split(self, ax=None, figsize=(36,4), ylims=[0.98,1.02], color="C0", **kw):
        i_split, _ = self.split_time()
        if ax is None:
            fig, ax = plt.subplots(ncols = len(i_split)-1, figsize=figsize, sharey=True)

        for i,(i0, i1) in enumerate(zip(i_split[:-1],i_split[1:])):
            ax[i].plot(self.time.value[i0:i1], self.flux[i0:i1], '.', color=color, **kw)
            ax[i].errorbar(self.time.value[i0:i1], self.flux[i0:i1], self.uncertainty[i0:i1], fmt='.', color=color,
                           alpha=0.1, **kw)
        ax[0].set_ylim(ylims[0], ylims[1])
        ax[0].set_ylabel("Flux")
        ax[-1].set_xlabel("Time [d]")
        ax[-1].legend()
        return ax

    def phasefold(self, period):
        new_lc = self._create_copy()
        t = self.time.value
        ph_d = ((t - min(t)) % period) / period
        ph_d[ph_d > 0.5] -= 1
        new_lc.timelike['phase'] = ph_d
        return new_lc



    from .remove_transit import (
        mask_existing_transit
    )

    from ..cleaning import (
        clean,
        mask_timelike_threshold,
        mask_bad_weather,
        mask_cosmics,
        mask_dust,
        get_clean_mask,
        get_clean_timelike,
        apply_masks,
        # clean_time,
        # clean_flux,
        # clean_uncertainty
    )
    from ..binning import (
        bin,
        split_time
    )
    from ..lightcurve_detrending import (
        gp_detrend,
    )
    from ..transit_detecting import (
        find_transits,
        bls,
    )
    from ..injection import (
        inject_transit,
        pool_inject_transit,
        create_lots_of_transit_params,
        inject_lots_of_transits,
    )

    from ..recovery import (
        was_injected_planet_recovered,
        full_injection_recovery,
        single_injection_recovery,
        split_lightcurve,
        was_planet_observed,
    )

    # from ..flare_finding import *
    # from ..read import *
    # from ..write import *
    # from recovery import *


