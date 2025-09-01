from ..imports import *
from astropy.timeseries import LombScargle
from scipy.stats import norm

def lombscargle_detrend(self, ls_binning=30*u.minute, plot=False, verbose=True, **ls_kw):
    binned_lc = self.bin(dt=ls_binning).remove_nans()
    detrended_lightcurve = self._create_copy()

    ls, periods, power, fap = binned_lc.lombscargle(plot=plot, **ls_kw)
    model = ls.model(self.time, 1/periods[0])

    if power[0] < fap:
        if 'nsigma' in ls_kw.keys():
            if verbose:
                print(f"None of the recovered periods are > the false-alarm probability of {fap} ({ls_kw['nsigma']} sigma)")
        else:
            if verbose:
                print(f"None of the recovered periods are > the false-alarm probability of {fap} (2 sigma)")

        detrended_lightcurve.metadata['ls_period'] = periods[0]
        detrended_lightcurve.metadata['ls_power'] = power[0]
        detrended_lightcurve.metadata['fap'] = fap
        detrended_lightcurve.metadata['recovered_period'] = False
    else:
        if verbose:
            print(f"The best recovered period={periods[0]}")
        detrended_lightcurve.metamethods['ls_model'] = ls.model
        detrended_lightcurve.metadata['ls_period'] = periods[0]
        detrended_lightcurve.metadata['ls_power'] = power[0]
        detrended_lightcurve.metadata['ls'] = ls
        detrended_lightcurve.metadata['fap'] = fap
        detrended_lightcurve.metadata['recovered_period'] = True
        detrended_lightcurve.timelike['ls_model'] = model
        detrended_lightcurve.timelike['original_flux'] = detrended_lightcurve.timelike['flux'] * 1
        detrended_lightcurve.timelike['flux'] = detrended_lightcurve.timelike['flux'] / model
        detrended_lightcurve._set_name(detrended_lightcurve.name + "_lsdetrend")

    if plot:
        axs = self.plot()
        i_split, t_split = self.split_time(split=self.split_by)
        for i, ax in enumerate(axs):
            ax.plot(self.time.value[i_split[i]:i_split[i + 1]],
                    model[i_split[i]:i_split[i + 1]], c='orange', zorder=10)

    return detrended_lightcurve

def lombscargle(self, plot=False, npeaks=3, minimum_frequency = 0.03/u.d, maximum_frequency = 12/u.d, nsigma=2,
                ylims=[0.99, 1.01], save_periodogram=None, **ls_kw):
    # minimum_frequency = 0.03  # P = 33.3d
    # maximum_frequency = 12  # P = 0.8333d = 2hrs

    ls = LombScargle(self.time, self.flux, dy=self.uncertainty)
    frequency, power = ls.autopower(minimum_frequency=minimum_frequency,
                                    maximum_frequency=maximum_frequency,
                                    **ls_kw)
    periods = 1 / frequency
    highest_peaks = np.argsort(power)[::-1]


    if plot:
        plt.plot(1 / frequency, power, 'k')
        plt.ylabel("LS Power")
        plt.xlabel("Period [d]")
        for i, p in enumerate(highest_peaks[:npeaks]):
            period = periods[p]
            plt.axvline(period.to_value('d'), color=f"C{i}", label=f"{i + 1} Most-Likely Period = {period:.2f}")
        probability = 1-(norm.cdf(nsigma) - norm.cdf(-nsigma))

        fap = ls.false_alarm_level(probability)
        plt.axhline(fap, color='r', linestyle='--', label=f"False Alarm Probability of {nsigma}$\sigma$")
        plt.legend()
        if save_periodogram is not None:
            plt.savefig(save_periodogram)

        fig, ax = plt.subplots(figsize=(6, 4*npeaks), nrows=npeaks, ncols=2, sharex='col', sharey=True, width_ratios=[2, 1])
        t_fit = np.linspace(self.time[0], self.time[-1], 1000)
        for peak in range(npeaks):
            y_fit = ls.model(t_fit, frequency[highest_peaks[peak]])
            if npeaks > 1:
                plt.sca(ax[peak, 0])
            else:
                plt.sca(ax[0])
            plt.plot(self.time.value, self.flux, 'k.')
            plt.plot(t_fit.value, y_fit)
            plt.ylim(ylims[0], ylims[1])
            plt.title(f"{peak+1} Most Likely P={periods[highest_peaks[peak]]:.2f}d")
            plt.ylabel("Relative Flux")

            if npeaks>1:
                plt.sca(ax[peak, 1])
            else:
                plt.sca(ax[1])

            ph = phase(t=self.time.value, period=periods[highest_peaks[peak]].value)
            ph_fit = phase(t=t_fit.value, period=periods[highest_peaks[peak]].value)
            plt.plot(ph, self.flux, 'k.')
            plt.plot(ph_fit, y_fit, '.')
            plt.title(f"Phase-Folded")

        plt.tight_layout()


    return ls, periods[highest_peaks], power[highest_peaks], fap

def phase(t, period):
    ph_d = ((t - np.min(t)) % period) / period
    ph_d[ph_d > 0.5] -= 1
    return ph_d