from ..imports import *
from astropy.timeseries import LombScargle

def lombscargle(self, plot=False, npeaks=3, minimum_frequency = 0.03/u.d, maximum_frequency = 12/u.d, nsigma=2, **ls_kw):
    # minimum_frequency = 0.03  # P = 33.3d
    # maximum_frequency = 12  # P = 0.8333d = 2hrs

    ls = LombScargle(self.time, self.flux, dy=self.uncertainty)
    frequency, power = ls.autopower(minimum_frequency=minimum_frequency,
                                    maximum_frequency=maximum_frequency,
                                    **ls_kw)
    # period = 1 / frequency[np.argmax(power)]
    highest_peaks = np.argsort(power)[::-1]


    if plot:
        plt.plot(1 / frequency, power, 'k')
        plt.ylabel("LS Power")
        plt.xlabel("Period [d]")
        for i, p in enumerate(highest_peaks[:npeaks]):
            period = (1 / frequency[p])
            plt.axvline(period.to_value('d'), color=f"C{i}", label=f"{i + 1} Most-Likely Period = {period:.2f}")
        two_sigma = 1 - 0.9544
        three_sigma = 1 - 0.9973
        probability = {2:two_sigma, 3:three_sigma}

        fap = ls.false_alarm_level(probability[nsigma])
        plt.axhline(fap,color='r', linestyle='--', label=f"False Alarm Probability of {nsigma}$\sigma$")
        plt.legend()


        fig, ax = plt.subplots(nrows=npeaks, sharex=True, sharey=True)
        t_fit = np.linspace(self.time[0], self.time[-1], 1000)
        for peak in range(npeaks):
            y_fit = ls.model(t_fit, frequency[highest_peaks[peak]])
            if npeaks>1:
                plt.sca(ax[peak])
            plt.plot(self.time.value, self.flux, 'k.')
            plt.plot(t_fit.value, y_fit)
            plt.ylim(0.98, 1.02)
            plt.title(f"{peak+1} Most Likely P={(1 / frequency[highest_peaks[peak]]):.2f}d")
        plt.tight_layout()


    return ls, 1/frequency[highest_peaks], power[highest_peaks]