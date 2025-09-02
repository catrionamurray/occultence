import numpy as np

from ..imports import *
from .flare_model import flare_model_mendoza2022, flare_model_davenport2014
import scipy.optimize as opt
from scipy.integrate import simpson

def detect_flares():
    return


def detect_flares_sclip(self, min_flare_duration=5 * u.minute, min_flare_separation=10*u.minute,
                        n_consecutive_points=3,
                        n_points=3,
                        global_sc_kw={'nsigma_upper': 5, 'nsigma_lower': 5},
                        local_sc_kw={'nsigma_upper': 3, 'nsigma_lower': np.inf, 'running_median_boxsize': 0.4},
                        tbuffer_before = 5*u.minute,
                        tbuffer_after = 20*u.minute,
                        plot=False):
    clip_lc = self._create_copy()
    clip_lc = clip_lc.global_and_local_sigma_clip(global_sc_kw=global_sc_kw, local_sc_kw=local_sc_kw)

    # isolate every sigma-clipped point
    outliers = np.where(np.isnan(clip_lc.flux))[0]
    regions = []

    if len(outliers) == 0:
        print("No flares found.")

    else:
        current_region = [outliers[0]]

        # loop over every detected outlier region
        for o in outliers[1:]:
            # check if the outlier regions are separated by more than the user-defined minimum flare separation
            if (clip_lc.time[o] - clip_lc.time[current_region[-1]]) > min_flare_separation:
                diffs = np.array([(b - a) for a, b in zip(current_region[:-1], current_region[1:])])
                if (len(current_region) > n_points) & (np.count_nonzero(diffs == 1) > n_consecutive_points):
                    t = clip_lc.time[current_region]
                    dur = np.max(t) - np.min(t)
                    if (dur > min_flare_duration) and (dur < self.split_by):
                        regions.append(current_region)
                current_region = [o]
            else:
                current_region.append(o)

        if len(regions) > 0:
            diffs = np.array([(b - a) for a, b in zip(current_region[:-1], current_region[1:])])
            if (len(current_region) >= n_points) and (current_region[-1] != regions[-1][-1]) and\
                    (np.count_nonzero(diffs == 1) >= n_consecutive_points):
                t = clip_lc.time[current_region]
                dur = np.max(t) - np.min(t)
                if (dur > min_flare_duration) and (dur < self.split_by):
                    regions.append(current_region)

            print(f"""
            There are {len(regions)} flare candidates with >{min_flare_duration} over {local_sc_kw['nsigma_upper']}$\sigma$
            """)

            flare_regions = [[r[0], r[-1]] for r in regions]

            # flares = targ_flares_removed.metadata['flares_detected']
            nfl = len(flare_regions)

            buffered_flare_regions = []
            for i, f in enumerate(flare_regions):
                # Add user-defined buffer before and after flare
                t_before_flare = self.time[f[0]] - tbuffer_before
                t_after_flare = self.time[f[1]] + tbuffer_after
                start_t = np.where(self.time > t_before_flare)[0][0]
                end_t = np.where(self.time < t_after_flare)[0][-1]
                buffered_flare_regions.append([start_t, end_t])

            if plot:
                ax = clip_lc.plot_all(figsize=(20, 8))
                self.extract(outliers).plot_all(ax=ax, color='orange', label=f">{local_sc_kw['nsigma_upper']}$\sigma$")
                for i, r in enumerate(regions):
                    if i == 0:
                        self.extract(r).plot_all(ax=ax, color='red',
                                                 label=f">{min_flare_duration} >{local_sc_kw['nsigma_upper']}$\sigma$")
                    else:
                        self.extract(r).plot_all(ax=ax, color='red')

                plt.ylim(0.98, np.nanmax(self.flux))

                fig, ax = plt.subplots(ncols=nfl, figsize=(nfl * 4, 3))

                for i, f in enumerate(flare_regions):
                    if nfl > 1:
                        plt.sca(ax[i])
                    else:
                        plt.sca(ax)

                    i_start = find_nearest(self.time.value, self.time.value[start_t] - 0.35)
                    i_end = find_nearest(self.time.value, self.time.value[end_t] + 0.35)
                    plt.plot(self.time.value[i_start:i_end], self.flux[i_start:i_end], 'k.')
                    plt.plot(self.time.value[start_t:end_t], self.flux[start_t:end_t], 'g.')
                    plt.plot(self.time.value[f[0]:f[1]], self.flux[f[0]:f[1]], 'r.')

            clip_lc.metadata['flares_detected'] = {'nflares': len(buffered_flare_regions),
                                                   'thresholds': {'duration': min_flare_duration,
                                                                 'sigma': local_sc_kw['nsigma_upper'],
                                                                 'separation': min_flare_separation},
                                                   'flare_regions': buffered_flare_regions,
                                                   'flare_frequency': len(buffered_flare_regions)/self.total_time_observed,
                                                   'global_sc_kw':global_sc_kw,
                                                   'local_sc_kw':local_sc_kw
                                          }
            clip_lc.timelike['flux'] = self.flux.copy()
            for r in buffered_flare_regions:
                clip_lc.timelike['flux'][r[0]:r[1]] = np.nan
            clip_lc._set_name(clip_lc.name + "_flaresremoved")
        else:
            print("No flares found.")
    return clip_lc


def fit_flare(t, f, unc, p0, bounds, t_new=None, method="mendoza2022", plot=False, ax=None):
    if method == "mendoza2022":
        mod = flare_model_mendoza2022
    elif method == "davenport2014":
        mod = flare_model_davenport2014
    w, _ = opt.curve_fit(mod, t, f, sigma=unc, p0=p0, bounds=bounds)
    # print("Estimated Parameters", w)

    yfit = mod(t, *w)

    if t_new is not None:
        ynew = mod(t_new, *w)
    else:
        ynew = None

    if plot:
        if ax is None:
            fig, ax = plt.subplots()
        plt.sca(ax)
        plt.errorbar(t, f, unc, fmt='k.')
        plt.plot(t, yfit)
        plt.plot(t_new, ynew, c="green", label="Optimized Fit")
        initial_guess = mod(t_new, p0[0], p0[1], p0[2])
        plt.plot(t_new, initial_guess, c='orange', label="Initial Guess")
        plt.legend()

    return yfit, ynew, w


def model_each_flare(self, flares_to_model, t, f, unc, n_before=3, n_after=5, method="mendoza2022",
                     tgap=(1 * u.hour).to_value('d'), plot=False):

    model_flares = self._create_copy()
    model_flares.metadata['flares_detected']['flares_to_model'] = flares_to_model
    model_flares.metadata['flares_detected']['flares_model_used'] = method

    n = np.count_nonzero(flares_to_model)
    if plot:
        fig, ax = plt.subplots(ncols=n, figsize=(3*n, 3))
        if n == 1:
            ax = [ax]
    else:
        ax = [None] * n

    count = 0
    tpeaks, fwhms, amps, integ = [], [], [], []
    for b, (f_s, f_e) in zip(flares_to_model, self.metadata['flares_detected']['flare_regions']):
        # REGRESSION ------------------------------------------------------------------
        if b > 0:

            start_t = f_s - n_before
            while t[f_s] - t[f_s - n_before] > tgap:
                n_before = n_before - 1
                start_t = f_s - n_before

            end_t = f_e + n_after
            while t[f_e + n_after] - t[f_e] > tgap:
                n_after = n_after - 1
                end_t = f_e + n_after

            if n_before == 0:
                yscale = 1
            else:
                yscale = np.nanmean(f[start_t:f_s - 1])

            xsamp, ysamp, errsamp = t[start_t:end_t], f[start_t:end_t]-yscale, unc[start_t:end_t]
            peak_t = xsamp[np.argmax(ysamp)]
            peak_f = np.max(ysamp)
            p0 = [peak_t, 0.001, peak_f]  # guessed params
            bounds = ([t[start_t], 0, 0], [t[end_t], 100, 100])

            t_new = np.linspace(t[start_t], t[end_t], 200)
            y_fit, y_new, w = fit_flare(xsamp, ysamp, errsamp, p0, bounds, t_new,method=method, plot=plot, ax=ax[count])
            tpeaks.append(w[0])
            fwhms.append(w[1])
            amps.append(w[2])

            area = simpson(y_new, x=t_new)
            integ.append(area)
            # print("area =", area)
            count=count+1

        else:
            tpeaks.append(np.nan)
            fwhms.append(np.nan)
            amps.append(np.nan)
            integ.append(np.nan)

    model_flares.metadata['flares_detected']['model_parameters'] = {'tpeak': tpeaks, 'fwhm': fwhms, 'amp': amps,
                                                                    'integral':integ}
    return model_flares

def manual_clean_flares(self, n_keep):
    self.metadata['flares_detected']['flare_regions'] = [self.metadata['flares_detected']['flare_regions'][n] for n in n_keep]
    self.metadata['flares_detected']['nflares'] = len(n_keep)
    self.metadata['flares_detected']['flare_frequency'] = self.metadata['flares_detected']['nflares']/self.total_time_observed
    return self