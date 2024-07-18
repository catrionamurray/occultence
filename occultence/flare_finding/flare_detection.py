from ..imports import *

def detect_flares():
    return


def detect_flares_sclip(self, min_flare_duration=5 * u.minute, min_flare_separation=10*u.minute,
                        global_sc_kw={'nsigma_upper': 5, 'nsigma_lower': 5},
                        local_sc_kw={'nsigma_upper': 3, 'nsigma_lower': np.inf, 'running_mean_boxsize': 0.4},
                        plot=False):
    clip_lc = self._create_copy()
    clip_lc = clip_lc.global_and_local_sigma_clip(global_sc_kw=global_sc_kw, local_sc_kw=local_sc_kw)

    # isolate every sigma-clipped point
    # isolate every sigma-clipped point
    outliers = np.where(np.isnan(clip_lc.flux))[0]
    current_region = [outliers[0]]
    regions = []

    for o in outliers[1:]:
        # print(o)
        if (clip_lc.time[o] - clip_lc.time[current_region[-1]]) > min_flare_separation:
            if len(current_region) > 1:
                t = clip_lc.time[current_region]
                dur = np.max(t) - np.min(t)
                if (dur > min_flare_duration) and (dur < self.split_by):
                    regions.append(current_region)
            current_region = [o]
        else:
            current_region.append(o)

    if (len(current_region) > 1) and (current_region[-1] != regions[-1][-1]):
        t = clip_lc.time[current_region]
        dur = np.max(t) - np.min(t)
        if (dur > min_flare_duration) and (dur < self.split_by):
            regions.append(current_region)

    print(
        f"There are {len(regions)} flare candidates with >{min_flare_duration} over {local_sc_kw['nsigma_upper']}$\sigma$")

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

    clip_lc.metadata['flares_detected'] = {'nflares': len(regions), 'threshold_duration': min_flare_duration,
                                  'flare_regions': [[r[0],r[-1]] for r in regions], 'threshold_sigma': local_sc_kw['nsigma_upper'],
                                  }
    clip_lc._set_name(clip_lc.name + "_flaresremoved")
    return clip_lc
