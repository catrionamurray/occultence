from ..imports import *

def find_transits(self, transit_durations=0.01, minimum_period=0.5, maximum_period=30, limitperiod=False,
                  obj='likelihood', oversample=100.0, minpower=5, plot=True, figsize=(12,4)):
    bls_lightcurve = self._create_copy()
    transit_pd = {"period": [], "depth": [], 'duration': [], 'epoch':[], 'epoch_start':[], 'epoch_end':[], 'snr': []}
    bls_f_model, transit_params, stats, BLS_obj = self.bls(transit_durations=transit_durations,
                                                           minimum_period= minimum_period,
                                                           maximum_period=maximum_period,
                                                           limitperiod=limitperiod,
                                                           obj=obj,
                                                           oversample=oversample,
                                                           minpower=minpower)

    if len(stats)>0:
        bls_transits = stats['per_transit_count']
        transits = np.where(np.array(bls_f_model) < np.nanmedian(bls_f_model))[0]

        trans_num = 0
        transit_found = True
        if len(bls_transits) > 0:
            for b in range(len(bls_transits)):
                if bls_transits[b] > 0:
                    mid_transit = stats['transit_times'][b]#.value
                    transit_start = mid_transit - (0.5 * transit_params[2])#.to_value('d'))
                    transit_end = mid_transit + (0.5 * transit_params[2])#.to_value('d'))

                    recovered_per = transit_params[0]#np.log10(transit_params[0].to_value('d'))
                    recovered_dur = transit_end - transit_start
                    recovered_depth = stats['depth'][0]
                    snr = (recovered_depth * np.sqrt(bls_transits[b])) / np.nanmedian(self.uncertainty)

                    transit_pd["period"].append(recovered_per)
                    transit_pd["depth"].append(recovered_depth)
                    transit_pd['duration'].append(recovered_dur)
                    transit_pd['epoch'].append(mid_transit)
                    transit_pd['epoch_start'].append(transit_start)
                    transit_pd['epoch_end'].append(transit_end)
                    transit_pd['snr'].append(snr)

                    if plot:
                        plt.figure(figsize=figsize)
                        plt.errorbar(self.time.value, self.flux, self.uncertainty, color='k', fmt='.')
                        plt.plot(self.time.value, bls_f_model, 'orange')
                        plt.axvline(transit_start.value)
                        plt.axvline(transit_end.value)
                        plt.plot(self.time.value[transits], self.flux[transits], 'b.')
                        plt.title("SNR = %0.2f, Period = %0.2f" % (snr, recovered_per.value))
                        plt.xlim(transit_start.value - 0.2, transit_start.value + 0.2)
                        plt.show()
                        plt.close()

                    trans_num = trans_num + 1
    else:
        transit_found=False
        transits = []

    bls_lightcurve.timelike['flux'] = self.timelike['flux'] / bls_f_model
    bls_lightcurve.timelike['BLS_model'] = bls_f_model
    bls_lightcurve.metadata['BLS_transits_found'] = transit_found
    bls_lightcurve.metadata['BLS_transits_ind'] = transits
    bls_lightcurve.metadata['BLS_transits_params'] = transit_pd

    bls_lightcurve._set_name(bls_lightcurve.name + "_bls")

    return  bls_lightcurve


def bls(self, transit_durations, minimum_period, maximum_period, limitperiod, obj, oversample, minpower):
    print("Running BLS Search")

    if limitperiod:
        maximum_period = min(30, max(self.time) - min(self.time))

    periods = np.linspace(minimum_period, maximum_period, num=10000)
    BLS_d = BoxLeastSquares(self.time, self.flux, dy=self.uncertainty)
    pg_d = BLS_d.power(periods, transit_durations, objective=obj, oversample=oversample)
    pers, power_d, epoch_d, depth_d, durs = pg_d.period, pg_d.power, pg_d.transit_time, pg_d.depth, pg_d.duration

    print("Number of Periods Checked: ", len(pers))
    max_power = np.argmax(power_d)

    if np.count_nonzero(np.isfinite(power_d)) > 0:

        if power_d[max_power] > minpower:
            best_per_d = pers[max_power]
            best_t0_d = epoch_d[max_power]
            best_dur_d = durs[max_power]
            best_depth_d = depth_d[max_power]
            stats = BLS_d.compute_stats(best_per_d, best_dur_d, best_t0_d)

            f_model = BLS_d.model(t_model=self.time, period=best_per_d, duration=best_dur_d, transit_time=best_t0_d)
            return f_model, [best_per_d, best_t0_d, best_dur_d, best_depth_d], stats, BLS_d
        else:
            print("No transits detected!")
            return np.ones(len(self.time)), [], [], []

    else:
        print("No transits detected!")
        return np.ones(len(self.time)), [], [], []



