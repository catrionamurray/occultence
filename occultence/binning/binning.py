from ..imports import *

def bin(self, dt, bin_func=np.nanmedian, **kw):
    # ts = TimeSeries(self.timelike)
    # bin_ts = aggregate_downsample(ts, time_bin_size=dt, aggregate_func=bin_func, **kw)
    # replaced the above line with
    # begin replaced
    time_values = self.timelike["time"].value

    cuts = np.where(time_values[1:] - time_values[:-1] > 0.5)[0]
    cuts = np.hstack((0, cuts + 1, len(time_values)))

    for i in range(len(cuts) - 1):
        t1 = self.timelike.copy()

        for key in t1.keys():
            t1[key] = self.timelike[key][cuts[i]:cuts[i + 1]]

        t1s = TimeSeries(t1)

        bin_ts1 = aggregate_downsample(t1s, time_bin_size=dt, aggregate_func=bin_func, **kw)

        if i == 0:
            bin_ts = bin_ts1.copy()
        else:
            bin_ts = astropy.table.vstack([bin_ts, bin_ts1])

            """for j in range(len(bin_ts1["time_bin_start"].value)):
                bin_ts.add_row({'time_bin_start': bin_ts1["time_bin_start"][j],
                     'time_bin_size': bin_ts1["time_bin_size"][j],
                     'flux': bin_ts1["flux"][j],
                     'uncertainty':bin_ts1["uncertainty"][j],
                     'original_time_index': bin_ts1["original_time_index"][j]})"""


    binned_lc = self._create_copy()
    binned_lc.timelike = {}
    for t in self.timelike:
        if t == 'time':
            binned_lc.timelike[t] = bin_ts['time_bin_start'] + 0.5 * bin_ts['time_bin_size']
            binned_lc.timelike['time_bin_start'] = bin_ts['time_bin_start']
            binned_lc.timelike['time_bin_size'] = bin_ts['time_bin_size']
        else:
            if type(bin_ts[t]) != astropy.time.core.Time:
                if type(bin_ts[t]) == astropy.table.column.MaskedColumn:
                    try:
                        binned_lc.timelike[t] = bin_ts[t].filled(np.nan).value
                    except TypeError:
                        binned_lc.timelike[t] = bin_ts[t].value
                else:
                    binned_lc.timelike[t] = bin_ts[t].value
            else:
                binned_lc.timelike[t] = bin_ts[t]

    binned_lc._set_name(binned_lc.name + "_bin")
    return binned_lc

def split_time(self, split=0.5 * u.d):
    t0 = self.time[0]
    prev_obs_night = 0
    obs_nights, obs_nights_indexes = [], []

    for i, t in enumerate(self.time[1:]):
        if (t - t0) >= split:
            obs_nights.append(self.time[prev_obs_night:i + 1])
            obs_nights_indexes.append(i + 1)
            prev_obs_night = i + 1
        t0 = t

    if len(obs_nights) == 0:
        obs_nights = [self.time]

    obs_nights_indexes.insert(0, 0)
    obs_nights_indexes.append(len(self.time))

    return obs_nights_indexes, obs_nights

