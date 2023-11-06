from ..imports import *

def bin(self, dt, bin_func=np.nanmedian, split=True, **kw):
    ts = TimeSeries(self.timelike)

    if split:
        bin_ts = []
        i_split, _ = self.split_time()
        for i0, i1 in zip(i_split[:-1], i_split[1:]):
            if i0 == 0:
                bin_ts = aggregate_downsample(ts[i0:i1], time_bin_size=dt, aggregate_func=bin_func, **kw)
            else:
                bts = aggregate_downsample(ts[i0:i1], time_bin_size=dt, aggregate_func=bin_func, **kw)
                for b in bts:
                    bin_ts.add_row(b)
    else:
        bin_ts = aggregate_downsample(ts, time_bin_size=dt, aggregate_func=bin_func, **kw)

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

    obs_nights_indexes.insert(0,0)
    obs_nights_indexes.append(len(self.time))

    return obs_nights_indexes, obs_nights

