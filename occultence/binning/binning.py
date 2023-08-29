from ..imports import *

def bin(self, dt, bin_func=np.nanmedian, **kw):
    ts = TimeSeries(self.timelike)
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