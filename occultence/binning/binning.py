from ..imports import *


# def bin(self, dt, bin_func=np.nanmedian, **kw):
#     # cuts = np.where(time_values[1:] - time_values[:-1] > 0.5)[0]
#     cuts, _ = self.split_time(split=self.split_by)
#     # cuts = np.hstack((0, cuts + 1, len(time_values)))
#
#     # Pre-allocate list to collect binned results (much faster than vstack in loop)
#     bin_results = []
#
#     for i in range(len(cuts) - 1):
#
#         # Create view of sliced data without copying the full dictionary
#         start_idx, end_idx = cuts[i], cuts[i + 1]
#
#         # Build sliced timelike dict directly
#         t1 = {key: self.timelike[key][start_idx:end_idx] for key in self.timelike.keys()}
#         t1s = TimeSeries(t1)
#
#         bin_ts1 = aggregate_downsample(t1s, time_bin_size=dt, aggregate_func=bin_func, **kw)
#         bin_results.append(bin_ts1)
#
#     # Single vstack operation instead of repeated vstacks
#     bin_ts = astropy.table.vstack(bin_results) if len(bin_results) > 1 else bin_results[0]
#
#     binned_lc = self._create_copy()
#     binned_lc.timelike = {}
#
#     # Process time columns
#     binned_lc.timelike['time'] = bin_ts['time_bin_start'] + 0.5 * bin_ts['time_bin_size']
#     binned_lc.timelike['time_bin_start'] = bin_ts['time_bin_start']
#     binned_lc.timelike['time_bin_size'] = bin_ts['time_bin_size']
#
#     for t in self.timelike:
#         if t == 'time':
#             continue
#
#         col = bin_ts[t]
#
#         # Avoid redundant type checks - check once
#         if isinstance(col, astropy.time.core.Time):
#             binned_lc.timelike[t] = col
#         elif isinstance(col, astropy.table.column.MaskedColumn):
#             try:
#                 binned_lc.timelike[t] = col.filled(np.nan).value
#             except TypeError:
#                 binned_lc.timelike[t] = col.value
#         else:
#             binned_lc.timelike[t] = col.value
#
#     binned_lc._set_name(binned_lc.name + "_bin")
#     return binned_lc


def bin(self, dt, bin_func=np.nanmedian, **kw):
    """Ultra-fast vectorized binning with unit support."""
    cuts, _ = self.split_time_indices_only(split=self.split_by)

    binned_lc = self._create_copy()
    binned_lc.timelike = {}

    all_results = []
    time_type_info = None  # Store time type info from first segment

    for i in range(len(cuts) - 1):
        start_idx, end_idx = cuts[i], cuts[i + 1]

        # Extract time values
        time_data = self.timelike['time'][start_idx:end_idx]

        # Handle units and extract numeric values
        if isinstance(time_data, astropy.time.core.Time):
            time_vals = time_data.value
            if i == 0:
                time_type_info = ('time', time_data.format, time_data.scale)
            dt_val = dt.to(astropy.units.day).value if hasattr(dt, 'to') else dt
        elif hasattr(time_data, 'unit'):
            time_vals = time_data.value
            if i == 0:
                time_type_info = ('quantity', time_data.unit)
            dt_val = dt.to(time_data.unit).value if hasattr(dt, 'to') else dt
        else:
            time_vals = time_data
            if i == 0:
                time_type_info = ('plain', None)
            dt_val = dt if not hasattr(dt, 'value') else dt.value

        # Create bins
        t_min, t_max = np.nanmin(time_vals), np.nanmax(time_vals)
        bin_edges = np.arange(t_min, t_max + dt_val, dt_val)
        bin_indices = np.digitize(time_vals, bin_edges) - 1

        # Create temporary dataframe with numeric values only
        df_dict = {'bin': bin_indices}
        col_units = {}  # Track units for each column

        for key in self.timelike.keys():
            if key == 'time':
                continue
            data = self.timelike[key][start_idx:end_idx]
            if isinstance(data, astropy.time.core.Time):
                continue  # Skip Time objects

            if hasattr(data, 'unit'):
                col_units[key] = data.unit
                data_vals = data.value
            elif isinstance(data, astropy.table.column.MaskedColumn):
                data_vals = data.filled(np.nan)
                if hasattr(data, 'value'):
                    data_vals = data_vals if not hasattr(data_vals, 'value') else data_vals.value
            else:
                data_vals = data.value if hasattr(data, 'value') else data

            df_dict[key] = data_vals

        df = pd.DataFrame(df_dict)

        # Group and aggregate
        if bin_func == np.nanmedian:
            grouped = df.groupby('bin').median()
        elif bin_func == np.nanmean or bin_func == np.mean:
            grouped = df.groupby('bin').mean()
        else:
            grouped = df.groupby('bin').agg(lambda x: bin_func(x.values))

        # Extract results
        valid_bins = grouped.index.values
        valid_bins = valid_bins[(valid_bins >= 0) & (valid_bins < len(bin_edges) - 1)]

        result = {
            'time_bin_start': bin_edges[valid_bins],
            'time_bin_size': np.full(len(valid_bins), dt_val),
            'time': bin_edges[valid_bins] + 0.5 * dt_val,
            '_col_units': col_units  # Store units metadata
        }

        for key in grouped.columns:
            result[key] = grouped.loc[valid_bins, key].values

        all_results.append(result)

    # Concatenate all segments and restore types
    col_units = all_results[0].get('_col_units', {})

    for key in ['time_bin_start', 'time_bin_size', 'time']:
        concat_vals = np.concatenate([r[key] for r in all_results])

        # Restore appropriate type
        if time_type_info[0] == 'time':
            if key == 'time_bin_size':
                binned_lc.timelike[key] = concat_vals * astropy.units.day
            else:
                binned_lc.timelike[key] = astropy.time.Time(concat_vals, format=self.time.format)
        elif time_type_info[0] == 'quantity':
            binned_lc.timelike[key] = concat_vals * time_type_info[1]
        else:
            binned_lc.timelike[key] = concat_vals

    # Concatenate data columns with units
    for key in all_results[0].keys():
        if key in ['time_bin_start', 'time_bin_size', 'time', '_col_units']:
            continue

        concat_vals = np.concatenate([r[key] for r in all_results])

        if key in col_units:
            binned_lc.timelike[key] = concat_vals * col_units[key]
        else:
            binned_lc.timelike[key] = concat_vals

    binned_lc._set_name(binned_lc.name + "_bin")
    return binned_lc

# def split_time(self, split=0.5 * u.d):
#     t0 = self.time[0]
#     prev_obs_night = 0
#     obs_nights, obs_nights_indexes = [], []
#
#     for i, t in enumerate(self.time[1:]):
#         if (t - t0) >= split:
#             obs_nights.append(self.time[prev_obs_night:i + 1])
#             obs_nights_indexes.append(i + 1)
#             prev_obs_night = i + 1
#         t0 = t
#
#     if len(obs_nights) == 0:
#         obs_nights = [self.time]
#
#     obs_nights_indexes.insert(0, 0)
#     obs_nights_indexes.append(len(self.time))
#
#     return obs_nights_indexes, obs_nights

def split_time(self, split=0.5 * u.d):
    """Vectorized version - typically 100-1000x faster."""

    # Handle units
    time_data = self.time
    if hasattr(split, 'to'):
        # split has units
        if isinstance(time_data, astropy.time.core.Time):
            # Convert Time to numeric (days)
            time_vals = time_data.value
            split_val = split.to(u.day).value
        elif hasattr(time_data, 'unit'):
            # time is a Quantity
            time_vals = time_data.value
            split_val = split.to(time_data.unit).value
        else:
            # time is plain array
            time_vals = time_data
            split_val = split.value
    else:
        # split is plain number
        if isinstance(time_data, astropy.time.core.Time):
            time_vals = time_data.value
        elif hasattr(time_data, 'value'):
            time_vals = time_data.value
        else:
            time_vals = time_data
        split_val = split

    # Vectorized gap detection: compute all differences at once
    time_diffs = np.diff(time_vals)

    # Find where gaps exceed split threshold
    gap_indices = np.where(time_diffs >= split_val)[0]

    # Build split indices: [0, gap1+1, gap2+1, ..., len(time)]
    obs_nights_indexes = np.concatenate([[0], gap_indices + 1, [len(time_vals)]])

    # Build obs_nights list by slicing
    obs_nights = []
    for i in range(len(obs_nights_indexes) - 1):
        start, end = obs_nights_indexes[i], obs_nights_indexes[i + 1]
        obs_nights.append(self.time[start:end])

    return obs_nights_indexes.tolist(), obs_nights


def split_time_indices_only(self, split=0.5 * u.d):
    """Even faster version if you only need indices (not the sliced arrays)."""

    # Handle units
    time_data = self.time
    if hasattr(split, 'to'):
        if isinstance(time_data, astropy.time.core.Time):
            time_vals = time_data.value
            split_val = split.to(u.day).value
        elif hasattr(time_data, 'unit'):
            time_vals = time_data.value
            split_val = split.to(time_data.unit).value
        else:
            time_vals = time_data
            split_val = split.value
    else:
        if isinstance(time_data, astropy.time.core.Time):
            time_vals = time_data.value
        elif hasattr(time_data, 'value'):
            time_vals = time_data.value
        else:
            time_vals = time_data
        split_val = split

    # Vectorized gap detection
    time_diffs = np.diff(time_vals)
    gap_indices = np.where(time_diffs >= split_val)[0]
    obs_nights_indexes = np.concatenate([[0], gap_indices + 1, [len(time_vals)]])

    # Return indices only (much faster if you don't need the sliced arrays)
    return obs_nights_indexes.tolist(), None


def extract(self, ind):
    new_lc = self._create_copy()
    for k, v in self.timelike.items():
        new_lc.timelike[k] = v[ind]
    if new_lc.telescope is not None:
        new_lc.telescope = new_lc.telescope[ind]
    if new_lc.filter is not None:
        new_lc.filter = new_lc.filter[ind]
    return new_lc
