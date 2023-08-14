from .imports import *

def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return idx

def running_box(x,y,boxsize,operation):

    if operation == "std":
        op = np.nanstd
    elif operation == 'median':
        op = np.nanmedian
    elif operation == 'mean':
        op = np.nanmean
    elif operation == "clipped_std":
        op = clipped_std
    else:
        print("No running function selected - choose std, median or mean.")
        return np.nan

    l = len(x)
    boxsize = boxsize / 2
    dy = np.zeros((np.size(x)))

    for i in range(0, l):
        s = x[i]
        box_s = find_nearest(x, s - boxsize)
        box_e = find_nearest(x, s + boxsize)

        # calculate [OPERATION] for the current window
        try:
            y_ext = y[box_s:box_e + 1]
            d = op(y_ext)
            if not np.isnan(d) and len(y_ext)>=3:
                dy[i] = d

        except Exception as e:
            print(e)

    dy = np.ma.masked_where(dy == 0, dy)
    dy = np.ma.masked_invalid(dy)
    return dy

def clipped_std(x, sigma=3):
    from astropy.stats import sigma_clip
    return np.ma.std(sigma_clip(x,sigma=sigma))

def clipped_mean(x, sigma=3):
    from astropy.stats import sigma_clip
    return np.ma.mean(sigma_clip(x,sigma=sigma))

def calculate_running_rms(x,y,boxsize):
    dy = running_box(x,y,boxsize,'std')
    return dy

def calculate_running_clipped_rms(x,y,boxsize):
    dy = running_box(x,y,boxsize,'clipped_std')
    return dy

def calculate_running_median(x,y,boxsize):
    dy = running_box(x,y,boxsize,'median')
    return dy

def calculate_running_mean(x,y,boxsize):
    dy = running_box(x,y,boxsize,'mean')
    return dy

def sort_on_time(time,*sortarrs):
    sorted_arrs = []
    for arr in sortarrs:
        sorted_arrs.append(np.array([j for i, j in sorted(zip(time, arr))]))
    sorted_time = np.array(sorted(time))
    return sorted_time, (*sorted_arrs,)

