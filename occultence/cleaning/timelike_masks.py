from ..imports import *

operator_dict = { ">": operator.gt,
                  "<": operator.lt,
                  ">=": operator.ge,
                  "<=": operator.le,
                  "==": operator.eq,
                  "!=": operator.ne}

def mask_timelike_threshold(self, timelike_key, threshold, op, verbose=False):
    """
    Create a mask for data above a certain threshold

    :param self: LightCurve object
    :param timelike_key: The name of the timelike array to mask
    :param threshold: The threshold value above (or equal) which to mask
    :return:
    """

    mask = np.zeros(self.ntime)
    # mask[operator_dict(self.timelike[timelike_key])(threshold)] = 1
    mask[operator_dict[op](self.timelike[timelike_key],threshold)] = 1

    self.masks[timelike_key] = mask

    if verbose:
        print(f"""{(100 * np.divide(float(np.count_nonzero(self.masks[timelike_key])),
                len(self.masks[timelike_key]))):.2f}% of data is flagged as {op} the threshold for {timelike_key} of 
                {threshold}
            """)