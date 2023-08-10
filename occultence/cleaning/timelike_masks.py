from ..imports import *

def mask_timelike_threshold(self, timelike_key, threshold, op):
    """
    Create a mask for data above a certain threshold

    :param self: LightCurve object
    :param timelike_key: The name of the timelike array to mask
    :param threshold: The threshold value above (or equal) which to mask
    :return:
    """

    mask = np.zeros(self.ntime)
    mask[operator_dict(self.timelike[timelike_key])(threshold)] = 1

    self.masks[timelike_key] = mask

    print(f"""{100 * np.divide(float(np.count_nonzero(self.masks['timelike_key'])),
                len(self.masks['timelike_key']))}% of data is flagged as {op} the threshold for {timelike_key} of 
                {threshold}"
            """)