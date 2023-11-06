from occultence.imports import *

def mask_existing_transit(self, period, t0, duration, buffer=1 * u.hour):
    """
    Mask existing transits within the time series

    :param self:
    :param period: Period of transiting planet [time]
    :param t0: Mid-transit epoch of transiting planet [time]
    :param duration: Transit duration of planet [time]
    :param buffer: Buffer to mask either side of planet ingress/egress [time]
    :return:
    """

    i = 1
    if f"transit_{i}" not in self.masks:
        self.masks[f'transit_{i}'] = np.zeros(self.ntime)
    else:
        while f"transit_{i}" in self.masks:
            i += 1

        self.masks[f'transit_{i}'] = np.zeros(self.ntime)

    num_transits_since_t0 = int((self.time[0].value * u.d - t0)/period)

    transit_t = t0 + (num_transits_since_t0 * period)

    while transit_t < self.time[-1].value * u.d:
        self.masks[f'transit_{i}'][(self.time.value * u.d > transit_t-(0.5*duration)-buffer) & \
                                   (self.time.value * u.d < transit_t+(0.5*duration)+buffer)] = 1
        transit_t += period
