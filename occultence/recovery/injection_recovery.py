from ..imports import *

def was_injected_planet_recovered(self, condition_on_depth=None, condition_on_overlap=None, condition_on_epoch=None,
                                  condition_on_period=None):
    """
    Returns a list of booleans whether each transit injected into the light curve was recovered by BLS based on user-
    defined conditions.
    :param self:
    :param condition_on_depth: Fraction from 0-1 of injected depth to be recovered (e.g. 0.5 means the recovered depth
    must be at least 0.5*injected depth)
    :param condition_on_overlap: Fraction from 0-1 of injected depth to be recovered (e.g. 0.5 means the recovered
    depth must be at least 50% of injected depth)
    :param condition_on_epoch: Time value. The recovered epoch must be within injected epoch +/- condition_on_epoch to
    be recovered.
    :param condition_on_period: Fraction from 0-1 of injected period that needs to be recovered (e.g. 0.1 means the
    recovered period must match the injected period to within 10%)
    :return:
    """
    recovered_all_planets = []
    recovered_all_transits = []

    injected_params = self.metadata['injected_planet']
    recovered_params = self.metadata['BLS_transits_params']

    for planet in range(len(injected_params['depth'])):
        for transit in range(len(recovered_params['depth'])):
            recovered = True

            if condition_on_depth is not None:
                if recovered_params['depth'][transit] < (condition_on_depth * injected_params['depth'][planet]):
                    recovered=False

            if condition_on_overlap is not None:
                inj_transit_start = injected_params['epoch'][planet] - (0.5*injected_params['duration'][planet])
                inj_transit_end = injected_params['epoch'][planet] + (0.5 * injected_params['duration'][planet])
                overlap = min(inj_transit_end, recovered_params['epoch_end'][transit]) - \
                          max(inj_transit_start, recovered_params['epoch_start'][transit])

                if overlap < (condition_on_overlap * injected_params['duration']):
                    recovered = False

            if condition_on_epoch is not None:
                if abs((recovered_params['epoch'][transit] - injected_params['epoch'][planet]).to_value('jd')) > \
                        condition_on_epoch.to_value('d'):
                    recovered = False

            if condition_on_period is not None:
                period_sway = condition_on_period * injected_params['period'][planet]
                if (recovered_params['period'][transit] < injected_params['period'][planet] - period_sway) or \
                        (recovered_params['period'][transit] > injected_params['period'][planet] + period_sway):
                    recovered=False

            recovered_all_transits.append(recovered)
        recovered_all_planets.append(recovered_all_transits)

    return recovered_all_planets

