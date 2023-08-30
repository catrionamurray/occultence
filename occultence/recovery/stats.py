from ..imports import *

def probability_of_finding_a_planet(occurrence_rate, sensitivity, r_p, a, R_star):
    probability_of_transit = ((R_star + r_p) / a).decompose()
    probability_of_detection = sensitivity
    probability_of_planet = occurrence_rate
    completeness = probability_of_transit * probability_of_detection
    prob = completeness * probability_of_planet
    return prob

def expected_number_of_planets(list_of_lightcurves, occurrence_rate, sensitivity, r_p, a):
    p = 0
    for lc in list_of_lightcurves:
        p += probability_of_finding_a_planet(occurrence_rate, sensitivity, r_p, a, lc.metadata['R_star'])
    return p

def calculate_completeness(x,y,z,x_bins,y_bins):
    # twod_plot(x, y, z_randomalignment, x_bins, y_bins, "Completeness", statistic='mean', nsig='%0.3f')
    ret = binned_statistic_2d(x, y, z, statistic="mean", bins=[x_bins, y_bins])
    completeness = ret.statistic.T
    completeness[np.isnan(completeness)] = 0
    return completeness



def something(R_star, R_planet, ):
    x = [10 ** i for i in temp_planets['logP'].values]
    y = temp_planets['r_p'].values
    z = temp_planets['recovered'].values
    sensitivity_i.append(calculate_completeness(x, y, z, x_bins, y_bins))
    z_randomalignment = z * (
                            (temp_planets['r_s'].values * sunrad) + (temp_planets['r_p'].values * earthrad)) / \
                                        temp_planets['a'].values
    completeness = calculate_completeness(x, y, z_randomalignment, x_bins, y_bins)
    H_i.append(completeness)
    if len(expected_planets) > 0:
        expected_planets = expected_planets + completeness
        expected_planets_r = expected_planets_r + (completeness * r)
    else:
        expected_planets = completeness
        expected_planets_r = completeness * r

    if len(no_detect_100) > 0:
        no_detect_100 = no_detect_100 * (1 - completeness)
        no_detect_50 = no_detect_50 * (1 - (0.5 * completeness))
        no_detect_10 = no_detect_10 * (1 - (0.1 * completeness))
    else:
        no_detect_100 = (1 - completeness)
        no_detect_50 = (1 - (0.5 * completeness))
        no_detect_10 = (1 - (0.1 * completeness))