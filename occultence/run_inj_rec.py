from occultence import *
import glob
import pickle as pkl


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('i')
    args = parser.parse_args()

    clean_kw = {'dust_removal': False, 'bad_weather_removal': False, 'cosmics_removal': True, 'cosmic_boxsize': 0.08,
                'cosmic_nsigma': 3,
                'threshold_removal': True, 'thresholds': {'airmass': 2.0, 'fwhm': 6, 'bkg': 1000},
                'threshold_operators': {"airmass": ">", "fwhm": ">", "bkg": ">"}}
    gp_kw = {'do_first_sigma_clip': True, 'do_second_sigma_clip': True,
             'running_mean_boxsize': 0.08, 'nsigma': 3, 'plot': False}
    bls_kw = {"minimum_period": 0.5, "maximum_period": 10,
              'transit_durations': np.linspace(0.01, 0.1, 10), 'plot': False}
    bls_bin = 7.5 * u.minute
    recovery_kw = {'condition_on_epoch': 1 * u.hour}
    gp_bin = 20 * u.minute
    plot = False
    verbose = False

    dirname = "/Users/catrionamurray/Library/CloudStorage/OneDrive-UCB-O365/SPECULOOS/Sp2049+3336/injection_recovery/3000_injected_planets"
    self = pickle.load(open(f"{dirname}/lc_without_planet.pkl", 'rb'))
    lcs_with_transits = pickle.load(open(f"{dirname}/lcs_injected_planets.pkl", 'rb'))
    planets = pd.read_csv(f"{dirname}/injected_planets_df.csv")

    clean_lcs, gp_lcs, bls_lcs = [], [], []

    for i in range((int(args.i)-1)*300, int(args.i)*300):
        print(f"{i + 1}/{int(args.i) * 300}...")
        clean_targ, gp_targ, bls_targ, planets = self.single_injection_recovery(lcs_with_transits[i],
                                                                                planets,
                                                                                i,
                                                                                clean_kw,
                                                                                gp_bin,
                                                                                gp_kw,
                                                                                bls_kw,
                                                                                bls_bin,
                                                                                recovery_kw,
                                                                                plot,
                                                                                verbose)
        clean_lcs.append(clean_targ)
        gp_lcs.append(gp_targ)
        bls_lcs.append(bls_targ)

    planets[(int(args.i)-1)*300:int(args.i)*300].to_csv(f"{dirname}/injected_planets_df_{int(args.i)}.csv", index=False)
    pickle.dump(clean_lcs, open(f"{dirname}/lcs_clean_{int(args.i)}.pkl", 'wb'))
    pickle.dump(gp_lcs, open(f"{dirname}/lcs_gp_{int(args.i)}.pkl", 'wb'))
    pickle.dump(bls_lcs, open(f"{dirname}/lc_bls_{int(args.i)}.pkl", 'wb'))