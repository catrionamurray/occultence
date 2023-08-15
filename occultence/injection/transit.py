from ..imports import *
from astropy.constants import G
from pytransit import QuadraticModel
import multiprocessing as mp

def semi_major_axis(per, M, R):
    """
    Calculate the semi-major axis for a given period, mass and radius.
    :param per: Orbital period in days [d]
    :param M: Stellar mass in solar masses [M_Sun]
    :param R: Stellar radius in solar radii [R_Sun]
    :return: a: semi-major axis
    """
    a = ((per)**2 * G * M/ (4 * math.pi**2))**(1./3) # in units of length
    a_radii = a / R # dimensionless

    return a_radii

def pytransit_model(time,
                    Rp=1 * u.R_earth,
                    ldc=[0.65, 0.28],
                    t0=0.0,
                    p=1.0 * u.d,
                    R=0.1 * u.R_sun,
                    M=0.1 * u.M_sun,
                    i=0.5 * math.pi
                    ):
    """
    Return pytransit model with given planet parameters and for a given series of times. [d]
    :param time: timeseries to generate the transit model.
    :param Rp: Planet radius in Earth radii [R_earth].
    :param ldc: Quadratic limb-darkening coefficients.
    :param t0: Epoch to inject transit [d].
    :param p: Orbital period of planet [d]
    :param R: Stellar radius [R_sun]
    :param M: Stellar mass [M_sun]
    :param i: Inclination [radians]
    :return: pytransit flux model with shape = time.
    """
    a_Rs = semi_major_axis(p,M,R).decompose()
    k = (Rp / R).decompose()
    tm = QuadraticModel()
    tm.set_data(time)
    flux = tm.evaluate(k=k, ldc=ldc, t0=time[0] + t0.to_value('d'), p=p.to_value('d'), a=a_Rs, i=i.to_value('radian'))
    return flux

def inject_transit(self, per, epoch, inc, rp, M, R, ld):
    """

    :param self:
    :param per:
    :param epoch:
    :param inc:
    :param rp:
    :param M:
    :param R:
    :param ld:
    :return:
    """
    model = pytransit_model(time=self.time.value, Rp=rp, ldc=ld, t0=epoch, p=per,R=R, M=M, i=inc)
    a_Rs = semi_major_axis(per, M, R).decompose()
    b = a_Rs * math.cos(inc.to_value('radian'))
    sini = math.sin(inc.to_value('radian'))
    k =(rp / R).decompose()
    transit_depth = k ** 2
    transit_duration = (per / math.pi) * math.asin((1 / a_Rs) * np.sqrt((1 + k) ** 2 - b ** 2) / sini)

    injected_lc = self._create_copy()
    injected_lc.timelike['flux'] = model * self.flux
    injected_lc.timelike['transit_model'] = model
    if 'injected_planet' in injected_lc.metadata:
        injected_lc.metadata['injected_planet']['period'].append(per)
        injected_lc.metadata['injected_planet']['epoch'].append(epoch)
        injected_lc.metadata['injected_planet']['inc'].append(inc)
        injected_lc.metadata['injected_planet']['rp'].append(rp)
        injected_lc.metadata['injected_planet']['ld'].append(ld)
        injected_lc.metadata['injected_planet']['depth'].append(transit_depth)
        injected_lc.metadata['injected_planet']['duration'].append(transit_duration)
    else:
        injected_lc.metadata['injected_planet'] = {'period':[per], 'epoch':[self.time.value[0]*u.d + epoch], 'inc':[inc],
                                                   'rp':[rp], 'ld':[ld], 'depth':[transit_depth],
                                                   'duration':[transit_duration]}
    return injected_lc

def pool_inject_transit(self, list_of_logpers, list_of_phase, list_of_cosi, list_of_rp, ld, ncores, R_star=None,
                        M_star=None):
    """

    :param self:
    :param list_of_logpers:
    :param list_of_phase:
    :param list_of_cosi:
    :param list_of_rp:
    :param ld:
    :param ncores:
    :param R_star:
    :param M_star:
    :return:
    """

    if R_star is None or M_star is None:
        if "R_star" in self.metadata and "M_star" in self.metadata:
            R_star = self.metadata['R_star']
            M_star = self.metadata['M_star']
        else:
            warnings.warn("To generate a planet distribution for this star, we need its radius and mass. Please store"
                          "the radius and mass in the .metadata dict or pass them explicitly to this function!")
            return

    print("Number of processors: ", mp.cpu_count(), ". Using ", str(ncores), "cores")
    pool = mp.Pool(ncores)
    lcs = pool.starmap(self.inject_transit, [(10**logp, phase * 10**logp, math.acos(cosi), rp, M_star, R_star, ld) for
                                             logp, phase, cosi, rp in zip(list_of_logpers, list_of_phase, list_of_cosi,
                                                                          list_of_rp)])

    pool.close()

    return lcs

# i_t_s = inj_transit_start
# i_t_e = inj_transit_end
# injected_transit_times = []
# while i_t_s < time[-1]:
#     injected_transit_times.append([i_t_s, i_t_e])
#     i_t_s = i_t_s + 10 ** logp
#     i_t_e = i_t_e + 10 ** logp


def generate_planet_distribution(nfake, m_s, r_s, per=[np.log10(0.5),np.log10(10)],phase=[0,1], cosi=[0,1],
                                 radius=[0.5, 6],mode_per="uniform", mode_phase="uniform",mode_cosi="uniform",
                                 mode_radius="uniform",**kwargs):
    """
    [ADAPTED FROM ALTAIPONY] Function to generate a planet distribution over period, phase, inclination and radius.
    :param nfake: Number of artificial planets to create.
    :param m_s: Stellar mass
    :param r_s: Stellar radius
    :param per: Range of orbital periods (min period, max period) [d]
    :param phase: Range of orbital phases at which to inject planets (min phase, max phase), determines the transit epoch.
    Typically this will be [0,1].
    :param cosi: Range of cos(inclinations) - THIS HAS BEEN REMOVED - replaced with range (0, 1/a)
    :param radius: Range of planet radii (min radius, max radius) [R_earth]
    :param mode_per: The type of distribution to use for period (options: uniform, normal, lognormal)
    :param mode_phase: The type of distribution to use for phase (options: uniform, normal, lognormal)
    :param mode_cosi: The type of distribution to use for cosi (options: uniform, normal, lognormal)
    :param mode_radius: The type of distribution to use for radius (options: uniform, normal, lognormal)
    :param kwargs:
    :return:
    """

    def generate_range(n, tup, **kwargs):
        x = (mod_random(n, **kwargs) * (tup[1] - tup[0])) + tup[0]
        return x

    def generate_normal_range(n,tup,**kwargs):
        x = np.random.normal(tup[0], (tup[1] - tup[0]) / 3., 1000)

        for i in range(n):
            j = x[i]
            while j<tup[0] or j>tup[1]:
                j = np.random.normal(tup[0], (tup[1] - tup[0]) / 2., size=1)
            x[i]=j

        return x

    def generate_lognormal_range(n, tup, **kwargs):
        sigma = 15
        mean=0
        x=(np.random.lognormal(mean=mean, sigma=sigma, size=n))

        for i in range(n):
            j = x[i]
            while j<tup[0] or j>tup[1]:
                j = np.random.lognormal(mean=mean, sigma=sigma, size=1)
            x[i]=j

        return x

    ms = []

    # loop over period, phase and radius and generate distributions:
    for z,mode in zip([per,phase,radius],[mode_per,mode_phase,mode_radius]):

        if mode == 'uniform':
            m = generate_range(nfake, z, **kwargs)

        if mode == 'lognormal':
            m = generate_lognormal_range(nfake, z, **kwargs)

        if mode == 'normal':
            m = generate_normal_range(nfake, z, **kwargs)

        ms.append(m)

    ms.append(np.zeros(len(ms[0])))
    mode = mode_cosi
    m_temp = []

    for mi in ms[0]:
        # calculate the semi-major axis to determine the inclination limits.
        ascale = semi_major_axis(10**mi,m_s,r_s)
        cosi_ind = [0,1/ascale]

        if mode == 'uniform':
            m = generate_range(1, cosi_ind, **kwargs)

        if mode == 'lognormal':
            m = generate_lognormal_range(1, cosi_ind, **kwargs)

        if mode == 'normal':
            m = generate_normal_range(1, cosi_ind, **kwargs)

        m_temp.append(float(m))
    ms[3] = np.array(m_temp)

    return ms

def create_lots_of_transit_params(self, nfake=1000, R_star=None, M_star=None, T_eff=None, SpT=None,
                                 minimum_planet_radius=0.5 * u.R_earth, maximum_planet_radius=3 * u.R_earth,
                                 minimum_period=0.5 * u.d, maximum_period=30 * u.d, store_planets=True, fname=None,):

    if R_star is None or M_star is None:
        if "R_star" in self.metadata and "M_star" in self.metadata:
            R_star = self.metadata['R_star']
            M_star = self.metadata['M_star']
        else:
            warnings.warn("To generate a planet distribution for this star, we need its radius and mass. Please store"
                          "the radius and mass in the .metadata dict or pass them explicitly to this function!")
            return
    else:
        self.metadata['R_star'] = R_star
        self.metadata['M_star'] = M_star

    if T_eff is None and "T_eff" in self.metadata:
        T_eff = self.metadata['T_eff']

    if SpT is None and "SpT" in self.metadata:
        T_eff = self.metadata['SpT']

    if store_planets:
        params = generate_planet_distribution(nfake, M_star, R_star, radius=[minimum_planet_radius, maximum_planet_radius],
                                              per=[minimum_period, maximum_period])
        params[0] = np.log10(params[0])
        planets = pd.DataFrame({'logP': params[0], 'phase': params[1], 'cosi': params[3], 'r_p': params[2],
                                'depth': np.zeros(len(params[0])), 'duration': np.zeros(len(params[0])),
                                'epoch': np.zeros(len(params[0])),
                                'recovered': np.zeros(len(params[0])), 'log_Prec': np.zeros(len(params[0])),
                                'rec_depth': np.zeros(len(params[0])), 'rec_duration': np.zeros(len(params[0])),
                                'rec_epoch': np.zeros(len(params[0])), 'run': np.zeros(len(params[0])),
                                'snr': np.zeros(len(params[0])), 'target': [self.name] * len(params[0]),
                                'r_s': [R_star] * len(params[0]), 'm_s': [M_star] * len(params[0]),
                                'teff': [T_eff] * len(params[0]),
                                'spt': [SpT] * len(params[0])})

        planets.to_csv(fname, index=False)
    else:
        if os.path.exists(fname):
            planets = pd.read_csv(fname)
        else:
            message = f""" Watch out! The filename given {fname} does not exist. To create a new planet file for this
                        target pass `store_models=True` to this function!
                                """
            cheerfully_suggest(message)
            return

    return planets

def inject_lots_of_transits(self, nfake=1000, R_star=None, M_star=None, T_eff=None, SpT=None, ld = [0.385, 0.304],
                            minimum_planet_radius=0.5 * u.R_earth, maximum_planet_radius=3 * u.R_earth,
                            minimum_period=0.5 * u.d, maximum_period=30 * u.d, store_planets=True, fname=None,
                            pool=True, ncores=1):

    planets = self.create_lots_of_transit_params(nfake=nfake, R_star=R_star, M_star=M_star, T_eff=T_eff, SpT=SpT,
                               minimum_planet_radius=minimum_planet_radius, maximum_planet_radius=maximum_planet_radius,
                               minimum_period=minimum_period, maximum_period=maximum_period,
                               store_planets=store_planets, fname=fname)

    if pool:
        lcs_with_transits = self.pool_inject_transit(planets['logP'], planets['phase'], planets['cosi'], planets['r_p'],
                                ld, ncores, R_star=R_star, M_star=M_star)

    return lcs_with_transits