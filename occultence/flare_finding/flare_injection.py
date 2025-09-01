from ..imports import *
from .flare_model import generate_fake_flare_distribution, generate_range, flare_model


def inject_flares(self, nfake, amp_range=[1e-4, 2], fwhm_range=[0.005, 0.012], mode='uniform', model='davenport2014'):
    flare_lc = self._create_copy()

    fwhms, amps = generate_fake_flare_distribution(nfake, ampl=amp_range, dur=fwhm_range, mode=mode)

    count = 0
    observed = []

    while count < nfake:
        t0s = generate_range(10*nfake, (self.time[0].value, self.time[-1].value))

        for t0 in t0s:
            mindiff = np.min(np.abs(self.time.value - t0))
            if mindiff < self.dt.value:
                observed.append(t0)
                count += 1
            if count == nfake:
                break

    flare_lc.metadata['flares_injected'] = {'nflares_injected': nfake, 'flare_params': {'t0': [],'fwhm': [],'amp': []}}
    for t0, fwhm, a in zip(observed, fwhms, amps):
        fm = flare_model(model, self.time.value, t0, fwhm, a)
        fm[np.isnan(fm)] = 0

        flare_lc.timelike['flux'] = flare_lc.timelike['flux'] + fm
        flare_lc.metadata['flares_injected']['flare_params']['t0'].append(t0)
        flare_lc.metadata['flares_injected']['flare_params']['fwhm'].append(fwhm)
        flare_lc.metadata['flares_injected']['flare_params']['amp'].append(a)

    flare_lc._set_name(flare_lc.name + "_injectedflares")

    return flare_lc



