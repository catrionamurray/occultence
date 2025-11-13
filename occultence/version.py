__version__ = "0.5.15"


def version():
    return __version__


# v0.2.1 - added the ability to mask existing transits in the timeseries
# v0.2.2 - added split_time
# v0.3.0 - added mcmc and lsq
# v0.3.2 - added celerite and emcee and corner as requirements
# v0.4.0 - fixed multiprocessing pooling
# v0.5.0 - added rotation + flare capabilities
# v0.5.1 - fixed bug in ...
# v0.5.2 - fixed bug in bad weather removal
# v0.5.3 - fixed bugs in flare detection and dust removal, small fixes
# v0.5.5 - add nightly normalization
# v0.5.7 - updated recovery criteria
# v0.5.8 - sped up binning
# v0.5.9 - skip running the full IR if the planet is not observed
# v0.5.10 - slightly speed up BLS
# v0.5.11 - add the ability to save plots
# v0.5.12 - add the ability to save LCs
# v0.5.13 - change BLS from .autopower to .power
# v0.5.14 - added "min_save" kw