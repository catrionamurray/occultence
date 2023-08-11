from ..imports import *
from ..lightcurve import *
from ..utils import sort_on_time


def import_diff_fits(fname: str,
                     ap: int,
                     tel: str,
                     targ: str,
                     # date: str or int,
                     filt: str,
                     pwv: bool = True,
                     check_saturation: bool = True,
                     outputfits_filename: str = "",
                     ):
    """

    :param fname: The diff.fits filename. This can either be a night or global diff.fits file.
    :param ap: The aperture integer corresponding to that diff.fits file.
    :param tel: The telescope name.
    :param targ: The target name.
    :param filt: The filter name.
    :param pwv: Boolean whether to extract the PWV-corrected flux (Default = True)
    :param check_saturation: Boolean whether to check each time point for saturation. (Default = True)
    :param outputfits_filename: If check_saturation==True then we need the output.fits filename to extract the raw flux.
    (Default = "")
    :return: LightCurve object
    """

    if not os.path.exists(fname):
        warnings.warn(f"Warning! The diff.fits filename provided ({fname}) does not seem to exist.")
        return

    # extract all useful information from the diff.fits file!
    with fitsio.FITS(fname) as infile:
        lcs = infile["LIGHTCURVE_" + str(ap)].read()
        pwv_lcs = infile["PWV_LIGHTCURVE_" + str(ap)].read()
        bwmask = infile['FLAGS']['BW_FLAG'].read()
        hdr = infile[0].read_header()
        imagelist = infile['imagelist']
        t = imagelist['jd-obs'].read()
        cat = infile['catalogue']
        alc_table = infile['ALC_' + str(ap)]
        alc = alc_table['alc'].read()
        teffs = cat['teff'].read()
        exp = imagelist['exposure'].read()
        ramove = imagelist['ra_move'].read()
        decmove = imagelist['ra_move'].read()
        ccd_temp = imagelist['ccd-temp'].read()
        airmass = imagelist['airmass'].read()
        bkg = imagelist['skylevel'].read()
        fwhm = imagelist['fwhm'].read()
        focus = imagelist['FOCUSPOS'].read()
        targ_teff = teffs[0]

        try:
            izmags = cat['izmag'].read()
            targ_mag = izmags[0]
        except:
            targ_mag = np.nan

        targ_gaia_id = hdr['GAIA_ID']

    # extract target lightcurve
    if pwv:
        f = pwv_lcs[:, 0]
    else:
        f = lcs[:, 0]

    # If we are to check the saturation then check that the output fits filepath is provided and exists.
    if check_saturation:
        if outputfits_filename != "":
            if os.path.exists(outputfits_filename):
                # Extract the peak flux
                peak, rawflux = get_peak_flux(outputfits_filename, targ_gaia_id, ap)
            else:
                message = f""" Warning! The outputfits file path provided ({outputfits_filename}) does not seem to exist """
                cheerfully_suggest(message)
                peak, rawflux = np.zeros(len(t)), np.zeros(len(t))
        else:
            message = f""" Warning! The outputfits file path has not been provided."""
            cheerfully_suggest(message)
            peak, rawflux = np.zeros(len(t)), np.zeros(len(t))

    assert len(peak) == len(t)  # make sure the peak flux length matches the length of the extracted time

    t, f, alc, exp, ccd_temp, airmass, bkg, fwhm, focus, ramove, decmove, lf, peak = sort_on_time(t, f, alc, exp,
                                                                                                  ccd_temp, airmass,
                                                                                                  bkg, fwhm, focus,
                                                                                                  ramove, decmove, peak)

    metadata = {'telescope': tel, 'filter': filt, 'dates': dates, 'gaia_dr2': targ_gaia_id, 't_eff': targ_teff,
                'I+z_mag': targ_mag}
    timelike = {'airmass': airmass, 'background': bkg, 'fwhm': fwhm, 'ramove': ramove, 'decmove': decmove,
                'focus': focus, 'artificial_lightcurve': alc, 'ccd_temperature': ccd_temp, 'bad_weather': bwmask,
                'peak_flux': peak}

    lc = LightCurve(name=targ, time=t, flux=f, uncertainty=e, timelike=timelike, metadata=metadata)
    return lc


def get_peak_flux(outfits, targ_gaia, ap):
    """
    import the outfits file and extract the 'peak' flux (in an individual pixel within the aperture) and the average
    flux in the aperture ('FLUX') for the target
    :param outfits: Output.fits filepath
    :param targ_gaia: The Gaia DR2 name of the target
    :param ap: The aperture number to extract
    :return: peak flux, average (non-normalized) flux
    """

    with fitsio.FITS(outfits) as infile:
        cat = infile['catalogue']
        gaia_id = cat['gaia_dr2_id'].read()
        peak = infile['peak'].read()
        flux = infile["FLUX_" + str(int(ap))].read()

    i = np.where(np.array(gaia_id) == targ_gaia)[0][0]
    return peak[i], flux[i]
