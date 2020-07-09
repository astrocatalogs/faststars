"""Import tasks for LAMOST."""
import os
from decimal import Decimal

import astropy.constants as con
import astropy.coordinates as coord
import astropy.units as un
from astrocats.catalog.spectrum import SPECTRUM
from astrocats.catalog.utils import jd_to_mjd, pbar, pretty_num
from astropy.io import fits
from astropy.time import Time as astrotime
from astroquery.sdss import SDSS
import urllib

from ..faststars import FASTSTARS


def do_sdss_spectra(catalog):
    """Import spectra from LAMOST."""
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
    keys = list(catalog.entries.keys())

    fureps = {'erg/cm2/s/A': 'erg/s/cm^2/Angstrom',
              '1E-17 erg/cm^2/s/Ang': 'erg/s/cm^2/Angstrom'}

    c_kms = con.c.cgs.value / 1.0e5
    cntsdss = 0
    for oname in pbar(keys, task_str):
        # Some events may be merged in cleanup process, skip them if
        # non-existent.

        if (FASTSTARS.RA not in catalog.entries[oname] or
                FASTSTARS.DEC not in catalog.entries[oname]):
            continue
        else:
            xid = SDSS.query_region(coord.SkyCoord(
                ra=catalog.entries[oname][FASTSTARS.RA][0]['value'],
                dec=catalog.entries[oname][FASTSTARS.DEC][0]['value'],
                unit=(un.hourangle, un.deg), frame='icrs'), spectro=True)
            # xid = SDSS.query_region(coord.SkyCoord(
            #        ra='14:34:06.17',
            #        dec='+56:30:47.24',
            #        unit=(un.hourangle, un.deg), frame='icrs'), spectro=True)
            if xid is None:
                continue
            while len(xid) > 1:
                notstar = xid['z'].argmax()
                xid.remove_row(notstar)
            #print(xid)

            # star = None
            # for row in tab:
            #    if (row['objType'] == 'Star' and
            #            row['Class'].lower() in ['star', 'unknown']):
            #        star = row
            #        break
            # if not star:
            #    continue

            try:
                name, source = catalog.new_entry(
                    oname, bibcode='2015ApJS..219...12A', srcname='SDSS',
                    url='http://www.sdss.org/')
            except Exception:
                catalog.log.warning(
                    '"{}" was not found, suggests merge occurred in cleanup '
                    'process.'.format(oname))
                continue

            ffile = ('spec-' + str(xid['specobjid'][0]) + '.fits.gz')

            # furl = 'http://dr3.lamost.org/sas/fits/' + vname + '/' + ffile

            datafile = os.path.join(
                catalog.get_current_task_repo(), 'SDSS', ffile)

            if not os.path.exists(datafile):
                # Download spectra
                try:
                    sp = SDSS.get_spectra(matches=xid)[0]
                except urllib.error.HTTPError:
                    catalog.log.warning(
                    '"{}" threw an HTTP 404, must be error upstream. Will likely go away on the next run.'.format(oname))
                continue

                # Identify star
                # assert sp[2].data['class'][0]=='STAR'

                # Write spectra
                # print(catalog.entries[oname][FASTSTARS.RA][0]['value'],catalog.entries[oname][FASTSTARS.DEC][0]['value'])
                sp.writeto(datafile, overwrite=True)
                # open(datafile, 'wb').write(fr.content)

            # Open spectra
            hdulist = fits.open(datafile)

            # sp contains a list of fits datafiles, identify main one
            i_primary = 0
            i_coadd = 1
            i_specobj = 2
            assert hdulist[i_primary].name == 'PRIMARY'
            assert hdulist[i_coadd].name == 'COADD'
            assert (hdulist[i_specobj].name ==
                    'SPECOBJ' or hdulist[i_specobj].name == 'SPALL')

            # xid = SDSS.query_region(coord.SkyCoord(
            #     ra='12:11:50.27',
            #     dec='+14:37:16.2',
            #     unit=(un.hourangle, un.deg), frame='icrs'), spectro=True)

            # from SPECOBJ
            # print('.'+hdulist[i_specobj].data['ELODIE_SPTYPE'][0]+'.')
            if hdulist[i_specobj].data['ELODIE_SPTYPE'][0] != 'unknown':
                ST, SCfull = hdulist[i_specobj].data[
                    'ELODIE_SPTYPE'][0][:2], hdulist[
                        i_specobj].data['ELODIE_SPTYPE'][0][2:]
                if len(SCfull) > 0:
                    if 'IV' in SCfull:
                        SC = 'sg'
                    elif 'III' in SCfull:
                        SC = 'g'
                    elif 'V' in SCfull:
                        SC = 'd'
                    elif 'I' in SCfull:
                        SC = 'Sg'
                    else:
                        SC = False
                    if SC is not False:
                        catalog.entries[name].add_quantity(
                            FASTSTARS.STELLAR_CLASS, SC, source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.SPECTRAL_TYPE, ST, source=source)

            if hdulist[i_specobj].data['Z'][0] != 0.0:
                catalog.entries[name].add_quantity(
                    FASTSTARS.REDSHIFT, str(hdulist[i_specobj].data['Z'][0]),
                    e_value=str(hdulist[i_specobj].data['Z_ERR'][0]),
                    source=source)
                catalog.entries[name].add_quantity(
                    FASTSTARS.VELOCITY,
                    pretty_num(
                        float(hdulist[i_specobj].data['Z'][0]) * c_kms, sig=5),
                    e_value=pretty_num(
                        float(hdulist[i_specobj].data['Z_ERR'][0] * c_kms),
                        sig=5),
                    source=source)

            for oi, obj in enumerate(hdulist[0].header):
                if any(x in ['.', '/'] for x in obj):
                    del (hdulist[0].header[oi])
            hdulist[0].verify('silentfix')
            hdrkeys = list(hdulist[0].header.keys())
            # print(hdrkeys)
            # for key in hdulist[0].header.keys():
            #     print(key, hdulist[0].header[key])
            if hdulist[0].header['SIMPLE']:
                if 'JD' in hdrkeys:
                    mjd = str(jd_to_mjd(Decimal(str(hdulist[0].header['JD']))))
                elif 'MJD' in hdrkeys:
                    mjd = str(hdulist[0].header['MJD'])
                elif 'DATE-OBS' in hdrkeys:
                    if 'T' in hdulist[0].header['DATE-OBS']:
                        dateobs = hdulist[0].header['DATE-OBS'].strip()
                    elif 'UTC-OBS' in hdrkeys:
                        dateobs = hdulist[0].header['DATE-OBS'].strip(
                        ) + 'T' + hdulist[0].header['UTC-OBS'].strip()
                    mjd = str(astrotime(dateobs, format='isot').mjd)
                else:
                    raise ValueError("Couldn't find JD/MJD for spectrum.")
                if hdulist[i_coadd].header['NAXIS'] == 2:
                    waves = [str(x)
                             for x in list(10**hdulist[i_coadd].data['loglam'])]
                    fluxes = [str(x)
                              for x in list(hdulist[i_coadd].data['flux'])]
                else:
                    print('Warning: Skipping FITS spectrum `{}`.'.format(
                        datafile))
                    continue
            else:
                raise ValueError('Non-simple FITS import not yet supported.')
            if 'BUNIT' in hdrkeys:
                fluxunit = hdulist[0].header['BUNIT']
                if fluxunit in fureps:
                    fluxunit = fureps[fluxunit]
                if fluxunit[:3] == '1E-17':
                    fluxes = [str(x * 1e-17)
                              for x in list(hdulist[i_coadd].data['flux'])]
            else:
                if max([float(x) for x in fluxes]) < 1.0e-5:
                    fluxunit = 'erg/s/cm^2/Angstrom'
                else:
                    fluxunit = 'Uncalibrated'
            specdict = {
                SPECTRUM.U_WAVELENGTHS: 'Angstrom',
                SPECTRUM.WAVELENGTHS: waves,
                SPECTRUM.TIME: mjd,
                SPECTRUM.U_TIME: 'MJD',
                SPECTRUM.FLUXES: fluxes,
                SPECTRUM.U_FLUXES: fluxunit,
                SPECTRUM.FILENAME: ffile,
                SPECTRUM.SOURCE: source
            }
            if 'TELESCOP' in hdrkeys:
                specdict[SPECTRUM.TELESCOPE] = hdulist[0].header['TELESCOP']
            if 'INSTRUME' in hdrkeys:
                specdict[SPECTRUM.INSTRUMENT] = hdulist[0].header['INSTRUME']
            if 'SITENAME' in hdrkeys:
                specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['SITENAME']
            elif 'OBSERVAT' in hdrkeys:
                specdict[SPECTRUM.OBSERVATORY] = hdulist[0].header['OBSERVAT']
            if 'OBSERVER' in hdrkeys:
                specdict[SPECTRUM.OBSERVER] = hdulist[0].header['OBSERVER']
            if 'AIRMASS' in hdrkeys:
                specdict[SPECTRUM.AIRMASS] = hdulist[0].header['AIRMASS']
            catalog.entries[name].add_spectrum(**specdict)
            cntsdss += 1
            hdulist.close()
            catalog.journal_entries()
    print('`{}` have SDSS spectra.'.format(
        cntsdss))
    return
