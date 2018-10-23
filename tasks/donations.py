"""Import tasks for data directly donated to the Open FastStars Catalog."""
import csv
import json
import os
from decimal import Decimal
from glob import glob
from math import floor, isnan

import numpy as np
from astrocats.catalog.photometry import PHOTOMETRY, set_pd_mag_from_counts
from astrocats.catalog.spectrum import SPECTRUM
from astrocats.catalog.utils import (get_sig_digits, is_number, jd_to_mjd,
                                     pbar, pbar_strings, pretty_num, rep_chars)
from astropy.io.ascii import read
from astropy.time import Time as astrotime

from ..faststars import FASTSTARS


def do_donated_spectra(catalog):
    """Import donated spectra."""
    task_str = catalog.get_current_task_str()
    fpath = os.path.join(catalog.get_current_task_repo(), 'Spectra')
    with open(os.path.join(fpath, 'meta.json'), 'r') as f:
        metadict = json.loads(f.read())

    donationscnt = 0
    oldname = ''
    for fname in pbar(metadict, task_str):
        name = metadict[fname]['name']
        name = catalog.get_preferred_name(name)
        if oldname and name != oldname:
            catalog.journal_entries()
        oldname = name
        sec_bibc = metadict[fname]['bibcode']
        name, source = catalog.new_entry(name, bibcode=sec_bibc)

        date = metadict[fname].get('date', '')
        year, month, day = date.split('/')
        sig = get_sig_digits(day) + 5
        day_fmt = str(floor(float(day))).zfill(2)
        time = astrotime(year + '-' + month + '-' + day_fmt).mjd
        time = time + float(day) - floor(float(day))
        time = pretty_num(time, sig=sig)

        with open(os.path.join(fpath, fname), 'r') as f:
            specdata = list(
                csv.reader(
                    f, delimiter=' ', skipinitialspace=True))
            specdata = list(filter(None, specdata))
            newspec = []
            oldval = ''
            for row in specdata:
                if row[0][0] == '#':
                    continue
                if row[1] == oldval:
                    continue
                newspec.append(row)
                oldval = row[1]
            specdata = newspec
        haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][
            2] != 'NaN'
        specdata = [list(i) for i in zip(*specdata)]

        wavelengths = specdata[0]
        fluxes = specdata[1]
        errors = ''
        if haserrors:
            errors = specdata[2]

        specdict = {
            SPECTRUM.U_WAVELENGTHS: 'Angstrom',
            SPECTRUM.U_TIME: 'MJD',
            SPECTRUM.TIME: time,
            SPECTRUM.WAVELENGTHS: wavelengths,
            SPECTRUM.FLUXES: fluxes,
            SPECTRUM.ERRORS: errors,
            SPECTRUM.SOURCE: source,
            SPECTRUM.FILENAME: fname
        }
        if 'instrument' in metadict[fname]:
            specdict[SPECTRUM.INSTRUMENT] = metadict[fname]['instrument']
        if 'telescope' in metadict[fname]:
            specdict[SPECTRUM.TELESCOPE] = metadict[fname]['telescope']
        if 'yunit' in metadict[fname]:
            specdict[SPECTRUM.U_FLUXES] = metadict[fname]['yunit']
            specdict[SPECTRUM.U_ERRORS] = metadict[fname]['yunit']
        else:
            if max([float(x) for x in fluxes]) < 1.0e-5:
                fluxunit = 'erg/s/cm^2/Angstrom'
            else:
                fluxunit = 'Uncalibrated'
            specdict[SPECTRUM.U_FLUXES] = fluxunit
            specdict[SPECTRUM.U_ERRORS] = fluxunit
        catalog.entries[name].add_spectrum(**specdict)
        donationscnt = donationscnt + 1
        if (catalog.args.travis and
                donationscnt % catalog.TRAVIS_QUERY_LIMIT == 0):
            break

    catalog.journal_entries()
    return
