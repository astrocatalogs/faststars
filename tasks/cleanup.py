"""Cleanup catalog before final write to disk."""
import re
import statistics
import warnings
from decimal import Decimal
from math import log10, pi, sqrt
import numpy as np

from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.utils import (get_sig_digits, is_number, pbar,
                                     pretty_num, tprint, uniq_cdl)
from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value

from ..constants import CLIGHT, KM
from ..faststars import FASTSTARS


def do_cleanup(catalog):
    """Cleanup catalog after importing all data."""
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
    keys = list(catalog.entries.keys())

    cleanupcnt = 0
    for oname in pbar(keys, task_str):
        # Some events may be merged in cleanup process, skip them if
        # non-existent.
        try:
            name = catalog.add_entry(oname)
        except Exception:
            catalog.log.warning(
                '"{}" was not found, suggests merge occurred in cleanup '
                'process.'.format(oname))
            continue

        # Set the preferred name, switching to that name if name changed.
        name = catalog.entries[name].set_preferred_name()

        aliases = catalog.entries[name].get_aliases()
        catalog.entries[name].set_first_max_light()
        
        # Clean discoverer field
        if FASTSTARS.DISCOVERER in catalog.entries[name]:
            if len(catalog.entries[name][FASTSTARS.DISCOVERER]) > 1:
                POSSIBLEDISCOVERER = [catalog.entries[name][FASTSTARS.DISCOVERER][i]['value'] for i in range(len(catalog.entries[name][FASTSTARS.DISCOVERER]))]
                POSSIBLEDISCOVERER_DATE = [int(DATE['value']) for DATE in catalog.entries[name][FASTSTARS.DISCOVER_DATE]]
                POSSIBLEDISCOVERER_DATE_SOURCES = [DATE['source'] for DATE in catalog.entries[name][FASTSTARS.DISCOVER_DATE]]
                EARLIESTSOURCE = POSSIBLEDISCOVERER_DATE_SOURCES[np.argmin(POSSIBLEDISCOVERER_DATE)]
                EARLIESTDISCOVER_DATE = catalog.entries[name][FASTSTARS.DISCOVER_DATE][np.argmin(POSSIBLEDISCOVERER_DATE)]

                # Deal with case where a star was 'discovered' multiple times in one year
                if ',' in EARLIESTSOURCE:
                    EARLIESTSOURCE = EARLIESTSOURCE.split(',')[0]

                for DISCOVERER in catalog.entries[name][FASTSTARS.DISCOVERER]:
                    for DISCOVERERSOURCE in DISCOVERER['source'].split(','):
                        if DISCOVERERSOURCE == EARLIESTSOURCE:
                            EARLIESTDISCOVERER = DISCOVERER
                
                for DISCOVERER in catalog.entries[name][FASTSTARS.DISCOVERER]:
                    for DISCOVERERSOURCE in DISCOVERER['source'].split(','):
                        if DISCOVERERSOURCE == EARLIESTSOURCE:
                            EARLIESTDISCOVERER = DISCOVERER
                catalog.entries[name][FASTSTARS.DISCOVERER] = [EARLIESTDISCOVERER]
                catalog.entries[name][FASTSTARS.DISCOVER_DATE] = [EARLIESTDISCOVER_DATE]

        # Convert all distances to kpc.
        if FASTSTARS.LUM_DIST in catalog.entries[name]:
            for li, ld in enumerate(catalog.entries[name][FASTSTARS.LUM_DIST]):
                if ld.get('u_value') != 'kpc':
                    if ld.get('u_value') == 'pc':
                        catalog.entries[name][FASTSTARS.LUM_DIST][li]['value'] = str(Decimal(
                            catalog.entries[name][FASTSTARS.LUM_DIST][li]['value']) * Decimal('0.001'))
                    elif ld.get('u_value') == 'Mpc':
                        catalog.entries[name][FASTSTARS.LUM_DIST][li]['value'] = str(Decimal(
                            catalog.entries[name][FASTSTARS.LUM_DIST][li]['value']) * Decimal('1000'))
                    else:
                        raise ValueError('unknown distance unit')
                    catalog.entries[name][FASTSTARS.LUM_DIST][li]['u_value'] = 'kpc'

        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]):
            prefixes = [
                'SDSS'
            ]
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:6])):
                        noprefix = alias.split(':')[-1].replace(
                            prefix, '').replace('.', '')
                        decsign = '+' if '+' in noprefix else '-'
                        noprefix = noprefix.replace('+', '|').replace('-', '|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + \
                            ('.' + rastr[6:] if len(rastr) > 6 else '')
                        dec = (
                            decsign + ':'.join(
                                [decstr[:2], decstr[2:4], decstr[4:6]]) +
                            ('.' + decstr[6:] if len(decstr) > 6 else ''))
                        if catalog.args.verbose:
                            tprint('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = catalog.entries[name].add_self_source()
                        catalog.entries[name].add_quantity(
                            FASTSTARS.RA, ra, source, derived=True)
                        catalog.entries[name].add_quantity(
                            FASTSTARS.DEC, dec, source, derived=True)
                        break
                if FASTSTARS.RA in catalog.entries[name]:
                    break

        if (FASTSTARS.MAX_ABS_MAG not in catalog.entries[name] and
                FASTSTARS.MAX_APP_MAG in catalog.entries[name] and
                FASTSTARS.LUM_DIST in catalog.entries[name]):
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in catalog.entries[name][FASTSTARS.LUM_DIST]:
                sig = get_sig_digits(ld[QUANTITY.VALUE])
                if sig > bestsig:
                    bestld = ld[QUANTITY.VALUE]
                    bestsrc = ld[QUANTITY.SOURCE]
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                source = catalog.entries[name].add_self_source()
                sources = uniq_cdl([source] + bestsrc.split(','))
                bestldz = z_at_value(cosmo.luminosity_distance,
                                     float(bestld) * un.Mpc)
                pnum = (
                    float(catalog.entries[name][FASTSTARS.MAX_APP_MAG][0][
                        QUANTITY.VALUE]) - 5.0 *
                    (log10(float(bestld) * 1.0e6) - 1.0
                     ) + 2.5 * log10(1.0 + bestldz))
                pnum = pretty_num(pnum, sig=bestsig + 1)
                catalog.entries[name].add_quantity(
                    FASTSTARS.MAX_ABS_MAG, pnum, sources, derived=True)

        catalog.entries[name].sanitize()
        catalog.journal_entries(bury=True, final=True, gz=True)

        cleanupcnt = cleanupcnt + 1
        if catalog.args.travis and cleanupcnt % 1000 == 0:
            break

    catalog.save_caches()

    return

