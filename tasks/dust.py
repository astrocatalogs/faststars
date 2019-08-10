"""Import tasks for the dust maps.
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from dustmaps.sfd import SFDQuery

from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from ..faststars import FASTSTARS
from ..utils import name_clean, check_dustmaps


def do_dust(catalog):
    task_str = catalog.get_current_task_str()

    # Set preferred names, calculate some columns based on imported data,
    # sanitize some fields
    keys = list(catalog.entries.keys())

    check_dustmaps(catalog.get_current_task_repo())

    sfd = SFDQuery()

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
        
        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]):
            continue
        else:
            Mname = name
            Mradec = str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value'])
            Mdist = '-1'
            c=coord(Mradec,unit=(un.hourangle, un.deg), frame='icrs')
            reddening = sfd(c)
            source = catalog.entries[name].add_source(bibcode='1998ApJ...500..525S')
            catalog.entries[name].add_quantity(FASTSTARS.EBV, str(reddening), source, upperlimit=True, derived=True)

    catalog.journal_entries()
    return
