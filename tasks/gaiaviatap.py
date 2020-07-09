"""Query against the Gaia database.
"""
import re

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.utils.tap.core import TapPlus
from astropy.table import vstack
import warnings

from astrocats.catalog.photometry import PHOTOMETRY
from astrocats.catalog.utils import is_number, pbar, single_spaces, uniq_cdl
from astrocats.catalog.correlation import CORRELATION
from ..faststars import FASTSTARS
from ..utils import name_clean, parallax_to_distance
mas_to_deg = 1e-3/3600.0

# astroquery.gaia is very verbose
import sys
from io import StringIO 

class NullIO(StringIO):
    def write(self, txt):
        pass

def silent(fn):
    """Decorator to silence functions."""
    def silent_fn(*args, **kwargs):
        saved_stdout = sys.stdout
        sys.stdout = NullIO()
        result = fn(*args, **kwargs)
        sys.stdout = saved_stdout
        return result
    return silent_fn
    
import os, sys
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout
    
#silentgaiaquery = silent(Gaia.query_object)
#v=Vizier(columns=['RA_ICRS','DE_ICRS','Plx','pmRA','pmDE','<Gmag>','e_RA_ICRS','e_DE_ICRS','e_Plx','e_pmRA','e_pmDE','RADEcor','RAPlxcor','RApmRAcor','RApmDEcor','DEPlxcor','DEpmRAcor','DEpmDEcor','PlxpmRAcor','PlxpmDEcor','pmRApmDEcor'])
#v=Vizier(columns=['**'])
#ari = TapPlus(url="http://gaia.ari.uni-heidelberg.de/tap")
ari = TapPlus(url="https://gaia.ari.uni-heidelberg.de/tap")

#silentgaiaquery = silent(v.query_region)
silentari = ari.launch_job

##### Gaia parallax offset
#gaiaparallaxoffset = -0.029


def do_gaiaviatap(catalog):
    # Disable warnings
    warnings.filterwarnings("ignore")
    #catalog.log.warning(
    #            'Some warnings are thrown by the Gaia module which do not affect the result of the Gaia queries.')
    
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())

    cntL = 0
    for oname in pbar(keys, task_str):
        try:
            name = catalog.add_entry(oname)
        except Exception:
            catalog.log.warning(
                '"{}" was not found, suggests merge occurred in cleanup '
                'process.'.format(oname))
            continue
        
        ### Check if has DR2 ID
        dr2id = ''
        for alias in catalog.entries[name][FASTSTARS.ALIAS]:
            if alias['value'][:8] == 'Gaia DR2' and len(alias['value'])>len(dr2id):
                dr2id = alias['value']
        if len(dr2id)<8 or FASTSTARS.PARALLAX not in catalog.entries[name]:
            #print(name)
            continue
        else:
            with HiddenPrints():
                source = catalog.entries[name].add_source(bibcode='2018arXiv180410121B')
                job = silentari("SELECT r_len FROM gaiadr2_complements.geometric_distance where source_id = "+dr2id[8:], dump_to_file=False, verbose=False)
                lengthscale = str(job.get_results()['col_0'][0]/1e3)
                catalog.entries[name].add_quantity(FASTSTARS.DISTANCE_PRIOR_LENGTH_SCALE,lengthscale, source, u_value='kpc')
                cntL += 1        
                        
    catalog.log.warning(
                '"{}" have a more accurate scale length from Bailer-Jones et al. (2018).'.format(cntL)) 
    catalog.journal_entries()
    
    # Reactivate warnings
    warnings.filterwarnings("default")
    return
