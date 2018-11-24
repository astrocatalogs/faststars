"""Query against the Gaia database.
"""
import re
import os

import numpy as np
from astropy.coordinates import SkyCoord as coord
import astropy.units as un
from astroquery.vizier import Vizier
from astropy.table import vstack
import warnings
import pandas as pd
from scipy.interpolate import interp1d

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
    
#silentgaiaquery = silent(Gaia.query_object)
#v=Vizier(columns=['RA_ICRS','DE_ICRS','Plx','pmRA','pmDE','<Gmag>','e_RA_ICRS','e_DE_ICRS','e_Plx','e_pmRA','e_pmDE','RADEcor','RAPlxcor','RApmRAcor','RApmDEcor','DEPlxcor','DEpmRAcor','DEpmDEcor','PlxpmRAcor','PlxpmDEcor','pmRApmDEcor'])
v=Vizier(columns=['**'])
#silentgaiaquery = silent(v.query_region)
silentgaiaregionquery = v.query_region
silentgaiaobjectquery = v.query_object



def do_gaiaviavizier(catalog):
    # Disable warnings
    warnings.filterwarnings("ignore")
    catalog.log.warning(
                'Some warnings are thrown by the Gaia module which do not affect the result of the Gaia queries.')
    
    task_str = catalog.get_current_task_str()
    keys = list(catalog.entries.keys())

    ##### Gaia parallax offset
    gaiaparallaxoffset = -0.029
    datafile = os.path.join(catalog.get_current_task_repo(), 'AUXILIARY', 'lindegren_error_inflation.csv')
    lei_data = pd.read_csv(datafile, header=None)
    lei_G = lei_data[0].values # Gaia magnitude
    lei_EIR = lei_data[1].values # ratio of external to internal error
    lei_EIR_interp = interp1d(lei_G,lei_EIR)
    def add_lindegren_inflation(AST_ERRORS,GMAG):
        _EIR_factor = lei_EIR_interp(float(GMAG))
        if isinstance(AST_ERRORS,list):
            return [str(float(AE)*_EIR_factor) for AE in AST_ERRORS]
        else:
            return str(float(AST_ERRORS)*_EIR_factor)

    cntgphot = 0
    cntgast = 0
    for oname in pbar(keys, task_str):
        try:
            name = catalog.add_entry(oname)
        except Exception:
            catalog.log.warning(
                '"{}" was not found, suggests merge occurred in cleanup '
                'process.'.format(oname))
            continue
        
        ### Check if has DR1 or DR2 ID
        gaiaid = ''
        try:
            for alias in catalog.entries[name][FASTSTARS.ALIAS]:
                if alias['value'][:8] == 'Gaia DR2' and len(alias['value'])>len(gaiaid):
                    gaiaid = alias['value']
                elif alias['value'][:8] == 'Gaia DR1' and len(alias['value'])>len(gaiaid):
                    gaiaid = alias['value']
        except KeyError:
            print(name)
        hasgaiaid = True if len(gaiaid)>0 else False
        
        if (FASTSTARS.RA not in catalog.entries[name] or
                FASTSTARS.DEC not in catalog.entries[name]) and not hasgaiaid:
            continue
        else:
            if hasgaiaid:
                result = silentgaiaobjectquery(gaiaid, radius=un.Quantity(10.0, un.arcsecond), catalog='I/345/gaia2')
                if len(result)<1:
                    catalog.log.warning(
                '"{}" should really have had a cross-match.'.format(name))
                    print(gaiaid)
            else:
                Mradec=str(catalog.entries[name][FASTSTARS.RA][0]['value'])+str(catalog.entries[name][FASTSTARS.DEC][0]['value'])
                #print(name,Mradec)
                c=coord(Mradec,unit=(un.hourangle, un.deg),frame='icrs')
                
                queryradius=0.
                if (FASTSTARS.PROPER_MOTION_RA in catalog.entries[name] and FASTSTARS.PROPER_MOTION_DEC in catalog.entries[name]):
                    queryradius = 1e-3*3.*10.*np.sqrt(float(catalog.entries[name][FASTSTARS.PROPER_MOTION_RA][0]['value'])**2.+float(catalog.entries[name][FASTSTARS.PROPER_MOTION_DEC][0]['value'])**2.)
                queryradius = max(queryradius,3.)
                
                if name in []:
                    result = silentgaiaregionquery(c, radius=un.Quantity(30.0, un.arcsecond), catalog='I/345/gaia2')
                else:
                    result = silentgaiaregionquery(c, radius=un.Quantity(queryradius, un.arcsecond), catalog='I/345/gaia2')
                        
            if len(result)>0:
                if len(result[0])>1:
                    indexclosest = np.argmin(result[0]['_r'])
                    result=result[0][indexclosest]
                else:
                    result=result[0][0]
                    
                # Utility function for accessing the Vizier table's columns
                def gtab(parstring):
                    return str(result[parstring])
                
                cntgphot += 1
                source = catalog.entries[name].add_source(bibcode='2016A&A...595A...2G')
                catalog.entries[name].add_photometry(
                        time=57023,
                        u_time='MJD',
                        telescope='Gaia',
                        band='G',
                        magnitude=gtab('Gmag'),
                        e_magnitude=gtab('e_Gmag'),
                        source=source)
                catalog.entries[name].add_photometry(
                        time=57023,
                        u_time='MJD',
                        telescope='Gaia',
                        band='GBP',
                        magnitude=gtab('BPmag'),
                        e_magnitude=gtab('e_BPmag'),
                        source=source)
                catalog.entries[name].add_photometry(
                        time=57023,
                        u_time='MJD',
                        telescope='Gaia',
                        band='GRP',
                        magnitude=gtab('RPmag'),
                        e_magnitude=gtab('e_RPmag'),
                        source=source)
                if not hasgaiaid: catalog.entries[name].add_quantity(FASTSTARS.ALIAS, gtab('DR2Name'), source=source)
                gra,gde = coord(ra=float(gtab('RA_ICRS'))*un.deg,dec=float(gtab('DE_ICRS'))*un.deg,frame='icrs').to_string('hmsdms', sep=':').split()
                if gtab('Plx') == '--':
                    catalog.entries[name].add_quantity(FASTSTARS.RA,gra, source, e_value=add_lindegren_inflation(gtab('e_RA_ICRS'),gtab('Gmag')), u_e_value='mas', correlations=[{CORRELATION.VALUE:gtab('RADEcor'), CORRELATION.QUANTITY:FASTSTARS.DEC, CORRELATION.KIND:'Pearson'}])
                    catalog.entries[name].add_quantity(FASTSTARS.DEC,gde, source, e_value=add_lindegren_inflation(gtab('e_DE_ICRS'),gtab('Gmag')), u_e_value='mas', correlations=[{CORRELATION.VALUE:gtab('RADEcor'), CORRELATION.QUANTITY:FASTSTARS.RA, CORRELATION.KIND:'Pearson'}])
                else:
                    #print(gtab('Plx'))
                    cntgast += 1
                    #catalog.log.warning(
                    #    '"{}" has Gaia astrometry.'.format(name))
                    
                    ast_keys = [FASTSTARS.RA,FASTSTARS.DEC,FASTSTARS.PARALLAX,FASTSTARS.PROPER_MOTION_RA,FASTSTARS.PROPER_MOTION_DEC]
                    ast_values = [gra,gde,str(float(gtab('Plx'))-gaiaparallaxoffset),gtab('pmRA'),gtab('pmDE')]
                    ast_errors = add_lindegren_inflation([gtab('e_RA_ICRS'),gtab('e_DE_ICRS'),gtab('e_Plx'),gtab('e_pmRA'),gtab('e_pmDE')],gtab('Gmag'))
                    ast_corr = [['1.0',gtab('RADEcor'),gtab('RAPlxcor'),gtab('RApmRAcor'),gtab('RApmDEcor')], [gtab('RADEcor'),'1.0',gtab('DEPlxcor'),gtab('DEpmRAcor'),gtab('DEpmDEcor')], [gtab('RAPlxcor'),gtab('DEPlxcor'),'1.0',gtab('PlxpmRAcor'),gtab('PlxpmDEcor')], [gtab('RApmRAcor'),gtab('DEpmRAcor'),gtab('PlxpmRAcor'),'1.0',gtab('pmRApmDEcor')], [gtab('RApmDEcor'),gtab('DEpmDEcor'),gtab('PlxpmDEcor'),gtab('pmRApmDEcor'),'1.0']]
                    ast_names = ['ra','dec','parallax','propermotionra','propermotiondec']
                    ast_units = ['hms','dms','mas','mas/yr','mas/yr']
                    ast_error_units = ['mas','mas','mas','mas/yr','mas/yr']
                    
                    #print(ast_values[0])
                    n_ast = len(ast_keys) # number of columns
                    for i in range(n_ast):
                        corr_dict = [{CORRELATION.VALUE:ast_corr[i][(i+j) % n_ast], CORRELATION.QUANTITY:ast_names[(i+j) % n_ast], CORRELATION.KIND:'Pearson'} for j in list(range(1,n_ast))]
                        catalog.entries[name].add_quantity(ast_keys[i], ast_values[i], source, e_value=ast_errors[i], u_value=ast_units[i], e_u_value=ast_error_units[i], correlations=corr_dict)
                    
                    if gtab('RV') != '--':
                        catalog.entries[name].add_quantity(FASTSTARS.VELOCITY, gtab('RV'), source, e_value=gtab('e_RV'), u_value='km/s')
                    ##### This has been moved to boundprobability.py
                    # Convert parallax to distance
                    #if (FASTSTARS.LUM_DIST in catalog.entries[name]):
                    #    catalog.log.warning(
                    #    '"{}" has distance from photmetric distance prior.'.format(name)) 
                    #    distance, distance_error = parallax_to_distance(name,result['Plx'][0],result['e_Plx'][0],float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['value']),float(catalog.entries[name][FASTSTARS.LUM_DIST][0]['e_value']))
                    #    catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                    #else:
                    #    catalog.log.warning(
                    #    '"{}" has distance from Astraatmadja & Bailer-Jones (2016) prior.'.format(name)) 
                    #    distance, distance_error = parallax_to_distance(name,result['Plx'][0],result['e_Plx'][0])
                    #    catalog.entries[name].add_quantity(FASTSTARS.LUM_DIST, str(distance), e_value=str(distance_error), u_value='kpc', source=source, derived=True)
                        
    catalog.log.warning(
                '"{}" have Gaia photometry and "{}" have Gaia astrometry.'.format(cntgphot,cntgast)) 
    catalog.journal_entries()
    
    # Reactivate warnings
    warnings.filterwarnings("default")
    return
