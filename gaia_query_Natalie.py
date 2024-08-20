import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.io import ascii
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.table import Table
from hammer import Pixel,XY,LonLat,Hammer
from pnOrientationUtils import calcMoments,calcMean,rose_plot,vectorDiagram,plotHammerProjection,linearOrderDiagram,plotEPA
from myUtils import plotLBMarks,raDecToLonLat
import subprocess

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    try:
        h, m, s = [float(i) for i in string.split(':')]
        return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)
    except:
        print('hmsToDeg: string = <'+string+'>')

def degToHMS(degrees):
    h = int(degrees / 15.)
    m = int((degrees - (h * 15.)) * 4.)
    s = (degrees - (m/4.) - (h*15.)) * 240.
    sStr = '%.3f' % (s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (h,m,sStr)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    try:
        d, m, s = [float(i) for i in string.split(':')]
    #    print('dmsToDeg: string = <'+string+'>: d = ',d,', m = ',m,', s = ',s)
        if string[0] == '-':
            d = 0. - d
            return 0. - (s / 3600. + m / 60. + d)
        return s / 3600. + m / 60. + d
    except:
        return float(string)

def degToDMS(degrees):
    d = int(degrees)
    m = int((degrees - d) * 60)
    s = (degrees - d - (m/60.)) * 3600.
    sStr = '%.3f' % abs(s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (d,abs(m),sStr)


def readGaiaMainTable(ra_deg, dec_deg, rad_deg, row_limit=10000000):
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
    Gaia.ROW_LIMIT = row_limit  # Ensure the default row limit.
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))
    r = j.get_results()
    return r

def showGaiaTables():
    tables = Gaia.load_tables(only_names=True)
    for table in tables:
        print(table.get_qualified_name())
    return tables

def readGaia_astrophysical_parameters():
    """ THIS ONE """
    table = Gaia.load_table('gaiadr3.astrophysical_parameters')
    return table

def readCDSTable(tableName,ReadMeName):
    table = ascii.read(tableName,
            readme=ReadMeName)
    return table

all_tables = showGaiaTables()
gaia_stars = readGaiaMainTable(244.9175, -49.2331, 1./3600., row_limit=1)
print('gaia = ',gaia_stars)
print('v_rad = ',gaia_stars['radial_velocity'])
print('Teff = ',gaia_stars['teff_val'])

gaia_stars.write('/Users/azuri/daten/uni/HKU/interns_projects/natalie/gaia_PMR5.csv',format='pandas.csv')
sourceID = gaia_stars['SOURCE_ID']


for table in all_tables:
    print('checking table <'+table+'>')
    query_astro_param = "SELECT * FROM %s WHERE source_id=%d;" % (table,sourceID)#gaiadr3.astrophysical_parameters
    job=Gaia.launch_job_async(query_astro_param)
    result = job.get_results()
    print('result = ',result)
    if len(result) > 0:
        result.write('/Users/azuri/daten/uni/HKU/interns_projects/natalie/gaia_%s_PMR5.csv' % (table),format='pandas.csv')
    if table == 'gaiadr3.astrophysical_parameters':
        print('teff_gspspec = ',result['teff_gspspec'])
        print('teff_gspphot = ',result['teff_gspphot'])
        print('teff_esphs = ',result['teff_esphs'])
        print('teff_espucd = ',result['teff_espucd'])
        print('teff_msc1 = ',result['teff_msc1'])
        print('teff_msc2 = ',result['teff_msc2'])
