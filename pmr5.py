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
from myUtils import plotLBMarks,raDecToLonLat,angularDistancePyAsl
import subprocess
from gaia_distance_from_parallax import get_gaia_distance
import lightkurve as lk

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
print('dir(gaia) = ',dir(gaia_stars))
print('gaia_stars.colnames = ',gaia_stars.colnames)
print('v_rad = ',gaia_stars['radial_velocity'])
print('Teff = ',gaia_stars['teff_val'])
print('parallax = ',gaia_stars['parallax'])
dist,fwhm_lo,fwhm_hi,post,quant,accept_rate = get_gaia_distance(gaia_stars,display=True)
print('dist = ',dist)
print('fwhm_lo = ',fwhm_lo)
print('fwhm_hi = ',fwhm_hi)
print('post = ',post)
print('quant = ',quant)
print('accept_rate = ',accept_rate)
dist = np.array([dist,fwhm_lo,fwhm_hi])
diam_arcsec = 32.
diam_deg = diam_arcsec / 3600.
print('diam_deg = ',diam_deg)
diam_rad = np.deg2rad(diam_deg)
print('diam_rad = ',diam_rad)
diam = 2. * dist * np.sin(diam_rad) * 3.086e+13
print('diam = ',diam,' km')
age_in_sec = diam / 169.
age_in_yrs = age_in_sec / 3600. / 24. / 365.25
print('age = ',age_in_yrs)

gaia_stars.write('/Users/azuri/daten/uni/HKU/interns_projects/natalie/gaia_PMR5.csv',format='pandas.csv')
sourceID = gaia_stars['SOURCE_ID']

print('gaia.keys = ',dir(gaia_stars))

import csvFree,csvData
nov = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/natalie/historicalNovaeTable3.csv')
novae = pd.read_csv('/Users/azuri/daten/uni/HKU/interns_projects/natalie/historicalNovaeTable3.csv')
print('novae.to_string()')

RAs = novae[:]['RA']
RAs = [float(RA) for RA in RAs]
print('RAs = ',len(RAs),': ',RAs)
DECs = novae[:]['DE']
DECs = [float(DEC) for DEC in DECs]
print('DECs = ',len(DECs),': ',DECs)
radii = novae[:]['Radius']
radii = [float(rad) for rad in radii]
print('radii = ',len(radii),': ',radii)

ham = Hammer()
Ls = []
Bs = []
Xs = []
Ys = []
dists = []
for i in range(len(RAs)):
    l,b = raDecToLonLat(RAs[i],DECs[i])
    Ls.append(l)
    Bs.append(b)
    xy = ham.lonLatToXY(l,b)
    Xs.append(xy.x)
    Ys.append(xy.y)
    dists.append(angularDistancePyAsl(244.9175,-49.2331, RAs[i], DECs[i]))
minDistIdx = np.argsort(np.array(dists))[0]
minDist = dists[minDistIdx]
year = novae['Year'][minDistIdx]
print('minDist = ',minDist,', year = ',year)

if False:
    xy = ham.lonLatToXY(333.9295,0.6863)
    #                                    print('xy = ',xy)
    x = xy.x
    y = xy.y
    fig = plt.figure(figsize=(25,10))
    plt.axis('off')
    plt.scatter(x,y,marker='o',s=10,c='r')
    #l,b = raDecToLonLat(0.,0.)
    xy = ham.lonLatToXY(0,0)
    plt.scatter(xy.x,xy.y,marker='d',s=20)
    #l,b = raDecToLonLat(180.,0.)
    xy = ham.lonLatToXY(180,0)
    plt.scatter(xy.x,xy.y,marker='d',s=20)
    for i in range(len(Xs)):
        plt.scatter(Xs[i],Ys[i],marker='o',s=radii[i]*100.,c='b')
    plotLBMarks(10)
    plt.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(axis='y',          # changes apply to the y-axis
                    which='both',      # both major and minor ticks are affected
                    left=False,      # ticks along the bottom edge are off
                    right=False,         # ticks along the top edge are off
                    labelleft=False) # labels along the bottom edge are off
    fig.tight_layout()
    fNameOut = '/Users/azuri/daten/uni/HKU/interns_projects/natalie/PMR5_in_Hammer_projection.pdf'
    fNameOutTemp = fNameOut[:-3]+'.tmp.pdf'
    fig.savefig(fNameOutTemp, bbox_inches='tight')
    subprocess.run(["gs","-sDEVICE=pdfwrite","-dCompatibilityLevel=1.4","-dPDFSETTINGS=/ebook","-dNOPAUSE", "-dQUIET", "-dBATCH", "-sOutputFile="+fNameOut, fNameOutTemp])
    subprocess.run(["rm",fNameOutTemp])
    plt.show()
    plt.close(fig)

retrieval_type = 'ALL'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
data_structure = 'INDIVIDUAL'     # Options are: 'INDIVIDUAL' or 'RAW'
data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
print(Gaia.load_data.__code__.co_argcount,': ',Gaia.load_data.__code__.co_varnames)
datalink = Gaia.load_data(ids=[sourceID],
                          #data_release=data_release,
                          retrieval_type=retrieval_type,
                          #data_structure=data_structure,
                          verbose=True
                          )

print('datalink = ',datalink)
dl_keys  = [inp for inp in datalink.keys()]
dl_keys.sort()
print(f'The following Datalink products have been downloaded:')
for dl_key in dl_keys:
    print(f' * {dl_key}')


search_result = lk.search_lightcurve('16:19:40.19 -49:13:59.00')
print('search_result = ',search_result)
STOP

print('dir(all_tables[0]) = ',dir(all_tables[0]))
for table in all_tables:
    print('checking table <',table.name,'>: ',table.name)
    try:
        query_astro_param = "SELECT * FROM %s WHERE source_id=%d;" % (table.name,sourceID)#gaiadr3.astrophysical_parameters//
        job=Gaia.launch_job_async(query_astro_param)
        result = job.get_results()
        print('result = ',len(result),': ',result)
        if len(result) > 0:
            result.write('/Users/azuri/daten/uni/HKU/interns_projects/natalie/gaia_%s_PMR5.csv' % (table.name),format='pandas.csv')
        if table == 'gaiadr3.astrophysical_parameters':
            print('teff_gspspec = ',result['teff_gspspec'])
            print('teff_gspphot = ',result['teff_gspphot'])
            print('teff_esphs = ',result['teff_esphs'])
            print('teff_espucd = ',result['teff_espucd'])
            print('teff_msc1 = ',result['teff_msc1'])
            print('teff_msc2 = ',result['teff_msc2'])
    except:
        pass
STOP

def table_to_dataframe(table):
    return table.to_pandas()

tableName = '/Users/azuri/daten/uni/HKU/interns_projects/julia/clusters.dat'
ReadMeName = '/Users/azuri/daten/uni/HKU/interns_projects/julia/ReadMe'
memberName = '/Users/azuri/daten/uni/HKU/interns_projects/julia/members.dat'
clusters = readCDSTable(tableName=tableName,ReadMeName=ReadMeName)
members = readCDSTable(tableName=memberName,ReadMeName=ReadMeName)

def createHSCTable(clusters, members):
    """ ALSO CHECK gaiadr3.astrophysical_parameters """
    query = "SELECT * FROM gaiadr3.gaia_source WHERE source_id=%d;"
    query_astro_param = "SELECT * FROM gaiadr3.astrophysical_parameters WHERE source_id=%d;"

    cluster_data_hsc = []
    cluster_data_lynga = []
    for cluster in clusters:
        if cluster['Name'] =='HSC_2686':
            print('found cluster. name = ',cluster['Name'])
            cluster_data_hsc.append(cluster)
        elif cluster['Name'] =='Lynga_3':
            print('found cluster. name = ',cluster['Name'])
            cluster_data_lynga.append(cluster)

    nStars = 0
    nStarsWithParams = 0
    member_data_hsc =[]
    member_data_lynga =[]
    gaia_data_hsc = None
    astro_param_data_hsc = None
    gaia_data_lynga = None
    astro_param_data_lynga = None
    for member in members:
        if member['Name'] == 'HSC_2686':
            nStars += 1
            print('found ',nStars,' in star belonging to cluster HSC_2686')
            member_data_hsc.append(member)
            job=Gaia.launch_job_async(query % (member['GaiaDR3']))
            result = job.get_results()
            print('type(result) = ',type(result))
            if gaia_data_hsc is None:
                gaia_data_hsc = result
            else:
                print('gaia_data_hsc.colnames = ',len(gaia_data_hsc.colnames),': ',gaia_data_hsc.colnames)
                print('result.colnames = ',len(result.colnames),': ',result.colnames)
                print('result.values() = ',type(result.values()),': ',result.values())
                gaia_data_hsc.add_row(result.values())
            gaia_data_hsc.write(tableName[:tableName.rfind('/')]+'/HSC_2686.csv',format='pandas.csv')
            print('gaia_data_hsc = ',type(gaia_data_hsc),': ',gaia_data_hsc)
            print('dir(gaia_data_hsc) = ',dir(gaia_data_hsc))
            #print('Gaia_data_hsc from the table= ',gaia_data_hsc)
            job=Gaia.launch_job_async(query_astro_param % (member['GaiaDR3']))
            result = job.get_results()
            if len(result) > 0:
                print('job.get_results() = ',result)
                nStarsWithParams += 1
                if astro_param_data_hsc is None:
                    astro_param_data_hsc = result
                else:
                    astro_param_data_hsc.add_row(result.values())
                astro_param_data_hsc.write(tableName[:tableName.rfind('/')]+'/astro_param_HSC_2686.csv',format='pandas.csv')
        elif member['Name'] == 'Lynga_3':
            nStars += 1
            print('found ',nStars,' in star belonging to cluster Lynga_3')
            member_data_lynga.append(member)
            job=Gaia.launch_job_async(query % (member['GaiaDR3']))
            result = job.get_results()
            if gaia_data_lynga is None:
                gaia_data_lynga = result
            else:
                gaia_data_lynga.add_row(result.values())
            gaia_data_lynga.write(tableName[:tableName.rfind('/')]+'/Lynga3.csv',format='pandas.csv')
            #print('gaia_data_lynga from the table = ',gaia_data_lynga)
            job=Gaia.launch_job_async(query_astro_param % (member['GaiaDR3']))
            result = job.get_results()
            if len(result) > 0:
                nStarsWithParams += 1
                if astro_param_data_lynga is None:
                    astro_param_data_lynga = result
                else:
                    astro_param_data_lynga.add_row(result.values())
                astro_param_data_lynga.write(tableName[:tableName.rfind('/')]+'/astro_param_Lynga3.csv',format='pandas.csv')
    print('Column names of hsc_data in gaia data table:', gaia_data_hsc[0].colnames)
    print('Column names of lynga_data in gaia data table:', gaia_data_lynga[0].colnames)
    print('nStarsWithParams = ',nStarsWithParams)

    print('cluster_data_hsc = ',cluster_data_hsc)
    table_cluster_data_hsc = Table(rows=cluster_data_hsc)
    print('table_cluster_data_hsc = ',table_cluster_data_hsc)
    print('member_data_hsc = ',member_data_hsc)
    table_member_data_hsc = Table(rows=member_data_hsc)
    print('table_member_data_hsc = ',table_member_data_hsc)
    print('gaia_data_hsc = ', gaia_data_hsc)
    table_gaia_data_hsc = Table(rows=gaia_data_hsc)
    print('table_gaia_data_hsc = ', gaia_data_hsc)

    print('cluster_data_lynga = ',cluster_data_lynga)
    table_cluster_data_lynga = Table(rows=cluster_data_lynga)
    print('table_cluster_data_lynga = ',table_cluster_data_lynga)
    print('member_data_lynga = ',member_data_lynga)
    table_member_data_lynga = Table(rows=member_data_lynga)
    print('table_member_data_lynga = ',table_member_data_lynga)
    print('gaia_data_lynga = ', gaia_data_lynga)
    table_gaia_data_lynga = Table(rows=gaia_data_lynga)
    print('table_gaia_data_lynga = ', gaia_data_lynga)

    print('astro_param_data_hsc = ', astro_param_data_hsc)
    table_astro_param_data_hsc = Table(rows=astro_param_data_hsc)
    print('table_astro_param_data_hsc = ', table_astro_param_data_hsc)
    print('astro_param_data_lynga = ', astro_param_data_lynga)
    table_astro_param_data_lynga = Table(rows=astro_param_data_lynga)
    print('table_astro_param_data_lynga = ', table_astro_param_data_lynga)

    return (table_cluster_data_hsc, table_member_data_hsc, table_gaia_data_hsc,
            table_cluster_data_lynga, table_member_data_lynga, table_gaia_data_lynga,
            table_astro_param_data_hsc, table_astro_param_data_lynga)

cluster_data_hsc, member_data_hsc, gaia_data_hsc, cluster_data_lynga, member_data_lynga, gaia_data_lynga, astro_params_hsc, astro_params_lynga = createHSCTable(clusters, members)

member_dfs_hsc = [table_to_dataframe(t) for t in member_data_hsc]
gaia_dfs_hsc = [table_to_dataframe(t) for t in gaia_data_hsc]
member_dfs_lynga = [table_to_dataframe(t) for t in member_data_lynga]
gaia_dfs_lynga = [table_to_dataframe(t) for t in gaia_data_lynga]
astro_params_dfs_hsc = [table_to_dataframe(t) for t in astro_params_hsc]
astro_params_dfs_lynga = [table_to_dataframe(t) for t in astro_params_lynga]

member_df_hsc = pd.concat(member_dfs_hsc, ignore_index=True)
gaia_df_hsc = pd.concat(gaia_dfs_hsc, ignore_index=True)
member_df_lynga = pd.concat(member_dfs_lynga, ignore_index=True)
gaia_df_lynga = pd.concat(gaia_dfs_lynga, ignore_index=True)
astro_params_df_hsc = pd.concat(astro_params_dfs_hsc, ignore_index=True)
astro_params_df_lynga = pd.concat(astro_params_dfs_lynga, ignore_index=True)

merged_df_hsc = pd.merge(member_df_hsc, gaia_df_hsc, on='source_id', how='left', suffixes=('_member', '_gaia'))
merged_df_lynga = pd.merge(member_df_lynga, gaia_df_lynga, on='source_id', how='left', suffixes=('_member', '_gaia'))
merged_df_astro_params = pd.merge(astro_params_df_hsc, astro_params_df_lynga, on='source_id', how='left', suffixes=('_member', '_gaia'))

print('Merged DataFrame of HSC:', merged_df_hsc)
print('Merged DataFrame of Lynga:', merged_df_lynga)
print('Merged DataFrame of astrophysical parameters:', merged_df_astro_params)

common_columns = set(member_df_hsc.columns).intersection(gaia_df_hsc.columns) - {'source_id'}
for col in common_columns:
    diff_col_name = f'{col}_diff'
    merged_df_hsc[diff_col_name] = merged_df_hsc[f'{col}_gaia'] - merged_df_hsc[f'{col}_member']
print("\nDifferences in Parameters:")
for col in common_columns:
    diff_col_name = f'{col}_diff'
    print(f"\nDifferences in {col}:")
    print(merged_df_hsc[['source_id', diff_col_name]])

print('clusters = ',clusters)
print('dir(clusters) = ',dir(clusters))
print('members = ',members)
print('dir(members) = ',dir(members))
print('Column names in clusters table:', clusters.colnames)
print('Column names in members table: ',members.colnames)


"""
Hertzsprung-Russel Diagram (absolute magnitude : Temperature)
"""

def calculate_absolute_magnitude(df):
    df['distance_pc'] = 1 / df['parallax']  # Convert parallax from mas to parsecs
    df['absolute_mag_g'] = df['phot_g_mean_mag'] - 5 * np.log10(df['distance_pc']) + 5
    return df

gaia_df_hsc = calculate_absolute_magnitude(gaia_df_hsc)
gaia_df_lynga = calculate_absolute_magnitude(gaia_df_lynga)

temp_hsc = astro_params_df_hsc['teff_gspphot']
mag_hsc = gaia_df_hsc['absolute_mag_g']
temp_lynga = astro_params_df_lynga['teff_gspphot']
mag_lynga = gaia_df_lynga['absolute_mag_g']

# Plot the HR Diagram
plt.figure(figsize=(10, 8))
plt.scatter(temp_hsc, mag_hsc, color='blue', label='HSC_2686')
plt.scatter(temp_lynga, mag_lynga, color='red', label='Lynga_3')
plt.gca().invert_yaxis()

plt.xlabel('Effective Temperature (K)')
plt.ylabel('Absolute Magnitude (G)')
plt.title('Hertzsprung-Russell Diagram')
plt.legend()
plt.show()
