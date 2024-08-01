import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.io import ascii
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

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


def readGaiaMainTable(ra_deg, dec_deg):
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
    Gaia.ROW_LIMIT = 1  # Ensure the default row limit.
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(3.0/3600., u.deg))
    r = j.get_results()
    return r

def showGaiaTables():
    tables = Gaia.load_tables(only_names=True)
    for table in tables:
        print(table.get_qualified_name())

def readGaia_astrophysical_parameters():
    """ THIS ONE """
    table = Gaia.load_table('gaiadr3.astrophysical_parameters')
    return table

def readCDSTable(tableName,ReadMeName):
    table = ascii.read(tableName,
            readme=ReadMeName)
    return table

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
