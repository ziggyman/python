import numpy as np
import matplotlib.pyplot as plt
import csvFree,csvData


#gaiaFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/4289684146940013184.csv'
gaiaFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/HD185806_gaia.csv'
#gaiaFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/hip16518_gaia.csv'
#gaiaFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/hip2599_gaia.csv'
#gaiaFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/hip11891_gaia.csv'
imFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/nebula_wise.png'

pm_cs = 9.26#mas/yr
age = 25750#years
movement = pm_cs * age
print('central star movement = ',movement,' mas')

gaiaData = csvFree.readCSVFile(gaiaFileName)
vrad = gaiaData.getData('radial_velocity',0)
evrad = gaiaData.getData('radial_velocity_error',0)
print('vrad = ',vrad,' +/- ',evrad)
delta = float(gaiaData.getData('dec',0))
pm_ra = float(gaiaData.getData('pmra',0))
pm_ra_error = float(gaiaData.getData('pmra_error',0))
print('pm_ra = ',pm_ra,' +/- ',pm_ra_error)
#print('pm_ra  * np.cos(np.radians(delta)) = ',pm_ra * np.cos(np.radians(delta)))
pm_dec = float(gaiaData.getData('pmdec',0))
pm_dec_error = float(gaiaData.getData('pmdec_error',0))
print('pm_dec = ',pm_dec,' +/- ',pm_dec_error)
pm = float(gaiaData.getData('pm',0))
print('pm = ',pm)
print('sqrt(pmra^2+pmdec^2) = ',np.sqrt((pm_ra*pm_ra)+(pm_dec * pm_dec)))

phi = np.degrees(np.arcsin(pm_ra / pm))
print('phi = ',phi)
alpha = np.degrees(np.arcsin(pm_dec / pm))
print('alpha = ',alpha)

img = plt.imread(imFileName)
fig,ax = plt.subplots()
ax.imshow(img)
print('ax.get_xlim() = ',ax.get_xlim())
print('ax.get_ylim() = ',ax.get_ylim())
center = [277.,415.]
print(center)
#ax.plot([center[0],center[0]-(pm_ra*10.)],[center[1],center[1]+(pm_dec*10.)],'-',linewidth=5.,color='y')
scale = 10.
ax.annotate("", xy=(center[0]-(pm_ra*scale), center[1]-(pm_dec*scale)), xytext=(center[0], center[1]), arrowprops=dict(width=3,color='green'),)#arrowstyle="->",
ax.axis('off')
plt.show()

spectrumFileName = '/Users/azuri/daten/uni/HKU/publications/cocoon/GDR3_CONTINUOUS_MEAN_SPECTRUM_COMBINED.csv'
from gaiaxpy import convert
#converted_spectra, sampling = convert(spectrumFileName)
#converted_spectra
import pandas as pd

def convert_to_nparray_double(strarray):
    '''convert a string representing an array of double into a numpy array of type float64.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        nparray = np.array(eval(strarray), dtype=np.float64)
    else:
        raise ValueError("the string array is not a list")
    return nparray
def convert_to_nparray_float(strarray):
    '''convert a string representing an array of float into a numpy array of type float32.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        nparray = np.array(eval(strarray), dtype=np.float32)
    else:
        raise ValueError("the string array is not a list")
    return nparray
def convert_to_nparray_short(strarray):
    '''convert a string representing an array of short into a numpy array of type int16.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        nparray = np.array(eval(strarray), dtype=np.int16)
    else:
        raise ValueError("the string array is not a list")
    return nparray
def convert_to_nparray_int(strarray):
    '''convert a string representing an array of integer into a numpy array of type int32.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        nparray = np.array(eval(strarray), dtype=np.int32)
    else:
        raise ValueError("the string array is not a list")
    return nparray
def convert_to_nparray_long(strarray):
    '''convert a string representing an array of long into a numpy array of type int64.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        nparray = np.array(eval(strarray), dtype=np.int64)
    else:
        raise ValueError("the string array is not a list")
    return nparray
def convert_to_nparray_bool(strarray):
    '''convert a string representing an array of boolean into a numpy array of type bool.
    '''
    if strarray.startswith('[') and strarray.endswith(']'):
        strarray = strarray.replace('f', 'False').replace('t', 'True')
        nparray = np.array(eval(strarray), dtype=bool)
    else:
        raise ValueError("the string array is not a list")
    return nparray

converter_map = {
    'source_id': convert_to_nparray_long,
    'solution_id': convert_to_nparray_long,
    'bp_basis_function_id': convert_to_nparray_int,
    'bp_degrees_of_freedom': convert_to_nparray_int,
    'bp_n_parameters': convert_to_nparray_int,
    'bp_n_measurements': convert_to_nparray_int,
    'bp_n_rejected_measurements': convert_to_nparray_int,
    'bp_standard_deviation': convert_to_nparray_float,
    'bp_chi_squared': convert_to_nparray_float,
    #'bp_coefficients': convert_to_nparray_double,
    #'bp_coefficient_errors': convert_to_nparray_float,
    #'bp_coefficient_correlations': convert_to_nparray_float,
    'bp_n_relevant_bases': convert_to_nparray_int,
    'bp_relative_shrinking': convert_to_nparray_float,
    'rp_basis_function_id': convert_to_nparray_int,
    'rp_degrees_of_freedom': convert_to_nparray_int,
    'rp_n_parameters': convert_to_nparray_int,
    'rp_n_measurements': convert_to_nparray_int,
    'rp_n_rejected_measurements': convert_to_nparray_int,
    'rp_standard_deviation': convert_to_nparray_float,
    'rp_chi_squared': convert_to_nparray_float,
    #'rp_coefficients': convert_to_nparray_double,
    #'rp_coefficient_errors': convert_to_nparray_float,
    #'rp_coefficient_correlations': convert_to_nparray_float,
    'rp_n_relevant_bases': convert_to_nparray_int,
    'rp_relative_shrinking': convert_to_nparray_float,
}

if False:
    df = pd.read_csv(spectrumFileName, converters=converter_map)
    print(df.head())
    #df = pd.read_csv(spectrumFileName) # The values in the DataFrame can be edited if the user wishes to do so.
    #print('df = ',df)
    #converted_spectra, sampling = convert(df)
    #converted_spectra

# velocities
velocities = [[-75.,20.],
              [-100.,40.],
              [-85.,20.],
              [-50.,5.],
              [-65.,-7],
              [-55.,10.],
              [-85.,30.],
              [-90.,-20.],
              [-50.,5.],
              [-80.,-30.],
              [-55.,15.],
              [-80.,-40.],
              [-60.,-20.],
             ]

minMaxVelocities = [[-100.,40.],
                    [-75.,35.],
                    [-100.,30.],
                    [-80.,5.],
                    [-90.,15.],
                    ]

expVelos = np.array([np.absolute(a[0]-a[1])/2. for a in velocities])
meanExpVelo = np.mean(expVelos)
sigExpVelo = np.std(expVelos)
print('expVelo = ',meanExpVelo,' +/- ',sigExpVelo)

minMaxExpVelos = np.array([np.absolute(a[0]-a[1])/2. for a in minMaxVelocities])
meanMinMaxExpVelo = np.mean(minMaxExpVelos)
sigMinMaxExpVelo = np.std(minMaxExpVelos)
print('expVeloMinMax = ',meanMinMaxExpVelo,' +/- ',sigMinMaxExpVelo)

systVelo = np.array([(a[0] + a[1]) / 2. for a in velocities])
meanSystVelo = np.mean(systVelo)
sigSystVelo = np.std(systVelo)
print('systVelo = ',meanSystVelo,' +/- ',sigSystVelo)

systMinMaxVelo = np.array([(a[0] + a[1]) / 2. for a in minMaxVelocities])
meanMinMaxSystVelo = np.mean(systMinMaxVelo)
sigMinMaxSystVelo = np.std(systMinMaxVelo)
print('systMinMaxVelo = ',meanMinMaxSystVelo,' +/- ',sigMinMaxSystVelo)

dist = 1832.94812
minDist = 1621.85229
maxDist = 2098.93457
print('dist = ',dist,' +',maxDist-dist,' -',dist-minDist)
majDiam = 193.
minDiam = 134.

physDiam = 2. * dist * np.tan(np.radians((majDiam+minDiam)/(4.*3600.)))
physDiamMax = 2. * dist * np.tan(np.radians((majDiam)/(2.* 3600.)))
physDiamMin = 2. * dist * np.tan(np.radians((minDiam)/(2. * 3600.)))
print('physDiam = ',physDiam,': [',physDiamMin,':',physDiamMax,']')

minDistPhysDiam = 2. * minDist * np.tan(np.radians((majDiam+minDiam)/(4.*3600.)))
minDistPhysDiamMax = 2. * minDist * np.tan(np.radians((majDiam)/(2.* 3600.)))
minDistPhysDiamMin = 2. * minDist * np.tan(np.radians((minDiam)/(2. * 3600.)))
print('minDistPhysDiam = ',minDistPhysDiam,': [',minDistPhysDiamMin,':',minDistPhysDiamMax,']')

maxDistPhysDiam = 2. * maxDist * np.tan(np.radians((majDiam+minDiam)/(4.*3600.)))
maxDistPhysDiamMax = 2. * maxDist * np.tan(np.radians((majDiam)/(2.* 3600.)))
maxDistPhysDiamMin = 2. * maxDist * np.tan(np.radians((minDiam)/(2. * 3600.)))
print('maxDistPhysDiam = ',maxDistPhysDiam,': [',maxDistPhysDiamMin,':',maxDistPhysDiamMax,']')

age = physDiam * 3.086e13 / meanMinMaxExpVelo / (3600.*24*365.25)
print('age = ',age)
minAge = minDistPhysDiamMin * 3.086e13 / (meanMinMaxExpVelo+sigMinMaxExpVelo) / (3600.*24*365.25)
maxAge = maxDistPhysDiamMax * 3.086e13 / (meanMinMaxExpVelo-sigMinMaxExpVelo) / (3600.*24*365.25)
print('minAge = ',minAge,', maxAge = ',maxAge)
