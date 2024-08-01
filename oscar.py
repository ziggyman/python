import pickle

# write list to binary file
def store_list(list_name,file_name):
    # store list in binary file so 'wb' mode
    with open(file_name, "wb") as fp:
        pickle.dump(list_name, fp)

# read list from pickle
def read_list(file_name):
    with open(file_name, "rb") as fp:
        return pickle.load(fp)

# read the data

data_table = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/data table combined")
object = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/object combined")

# gaussian fit
if False:
#def gauss_fit(x, y):
    from scipy.optimize import curve_fit
    '''
    def Gauss(x, A, mu, sigma, background):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background if background is not None else return A * np.exp(- (x - mu)**2 / (2 * sigma**2))
    '''
    def Gauss(x, A, mu, sigma, background):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background

    # Initial guess for parameters A, mu, sigma
    initial_guess = [max(y), np.mean(x), np.std(x), 0]

    parameters, covariance = curve_fit(Gauss, x, y, p0=initial_guess, maxfev=5000)
    amp, mean, stddev, background = parameters
    return amp, mean, stddev, background

# gaussian fit
def gauss_fit(x, y, background = None):
    from scipy.optimize import curve_fit
    '''
    def Gauss(x, A, mu, sigma, background):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background if background is not None else return A * np.exp(- (x - mu)**2 / (2 * sigma**2))
    '''
    def Gauss(x, A, mu, sigma, background=None):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background if background is not None else A * np.exp(- (x - mu)**2 / (2 * sigma**2))

    # Initial guess for parameters A, mu, sigma
    if background is None:
        initial_guess = [max(y), np.mean(x), np.std(x)]
    else:
        initial_guess = [max(y), np.mean(x), np.std(x), background]

    parameters, covariance = curve_fit(Gauss, x, y, p0=initial_guess, maxfev=5000)
    #amp, mean, stddev, background = parameters
    return parameters#amp, mean, stddev, background

def gauss_fit_min(x, y):
    from scipy.optimize import curve_fit
    '''
    def Gauss(x, A, mu, sigma, background):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background if background is not None else return A * np.exp(- (x - mu)**2 / (2 * sigma**2))
    '''
    def Gauss(x, A, mu, sigma, background):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background

    # Initial guess for parameters A, mu, sigma
    initial_guess = [np.min(y) - np.max(y), np.mean(x), np.std(x), np.max(y)]
    print('initial guess: amplitude = ',initial_guess[0])
    print('initial guess: center = ',initial_guess[1])
    print('initial guess: sigma = ',initial_guess[2])
    print('initial guess: background = ',initial_guess[3])

    parameters, covariance = curve_fit(Gauss, x, y, p0=initial_guess, maxfev=5000)
    amp, mean, stddev, background = parameters
    print('fitted values: ',amp,mean,stddev,background)
    return amp, mean, stddev, background

# generating the gaussian function

def gauss_gen(A, mu, sigma, x, background=None):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2)) + background if background is not None else A * np.exp(- (x - mu)**2 / (2 * sigma**2))


# fit gaussian

def gauss_fit_gen(data_table,wLen_test):

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    '''

    i = 10 # object
    j = 25 # row

    wLen_test = 6562.8



    # Ensure data_table is properly defined
    temp_1 = np.abs(data_table[i][0] - (wLen_test - 10))
    temp_2 = np.abs(data_table[i][0] - (wLen_test + 10))

    low_index = np.argmin(temp_1)
    high_index = np.argmin(temp_2)


    '''

    wLen_obs = data_table[i][0][low_index:high_index]
    flux_obs = data_table[i][j][low_index:high_index]

    amp, mean, stddev = gauss_fit(wLen_obs, flux_obs, background = None)

    print("Amplitude:", amp)
    print("Mean:", mean)
    print("Standard Deviation:", stddev)

    x_fit = np.linspace(min(wLen_obs), max(wLen_obs), 1000)
    flux_fit = gauss_gen(amp, mean, stddev, x_fit, background=None)

    flux_lab = gauss_gen(amp,wLen_test,stddev,x_fit)

    plt.plot(wLen_obs, flux_obs, label='Real Data')
    plt.plot(x_fit, flux_fit, label='Fitted Gaussian', linestyle='--')
    plt.plot(x_fit, flux_lab, label="Lab Gaussian", linestyle="-.")
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

    return (amp,mean,stddev,wLen_obs)


# doppler correction of wavelength where wavelength is an list

def dop_correction (speed, wavelength):
    c = 3e8
    return [(1 + (speed/c)) * element for element in wavelength]

def chi_sq(spec_obs, wLen_obs, amp, mean, stddev):

    import numpy as np
    import matplotlib.pyplot as plt

    total = []
    c = 3e8  # speed of light in m/s

    # Loop over speed range in km/s
    for i in range(-100, 101):
        speed = i * 1e3  # convert to m/s
        dop_wLen = dop_correction(speed, wLen_obs)
        dop_mean = (1 + speed / c) * mean
        slit_real = gauss_gen(amp, dop_mean, stddev, wLen_obs)#,background)

        # Interpolate the observed spectrum to align with Doppler-shifted wavelengths
        slit_obs_interp = np.interp(dop_wLen, wLen_obs, spec_obs[low_index:high_index])

        # Calculate chi-square value
        chi_square = np.sum(((slit_obs_interp - slit_real) ** 2) )
        total.append(chi_square)
        v_range = np.linspace(-100, 100, 201)

    # Plotting for visualization
    '''
    plt_index = np.argmin(total)
    plt_high_index = plt_index + 40
    plt_low_index = plt_index - 40
    '''

    # Fit a Gaussian to the chi-square values to find the minimum
    amp_fit, mean_fit, stddev_fit,background = gauss_fit_min(v_range, total)
    fitting_y = gauss_gen(amp_fit, mean_fit, stddev_fit, v_range,background)
    index_fit = np.argmin(fitting_y)
    speed = v_range[index_fit]

    plt.plot(v_range, total, label='Chi-square values')
    plt.plot(v_range, fitting_y, label='Fitted Gaussian', linestyle='--')
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Chi-square')
    plt.legend()
    plt.show()

    return speed

# overall program of finding the radio velocity


import pandas as pd
import numpy as np

df = pd.read_csv ("/Users/azuri/daten/uni/HKU/interns_projects/oscar/element_wLen.csv")

element = df["Element"].values
wLen_ele = df["wLen"].values

radial_vel = np.empty([len(data_table),43,len(wLen_ele)])       # def the storing variable where x is the element, y is the rows, z is the object
'''

for i in range(len(data_table)):  # object
    for j in range(len(data_table[i])):   # row

'''
i = 1
j = 1


for m in range (len(wLen_ele)):

    wLen_test = wLen_ele[m]

    temp_1 = np.abs(data_table[i][0] - (wLen_test - 10))
    temp_2 = np.abs(data_table[i][0] - (wLen_test + 10))

    low_index = np.argmin(temp_1)
    high_index = np.argmin(temp_2)

    # checking if there is any signal in the range

    # finding the noise
    # noise range
    noise_data = data_table[i][j][0:600]
    noise = np.std(noise_data)

    found = True
    n = low_index
    while n < (high_index + 1) and found:

        if data_table[i][j][n] > (10 * noise):
            found = False
        else:
            n = n + 1

    if found == False:
        amp,mean,stddev,wLen_obs = gauss_fit_gen(data_table,wLen_test)
        speed = chi_sq(data_table[i][j],wLen_obs,amp,wLen_ele[m],stddev)

    radial_vel[i,j,m]
