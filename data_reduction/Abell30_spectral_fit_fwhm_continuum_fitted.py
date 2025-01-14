#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:19:21 2024

@author: kam
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy import modeling

# Open the .fits file
hdul = fits.open('/Users/azuri/daten/uni/HKU/Kamila/Knot_Spectra/spectrumB_j4centre.fits')
print('len(hdul) = ',len(hdul))

# Assuming the spectrum is in the first extension, you can access it like this
spectrum_data = hdul[1].data
print('spectrum_data = ',spectrum_data)
wLen = np.array([x[0] for x in spectrum_data])
flux = np.array([x[1] for x in spectrum_data])

#STOP
# Find the peaks in the spectrum
peaks, _ = find_peaks(flux, prominence=0.1)  # Adjust prominence as needed
print('peaks = ',len(peaks),': ',peaks)
print('flux[peak] = ',flux[peaks])

# Define the composite model with multiple Gaussian components and a linear model for the continuum
composite_init = None
for peak, w in zip(peaks, wLen[peaks]):
    if composite_init is None:
        composite_init = modeling.models.Gaussian1D(amplitude=flux[peak], mean=w, stddev=1)
    else:
        gaussian = modeling.models.Gaussian1D(amplitude=flux[peak], mean=w, stddev=1)
        composite_init += gaussian
composite_init += modeling.models.Polynomial1D(degree=5)  # Linear model for the continuum

# Initialize a fitter
fit_g = modeling.fitting.LevMarLSQFitter()

# Fit the composite model to the data
fitted_model = fit_g(composite_init, wLen, flux)

# Extract the mean wavelength, fitted flux, and FWHM from the individual Gaussian components
mean_wavelength = [component.mean.value for component in fitted_model[:-1]]  # Exclude the last component (linear model)
fitted_flux = [component.amplitude.value for component in fitted_model[:-1]]
fwhm_values = [2.355 * component.stddev.value for component in fitted_model[:-1]]  # Calculate FWHM for each peak

# Plot the spectrum and the fitted composite model
plt.plot(wLen, flux, label='Spectrum')
plt.plot(wLen, fitted_model(wLen), label='Gaussian Fit')  # Plot the fitted composite model

plt.xlabel('Wavelength/Angstrom')
plt.ylabel('Flux')
plt.title('Spectrum with Gaussian Fit')

# Set the x-axis limits to zoom in at a specific range of wavelengths
#plt.xlim(4900, 5010)  # Adjust the wavelength range as needed
#plt.ylim(-10,50)

plt.legend()
plt.show()

# Save the parameters, errors, and mean wavelength to a .txt file
with open('/Users/azuri/daten/uni/HKU/Kamila/fitting_results.txt', 'w') as file:
    file.write('Mean Wavelength (Angstrom)  Fitted Flux  FWHM\n')
    for w, f, fwhm in zip(mean_wavelength, fitted_flux, fwhm_values):
        file.write(f'{w}  {f}  {fwhm}\n')
