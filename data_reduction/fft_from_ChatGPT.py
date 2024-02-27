import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Replace 'your_spectrum.fits' with the path to your FITS file
file_path = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-06/SCIENCE_BMP1457-5413_dbs00051r_otzxfifEcdF.fits'
with fits.open(file_path) as hdulist:
    # Get the header and data from the primary HDU (Header Data Unit)
    header = hdulist[0].header
    data = hdulist[0].data

# Get the number of pixels, starting wavelength, and dispersion from the header
n_pixels = header['NAXIS1']
start_wavelength = header['CRVAL1']
dispersion = header['CDELT1']

# Create an array of wavelength values
wavelength = np.arange(n_pixels) * dispersion + start_wavelength

plt.figure()
plt.plot(wavelength, data)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux')
plt.title('Spectrum')
plt.show()

# Perform the FFT
spectrum_fft = np.fft.fft(data)

# Calculate the frequencies corresponding to the FFT
sampling_rate = 1.0 / dispersion
fft_freqs = np.fft.fftfreq(n_pixels, d=1.0 / sampling_rate)

# Calculate the power spectrum (magnitude squared)
power_spectrum = np.abs(spectrum_fft)**2

# Shift the zero-frequency component to the center of the spectrum
power_spectrum_shifted = np.fft.fftshift(power_spectrum)
fft_freqs_shifted = np.fft.fftshift(fft_freqs)

# Plot the power spectrum
plt.figure()
plt.plot(fft_freqs_shifted, power_spectrum_shifted)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power Spectrum')
plt.show()
