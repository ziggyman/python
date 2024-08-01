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


1. Import the required libraries in your Python script:

<pre style="background-color: rgb(43, 43, 43);margin-right: 15px;"><div class="pre-code-area"><code class="language-javascript" style="white-space: pre-wrap;">python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import convolve
</code></div></pre>

2. Open the FITS file, extract the data, create the wavelength array, and perform the FFT, as shown in the previous examples.

3. Create the top-hat function and convolve it with the power spectrum:

<pre style="background-color: rgb(43, 43, 43);margin-right: 15px;"><div class="pre-code-area"><code class="language-javascript" style="white-space: pre-wrap;">python
# Define the width of the top-hat function
top_hat_width = 0.5

# Calculate the number of points for the top-hat function
num_top_hat_points = int(top_hat_width * sampling_rate)

# Create the top-hat function
top_hat = np.ones(num_top_hat_points) / num_top_hat_points

# Convolve the power spectrum with the top-hat function
convolved_power_spectrum = convolve(power_spectrum_shifted, top_hat, mode='same')
</code></div></pre>

4. Perform the inverse FFT on the convolved power spectrum:

<pre style="background-color: rgb(43, 43, 43);margin-right: 15px;"><div class="pre-code-area"><code class="language-javascript" style="white-space: pre-wrap;">python
# Shift the zero-frequency component back to the edges
convolved_power_spectrum_unshifted = np.fft.ifftshift(convolved_power_spectrum)

# Perform the inverse FFT
convolved_spectrum_fft_inv = np.fft.ifft(convolved_power_spectrum_unshifted)
</code></div></pre>

5. Plot the original spectrum and the convolved spectrum:

<pre style="background-color: rgb(43, 43, 43);margin-right: 15px;"><div class="pre-code-area"><code class="language-javascript" style="white-space: pre-wrap;">python
# Plot the original spectrum
plt.figure()
plt.plot(wavelength, data, label='Original Spectrum')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux')
plt.title('Spectrum')

# Plot the convolved spectrum
plt.plot(wavelength, np.abs(convolved_spectrum_fft_inv), label='Convolved Spectrum')
plt.legend()
plt.show()
</code></div></pre>

This will display the plot of the original spectrum and the convolved spectrum after applying the top-hat function and performing the inverse Fast Fourier Transform. Note that we used the absolute value of the inverse FFT result, as the output may have a small imaginary component due to numerical errors in the FFT computation.
