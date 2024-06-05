import math

def calculate_limiting_magnitude1(aperture_diameter, pixel_size, quantum_efficiency, read_noise, dark_current, gain):
    #h = 6.62607015 * 10**-34  # Planck's constant (Js)
    #c = 299792458  # Speed of light (m/s)
    exposure_time = 100  # Exposure time in seconds

    # Calculate the total number of photons collected
    aperture_radius = aperture_diameter / 2
    aperture_area = math.pi * (aperture_radius**2)
    pixel_area = pixel_size**2
    total_photons = aperture_area * quantum_efficiency * exposure_time * pixel_area * c / h

    # Calculate the signal and noise components
    signal = total_photons
    dark_current_noise = dark_current * exposure_time
    read_noise_squared = read_noise**2
    shot_noise_squared = signal
    noise_squared = dark_current_noise + read_noise_squared + shot_noise_squared

    # Calculate the signal-to-noise ratio (SNR)
    snr = signal / math.sqrt(noise_squared)

    # Calculate the limiting magnitude
    limiting_magnitude = -2.5 * math.log10(snr)

    return limiting_magnitude

# Example usage
aperture_diameter = 80  # Aperture diameter in mm
pixel_size = 4.5  # Pixel size in micrometers
quantum_efficiency = 0.8  # Quantum efficiency (between 0 and 1)
read_noise = 10  # Read noise in electrons
dark_current = 0.1  # Dark current in electrons per second
gain = 1  # Gain in electrons per ADU

limiting_mag = calculate_limiting_magnitude1(aperture_diameter, pixel_size, quantum_efficiency, read_noise, dark_current, gain)
print("Limiting Magnitude:", limiting_mag)

def calculate_limiting_magnitude2(focal_length, entrance_pupil_diameter, f_number, spectral_range, optical_design_mtf, distortion, average_transmission, window_diameter, maximum_field_of_view, sunshade_length, airy_disk_diameter):
    exposure_time = 100  # Exposure time in seconds

    # Calculate the effective aperture diameter
    effective_aperture_diameter = entrance_pupil_diameter / f_number

    # Calculate the total number of photons collected
    field_of_view_area = math.radians(maximum_field_of_view) * focal_length**2 / f_number
    pixel_area = airy_disk_diameter**2
    total_photons = field_of_view_area * spectral_range * optical_design_mtf * distortion * average_transmission * exposure_time * pixel_area

    # Calculate the signal and noise components
    signal = total_photons
    noise_squared = signal

    # Calculate the signal-to-noise ratio (SNR)
    snr = signal / math.sqrt(noise_squared)

    # Calculate the limiting magnitude
    limiting_magnitude = -2.5 * math.log10(snr)

    return limiting_magnitude

# Camera parameters
focal_length = 65.3  # Focal length in mm
entrance_pupil_diameter = 14  # Entrance pupil diameter in mm
f_number = 4.66
spectral_range = 800 - 450  # Spectral range in nm
optical_design_mtf = 0.35  # MTF at 110 lp/mm
distortion = 0.1  # Distortion in percentage
average_transmission = 0.75
window_diameter = 34  # Window diameter in mm
maximum_field_of_view = 15  # Maximum field of view in degrees
sunshade_length = 50  # Sunshade length in mm
airy_disk_diameter = 9  # Airy disk diameter in micrometers

limiting_mag = calculate_limiting_magnitude2(focal_length, entrance_pupil_diameter, f_number, spectral_range, optical_design_mtf, distortion, average_transmission, window_diameter, maximum_field_of_view, sunshade_length, airy_disk_diameter)
print("Limiting Magnitude:", limiting_mag)
