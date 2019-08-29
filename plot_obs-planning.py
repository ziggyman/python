# -*- coding: utf-8 -*-
"""
===================================================================
Determining and plotting the altitude/azimuth of a celestial object
===================================================================

This example demonstrates coordinate transformations and the creation of
visibility curves to assist with observing run planning.

In this example, we make a `~astropy.coordinates.SkyCoord` instance for target.
The altitude-azimuth coordinates are then found using
`astropy.coordinates.EarthLocation` and `astropy.time.Time` objects.

This example is meant to demonstrate the capabilities of the
`astropy.coordinates` package. For more convenient and/or complex observation
planning, consider the `astroplan <https://astroplan.readthedocs.org/>`_
package.

-------------------

*By: Erik Tollerud, Kelle Cruz*

*License: BSD*

-------------------

"""

##############################################################################
# Let's suppose you are planning to visit picturesque Bear Mountain State Park
# in New York, USA. You're bringing your telescope with you (of course), and
# someone told you target is a great target to observe there. You happen to know
# you're free at 11:00 pm local time, and you want to know if it will be up.
# Astropy can answer that.
#
# Make print work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:

import numpy as np
import matplotlib.pyplot as plt
from astroplan import Observer
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)


##############################################################################
# Import the packages necessary for finding coordinates and making
# coordinate transformations

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

##############################################################################
# `astropy.coordinates.SkyCoord.from_name` uses Simbad to resolve object
# names and retrieve coordinates.
#
# Get the coordinates of target:

#target = SkyCoord.from_name('target')
#observatoryLocation = EarthLocation(lat=41.3*u.deg, lon=-74*u.deg, height=390*u.m)
#observatoryName: must be resolvable by Observer.at_site("Siding Spring Observatory")
#utcoffset = -4*u.hour
#date = '2019-09-05'

def plot_target(target, observatoryLocation, observatoryName, utcoffset, date):
    midnight = Time(date+' 00:00:00') - utcoffset
    observatory = observatoryLocation
    obs = Observer.at_site(observatoryName)#, timezone='Eastern Standard Time')
    observationStartTime = obs.twilight_evening_astronomical(midnight)#Time('2019-9-5 19:10:00') - utcoffset
    observationEndTime = obs.twilight_morning_astronomical(midnight)#Time('2019-9-6 4:55:00') - utcoffset
    observationStartTime.format = 'iso'
    observationEndTime.format = 'iso'
    ##############################################################################
    # Use `astropy.coordinates.EarthLocation` to provide the location of Bear
    # Mountain and set the time to 11pm EDT on 2012 July 12:


    ##############################################################################
    # `astropy.coordinates.EarthLocation.get_site_names` and
    # `~astropy.coordinates.EarthLocation.get_site_names` can be used to get
    # locations of major observatories.
    #
    # Use `astropy.coordinates` to find the Alt, Az coordinates of target at as
    # observed from Bear Mountain at 11pm on 2012 July 12.

    #time = Time('2012-7-12 23:00:00') - utcoffset
    #targetaltaz = target.transform_to(AltAz(obstime=time,location=observatory))
    #print("target's Altitude = {0.alt:.2}".format(targetaltaz))

    ##############################################################################
    # This is helpful since it turns out target is barely above the horizon at this
    # time. It's more informative to find target's airmass over the course of
    # the night.
    #
    # Find the alt,az coordinates of target at 100 times evenly spaced between 10pm
    # and 7am EDT:

    time = observationStartTime
    times = [time]
    while time < observationEndTime:
        time += 1*u.minute
        print('time = ',time)
        times.append(time)
    delta_midnight = np.linspace(-10, 10, 100)*u.hour
    frameNight = AltAz(obstime=midnight+delta_midnight,
                              location=observatory)
    targetaltazsNight = target.transform_to(frameNight)

    ##############################################################################
    # convert alt, az to airmass with `~astropy.coordinates.AltAz.secz` attribute:

    targetairmasssNight = targetaltazsNight.secz

    ##############################################################################
    # Plot the airmass as a function of time:

    plt.plot(delta_midnight, targetairmasssNight)
    plt.xlim(-2, 10)
    plt.ylim(1, 4)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Airmass [Sec(z)]')
    plt.show()

    ##############################################################################
    # Use  `~astropy.coordinates.get_sun` to find the location of the Sun at 1000
    # evenly spaced times between noon on July 12 and noon on July 13:

    from astropy.coordinates import get_sun
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=observatory)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)


    ##############################################################################
    # Do the same with `~astropy.coordinates.get_moon` to find when the moon is
    # up. Be aware that this will need to download a 10MB file from the internet
    # to get a precise location of the moon.

    from astropy.coordinates import get_moon
    moon_July12_to_13 = get_moon(times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)

    ##############################################################################
    # Find the alt,az coordinates of target at those same times:

    targetaltazs_July12_to_13 = target.transform_to(frame_July12_to_13)

    ##############################################################################
    # Make a beautiful figure illustrating nighttime and the altitudes of target and
    # the Sun over that time:

    plt.plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
    plt.scatter(delta_midnight, targetaltazs_July12_to_13.alt,
                c=targetaltazs_July12_to_13.az, label='target', lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_July12_to_13.alt < -18*u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(-12, 12)
    plt.xticks(np.arange(13)*2 -12)
    plt.ylim(0, 90)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.show()

observatoryName = "Siding Spring Observatory"
observatoryLocation = EarthLocation(lat=-31.2749*u.deg, lon=149.0685*u.deg, height=1165*u.m)
utcoffset = 10*u.hour
date = '2019-09-05'
targetCoord = SkyCoord(ra=293.95888*u.deg, dec=16.47986*u.deg, frame='icrs')

plot_target(targetCoord, observatoryLocation, observatoryName, utcoffset, date)
