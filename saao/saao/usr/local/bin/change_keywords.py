#!/usr/bin/python

"""A utility for users to change header keywords, and to add new keywords ."""

__version__ = '1.0'

import os
import argparse

from astropy.io import fits
from astropy.time import Time

try:
    input = raw_input
except NameError:
    pass


def get_input(prompt, validator=lambda r: r):
    """Tries to get user input until the validator function succeeds."""
    while True:
        response = input(prompt)
        if validator(response):
            break
    return response


def convert(cube, forward=True):
    """Convert between different versions of FITS keywords."""
    basename, extension = os.path.splitext(cube.filename())
    hdu = cube[0]

    print("Converting '{0}'.".format(cube.filename()))

    for old, new in KEYWORDS:
        try:
            if forward:
                hdu.header.rename_keyword(old, new)
            else:
                hdu.header.rename_keyword(new, old)
        except ValueError as e:
            print(e)

    filename = '{0}_converted.fits'.format(basename)
    while os.path.exists(filename):
        overwrite = get_input(
            "The file '{0}' already exists. Would you like to overwrite it? [y/n] ".format(filename),
            validator=lambda r: r.lower() in ('y', 'n')).lower()

        if overwrite == 'y':
            overwrite = get_input(
                "Are you sure that you want to overwrite '{0}'? [y/n] ".format(filename),
                validator=lambda r: r.lower() in ('y', 'n')).lower()

        if overwrite == 'n':
            filename = get_input('Please enter a new filename: ')
        else:
            break
    hdu.writeto(filename, clobber=True)

    print("Keywords successfully updated and saved to '{0}'.".format(filename))
    print('')

def add_hjd(cube):
    """Add Heliocentric Julian Date """
    basename, extension = os.path.splitext(cube.filename())
    hdu = cube[0]
    print("Converting '{0}'.".format(cube.filename()))

    date = hdu.header['DATE-OBS']
    ut = hdu.header['UT-START']
    datetime = "{}T{}".format(date, ut)
    print("Datetime = {}".format(datetime))
    t = Time(datetime, format='isot', scale='utc')
    hjd = t.mjd + 0.5
    print("HJD-OBS = {}".format(hjd))
    hdu.header["HJD-OBS"] = (hjd, 'The heliocentric Julian date')

    filename = '{0}_converted.fits'.format(basename)
    while os.path.exists(filename):
        overwrite = get_input(
            "The file '{0}' already exists. Would you like to overwrite it? [y/n] ".format(filename),
            validator=lambda r: r.lower() in ('y', 'n')).lower()

        if overwrite == 'y':
            overwrite = get_input(
                "Are you sure that you want to overwrite '{0}'? [y/n] ".format(filename),
                validator=lambda r: r.lower() in ('y', 'n')).lower()

        if overwrite == 'n':
            filename = get_input('Please enter a new filename: ')
        else:
            break
    hdu.writeto(filename, clobber=True)

    print("Keywords successfully updated and saved to '{0}'.".format(filename))
    print('')


def process_arguments():
    parser = argparse.ArgumentParser(description=convert.__doc__)
    parser.add_argument('files', nargs='+', help='Files to convert.')
    parser.add_argument('-r', '--revert', action='store_false', help='Revert the headers of an already converted file.')
    parser.add_argument('keyword_file', nargs='?', help='A file containing keywords to convert.')
    parser.add_argument('-m', '--hjd', action='store_true', help='Add the Heliocentric Julian Date with HJD keyword')
    arguments = parser.parse_args()
    return arguments

if __name__ == '__main__':

    arguments = process_arguments()
    for path in arguments.files:
        try:
            cube = fits.open(path, do_not_scale_image_data=True)
        except IOError as e:
            print(e)
        else:
            if arguments.hjd:
                add_hjd(cube)
                cube.close()
                exit(0)
            if arguments.keyword_file:
                with open(arguments.keyword_file) as kwf:
                    KEYWORDS = []
                    for line in kwf:
                        old, new = line.strip().split(',')
                        KEYWORDS.append((old, new))
            convert(cube, forward=arguments.revert)
            cube.close()
