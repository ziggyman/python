import ephem
from astroplan import Observer
#from astroplan import download_IERS_A
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

#longitude = 20.8105
#utcoffset = 2*u.hour
#localTime = Time('2019-9-7 22:00:00')
def ltToLST(longitude, utcoffset, localTime):
    #Sidereal Time and Julian Date Calculator
    #Revision history: Justine Haupt, v1.0 (11/23/17)
    UTC = localTime - utcoffset
    print('UTC = ',UTC)

    #Only valid for dates between 1901 and 2099. Accurate to within 1.1s.

    #References:
    #http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    #http://aa.usno.navy.mil/faq/docs/GAST.php and

    lon = longitude
    #Calculate longitude in DegHHMM format for edification of user:
    lonMin = (lon - int(lon))*60
    lonSec = (lonMin - int(lonMin))*60
    lonMin = int(lonMin)
    lonSec = int(lonSec)

    #split TD into individual variables for month, day, etc. and convert to floats:
    MM = float(UTC.strftime("%m"))
    DD = float(UTC.strftime("%d"))
    YY = float(UTC.strftime("%y"))
    YY = YY+2000
    hh = float(UTC.strftime("%H"))
    mm = float(UTC.strftime("%M"))

    #convert mm to fractional time:
    mm = mm/60

    #reformat UTC time as fractional hours:
    UT = hh+mm

    #calculate the Julian date:
    JD = (367*YY) - int((7*(YY+int((MM+9)/12)))/4) + int((275*MM)/9) + DD + 1721013.5 + (UT/24)
    print('\nJulian Date: JD%s' %(JD))

    #calculate the Greenwhich mean sidereal time:
    GMST = 18.697374558 + 24.06570982441908*(JD - 2451545)
    GMST = GMST % 24    #use modulo operator to convert to 24 hours
    GMSTmm = (GMST - int(GMST))*60          #convert fraction hours to minutes
    GMSTss = (GMSTmm - int(GMSTmm))*60      #convert fractional minutes to seconds
    GMSThh = int(GMST)
    GMSTmm = int(GMSTmm)
    GMSTss = int(GMSTss)
    print('\nGreenwhich Mean Sidereal Time: %s:%s:%s' %(GMSThh, GMSTmm, GMSTss))

    #Convert to the local sidereal time by adding the longitude (in hours) from the GMST.
    #(Hours = Degrees/15, Degrees = Hours*15)
    lon = lon/15      #Convert longitude to hours
    LST = GMST+lon     #Fraction LST. If negative we want to add 24...
    print('LST = ',LST)
    if LST < 0:
        LST = LST +24
        print('a LST = ',LST)
    if LST >= 24:
        LST = LST - 24
        print('b LST = ',LST)
    LSTmm = (LST - int(LST))*60          #convert fraction hours to minutes
    LSTss = (LSTmm - int(LSTmm))*60      #convert fractional minutes to seconds
    LSThh = int(LST)
    LSTmm = int(LSTmm)
    LSTss = int(LSTss)

    print('\nLocal Sidereal Time %s:%s:%s \n\n' %(LSThh, LSTmm, LSTss))
    year = UTC.strftime("%Y")
    month = UTC.strftime("%m")
    day = UTC.strftime("%d")
    timeStr = '%s-%s-%s %s:%s:%s' % (year, month, day, LSThh, LSTmm, LSTss)
    print('timeStr = <'+timeStr+'>')
    returnValue = Time(timeStr)
    return returnValue

#longitude = 20.8105
#utcoffset = 2*u.hour
#localTime = Time('2019-9-7 22:00:00')

#print(ltToLST(longitude, utcoffset, localTime))