from astropy.coordinates import SkyCoord
import astropy.units as u

inputFileName = 'test.txt'
#NameObjet HH:MM:SS.ss DD:MM:SS.ss

outputFileName = 'testOut.txt'
#NameObjet dd.ddddd dd.ddddd

def raDecToLonLat(ra, dec):
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    return [float(c.galactic.l.value),float(c.galactic.b.value)]

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip().split(' ') for line in lines]
    return linesOut

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    h, m, s = [float(i) for i in string.split(':')]
    return (s / 240.) + (m / 4.) + (h * 15.)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    d, m, s = [float(i) for i in string.split(':')]
    if d < 0.:
        d = 0. - d
        return 0. - (s / 3600. + m / 60. + d)
    return s / 3600. + m / 60. + d

if __name__ == '__main__':
    inputLines = readFileToArr(inputFileName)
    with open(outputFileName,'w') as f:
        for line in inputLines:
            print(line)
            l, b = raDecToLonLat(hmsToDeg(line[1]), dmsToDeg(line[2]))
            f.write('%s %.5f %.5f\n' % (line[0],l,b))
