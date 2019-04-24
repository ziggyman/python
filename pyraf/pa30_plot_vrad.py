#from myUtils import getDate, findClosestDate,...
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import re

#from myUtils import getWavelength

speedOfLight = 299792.458# km/s

tab = [
{'column': 1135, 'centerDistance' : -85.12, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6724.81, 'SII6731b' : 6737.28, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1144, 'centerDistance' : -82.84, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6722.17, 'SII6731b' : 6734.94, 'unidentified1' : 0.0, 'ArIII7136' : 7127.61, 'unidentified2' : 7140.67},
{'column': 1175, 'centerDistance' : -74.95, 'SII6716a': 6708.83, 'SII6731a' : 0.0, 'SII6716b' : 6719.52, 'SII6731b' : 6732.37, 'unidentified1' : 7007.79, 'ArIII7136' : 7119.82, 'unidentified2' : 7137.82},
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': 1323, 'centerDistance' : -37.31, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6736.7, 'SII6731b' : 0.0, 'unidentified1' : 6999.17, 'ArIII7136' : 7120.80, 'unidentified2' : 7160.06},
{'column': 1324, 'centerDistance' : -37.06, 'SII6716a': 6697.93, 'SII6731a' : 6712.08, 'SII6716b' : 6738.15, 'SII6731b' : 0.0, 'unidentified1' : 7001.46, 'ArIII7136' : 7115.85, 'unidentified2' : 7160.63},
{'column': 1325, 'centerDistance' : -36.79, 'SII6716a': 6698.82, 'SII6731a' : 6711.64, 'SII6716b' : 6736.49, 'SII6731b' : 0.0, 'unidentified1' : 6999.32, 'ArIII7136' : 7115.25, 'unidentified2' : 7160.82},
{'column': 1326, 'centerDistance' : -36.54, 'SII6716a': 6699.99, 'SII6731a' : 6713.83, 'SII6716b' : 6736.79, 'SII6731b' : 6753.37, 'unidentified1' : 7000.49, 'ArIII7136' : 0.0, 'unidentified2' : 7158.87},
{'column': 1327, 'centerDistance' : -36.28, 'SII6716a': 6700.85, 'SII6731a' : 6711.92, 'SII6716b' : 6738.88, 'SII6731b' : 6753.73, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1332, 'centerDistance' : -35.02, 'SII6716a': 6703.25, 'SII6731a' : 6711.63, 'SII6716b' : 6739.34, 'SII6731b' : 6752.67, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1357, 'centerDistance' : -28.65, 'SII6716a': 6696.05, 'SII6731a' : 6713.10, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': 1387, 'centerDistance' : -21.02, 'SII6716a': 6696.69, 'SII6731a' : 6711.05, 'SII6716b' : 6735.63, 'SII6731b' : 6749.30, 'unidentified1' : 6995.58, 'ArIII7136' : 7115.21, 'unidentified2' : 0},
{'column': 1388, 'centerDistance' : -20.77, 'SII6716a': 6697.10, 'SII6731a' : 6710.98, 'SII6716b' : 6735.32, 'SII6731b' : 6750.59, 'unidentified1' : 6997.59, 'ArIII7136' : 7114.89, 'unidentified2' : 0},
{'column': 1389, 'centerDistance' : -20.51, 'SII6716a': 6697.70, 'SII6731a' : 6711.08, 'SII6716b' : 6736.75, 'SII6731b' : 6750.92, 'unidentified1' : 6997.14, 'ArIII7136' : 7114.21, 'unidentified2' : 0},
{'column': 1390, 'centerDistance' : -20.26, 'SII6716a': 6697.63, 'SII6731a' : 6711.51, 'SII6716b' : 6736.30, 'SII6731b' : 6749.04, 'unidentified1' : 6996.34, 'ArIII7136' : 7114.80, 'unidentified2' : 0},
{'column': 1391, 'centerDistance' : -20.01, 'SII6716a': 6697.36, 'SII6731a' : 6712.50, 'SII6716b' : 6733.98, 'SII6731b' : 6749.16, 'unidentified1' : 6994.28, 'ArIII7136' : 7115.17, 'unidentified2' : 0},
{'column': 1392, 'centerDistance' : -19.75, 'SII6716a': 6697.96, 'SII6731a' : 6712.36, 'SII6716b' : 6734.90, 'SII6731b' : 6748.21, 'unidentified1' : 6994.75, 'ArIII7136' : 7114.58, 'unidentified2' : 7155.94},
{'column': 1393, 'centerDistance' : -19.50, 'SII6716a': 6697.19, 'SII6731a' : 6712.32, 'SII6716b' : 6734.62, 'SII6731b' : 6748.58, 'unidentified1' : 6996.99, 'ArIII7136' : 7115.46, 'unidentified2' : 7155.13},
{'column': 1394, 'centerDistance' : -19.25, 'SII6716a': 6697.36, 'SII6731a' : 6713.20, 'SII6716b' : 6733.76, 'SII6731b' : 6746.37, 'unidentified1' : 6995.50, 'ArIII7136' : 7116.50, 'unidentified2' : 7154.88},
{'column': 1395, 'centerDistance' : -18.99, 'SII6716a': 6698.82, 'SII6731a' : 6714.39, 'SII6716b' : 6733.44, 'SII6731b' : 6746.61, 'unidentified1' : 6995.44, 'ArIII7136' : 7114.62, 'unidentified2' : 7154.55},
{'column': 1396, 'centerDistance' : -18.74, 'SII6716a': 6699.46, 'SII6731a' : 6713.44, 'SII6716b' : 6733.95, 'SII6731b' : 6747.31, 'unidentified1' : 6998.17, 'ArIII7136' : 7114.56, 'unidentified2' : 7157.68},
{'column': 1397, 'centerDistance' : -18.48, 'SII6716a': 6699.42, 'SII6731a' : 6713.19, 'SII6716b' : 6733.49, 'SII6731b' : 6746.34, 'unidentified1' : 7001.56, 'ArIII7136' : 7116.51, 'unidentified2' : 0},
{'column': 1398, 'centerDistance' : -18.23, 'SII6716a': 6700.47, 'SII6731a' : 6711.84, 'SII6716b' : 6732.95, 'SII6731b' : 6747.07, 'unidentified1' : 6997.28, 'ArIII7136' : 7114.01, 'unidentified2' : 7152.35},
{'column': 1399, 'centerDistance' : -17.98, 'SII6716a': 6700.86, 'SII6731a' : 6716.18, 'SII6716b' : 6732.18, 'SII6731b' : 6746.45, 'unidentified1' : 6997.80, 'ArIII7136' : 7113.63, 'unidentified2' : 7154.34},
{'column': 1400, 'centerDistance' : -17.72, 'SII6716a': 6700.96, 'SII6731a' : 6716.04, 'SII6716b' : 6732.86, 'SII6731b' : 6745.74, 'unidentified1' : 6994.81, 'ArIII7136' : 0.0, 'unidentified2' : 7152.56},
{'column': 1416, 'centerDistance' : -13.65, 'SII6716a': 6696.77, 'SII6731a' : 6711.62, 'SII6716b' : 6735.77, 'SII6731b' : 6744.59, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1434, 'centerDistance' : -9.07, 'SII6716a': 6697.23, 'SII6731a' : 6712.33, 'SII6716b' : 6737.51, 'SII6731b' : 6743.87, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1476, 'centerDistance' : 1.61, 'SII6716a': 6701.552, 'SII6731a' : 6709.95, 'SII6716b' : 6739.47, 'SII6731b' : 6750.99, 'unidentified1' : 0.0, 'ArIII7136' : 7111.97, 'unidentified2' : 0.0},
{'column': 1486, 'centerDistance' : 4.15, 'SII6716a': 6698.65, 'SII6731a' : 6711.70, 'SII6716b' : 6737.38, 'SII6731b' : 6753.86, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1491, 'centerDistance' : 5.43, 'SII6716a': 6698.46, 'SII6731a' : 6711.77, 'SII6716b' : 6739.16, 'SII6731b' : 6754.07, 'unidentified1' : 0.0, 'ArIII7136' : 7110.81, 'unidentified2' : 0.0},
{'column': 1501, 'centerDistance' : 7.97, 'SII6716a': 6696.42, 'SII6731a' : 6710.25, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7116.58, 'unidentified2' : 0.0},
{'column': 1511, 'centerDistance' : 10.51, 'SII6716a': 6697.56, 'SII6731a' : 6711.89, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7110.14, 'unidentified2' : 0.0},
{'column': 1521, 'centerDistance' : 13.05, 'SII6716a': 6697.23, 'SII6731a' : 6713.88, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1531, 'centerDistance' : 15.60, 'SII6716a': 6699.35, 'SII6731a' : 6713.37, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7116.99, 'unidentified2' : 0.0},
{'column': 1541, 'centerDistance' : 18.14, 'SII6716a': 6697.23, 'SII6731a' : 6711.86, 'SII6716b' : 6730.73, 'SII6731b' : 6745.77, 'unidentified1' : 0.0, 'ArIII7136' : 7111.45, 'unidentified2' : 0.0},
{'column': 1551, 'centerDistance' : 20.68, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6734.87, 'SII6731b' : 6747.23, 'unidentified1' : 0.0, 'ArIII7136' : 7116.50, 'unidentified2' : 0.0},
{'column': 1561, 'centerDistance' : 23.23, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6735.39, 'SII6731b' : 6748.11, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1571, 'centerDistance' : 25.77, 'SII6716a': 6698.32, 'SII6731a' : 6713.72, 'SII6716b' : 6739.63, 'SII6731b' : 6750.23, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1581, 'centerDistance' : 28.32, 'SII6716a': 6698.58, 'SII6731a' : 6712.00, 'SII6716b' : 6743.03, 'SII6731b' : 6754.47, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1591, 'centerDistance' : 30.86, 'SII6716a': 6698.72, 'SII6731a' : 6712.02, 'SII6716b' : 6743.41, 'SII6731b' : 6756.56, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1601, 'centerDistance' : 33.40, 'SII6716a': 6699.62, 'SII6731a' : 6711.57, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1611, 'centerDistance' : 35.95, 'SII6716a': 6699.85, 'SII6731a' : 6710.85, 'SII6716b' : 6740.03, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7116.70, 'unidentified2' : 0.0},
{'column': 1621, 'centerDistance' : 38.49, 'SII6716a': 6697.23, 'SII6731a' : 6710.71, 'SII6716b' : 6736.62, 'SII6731b' : 6748.11, 'unidentified1' : 0.0, 'ArIII7136' : 7114.45, 'unidentified2' : 0.0},
{'column': 1631, 'centerDistance' : 41.04, 'SII6716a': 6697.23, 'SII6731a' : 6712.07, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7112.87, 'unidentified2' : 0.0},
{'column': 1641, 'centerDistance' : 43.58, 'SII6716a': 6701.47, 'SII6731a' : 6712.07, 'SII6716b' : 6734.85, 'SII6731b' : 6748.40, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1651, 'centerDistance' : 46.12, 'SII6716a': 6699.04, 'SII6731a' : 6712.07, 'SII6716b' : 6730.49, 'SII6731b' : 6745.99, 'unidentified1' : 0.0, 'ArIII7136' : 7116.99, 'unidentified2' : 0.0},
{'column': 1661, 'centerDistance' : 48.67, 'SII6716a': 6698.51, 'SII6731a' : 6711.99, 'SII6716b' : 6731.58, 'SII6731b' : 6745.67, 'unidentified1' : 0.0, 'ArIII7136' : 7116.19, 'unidentified2' : 0.0},
{'column': 1671, 'centerDistance' : 51.21, 'SII6716a': 6697.60, 'SII6731a' : 6711.52, 'SII6716b' : 6732.70, 'SII6731b' : 6746.80, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1681, 'centerDistance' : 53.76, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6732.81, 'SII6731b' : 6745.74, 'unidentified1' : 0.0, 'ArIII7136' : 7113.97, 'unidentified2' : 0.0},
{'column': 1691, 'centerDistance' : 56.30, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6724.79, 'SII6731b' : 6739.63, 'unidentified1' : 0.0, 'ArIII7136' : 7112.75, 'unidentified2' : 0.0},
{'column': 1701, 'centerDistance' : 58.84, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6717.14, 'SII6731b' : 6736.10, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1711, 'centerDistance' : 61.39, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6719.39, 'SII6731b' : 6735.39, 'unidentified1' : 0.0, 'ArIII7136' : 7112.09, 'unidentified2' : 0.0},
{'column': 1721, 'centerDistance' : 63.93, 'SII6716a': 6703.04, 'SII6731a' : 6711.05, 'SII6716b' : 6718.29, 'SII6731b' : 6735.05, 'unidentified1' : 0.0, 'ArIII7136' : 7111.55, 'unidentified2' : 0.0},
{'column': 1731, 'centerDistance' : 66.48, 'SII6716a': 6704.06, 'SII6731a' : 6712.07, 'SII6716b' : 6718.31, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7113.60, 'unidentified2' : 0.0},
{'column': 1741, 'centerDistance' : 69.02, 'SII6716a': 6699.44, 'SII6731a' : 6713.11, 'SII6716b' : 6726.47, 'SII6731b' : 6740.91, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1751, 'centerDistance' : 71.56, 'SII6716a': 6701.43, 'SII6731a' : 6716.31, 'SII6716b' : 6726.65, 'SII6731b' : 6740.05, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1761, 'centerDistance' : 74.11, 'SII6716a': 6713.54, 'SII6731a' : 6729.03, 'SII6716b' : 6731.15, 'SII6731b' : 6743.87, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1771, 'centerDistance' : 76.65, 'SII6716a': 6703.59, 'SII6731a' : 6714.19, 'SII6716b' : 6729.03, 'SII6731b' : 6741.75, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1781, 'centerDistance' : 79.19, 'SII6716a': 6714.19, 'SII6731a' : 6729.03, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1791, 'centerDistance' : 81.73, 'SII6716a': 6712.19, 'SII6731a' : 0.0, 'SII6716b' : 6721.95, 'SII6731b' : 6734.55, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1801, 'centerDistance' : 84.27, 'SII6716a': 6712.07, 'SII6731a' : 0.0, 'SII6716b' : 6723.62, 'SII6731b' : 6734.95, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1811, 'centerDistance' : 86.82, 'SII6716a': 6711.56, 'SII6731a' : 6726.88, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7136.44, 'unidentified2' : 0.0},
{'column': 1821, 'centerDistance' : 89.36, 'SII6716a': 6712.39, 'SII6731a' : 6729.49, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7135.08, 'unidentified2' : 0.0},
{'column': 1831, 'centerDistance' : 91.91, 'SII6716a': 6712.05, 'SII6731a' : 6725.03, 'SII6716b' : 6720.55, 'SII6731b' : 6733.52, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1841, 'centerDistance' : 94.45, 'SII6716a': 6713.51, 'SII6731a' : 6729.03, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7134.07, 'unidentified2' : 0.0},
{'column': 1851, 'centerDistance' : 97.00, 'SII6716a': 6714.53, 'SII6731a' : 6726.91, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 0.0, 'unidentified2' : 0.0},
{'column': 1861, 'centerDistance' : 99.53, 'SII6716a': 0.0, 'SII6731a' : 6726.91, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7136.07, 'unidentified2' : 0.0},
{'column': 1871, 'centerDistance' : 102.0, 'SII6716a': 0.0, 'SII6731a' : 6732.85, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7136.26, 'unidentified2' : 0.0},
{'column': 1881, 'centerDistance' : 104.5, 'SII6716a': 6714.19, 'SII6731a' : 6729.03, 'SII6716b' : 0.0, 'SII6731b' : 0.0, 'unidentified1' : 0.0, 'ArIII7136' : 7136.07, 'unidentified2' : 0.0},
]

#plt.plot([0.0, 0.0], [0.0, 9000], 'g')
dist = [a['centerDistance'] for a in tab]
print('dist = ',dist)
SII6716a = [a['SII6716a'] for a in tab]
SII6731a = [a['SII6731a'] for a in tab]
SII6716b = [a['SII6716b'] for a in tab]
SII6731b = [a['SII6731b'] for a in tab]
unidentified1 = [a['unidentified1'] for a in tab]
ArIII7136 = [a['ArIII7136'] for a in tab]

lambda0SII6716 = 6716.44
vradSII6716a = (np.array(SII6716a) - lambda0SII6716) * speedOfLight / lambda0SII6716
vradSII6716b = (np.array(SII6716b) - lambda0SII6716) * speedOfLight / lambda0SII6716

lambda0SII6731 = 6730.815
vradSII6731a = (np.array(SII6731a) - lambda0SII6731) * speedOfLight / lambda0SII6731
vradSII6731b = (np.array(SII6731b) - lambda0SII6731) * speedOfLight / lambda0SII6731

lambda0ArIII7136 = 7135.8
vradArIII7136 = (np.array(ArIII7136) - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136

print('minimum velocity = ',min([min(vradSII6716a), min(vradSII6716b), min(vradSII6731a), min(vradSII6731b), min(vradArIII7136)]))
print('maximum velocity = ',max([max(vradSII6716a), max(vradSII6716b), max(vradSII6731a), max(vradSII6731b), max(vradArIII7136)]))

plt.plot(dist, vradSII6716a, 'b*', label='[SII] 6716')
plt.plot(dist, vradSII6716b, 'b*')
plt.plot(dist, vradSII6731a, 'r*', label='[SII] 6731')
plt.plot(dist, vradSII6731b, 'r*')
plt.plot(dist, vradArIII7136, 'g*', label='[ArIII] 7136')

plt.ylim(-1300., 1300.)
plt.xlim(-130., 130.)
#plt.xlim(-50., -15.)
plt.xlabel('center distance [arc sec]')
plt.ylabel('radial velocity [km/s]')
plt.legend()
plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_GTC.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


def getWavelength(header, axis=2):
    nPix = int(header['NAXIS'+str(axis)])
    print('getWavelength: nPix = ',nPix)
    crPix = int(header['CRPIX'+str(axis)])
    print('getWavelength: crPix = ',crPix)
    crVal = float(header['CRVAL'+str(axis)])
    print('getWavelength: crVal = ',crVal)
    cDelt = float(header['CDELT'+str(axis)])
    print('getWavelength: cDelt = ',cDelt)
  #    lam = np.ndarray(nPix, dtype=np.float32)
  #    lam[0] =
  #    for i in np.arange(1,nPix):
    lam = np.arange(crVal, crVal + ((nPix-0.5)*cDelt), cDelt, dtype=np.float32)
    print('getWavelength: lam = ',len(lam),': ',lam)
    return lam

def lambdaToY(lam, lambdaArr):
  return (lam-lambdaArr[0]) * len(lambdaArr) / (lambdaArr[len(lambdaArr)-1] - lambdaArr[0]) - 0.5

def getRaDecXY(string):
  strs = re.sub( '\s+', ' ', string ).strip()
  #  print('getRaDecXY: strs = ',strs)
  strs = strs.rstrip().split(' ')
  #print('getRaDecXY: strs = ',strs)
  return [strs[0], strs[1], float(strs[len(strs)-2]), float(strs[len(strs)-1])]

def getArcsecDistance(fitsName, x1, y1, x2, y2):
  result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
  result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
  #print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
  #print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
  raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
  raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
  #print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
  #print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
  mm1 = SkyCoord(ra=raHMS1, dec=decHMS1, unit=(u.hourangle, u.deg))
  mm2 = SkyCoord(ra=raHMS2, dec=decHMS2, unit=(u.hourangle, u.deg))
  return mm1.separation(mm2).arcsecond

# find x in Rows from [xFrom, yFrom] for column yTo, do it the crude way...
def findXForArcSecDistanceFrom(xFrom, yFrom, yTo, dist, xRange, fitsName):
  #print('findXForArcSecDistanceFrom: xFrom = ',xFrom,', yFrom =',yFrom,', yTo = ',yTo,', dist = ',dist,', xRange = ',xRange,', fitsName = ',fitsName)
  xRangeTemp = xRange
  distTemp = dist
  if dist < 0:
    distTemp = 0. - dist
    xRangeTemp = [xRange[0], xFrom]
  else:
    if xFrom > xRange[0]:
      xRangeTemp = [xFrom, xRange[1]]

#print('findXForArcSecDistanceFrom: dist = ',dist,', distTemp = ',distTemp)
  hdulist = pyfits.open(fitsName)
  nCols = hdulist[0].data.shape[1]
#print('findXForArcSecDistanceFrom: nCols = ',nCols)
  previousDist = 100000000.
  for x in np.arange(xRangeTemp[0], xRangeTemp[1]):
    #   print('findXForArcSecDistanceFrom: xFrom = ',xFrom,', yFrom =',yFrom,', yTo = ',yTo,', distTemp = ',distTemp,', xRange = ',xRange,', fitsName = ',fitsName)
    tempDist = getArcsecDistance(fitsName, xFrom, yFrom, x, yTo)
    #print('findXForArcSecDistanceFrom): distTemp = ',distTemp,': x = ',x,': tempDist = ',tempDist)
    if (((dist > 0.) and ((previousDist < distTemp) and (tempDist >= distTemp)))
        or ((dist < 0.) and ((previousDist > distTemp) and (tempDist <= distTemp)))
        or ((distTemp == 0) and (tempDist == 0.))):
      #print('findXForArcSecDistanceFrom: returning x = ',x)
      return x

    previousDist = tempDist
  raise('findXForArcSecDistanceFrom: ERROR: could not determine x')

def findYForWLen(wLenArr, wLen):
  previousWLen = 10000000.
  pos = 0
  for w in wLenArr:
    #    print('findYForWLen: previousWLen = ',previousWLen,', w = ',w,', wLen = ',wLen)
    if (previousWLen < wLen) and (w >= wLen):
      #  print('findYForWLen: returning pos = ',pos)
      return pos
    pos += 1
    previousWLen = w
  raise('findYForWLen: ERROR: wLen ',wLen,' not found')

minRow = 1573
maxRow = 1613
minCol = 1000
maxCol = 1920
starCol = 1463
xCenter = starCol - minCol
yCenter = 994

twoDSpecFile = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_mSky_obs_not_smoothed_minComb.fits'
imFile = '/Users/azuri/daten/uni/HKU/Pa30/gtc_object_wcsIm_sky_0000956342.fits'
hdulist = pyfits.open(twoDSpecFile)
header = hdulist[0].header
twoDSpec = hdulist[0].data[minRow:maxRow,minCol:maxCol]
print('twoDSpec.shape = ',twoDSpec.shape)

wavelength = getWavelength(header)[minRow:maxRow]
print('wavelength = ',wavelength.shape,': ',wavelength)

#remove background
upper = hdulist[0].data[maxRow:maxRow+20,minCol:maxCol]
lower = hdulist[0].data[minRow-20:minRow,minCol:maxCol]
for col in range(twoDSpec.shape[1]):
    twoDSpec[:,col] = twoDSpec[:,col] - np.mean([np.mean(upper[:,col]), np.mean(lower[:,col])])

rows = range(twoDSpec.shape[0])
cols = range(twoDSpec.shape[1])

xlim = [-0.5,len(cols)-0.5]
ylim = [-0.5,len(rows)-0.5]

# show 2D spectrum and reverse y (because of how matrixes are stored / interpreted)
twoDSpecPlot = np.ndarray(twoDSpec.shape)
for row in range(twoDSpec.shape[0]):
  twoDSpecPlot[row,:] = twoDSpec[twoDSpec.shape[0]-1-row,:]
plt.imshow(twoDSpecPlot, cmap='Greys',vmin=0., vmax=2.0e-19, extent=(xlim[0],xlim[1],ylim[0],ylim[1]), aspect='auto')

# mark CS with X
plt.plot([len(cols)/2.],[len(rows)/2.],'rx',markersize=10)
limits = plt.axis()

# plot measured lines with vrad as color
cm = plt.cm.get_cmap('rainbow')

i = 0
sc = None
for lam in [SII6716a,SII6716b,SII6731a,SII6731b]:
    plotVRadY = lambdaToY(lam, wavelength)
    #    plotVRadY = len(wavelength) - lambdaToY(lam, wavelength)
    print('plotVRadY = ',plotVRadY)
    plotVRadX = np.array([a['column'] for a in tab]) - minCol
    print('plotVRadX = ',plotVRadX)
      #    if i < 2:
      #  marker = 'r+'
      #else:
      #  marker = 'b+'
          #    plt.plot(plotVRadX,
          #   plotVRadY,
          #   marker)
    marker = None
    label = None
    if i == 0:
        vrad = vradSII6716a
        marker = 'o'
        label = '[SII] 6716'
    elif i == 1:
        vrad = vradSII6716b
        marker = 'o'
    #        label = '[SII] 6716'
    elif i == 2:
        vrad = vradSII6731a
        marker = '^'
        label = '[SII] 6731'
    elif i == 3:
        vrad = vradSII6731b
        marker = '^'
#    label = '[SII] 6731'

    sc = plt.scatter(plotVRadX,
                     plotVRadY,
                     c=vrad,
                     vmin=-1200.,
                     vmax=1200.,
                     s=15,
                     cmap=cm,
                     marker = marker,
                     label = label
                    )
    i += 1
plt.xlim = xlim
plt.ylim = ylim
plt.axis(limits)
plt.autoscale(False)
plt.colorbar(sc,label='radial velocity [km/s]')
plt.legend()

# change x-Axis
xTicksDist = [-100., -80., -60., -40., -20., 0., 20., 40., 60.,80.,100.]
xTicksCols = []
xRange = [minCol, maxCol]
for xTick in xTicksDist:
  if xTick != xTicksDist[0]:
    xRange = [xTicksCols[len(xTicksCols)-1], xRange[1]]
  xTicksCols.append(findXForArcSecDistanceFrom(xCenter+minCol, yCenter, yCenter, xTick, xRange, imFile))
print('xTicksCols = ',xTicksCols)

xTicks = []
for x in xTicksCols:
  xTicks.append(x - minCol)
print('xTicks = ',xTicks)

xTicksDistStr =()
for x in xTicksDist:
  xTicksDistStr = xTicksDistStr + ('%.0d' % (x),)
print('findXForArcSecDistanceFrom: xTicksDistStr = ',xTicksDistStr)

locks, labels = plt.xticks()
print('findXForArcSecDistanceFrom: locks = ',locks,', labels[0] = ',labels[0])

plt.xticks(xTicks,
           xTicksDistStr)

plt.xlabel('center distance [arcsec]')

yTicksWLen = ('%.0d' % ((int(getWavelength(header)[minRow] / 10.) + 1.) * 10.),)
yTicksRow = [findYForWLen(wavelength, int(yTicksWLen[len(yTicksWLen)-1]))]
while True:
  if int(yTicksWLen[len(yTicksWLen)-1]) + 10. > getWavelength(header)[maxRow]:
    break
  yTicksWLen = yTicksWLen + ('%.0d' % (int(yTicksWLen[len(yTicksWLen)-1]) + 10.),)
  yTicksRow.append(findYForWLen(wavelength, int(yTicksWLen[len(yTicksWLen)-1])))
print('findXForArcSecDistanceFrom: yTicksRow = ',yTicksRow)
print('findXForArcSecDistanceFrom: yTicksWLen = ',yTicksWLen)
plt.yticks(yTicksRow, yTicksWLen)
plt.ylabel('wavelength [$\mathrm{\AA}$]')

# insert magnification of most striking features
a = plt.axes([.135, .67, .15, .18])
magMinRow = 4
magMaxRow = twoDSpecPlot.shape[0] - 4
magMinCol = 300
magMaxCol = 400
plt.imshow(twoDSpecPlot[magMinRow:magMaxRow,magMinCol:magMaxCol], cmap='Greys',vmin=0., vmax=2.0e-19, extent=(xlim[0],xlim[1],ylim[0],ylim[1]), aspect='auto')
#plt.plot([1.,2.], [1.,2.], 'b+')
#plt.title('Probability')
plt.xticks([])
plt.yticks([])


plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_map_on_2dspec.eps', format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
