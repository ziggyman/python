execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

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

lambda0SII6716 = 6716.0
vradSII6716a = (np.array(SII6716a) - lambda0SII6716) * speedOfLight / lambda0SII6716
vradSII6716b = (np.array(SII6716b) - lambda0SII6716) * speedOfLight / lambda0SII6716

lambda0SII6731 = 6731.0
vradSII6731a = (np.array(SII6731a) - lambda0SII6731) * speedOfLight / lambda0SII6731
vradSII6731b = (np.array(SII6731b) - lambda0SII6731) * speedOfLight / lambda0SII6731

lambda0ArIII7136 = 7136.0
vradArIII7136 = (np.array(ArIII7136) - lambda0ArIII7136) * speedOfLight / lambda0ArIII7136

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
plt.show()
