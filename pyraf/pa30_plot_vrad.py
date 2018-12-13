execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

speedOfLight = 299792.458# km/s

tab = [
{'column': 1323, 'centerDistance' : -37.31, 'SII6716a': 0.0, 'SII6731a' : 0.0, 'SII6716b' : 6736.7, 'SII6731b' : 0.0, 'unidentified1' : 6999.17, 'ArIII7136' : 7120.80, 'unidentified2' : 7160.06},
{'column': 1324, 'centerDistance' : -37.06, 'SII6716a': 6697.93, 'SII6731a' : 6712.08, 'SII6716b' : 6738.15, 'SII6731b' : 0.0, 'unidentified1' : 7001.46, 'ArIII7136' : 7115.85, 'unidentified2' : 7160.63},
{'column': 1325, 'centerDistance' : -36.79, 'SII6716a': 6698.82, 'SII6731a' : 6711.64, 'SII6716b' : 6736.49, 'SII6731b' : 0.0, 'unidentified1' : 6999.32, 'ArIII7136' : 7115.25, 'unidentified2' : 7160.82},
{'column': 1326, 'centerDistance' : -36.54, 'SII6716a': 6699.99, 'SII6731a' : 6713.83, 'SII6716b' : 6736.79, 'SII6731b' : 6753.37, 'unidentified1' : 7000.49, 'ArIII7136' : 0.0, 'unidentified2' : 7158.87},
{'column': 1327, 'centerDistance' : -36.28, 'SII6716a': 6702.50, 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
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
{'column': 1400, 'centerDistance' : -17.72, 'SII6716a': 6700.96, 'SII6731a' : 6716.04, 'SII6716b' : 6732.86, 'SII6731b' : 6745.74, 'unidentified1' : 6994.81, 'ArIII7136' : 0, 'unidentified2' : 7152.56},
{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
#{'column': , 'centerDistance' : , 'SII6716a': , 'SII6731a' : , 'SII6716b' : , 'SII6731b' : , 'unidentified1' : , 'ArIII7136' : , 'unidentified2' : },
]

plt.plot([0.0, 0.0], [0.0, 9000], 'g')
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

plt.plot(dist, vradSII6716a, 'b*')
plt.plot(dist, vradSII6716b, 'bx')
plt.plot(dist, vradSII6731a, 'r*')
plt.plot(dist, vradSII6731b, 'rx')
plt.plot(dist, vradArIII7136, 'g*')

plt.ylim(-1100., 1100.)
plt.xlabel('center distancd in arc sec')
plt.ylabel('radial velocity')
plt.show()