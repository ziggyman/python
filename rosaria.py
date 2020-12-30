import numpy as np
import csv
import matplotlib.pyplot as plt

fName = '/Users/azuri/Downloads/rosaria.csv'

h = []
d = []
nRows = 0
with open(fName,'r') as f:
    data = csv.DictReader(f)
    for row in data:
        h.append(float(row['height']))#km
        d.append(float(row['density']) * np.power(10.,15))#g/km^3
        nRows += 1
print('nRows = ',nRows)
print('height = ',h)
print('density = ',d)

h = np.array(h)#km
d = np.array(d)

mu = 398600.4418#km^3 / s^2
re = 6371.#km
cd = 2.2
m = 8000.#g
area = 3. * np.power(10.,-8)#km^2

a = h + re#km
b = re + 80.#km

v = np.sqrt(mu/a)
dvf = v * (np.sqrt(2. * b / (a + b)) - 1.)
dvdrag = (d * v * v * cd * area / (2. * m)) * 6307200.

fig = plt.figure(figsize=(10, 6))
plt.plot(h,dvf)
plt.xlabel('height [km]')
plt.ylabel('$\Delta V_f$ [km/s]')
fig.tight_layout()
fig.savefig('/Users/azuri/Downloads/rosaria_h_vs_dvf.pdf', bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10, 6))
plt.plot(h,dvdrag * 1000.)
plt.xlabel('height [km]')
plt.ylabel('$\Delta V_{drag}$ [m/s]')
fig.tight_layout()
fig.savefig('/Users/azuri/Downloads/rosaria_h_vs_dvdrag.pdf', bbox_inches='tight')
plt.show()
