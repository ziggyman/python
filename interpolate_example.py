from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
print('xx = ',xx.shape,': ',xx)
print('yy = ',yy.shape,': ',yy)

z = np.sin(xx**2+yy**2)
print('z = ',z.shape,': ',z)
f = interpolate.interp2d(x, y, z, kind='cubic')
xnew = np.arange(-5.01, 5.01, 1e-2)
ynew = np.arange(-5.01, 5.01, 1e-2)
znew = f(xnew, ynew)
plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
plt.show()