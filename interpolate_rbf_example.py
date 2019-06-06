import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#%matplotlib inline
import random
# set seed to reproducible
random.seed(1)
data_size = 51
max_value_range = 132651
x = np.array([random.random()*max_value_range for p in range(0,data_size)])
y = np.array([random.random()*max_value_range for p in range(0,data_size)])
z = 2*x*x*x + np.sqrt(y)*y + random.random()
fig = plt.figure(figsize=(10,6))
ax = axes3d.Axes3D(fig)
ax.scatter3D(x,y,z, c='r')

x_grid = np.linspace(0, 132651, 1000*len(x))
y_grid = np.linspace(0, 132651, 1000*len(y))
B1, B2 = np.meshgrid(x_grid, y_grid, indexing='xy')
Z = np.zeros((x.size, z.size))

import scipy as sp
import scipy.interpolate
spline = sp.interpolate.Rbf(x,y,z,function='thin_plate',smooth=5, episilon=5)

Z = spline(B1,B2)
fig = plt.figure(figsize=(10,6))
ax = axes3d.Axes3D(fig)
ax.plot_wireframe(B1, B2, Z)
ax.plot_surface(B1, B2, Z,alpha=0.2)
ax.scatter3D(x,y,z, c='r')

predict_data_size = 132000
x_predict = np.array([random.random()*max_value_range for p in range(0,predict_data_size)])
y_predict = np.array([random.random()*max_value_range for p in range(0,predict_data_size)])
z_predict = spline(x_predict, y_predict)
fig = plt.figure(figsize=(10,6))
ax = axes3d.Axes3D(fig)
ax.plot_wireframe(B1, B2, Z)
ax.plot_surface(B1, B2, Z,alpha=0.2)
ax.scatter3D(x_predict,y_predict,z_predict, c='r')
