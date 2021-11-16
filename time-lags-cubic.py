import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams['font.family'] = 'serif'
from scipy.interpolate import interp1d

df1
x = df1.iloc[:,0]
y = df1.iloc[:,2]
f = interp1d(x, y,fill_value='extrapolate')
f2 = interp1d(x, y, kind='cubic', bounds_error=False, fill_value='extrapolate')

x_new = np.linspace(0.00006,0.006, num=1000, endpoint=True)
plt.plot(x, y, 'o', x_new, f(x_new), '--', x_new, f2(x_new), '-', linewidth='1.5')
plt.xlabel(r'Frequency $(Hz)$', fontsize=14)
plt.ylabel(r'Time Lag $(s)$', fontsize=14)
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.xscale('log')
plt.xlim(5e-5,1e-2)
line = ([1e-5, 1e-2], [0,0])
plt.plot(line, linestyle='--', linewidth='0.3', color = 'black')
#plt.show()
plt.savefig('1H0707-comb-data-cubicspline.png')
plt.close()
