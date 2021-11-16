import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams['font.family'] = 'serif'

# Import csv and name cloumn headers
df1cols=['Freq','Freq_err', 'Lag','Lag_err']
df1 = pd.read_csv('/home/steff075/phd_ecm/ecm-for-python/1H0707-495-combined-data-lag.csv', names=df1cols, header=0)
df1

# gca stands for 'get current axis'
ax = plt.gca()
df1.plot(kind='line',x='Freq',y='Lag', xerr='Freq_err',  yerr='Lag_err', fmt='o', c='steelblue', ecolor='grey', linestyle='-', capsize=3, ax=ax)
plt.xlim(5e-5,2e-2)
plt.xscale('log')
figsize=(8,6)
plt.xlabel(r'Frequency $(Hz)$', fontsize=14)
plt.ylabel(r'Time Lag $(s)$', fontsize=14)
line = ([-1e-5, 1.0], [0,0])
plt.plot(line, linestyle='--', linewidth='0.3', color = 'black')
plt.show()
#plt.savefig('1H0707-comb-data.png')
plt.close()

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
