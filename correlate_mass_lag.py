# Pearson correlation coefficient is a statistic that measures linear correlation between two variables X and Y. It has a value between +1 and −1. A value of +1 is total positive linear correlation, 0 is no linear correlation, and −1 is total negative linear correlation.

import pandas as pd 
from scipy.stats import spearmanr
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use('classic')
plt.rcParams['font.family'] = 'serif'

# Import the spreadsheet
df = pd.read_csv("results-Final-spectra.csv") 
#df = pd.read_csv("results-Final-spectra-COMBINED.csv") 
#df = pd.read_csv("results-Final-spectra-HIFLUX.csv") 
#df = pd.read_csv("results-Final-spectra-LOFLUX.csv") 

# Define columns
# Define columns
M = df.iloc[:, 1] # Log M
M_hi_err = df.iloc[:, 2]
M_lo_err = df.iloc[:, 3]
z = df.iloc[:, 4] 
Edd_Frac = df.iloc[:, 5] 
L_edd = df.iloc[:, 6] 
L_bol = df.iloc[:, 7] 
Edd_ratio = df.iloc[:, 9]
Lag = df.iloc[:, 11] 
Lag_err = df.iloc[:, 12] 
Lag_freq = df.iloc[:, 13]
cts = df.iloc[:, 15]
L03_10 = df.iloc[:, 17] # Luminosity 0.3-10 keV
L2_10 = df.iloc[:, 19] # Luminosity 2-10 keV
Ref_F = df.iloc[:, 24] # Relflection flux
Plaw_F = df.iloc[:, 27] # Powerlaw flux
Gamma_Plaw = df.iloc[:, 30]
Gamma_Plaw_hi_err = df.iloc[:, 31]
Gamma_Plaw_lo_err = df.iloc[:, 32]
Gamma_Relx = df.iloc[:, 33]
xi = df.iloc[:, 36] # ionisation
Fe = df.iloc[:, 39]
RF = df.iloc[:, 42]
inc = df.iloc[:, 45]
Cvr = df.iloc[:, 48] # Covering fraction
nH = df.iloc[:, 51]


print('Cvr - Gamma_Plaw', stats.spearmanr(Cvr,Gamma_Plaw, nan_policy='omit'))
print('Cvr - Gamma_Relx', stats.spearmanr(Cvr,Gamma_Relx, nan_policy='omit'))
print('Cvr - nH', stats.spearmanr(Cvr,nH, nan_policy='omit'))
print('Cvr - xi', stats.spearmanr(Cvr,xi, nan_policy='omit'))
print('Cvr - Fe', stats.spearmanr(Cvr,Fe, nan_policy='omit'))
print('Cvr - Edd_Frac', stats.spearmanr(Cvr,Edd_Frac, nan_policy='omit'))
print('Cvr - Edd_ratio', stats.spearmanr(Cvr,Edd_ratio, nan_policy='omit'))


xdata = 10**M
xerr_hi = df.iloc[:, 2]
xerr_lo = df.iloc[:, 3] 
asymmetric_xerror = np.array(list(zip(xerr_lo, xerr_hi))).T

ydata = Lag
yerr = Lag_err


mean_M = np.mean(xdata)
sigma_M = np.std(xdata)
print('mean_M: %0.3f' % mean_M)
print('sigma_M: %0.3f' % sigma_M)

mean_Lag = np.mean(ydata)
sigma_Lag = np.std(ydata)
print('mean_Lag: %0.3f' % mean_Lag)
print('sigma_Lag: %0.3f' % sigma_Lag)


##############################################
# Plot the data and polyfit
##############################################

logx, logy, log_sig_lag = np.log(xdata), np.log(ydata), np.log(sigma_Lag)
p = np.polyfit(logx, logy, 1)
y_fit = np.exp(np.polyval(p, logx))

print('log_sig_lag: %0.3f' % log_sig_lag)

# define values for upper and lower sigma regions
p_sig_hi = np.polyfit(logx, logy + log_sig_lag/2, 1)
p_sig_lo = np.polyfit(logx, logy - log_sig_lag/2, 1)
p_sig_fit_hi  = np.exp(np.polyval(p_sig_hi, logx))
p_sig_fit_lo  = np.exp(np.polyval(p_sig_lo, logx))

# Plot the data and polyfit with 1 sigma regions
plt.errorbar(xdata, ydata, yerr=yerr, xerr=asymmetric_xerror, markersize=4.5,fmt='o', c='steelblue', ecolor = 'lightgrey')
plt.plot(xdata, y_fit,'r-', lw=0.5)
plt.plot(xdata, p_sig_fit_hi, '-', c = 'lightgrey', lw=0.5)
plt.plot(xdata, p_sig_fit_lo, '-', c = 'lightgrey', lw=0.5)
#plt.fill_between(xdata, y_fit-p_sig_fit_lo, y_fit+p_sig_fit_hi, alpha=.25)
plt.xlabel(r'$M/M_{\odot}$', fontsize=16)
plt.ylabel(r'Time lag (s)', fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e5,2e9)
plt.ylim(5e-2,2e4)
#plt.show(block=False)
#plt.show()
plt.savefig('Mass_lag_log_polyfit_ALL.png')
plt.savefig('Mass_lag_log_polyfit_ALL.pdf')
plt.close()



