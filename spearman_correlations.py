import pandas as pd 
from scipy.stats import spearmanr
from scipy import stats

# Import the spreadsheet
#df = pd.read_csv("1H0707_relxill_RF_ECM_TM.csv") 
df = pd.read_csv("1H0707_relxill_RF_ECM_SIM.csv") 
#df = pd.read_csv("IRAS_relxill_RF_ECM_SIM.csv") 
#df = pd.read_csv("1H_and_IRAS_relxill_RF_ECM_SIM.csv")
 
# Define columns
L_edd = df.iloc[:, 1] 
L_bol = df.iloc[:, 2] 
Edd_frac = df.iloc[:, 4]
Lag = df.iloc[:, 7] 
Lag_err = df.iloc[:, 8] 
Lag_freq = df.iloc[:, 9]
Lag_freq_err = df.iloc[:, 10]
cts = df.iloc[:, 11]
L03_10 = df.iloc[:, 13] # Luminosity 0.3-10 keV
L2_10 = df.iloc[:, 14] # Luminosity 2-10 keV
Ref_F = df.iloc[:, 17] # Relflection flux
Ref_F_hi = df.iloc[:, 18]
Ref_F_lo = df.iloc[:, 19]
Plaw_F = df.iloc[:, 20]
Plaw_F_hi = df.iloc[:, 21] # Powerlaw flux
Plaw_F_lo = df.iloc[:, 22] # Powerlaw flux

Gamma_Plaw = df.iloc[:, 23] # Powerlaw 
Gamma_Plaw_hi = df.iloc[:, 24]
Gamma_Plaw_lo = df.iloc[:, 25]

Gamma_Relx = df.iloc[:, 26]
Gamma_Relx_hi = df.iloc[:, 27]
Gamma_Relx_lo = abs(df.iloc[:, 28])

xi = df.iloc[:, 29] # ionisation
xi_ho = df.iloc[:, 30]
xi_lo = abs(df.iloc[:, 31])

Fe = df.iloc[:, 32]
Fe_hi = df.iloc[:, 33]
Fe_lo = abs(df.iloc[:, 34])

RF = df.iloc[:, 35]
RF_hi = df.iloc[:, 36]
RF_lo = abs(df.iloc[:, 37])

inc = df.iloc[:, 38]
inc_hi = df.iloc[:, 39]
inc_lo = abs(df.iloc[:, 40])

Cvr = df.iloc[:, 41] # Covering fraction
Cvr_hi = df.iloc[:, 42] # Covering fraction
Cvr_lo = abs(df.iloc[:, 43]) # Covering fraction

nH = df.iloc[:, 44]
nH_hi = df.iloc[:, 45]
nH_lo = abs(df.iloc[:, 46])


#ecm Results
h2_ecm = df.iloc[:, 49] 
h2_ecm_hi = df.iloc[:, 50]
h2_ecm_lo = abs(df.iloc[:, 51])

g1_ecm = df.iloc[:, 52]
g1_ecm_hi = df.iloc[:, 53]
g1_ecm_lo = abs(df.iloc[:, 54])

g2_ecm = df.iloc[:, 55]
g2_ecm_hi = df.iloc[:, 56]
g2_ecm_lo = abs(df.iloc[:, 57])

xi_ecm = df.iloc[:, 58]
xi_ecm_hi = df.iloc[:, 59]
xi_ecm_lo = abs(df.iloc[:, 60])

b_ecm = df.iloc[:, 61]
b_ecm_hi = df.iloc[:, 62]
b_ecm_lo= abs(df.iloc[:, 63])

tmax_ecm = df.iloc[:, 64]
tmax_ecm_hi = df.iloc[:, 65]
tmax_ecm_lo = abs(df.iloc[:, 66])

tshift_ecm = df.iloc[:, 67]
tshift_ecm_hi = df.iloc[:, 68]
tshift_ecm_lo = abs(df.iloc[:, 69])

g_delta = df.iloc[:, 70]
g_ratio = df.iloc[:, 71]
h2_delta = df.iloc[:, 72]


#M = df.iloc[:, 70]
#M_hi = df.iloc[:, 71]
#M_lo = abs(df.iloc[:, 72])

# To find the correlation among the columns using Spearman method 
print('g2_ecm - L_bol', stats.spearmanr(g2_ecm, L_bol))
print('g2_ecm - Edd_frac', stats.spearmanr(g2_ecm, Edd_frac))
print('g2_ecm - Ref_F', stats.spearmanr(g2_ecm, Ref_F))
print('g2_ecm - Plaw_F', stats.spearmanr(g2_ecm, Plaw_F))
print('g2_ecm - Gamma_Plaw', stats.spearmanr(g2_ecm, Gamma_Plaw))
print('g2_ecm - Gamma_Relx', stats.spearmanr(g2_ecm, Gamma_Relx))
print('g2_ecm - b_ecm', stats.spearmanr(g2_ecm, b_ecm))
print('g2_ecm - tmax_ecm', stats.spearmanr(g2_ecm, tmax_ecm))
print('g2_ecm - tshift_ecm', stats.spearmanr(g2_ecm, tshift_ecm))
print('g2_ecm - xi', stats.spearmanr(g2_ecm, xi))
print('g2_ecm - Fe', stats.spearmanr(g2_ecm, Fe))
print('g2_ecm - RF', stats.spearmanr(g2_ecm,RF ))
print('g2_ecm - inc', stats.spearmanr(g2_ecm, inc))
