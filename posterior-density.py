import numpy as np
import matplotlib.pyplot as plt

#Import data
data = np.loadtxt('1H0707-495_results_single_obs.csv', delimiter=',', skiprows=1)

mass = data[:,101]

ndim = 15 # define number of dimensions of the data

cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov, cov)

nwalkers = 100

# Guess sarting point for each of the 100 walkers.
p0 = np.random.rand(nwalkers, ndim)

# Define function that returns the density p(x⃗ ) for specific values of x⃗ , μ⃗  and Σ−1. In fact, emcee actually requires the logarithm of p.

def log_prob(x, mu, cov):
    diff = x - mu
    return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))

# Import EMCEE and define the Ensembler
import emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[mass, cov])

# Define the probability function (log_prob) for this data
log_prob(p0[0], mass, cov)

print("log_prob: {0:.3f}".format(log_prob(p0[0], mass, cov)))

# Setup some initial 'burn- in' steps for the MCMC chain to allow walkers to explore parameter space to settle into the maximum of the density
print("Running burn-in...")
state = sampler.run_mcmc(p0, 100) # Save the final position of the walkers after N steps as 'state'

sampler.reset() # Reset the sampler before running MCMC

# Run MCMC for N steps starting at the final 'burn-in' locations
print("Running MCMC chain...")
sampler.run_mcmc(state, 50000);

mass_mean = np.mean(mass)
mass_sigma = np.std(mass)

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
print("Mean autocorrelation time: {0:.3f}steps".format(np.mean(sampler.get_autocorr_time())))
print("Mass mean: {0:.3f}".format(np.mean(mass)))
print("Mass sigma: {0:.3f}".format(np.std(mass)))

conf68 = stats.norm.interval(0.68, loc=mass_mean, scale=mass_sigma/np.sqrt(len(mass)))
conf95 = stats.norm.interval(0.95, loc=mass_mean, scale=mass_sigma/np.sqrt(len(mass)))

#conf68_lines = (conf68[:,0], mass_median, conf68[:,1])

samples = sampler.get_chain(flat=True)
plt.hist(samples[:, 0], 100, color="k", histtype="step")
plt.axvline(mass_mean, lw =0.5, c="k") # Mean
plt.axvline(mass_mean - mass_sigma, lw =0.8, linestyle='dashed', c="grey") # sigma_lo
plt.axvline(mass_mean + mass_sigma, lw =0.8, linestyle='dashed', c="grey") # sigma_hi
plt.axvspan(mass_mean - mass_sigma, mass_mean + mass_sigma, alpha=0.2, color="b") # Fill sigma range
plt.xlim(1,10)
plt.xlabel(r"$M_{BH}~logM_{\odot}$", fontsize=16)
plt.ylabel(r"Posterior Density", fontsize=16)
plt.gca().set_yticks([]);
#plt.show(block=False)

plt.savefig('1H0707-mass-density-histogram-fmax.png')
plt.savefig('1H0707-mass-density-histogram-fmax.pdf')
