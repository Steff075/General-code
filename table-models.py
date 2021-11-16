# Run this in Python 3

# Updated to remove some of the gaps in the tables
# It is OK for zeros when gamma1 < gamma2 (by design gamma1 should be larger than gamma2)
# Also updated code to check for NaN in the lag files

# We can use this to check to see if a file exists
import os.path
import sys

# HEASP has to be on the path - this is done by initialising HEASOFT first
from heasp import *

# Need to read FITS files to get "spectra"
from astropy.io import fits

# We need numpy for some arrays
import numpy as np

# path to files containing the lag data
path_to_rev = "/data/typhon1/reverb"

# need to know grid size because we don't read in the grid until we find the first valid file
# the code checks that this value is correct in any case
n_e = 491

# model parameters
rev_par_names = ["h1", "h2", "i", "gamma1", "gamma2", "A", "xi", "p", "r1", "r2", "b", "q2", "tmax", "tshift"]
# tabulated values
rev_par_values = [
	[2.0], 											# h1
	[3.0, 4.0, 5.0, 6.0, 7.0, 8.0 ,9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0],	# h2
	[53.0],											# i
	[1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3],	# gamma1 (> gamma2)
	[1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2], # gamma2
	[1.0],											# A
	[0.0, 1.0, 2.0, 3.0],									# xi
	[0.0],											# p
	[1.0],											# r1
	[1.0],											# r2
	[1.0, 2.0, 3.0],									# b
	[0.5],											# q2
	[50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0], # tmax
	[10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0]	# tshift
	]
# initial values (gamma2 < gamma1), delta
rev_par_initial = [2.0, 6.0, 53.0,  2.5,  2.6,  1.0,  1.0,  0.0,  1.0,  1.0,  2.0,  0.5, 100.0,  70.0]
rev_par_delta =   [0.5, 0.5, -1.0,  0.1,  0.1, -1.0, 0.25, -1.0, -1.0, -1.0,  0.5, -1.0, -1.0, -1.0]

# Define reverberation table model
rev_tab = table()

# set table descriptors and the energy array
rev_tab.setModelName("reverb")
rev_tab.setModelUnits(" ")
rev_tab.setisRedshift(False)
rev_tab.setisAdditive(True)
rev_tab.setisError(False)

# read in the energy grid
# it is obviously essential that all models are computed on the same grid
# we'll read in a representative file to get the grid
e_grid_file = 'lag-files/i53/a0.998/h02.0+06.0/g2.9+2.7/A01.0/x2.0/p00.0/b02.0/q0.5/tmax0300/tshift0080/lag.fits'
hdulist = fits.open(path_to_rev + "/" + e_grid_file)
tabdata = hdulist[1].data
# put the time lags into the table model
if (n_e != len(tabdata)):
	sys.exit("Check data table length")
cur_e_grid = []
for j in range(n_e):
	cur_e_grid.append(tabdata[j][0])
rev_tab.setEnergies(np.array(cur_e_grid))

# energy bin widths; repeat last value even though it is not really needed
e_bin_width = []
for j in range(n_e-1):
	e_bin_width.append(cur_e_grid[j+1]-cur_e_grid[j])
e_bin_width.append(e_bin_width[-1])

# all interpolated parameters
rev_tab.setNumIntParams(len(rev_par_names))
rev_tab.setNumAddParams(0)
print("Number of parameters = %i" % len(rev_par_names))

# add the list of parameters to table model
for i in range(len(rev_par_names)):
	print( "par[%i] = %s with %i values" % (i, rev_par_names[i], len(rev_par_values[i])) )
	cur_par = tableParameter()
	cur_par.setName(rev_par_names[i])
	cur_par.setInterpolationMethod(0)	# 0 - linear; 1 - logarithmic
	cur_par.setInitialValue(rev_par_initial[i])
	cur_par.setDelta(rev_par_delta[i])
	cur_par.setMinimum(rev_par_values[i][0])
	cur_par.setBottom(rev_par_values[i][0])
	cur_par.setTop(rev_par_values[i][-1])
	cur_par.setMaximum(rev_par_values[i][-1])
	cur_par_tab_vals = []
	for j in range(len(rev_par_values[i])):
		cur_par_tab_vals.append(rev_par_values[i][j])
	cur_par.setTabulatedValues(np.array(cur_par_tab_vals))

	# and push it onto the vector of parameters
	rev_tab.pushParameter(cur_par)

# add spectra to table model
# the first time we do this we will also add the energy grid
index = []
for i in range(len(rev_par_names)):
	index.append(0)

carry = 0

while (carry==0):
	#print "%i %i %i %i %i %i %i %i" % (index[0], index[1], index[2], index[3], index[4], index[5], index[6], index[7])

	cur_file = "lag-files/i{:02d}/a0.998/h{:04.1f}+{:04.1f}/g{:03.1f}+{:03.1f}/A{:04.1f}/x{:03.1f}/p{:04.1f}/b{:04.1f}/q{:03.1f}/tmax{:04.0f}/tshift{:04.0f}/lag.fits" . format( int(rev_par_values[2][index[2]]), rev_par_values[0][index[0]], rev_par_values[1][index[1]], rev_par_values[3][index[3]], rev_par_values[4][index[4]], rev_par_values[5][index[5]], rev_par_values[6][index[6]], rev_par_values[7][index[7]], rev_par_values[10][index[10]], rev_par_values[11][index[11]], rev_par_values[12][index[12]], rev_par_values[13][index[13]] )

	print("  file = %s" % cur_file)

	# clear current "spectrum" and add parameter values
	cur_spec = tableSpectrum()

	# make list of parameter rev_par_values
	cur_spec_par_vals = []
	for i in range(len(rev_par_names)):
		cur_spec_par_vals.append(rev_par_values[i][index[i]])
	cur_spec.setParameterValues(np.array(cur_spec_par_vals))

	# if file is non-existent see if we can find a replacement with slightly smaller t_max or t_shift
	if not os.path.isfile(path_to_rev + "/" + cur_file):
		# try to find a replacement if gamma1 > gamma2
		# gamma1 is parameter [3]
		# gamma2 is parameter [4]
		if rev_par_values[3][index[3]] > rev_par_values[4][index[4]]:

			# t_max is parameter [12]
			# t_shift is parameter [13]

			# first try the previous t_shift
			try_cur_file = "lag-files/i{:02d}/a0.998/h{:04.1f}+{:04.1f}/g{:03.1f}+{:03.1f}/A{:04.1f}/x{:03.1f}/p{:04.1f}/b{:04.1f}/q{:03.1f}/tmax{:04.0f}/tshift{:04.0f}/lag.fits" . format( int(rev_par_values[2][index[2]]), rev_par_values[0][index[0]], rev_par_values[1][index[1]], rev_par_values[3][index[3]], rev_par_values[4][index[4]], rev_par_values[5][index[5]], rev_par_values[6][index[6]], rev_par_values[7][index[7]], rev_par_values[10][index[10]], rev_par_values[11][index[11]], rev_par_values[12][index[12]], rev_par_values[13][index[13]-1] )
			if os.path.isfile(path_to_rev + "/" + try_cur_file):
				# success - we have found a replacement file
				print("  found replacement t_shift")
				cur_file = try_cur_file
			else:
				# let's try the prevoius t_max instead
				try_cur_file = "lag-files/i{:02d}/a0.998/h{:04.1f}+{:04.1f}/g{:03.1f}+{:03.1f}/A{:04.1f}/x{:03.1f}/p{:04.1f}/b{:04.1f}/q{:03.1f}/tmax{:04.0f}/tshift{:04.0f}/lag.fits" . format( int(rev_par_values[2][index[2]]), rev_par_values[0][index[0]], rev_par_values[1][index[1]], rev_par_values[3][index[3]], rev_par_values[4][index[4]], rev_par_values[5][index[5]], rev_par_values[6][index[6]], rev_par_values[7][index[7]], rev_par_values[10][index[10]], rev_par_values[11][index[11]], rev_par_values[12][index[12]-1], rev_par_values[13][index[13]] )
				if os.path.isfile(path_to_rev + "/" + try_cur_file):
					# success - we have found a replacement file
					print("  found replacement t_max")
					cur_file = try_cur_file
				else:
					# try changing both
					try_cur_file = "lag-files/i{:02d}/a0.998/h{:04.1f}+{:04.1f}/g{:03.1f}+{:03.1f}/A{:04.1f}/x{:03.1f}/p{:04.1f}/b{:04.1f}/q{:03.1f}/tmax{:04.0f}/tshift{:04.0f}/lag.fits" . format( int(rev_par_values[2][index[2]]), rev_par_values[0][index[0]], rev_par_values[1][index[1]], rev_par_values[3][index[3]], rev_par_values[4][index[4]], rev_par_values[5][index[5]], rev_par_values[6][index[6]], rev_par_values[7][index[7]], rev_par_values[10][index[10]], rev_par_values[11][index[11]], rev_par_values[12][index[12]-1], rev_par_values[13][index[13]-1] )
					if os.path.isfile(path_to_rev + "/" + try_cur_file):
						# success - we have found a replacement file
						print("  found replacement t_shift and t_max")
						cur_file = try_cur_file
			# we've failed to find a replacement
			#input("  no replacement")
	# if we don't find a replacement we'll just fill in with zeros

	# see if file exists (non-existent file might mean we have an invalid choice of parameters)
	if os.path.isfile(path_to_rev + "/" + cur_file):
		#print("  exists")
		hdulist = fits.open(path_to_rev + "/" + cur_file)
		tabdata = hdulist[1].data
		# put the time lags into the table model
		cur_spec_flux = []
		# note that because we have bin hi and lo there are one fewer bins than there are energies
		for j in range(n_e-1):
			# table models contain quantities integrated over the bin width
			cur_spec_flux.append(tabdata[j][1] * e_bin_width[j])
		cur_spec.setFlux(np.array(cur_spec_flux))
	else:
		#print("  does not exist")
		if rev_par_values[3][index[3]] > rev_par_values[4][index[4]]:
			# we don't expect lag files to exist for gamma1 <= gamma2
			print("  does not exist moving to next file")
		# put zeros into the model if file does not exist (parameter choice is invalid)
		zero_flux = []
		for i in range(n_e-1):
			zero_flux.append(0.0)
		cur_spec.setFlux(np.array(zero_flux))
	rev_tab.pushSpectrum(cur_spec)

	# go to next file in the list
	carry = 1
	for i in range(len(rev_par_names)-1,-1,-1):
		if (carry == 1):
			index[i] = index[i] + 1
			carry = 0
			if (index[i] >= len(rev_par_values[i])):
				carry = 1
				index[i] = 0

# save table model
rev_tab.write("tab_rev.mod")
