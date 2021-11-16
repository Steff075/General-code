#! /usr/bin/python3

# Look at all the best fits and pick the one with the smallest statistic

import glob
import os
import string

files = glob.glob('best_fit_??.txt')
best_stat = 1.0E10;
best_file = "";
for file in files:
	with open(file) as f:
		cur_stat = float(next(f))
		if (cur_stat < best_stat):
			best_stat = cur_stat
			best_file = file

print("Best fit is", best_stat, "in file", best_file)
os.system("cp " + best_file + " best_fit.txt")
best_file.replace(".txt", ".par")
os.system("cp " + best_file + " best_fit.par")
best_file.replace(".par", ".ps")
os.system("cp " + best_file + " best_fit.ps")
