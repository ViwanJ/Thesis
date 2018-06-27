#!/usr/bin/env python

## This is a script to average any analysis over different trajectory
import numpy as np
import pandas as pd
import os, sys
import MDAnalysis
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

def load_data(input_filename, Data_set) :
	tmp = {}
	global path
	for f in Data_set.keys() :
		tmp[f] = do_stat([np.genfromtxt(path[x]+input_filename) for x in Data_set[f]])
	return tmp

class do_stat :
	def __init__(self, datas) :
		self.av = np.average(np.vstack(datas), axis = 0)
		self.std = np.std(np.vstack(datas), axis = 0)
		
class diff :
	def __init__(self, sel, input1 , input2):
		self.av = sel[input1].av - sel[input2].av
		self.std = np.sqrt(sel[input1].std**2 + sel[input2].std**2)
		self.x = np.arange(len(self.av)) + 1

	def plot(self, color='Blue') :
		plt.plot(self.x, self.av, color = color)
		plt.fill_between(self.x, self.av - self.std, self.av + self.std, alpha = 0.1, color = color)
		plt.xlim([1,self.x[-1]+5])
		plt.ticklabel_format(style = 'sci', scilimits = (0,4), axis = 'y')
		plt.xlabel('resid')
	
	def plot_save(self, filename) :
		plt.savefig('./Analysis_OverTJ/'+filename + '.svg', dpi = 300)
		plt.close()

	def set_bfactor(self, gro_path, outname) :
		u = MDAnalysis.Universe(gro_path)
		Protein = u.select_atoms('protein')
		for i in Protein.residues.resids :
			Protein.select_atoms('resid ' + str(i)).bfactors = [abs(self.av[i-1])/max(self.av)] * len(Protein.select_atoms('resid ' + str(i)))
		Protein.write('./Analysis_OverTJ/' + outname + '.pdb')		

if __name__ == '__main__' :
	if not os.path.exists('./Analysis_OverTJ'):
		os.makedirs('Analysis_OverTJ')
