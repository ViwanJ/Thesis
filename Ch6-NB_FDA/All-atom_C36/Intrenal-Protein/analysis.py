#!/usr/bin/env python

import sys,os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import MDAnalysis
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


class PearsonR_Pair :
	'''
	Example : tmp = PearsonR_Pair(u)
		  tmp.load_file_data('pair_SC.dat', (11,521,521))
		  tmp.pearsonr_coef()
	  	  d_resid, d_resname = Sim_to_Xtal(u, 73, 261)
	'''
	def __init__(self, u) :
		ref = u.select_atoms('name CA').positions
		mobile = [u.select_atoms('name CA').positions for t in u.trajectory]
		self.rmsd = [rmsd(ref, mobile[i] ,superposition=True) for i in range(len(mobile))]
	
	def load_file_data(self,data, shape) :  ## shape = (nframe, nresidue, nresidue)
		self.data = np.loadtxt(data).reshape(shape)

	def load_data(self,data) :
		self.data = data
		
	def pearsonr_coef(self) :
		foo = np.array([np.array([stats.pearsonr(self.rmsd, self.data[:,i,j]) for j in range(len(self.data[0]))]) for i in range(len(self.data[0]))])
		self.coef = foo[:,:,0]
		self.p_value = foo[:,:,1]
		
	def hist(self):
		plt.hist(np.concatenate(np.nan_to_num(self.coef)), bins = 40)
		
def sel_pairs(data, criteria = 0.80) :
	con = abs(data) >= criteria 
	return [pairs for pairs in zip(np.where(con)[0]+1, np.where(con)[1] +1) ] ## +1 convert matrix index (0) to resid (1)

def add_vmd_labels(u, pairs_lists, output_file) :
	if len(pairs_lists) != 0 :	
		foo = u.select_atoms('name CA')
		Res_to_CaIndex = {r:i for r, i in zip(foo.resids, foo.indices)}
		residues_lists = set(np.concatenate(pairs_lists))
		resids = ''
		for res in  residues_lists :
			resids = resids + str(int(res)) + ' '
		with file(output_file, 'w') as outfile:
			outfile.write("atomselect macro sel_resids {resid " + resids + " and protein}\n")  
			outfile.write("mol representation Licorice 0.300000 10.000000 10.000000\n")
			outfile.write("mol selection sel_resids \nmol addrep top \n color Labels Bonds red \n")
			for res_i, res_j in pairs_lists :
				outfile.write("label add Bonds 0/" + str(Res_to_CaIndex[res_i]) + " 0/" + str(Res_to_CaIndex[res_j]) + "\n")
			outfile.write("label show Bonds all\n")	
	else :
		print "\n Warning !! no objects in your selected list. \n"

def Sim_to_Xtal(u, start, nres) :
	Protein = u.select_atoms('name CA')
	return {i : i + (start-1) - (i/(nres+1) * 261) for i in Protein.resids}, {i : Protein.resnames[i-1]  for i in Protein.resids}

def Array2DataFrame(array2D, d_resid , d_resname ) :
	df = pd.DataFrame(array2D).unstack().dropna()
	df.index.names = ['res_i', 'res_j']
	df = df.reset_index()
	df['res_i'] = df['res_i'] + 1
	df['res_j'] = df['res_j'] + 1
	df['name_i'] = [d_resname[i] for i in df['res_i']]
	df['name_j'] = [d_resname[i] for i in df['res_j']]
	df['Xid_i'] = [d_resid[i] for i in df['res_i']]
	df['Xid_j'] = [d_resid[i] for i in df['res_j']]
	df['abs_value'] = abs(df[0])
	return df





if __name__ == '__main__' :
	import argparse, ast 
	parser = argparse.ArgumentParser(description='ipython -i ~/FDA/C36-Pro-Pro-ed3/extension/analysis.py  -- -f ../Pro_lip.gro -x ../Pro_lip.xtc -d pair_SC.dat -sh [201,521,521]', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-f', '--gro', help='input gro file', required=True)
	parser.add_argument('-x', '--xtc', help='input xtc files', required=True)
	parser.add_argument('-d', '--data', help='data file', required=True)
	parser.add_argument('-sh', '--shape', help='data file shape [nframe, nresidue, nresidue]', required=True)
	parser.add_argument("-cri", "--criteria", help="Correlation coefficient cutoff", default = 0.7)
	args = parser.parse_args()

	u = MDAnalysis.Universe(args.gro, args.xtc)
	bar = PearsonR_Pair(u)
	plt.plot(bar.rmsd)
	plt.show()
	bar.load_file_data(args.data, ast.literal_eval(args.shape))
	print bar.data.shape
	bar.pearsonr_coef()
	bar.hist()
	plt.show()
	
	d_resid, d_resname = Sim_to_Xtal(u, 73,261)
	df = Array2DataFrame(bar.coef, d_resid=d_resid, d_resname= d_resname)

	SumonXid = df.groupby(['res_i', 'Xid_i', 'name_i'])['abs_value'].sum()
	SumonXid.to_csv('SumAbsvalueOnRes_i.csv')
	
	Protein = u.select_atoms('protein')
	for i in Protein.residues.resids : 
	    Protein.select_atoms('resid ' + str(i)).bfactors = [SumonXid[i-1]] * len(Protein.select_atoms('resid ' + str(i)))
	Protein.write('./Fri_Absvalue.pdb')

	
	'''
	df[df['Xid_i'] == 276][df['abs_value']> 0.3]
	add_vmd_labels(u, sel_pairs(bar.coef, args.criteria), './re_checkWithCriterai' + str(args.criteria)+ '.tcl')

	d_resid, d_resname = analysis.Sim_to_Xtal(u, 73,261)
	df = Array2DataFrame(bar.coef, d_resid=d_resid, d_resname= d_resname)
	foo = df.groupby(['res_i', 'Xid_i', 'name_i'])['abs_value'].max()
	foo = foo.sort_values(ascending=False).reset_index()
	tmp =  df.loc[df['abs_value'] >= criteria]
	tmp2 = pd.DataFrame(tmp.groupby(['res_i'])['abs_value'].count())
	bfac = np.zeros(len(Protein.residues.resids))

	for i, j in tmp2.iterrows() :
	    bfac[i-1] = j['abs_value']

	
	## networkx##
	g = nx.Graph()
	for i, j in tmp.iterrows() :
		g.add_edge(j['res_i'], j['res_j'], weight=j['abs_value'])
	N_node = len(g.nodes())
	Sn = sorted(g.nodes())
	color = [i[1] for i in g.degree_iter()]
	pos=nx.shell_layout(g, nlist = [Sn[:N_node/9],Sn[N_node/9:(N_node*3)/9],Sn[N_node*3/9:(N_node*6)/9],Sn[(N_node*6)/9:]] )
	#pos=nx.spring_layout(g, weight = 'weight')
	edges = g.edges()
	weights = [g[u][v]['weight'] for u,v in edges]
	plt.figure(figsize=(10,10))
	nx.draw(g,with_labels=True,node_size=800, node_color= color , pos=pos, alpha = 0.5,  cmap=plt.cm.Blues, edge_color = np.random.rand(len(edges)),edge_cmap = plt.cm.gist_rainbow)
	plt.savefig('./FDA_Pro_Pro/Node_Map' + '.png', dpi = 200)
	plt.cla()

	plt.figure(figsize=(10,10))
	nx.draw(g,with_labels=False,node_size=800, node_color= color , pos=pos, alpha = 0.5,  cmap=plt.cm.Blues, edge_color = np.random.rand(len(edges)),edge_cmap = plt.cm.gist_rainbow) 
	nx.draw_networkx_labels(g,pos=pos,labels = {i : str(d_amino_1letter[d_resname[i]])+str(d_resid[i]) for i in g.nodes()}, font_size=7)
	plt.savefig('./FDA_Pro_Pro/Node_MapXtal' + '.png', dpi = 200)
	plt.cla()
	plt.close()
	'''

##	d_resid, d_resname = Sim_to_Xtal(u, 73, 261)
