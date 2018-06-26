#!/usr/bin/env python

import sys,os
import pandas as pd
import matplotlib
import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import timeit, math, json, argparse, ast
from matplotlib import colors
from matplotlib import rc, rcParams

import analysis

rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'

start = timeit.default_timer()

###############################################################
#                                                             # 
#          Set up parameters and locate all input files       #
#                                                             #
###############################################################
### Read coordinate from gro file, atom typr and charge from itp file, sigma and epsilon from charmm36

parser = argparse.ArgumentParser(description='FDA--Protein interaction use as : ipython -i FDA_Pro-Pro_ed3 -- -f trj.gro -x trj.xtc -ff ff_charmm.itp -pt TREK2_charmm.itp -v', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
### required parameters
parser.add_argument('-f', '--gro', help='input gro file', required=True)
parser.add_argument('-x', '--xtc', help='input xtc files', required=True)
parser.add_argument('-ff', '--forcefield', help='input forcefile file', required=True)
parser.add_argument('-pt', '--portein_topol', help='input protein topology file', required=True)
### optional parameters
parser.add_argument('-r', '--r_cutoff', help='distance cutoff in Angstrom', default = 6)
parser.add_argument('-BB', '--backbone', help='Backbone atom name CA for Charmm36, BB for CG-martini', default = "CA", choices=["CA", "BB"], nargs='?')
parser.add_argument('-e_r', help='dielectric part in Coulomb potentail all-atom=1, CG-martini=15', default = 1, nargs='?')
parser.add_argument('-b', help='First frame to read from trajectory', default = 0, nargs='?')
parser.add_argument('-e', help='Last frame to read from trajectory', default = -1, nargs='?')
parser.add_argument('-pbc', help='fix PBC problem, this will slow down calulation time', default = 'Yes', choices=["Yes", "No"], nargs='?')
parser.add_argument('-saveall', help='save all data in atom-level, this will require a lot of space', default = 'No', choices=["Yes", "No"], nargs='?')
parser.add_argument('-savestep', help='frequency for saving data, 1 is every frame ', default = 1)
parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
parser.add_argument("-cri", "--criteria", help="Correlation coefficient cutoff", default = 0.7, nargs='?')


args = parser.parse_args()
print "\n\n\n"
print args, "\n\n"

gro = args.gro
xtc = args.xtc
ff = args.forcefield
protein_itp = args.portein_topol
BB = args.backbone 
e_r = args.e_r
fix_pbc = args.pbc
r_cutoff = args.r_cutoff
save_all_data = args.saveall
save_step = args.savestep
criteria = args.criteria
##### Create directories for files ####

if not os.path.exists('./FDA_Pro_Pro') :
	os.makedirs('./FDA_Pro_Pro')
'''	
if not os.path.exists('./FDA_Pro_Pro/Force_Stress'):
	os.makedirs('./FDA_Pro_Pro/Force_Stress')

if not os.path.exists('./FDA_Pro_Pro/Force_Network'):
	os.makedirs('./FDA_Pro_Pro/Force_Network')

if not os.path.exists('./FDA_Pro_Pro/Figure'):
	os.makedirs('./FDA_Pro_Pro/Figure')
'''

d_amino_1letter = { 'GLY' : 'G', 'ALA' : 'A', 'VAL' : 'V', 'LEU' : 'L', 'MET' : 'M', 'ILE' : 'I', 'SER' : 'S',
          'THR' : 'T', 'CYS' : 'C', 'PRO' : 'P', 'ASP' : 'D', 'GLN' : 'Q', 'PHE' : 'F', 'TYR' : 'Y',
          'TRP' : 'W', 'LYS' : 'K', 'HIS' : 'H', 'ASN' : 'N', 'GLU' : 'E', 'ARG' : 'R'} 
          
######## Read input files #########
itpProtein = pd.read_table(protein_itp, sep=r"\s*", engine='python')
ff = pd.read_table(ff, sep=r"\s*", engine='python')

################## Build up Force Filed Martix (ff_M) #########################
natom_build = len(itpProtein)
ntype_build = len(ff)
## dictionary for Charmm36  
ff_charmm36 = {j['name']:i for i,j in ff.iterrows()}

#### Check if all typr defined in dictionary
#for i in range(0, natom_build) :
#	print ff_charmm36[itp.type[i]]

#### Calculation combination rule (type2)
ff_M = np.zeros([len(ff_charmm36),len(ff_charmm36),2]) # [i,j,0=c6, 1=c12] 

for i in range(0, ntype_build):
	for j in range(0, ntype_build):
		epsilon = np.sqrt(ff.epsilon[i] * ff.epsilon[j])
		sigma = (ff.sigma[i] + ff.sigma[j])/float(2)
		c6 = 4 * epsilon * (sigma**6)
		c12 = 4 * epsilon * (sigma**12)
		ff_M[ ff_charmm36[ ff.name[i] ], ff_charmm36[ ff.name[j] ], 0 ] = c6
		ff_M[ ff_charmm36[ ff.name[i] ], ff_charmm36[ ff.name[j] ], 1 ] = c12

#### Check the calculation ####
#for i in range(0,ntype_build):
#	print ff_M[i,i]

####################################################################

#########################################
#            Load trajectory            #
#########################################
u = MDAnalysis.Universe(gro,xtc)
Natoms = len(u.atoms) ### total number of atoms (or particles) 
nframes = len(u.trajectory)
Protein = u.select_atoms('protein')
Nresidues = len(u.residues)
#### Add data in u.types and u.charges for all atomes in universe 
u.atoms[Protein.indices].types = itpProtein['type'].tolist()
u.atoms[Protein.indices].charges = itpProtein['charge'].tolist()

itpProtein['Part'] = 'SC'
BB_list = ['HT1','HT2','N','HN', 'CA', 'HA', 'C', 'O', 'OT1', 'OT2']
for i in BB_list :
    itpProtein.loc[itpProtein.atom == i, 'Part'] = 'BB'


#Frirj_av = np.zeros([len(u.residues),len(u.residues)])
print "number of frame is", nframes, "\n"
print "program is ruuning \n"
check1 = timeit.default_timer()
print "time to setup " + str(check1 - start )

###build collective matrix for |force| overtime 

F_norm_t, Pairs_t, Frirj_t, Fri_t, Frirj_SC_t = [], [], [], [], []

#Vi_t = np.zeros([nframes,len(u.residues)])

######################################################
##              Loop over time                      ##
######################################################
# Run on each new frame:
start = timeit.default_timer()
t_step = 0
### strat loop
if int(args.e) == -1 :
	args.e = len(u.trajectory)
elif int (args.e) > u.trajectory[-1].frame :
	print "\n\n Warning!! The last frame defined by user is longer than trajectory frame, the last frame will be last trajectory frame \n\n"
	args.e = u.trajectory[-1].frame + 1
else :
	args.e = int(args.e) + 1
for ts in u.trajectory[int(args.b):args.e] :
	if args.verbose :
		print ts 
### Build up Neighbour list Metric and all parameter that required #############
	check4 = timeit.default_timer()
	pairs, Charge_ij, C6_ij, C12_ij = [], [], [], []
	RattimeT = np.empty((0,3))
	for k in Protein.residues.resids : #index
		b = Protein.select_atoms("resid " + str(k))
		a = Protein.select_atoms("not resid " + str(k) + " and around " + str(r_cutoff) + " resid " + str(k)) ## bynum is starting from 1, index from 0, SO bynum = index+1
		if (len(a) != 0) :
			for n in range(len(b)) :
				r = a.positions - b.positions[n]
				RattimeT = np.vstack((RattimeT,r)) # for data less than 2 millions points
				for j in a.indices:
					i = b.indices[n]
					c6 = ff_M[ff_charmm36[u.atoms[i].type],ff_charmm36[u.atoms[j].type],0]
					c12 = ff_M[ff_charmm36[u.atoms[i].type],ff_charmm36[u.atoms[j].type],1]
					qq = u.atoms[i].charge * u.atoms[j].charge  ## Charge * Charge
					pairs.append([i,j]) ### atom indices
					Charge_ij.append(qq)
					C6_ij.append(c6)
					C12_ij.append(c12)
	C6_ij = np.array(C6_ij) 
	C12_ij = np.array(C12_ij)
	Charge_ij = np.array(Charge_ij)
	pairs = np.array(pairs)
	N_pairs = len(pairs)
	check2 = timeit.default_timer()
	if args.verbose :
		print " #pairs is  " + str(N_pairs) + "  pairs and #distances is " + str(len(RattimeT))
		print "max distance is " + str(RattimeT.max())
		print "time to loop all indices " + str(check2 - check4 )
		print "Last protein index is " + str(i)

####################################################################################
	RattimeT_nm = RattimeT/float(10)  ## A to nm
	r_norm = np.linalg.norm(RattimeT_nm, axis =1)
	if r_norm.max() > ts.dimensions[0]/20.0 : ##A to nm, distance larger than half of simulation box (2*10 = 20)
		con =  r_norm > ts.dimensions[0]/20.0 ## A to nm
		print 'There is PBC problem need to fix your trajectory \nuse -pbc option or reduce r_cuttoff \nnumber of pairs that need to fix is ' + str(len(RattimeT_nm[con]))
#		print RattimeT_nm[con], pairs[con]
		if fix_pbc == 'Yes' :
			for l in np.where(con)[0] :
#				print l
				tmpx = np.array([RattimeT_nm[l][0], RattimeT_nm[l][0]+ts.dimensions[0]/10.0, RattimeT_nm[l][0]-ts.dimensions[0]/10.0])
				tmpy = np.array([RattimeT_nm[l][1], RattimeT_nm[l][1]+ts.dimensions[1]/10.0, RattimeT_nm[l][1]-ts.dimensions[1]/10.0])
				RattimeT_nm[l][0] = tmpx[np.abs(tmpx).argmin()]
				RattimeT_nm[l][1] = tmpy[np.abs(tmpy).argmin()]
### re compute r_norm
			r_norm = np.linalg.norm(RattimeT_nm, axis =1)
### force from atom j on atom i
	FLJ_ij_C =  ((12*C12_ij/((r_norm)**13) - 6*C6_ij/((r_norm)**7)) ) 
	FLJ_ij = FLJ_ij_C.reshape(N_pairs,1) * (RattimeT_nm/r_norm.reshape(N_pairs,1))
	FC_ij_C = ((138.935485*Charge_ij/(r_norm**2)*e_r)) 
	FC_ij = FC_ij_C.reshape(N_pairs,1) * (RattimeT_nm/r_norm.reshape(N_pairs,1))
	F_norm = np.linalg.norm(FLJ_ij+FC_ij, axis = 1)
	F_direction = np.sign(FLJ_ij_C + FC_ij_C)

### Potential ####
#	VLJ_ij = C12_ij/(r_norm)**12) - C6_ij/(r_norm)**6)
#	VC_ij = (138.935485*Charge_ij/(r_norm)*e_r)
#	V_ij = VLJ_ij + VC_ij
######### sum pf norm(forces) RESIDUE j to residue i so called 'stress'  
### Sum of norm(forces) from RESIDUE j on residue i
	df = pd.DataFrame({"B_index": pairs[:,0] , "A_index": pairs[:,1] , "B_resid": u.atoms[pairs[:,0]].resids, "A_resid": u.atoms[pairs[:,1]].resids , "B_atom": u.atoms[pairs[:,0]].names , "A_atom": u.atoms[pairs[:,1]].names , "B_resname" : u.atoms[pairs[:,0]].resnames, "A_resname" : u.atoms[pairs[:,1]].resnames, "F": F_norm, "F_direction" : F_direction,"r": r_norm, "frame": [int(ts.frame)]*len(pairs)})
	df['B_Part'] = np.array(itpProtein['Part'])[pairs[:,0]]
	df['A_Part'] = np.array(itpProtein['Part'])[pairs[:,1]]
	
	### save data every frame
	Frirj = df.groupby(["B_resid", "A_resid"])['F'].sum().unstack() 
	Fri = df.groupby(["B_resid"])['F'].sum()
	Frirj_SC = df[df['B_Part'] == 'SC'][df['A_Part'] == 'SC']
	Frirj_SC = Frirj_SC.groupby(["B_resid", "A_resid"])['F'].sum().unstack() 
	Frirj_t.append(Frirj.fillna(0).as_matrix())
	Frirj_SC_t.append(Frirj_SC.fillna(0).as_matrix())
	Fri_t.append(Fri)

	t_step = t_step + 1

Frirj_SC_t = np.array(Frirj_SC_t)
Frirj_t = np.array(Frirj_t)
Fri_t = np.array(Fri_t)
        
ana = analysis.PearsonR_Pair(u)
ana.load_data(Frirj_SC_t)
ana.pearsonr_coef()
analysis.add_vmd_labels(u, analysis.sel_pairs(ana.coef, criteria), './FDA_Pro_Pro/vmd_full.tcl')
d_resid, d_resname = analysis.Sim_to_Xtal(u, 73,261)
df = analysis.Array2DataFrame(ana.coef, d_resid=d_resid, d_resname= d_resname)

foo = df.groupby(['res_i', 'Xid_i', 'name_i'])['abs_value'].max()
foo = foo.sort_values(ascending=False).reset_index()
tmp =  df.loc[df['abs_value'] >= criteria]
tmp2 = pd.DataFrame(tmp.groupby(['res_i'])['abs_value'].count())
bfac = np.zeros(len(Protein.residues.resids))

###################################################################
##                           Output Format                       ##
###################################################################
### Save to text ###
with file('./FDA_Pro_Pro/pair_SC.dat', 'w') as outfile:
    outfile.write('# Array shape: {0}\n'.format(Frirj_SC_t.shape))
    for data_slice in Frirj_SC_t:
        np.savetxt(outfile, data_slice )
        outfile.write('# New slice\n')

with file('./FDA_Pro_Pro/pair_all.dat', 'w') as outfile:
    outfile.write('# Array shape: {0}\n'.format(Frirj_t.shape))
    for data_slice in Frirj_t:
        np.savetxt(outfile, data_slice )
        outfile.write('# New slice\n')

np.savetxt('./FDA_Pro_Pro/Fri_t.dat',Fri_t)
        
## PDB ##
for i, j in tmp2.iterrows() :
    bfac[i-1] = j['abs_value']

for i in Protein.residues.resids : 
	Protein.select_atoms('resid ' + str(i)).bfactors = [bfac[i-1]] * len(Protein.select_atoms('resid ' + str(i)))
Protein.write('./FDA_Pro_Pro/Fri.pdb')

## networkx##
g = nx.Graph()
for i, j in tmp.iterrows() :
	g.add_edge(j['res_i'], j['res_j'], weight=j['abs_value'])
N_node = len(g.nodes())
Sn = sorted(g.nodes())
color = [i[1] for i in g.degree_iter()]
if N_node > 8 :
	pos=nx.shell_layout(g, nlist = [Sn[:N_node/9],Sn[N_node/9:(N_node*3)/9],Sn[N_node*3/9:(N_node*6)/9],Sn[(N_node*6)/9:]] )
else :
	pos=nx.spring_layout(g, weight = 'weight')
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

stop = timeit.default_timer()
t_run = stop - start
t_run = str(t_run)
print "               Program running time " + t_run + " s"
