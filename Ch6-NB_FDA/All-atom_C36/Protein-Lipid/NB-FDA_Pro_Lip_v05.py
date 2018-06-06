#!/usr/bin/env python

'''
This is version 0.5 of NB-FDA Protein-Lipid interactions 
written by : Viwan Jarerattanachat 
Upadted :May 2018
'''

import sys,os
import pandas as pd
import MDAnalysis
import MDAnalysis.analysis.leaflet  
import numpy as np
import matplotlib.pyplot as plt
import timeit, math, json, argparse, ast
from matplotlib import rc, rcParams
rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'
start = timeit.default_timer()
###############################################################
#                                                             # 
#          Set up parameters and locate all input files       #
#                                                             #
###############################################################
### Read coordinate from gro file, atom typr and charge from itp file, sigma and epsilon from charmm36

parser = argparse.ArgumentParser(description='FDA--Protein-Lpid interaction use as : python FDAed5_Protein_Lipid.py -f trj.gro -x trj.xtc -ff ff_charmm.itp -pt TREK2_charmm.itp -lt POPC_charmm36.itp -lm ./POPC.C36.map.json -v', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
### required parameters
parser.add_argument('-f', '--gro', help='input gro file', required=True)
parser.add_argument('-x', '--xtc', help='input xtc files', required=True)
parser.add_argument('-ff', '--forcefield', help='input forcefile file', required=True)
parser.add_argument('-pt', '--portein_topol', help='input protein topology file', required=True)
parser.add_argument('-lt', '--lipid_topol', help='input lipid topology file', required=True)
parser.add_argument('-lm', '--lipid_map', help='input lipid mapping file', required=True)
### optional parameters
parser.add_argument('-r', '--r_cutoff', help='distance cutoff in Angstrom', default = 6)
parser.add_argument('-BB', '--backbone', help='Backbone atom name CA for Charmm36, BB for CG-martini', default = "CA", choices=["CA", "BB"], nargs='?')
parser.add_argument('-phosphate', help='phosphate group atom name', default = "P", choices=["P", "PO4"], nargs='?')
parser.add_argument('-e_r', help='dielectric part in Coulomb potentail all-atom=1, CG-martini=15', default = 1, nargs='?')
parser.add_argument('-b', help='First frame to read from trajectory', default = 0, nargs='?')
parser.add_argument('-e', help='Last frame to read from trajectory', default = -1, nargs='?')
parser.add_argument('-pbc', help='fix PBC problem, this will slow down calulation time', default = 'Yes', choices=["Yes", "No"], nargs='?')
parser.add_argument('-saveall', help='save all data in atom-level, this will require a lot of space', default = 'No', choices=["Yes", "No"], nargs='?')
parser.add_argument('-savestep', help='frequency for saving data, 1 is every frame ', default = 1)
parser.add_argument('-non_prot', help='list of non protein in the gro file, excluding water. List need to be in the same order as in gro file', default = ['POPC', 'POPG', 'POPE', 'POPS', 'CHOL'])
parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")

args = parser.parse_args()
print "\n\n\n"
print args, "\n\n"

gro = args.gro
xtc = args.xtc
ff = args.forcefield
lipid_itp = args.lipid_topol
protein_itp = args.portein_topol
BB = args.backbone 
e_r = args.e_r
fix_pbc = args.pbc
r_cutoff = args.r_cutoff
save_all_data = args.saveall
save_step = args.savestep
nonProtein_type = args.non_prot
P = args.phosphate

'''
gro = "./trj.gro"  ## must have protein+lipids
protein_itp = "./TREK2_charmm.itp"
lipid_itp ="./POPC_charmm36.itp"
ff = "./ff_charmm.itp"
xtc = "./trj.xtc"
nonProtein_type = ['POPC', 'POPG', 'POPE', 'POPS', 'CHOL'] ## need to be in order that match .gro file
r_cutoff = 6 ## A
BB = "CA"
e_r = 1
save_all_data = 'Yes'
save_step = 1
P = 'P'
'''

###############################################################
#                                                             # 
#               	Read all input files                  #
#                                                             #
###############################################################
#helix = np.genfromtxt(helix_file, usecols=(0))
itpProtein = pd.read_table(protein_itp, sep=r"\s*", engine='python')
itpLip = pd.read_table(lipid_itp, sep=r"\s*" , engine='python')
ff = pd.read_table(ff, sep=r"\s*", engine='python')

###############################################################
#                                                             # 
#                Output files control options                 #
#                                                             #
###############################################################

##### Create directories for files ####
if not os.path.exists('./FDA_Pro_Lip'):
	os.makedirs('FDA_Pro_Lip')

if not os.path.exists('./FDA_Pro_Lip/PNG'):
	os.makedirs('./FDA_Pro_Lip/PNG')

if not os.path.exists('./FDA_Pro_Lip/t_step'):
	os.makedirs('./FDA_Pro_Lip/t_step')
###############################################################
#                                                             # 
#                     Defined dictionaries                    #
#                                                             #
###############################################################
d_polar = { 'GLY' : 0 , 'ALA' : 0, 'VAL' : 0, 'LEU' : 0, 'MET' : 0, 'ILE' : 0, 'SER' : 1,
          'THR' : 1, 'CYS' : 0, 'PRO' : 0, 'ASP' : 3, 'GLN' : 1, 'PHE' : 0, 'TYR' : 1,
          'TRP' : 0, 'LYS' : 2, 'HIS' : 2, 'ASN' : 1, 'GLU' : 3, 'ARG' : 2}  ## 0 nonpolar 1 polar 2 basic polar 3 acidic polar

d_num2polar = {0: 'nonpolar', 1 : 'polar', 2 : 'basic polar', 3 : 'acidic polar'}

d_amino2num = { 'GLY' : 8, 'ALA' : 7, 'VAL' : 3, 'LEU' : 2, 'MET' : 12, 'ILE' : 1, 'SER' : 15,
          'THR' : 14, 'CYS' : 13, 'PRO' : 6, 'ASP' : 20, 'GLN' : 11, 'PHE' : 4, 'TYR' : 9,
          'TRP' : 5, 'LYS' : 18, 'HIS' : 10, 'ASN' : 16, 'GLU' : 19, 'ARG' : 17} 

d_num2amino = {d_amino2num[x]:x for x in d_amino2num.keys()}

#### Colour options for amino acids or polrity
d_color_polar = {0: 'Gray', 1: 'green', 2 : 'Blue', 3 : 'Red'}  ## 0 nonpolar 1 polar 2 basic polar 3 acidic polar

d_color_amino ={0:'Red' ,1:'Blue' ,2:'Green', 3:'Pink', 4:'Brown' , 5:'Orange', 6:'Gray', 7: 'Olive',
		8: 'Cyan', 9 : 'Plum', 10: 'magenta', 11:'lime', 12:'yellow', 13:'indigo', 14:'salmon', 15:'crimson',
		16: 'orchid', 17: 'gold', 18: 'silver', 19:'purple', 20:'skyblue' }

##### Build Lipid information #########
d_Lip_part = {'NC3': 'Head', 'PO4' : 'Head', 
              'GL' : 'GL',
              'C1A' : 'Tail', 'D2A' : 'Tail', 'C2A' : 'Tail', 'C3A' : 'Tail', 'C4A' : 'Tail', 'C5A' : 'Tail',
              'C1B' : 'Tail', 'C2B' : 'Tail', 'C3B' : 'Tail', 'C4B' : 'Tail','C5B' : 'Tail'}

with open(args.lipid_map) as POPC_map :
     tmp =  POPC_map.read()
tmp2 = ast.literal_eval(tmp.replace('\n',''))
d_POPC_map = {}
for x in tmp2.keys() :
    d_POPC_map[tmp2[x]] = range(x[0],x[1]+1)
del tmp, tmp2
POPC_map.close() 

POPC_type = np.chararray(len(itpLip), itemsize = 5)
for i in d_POPC_map.keys() :
    POPC_type[np.array(d_POPC_map[i])-1] = i
itpLip['part'] = POPC_type   
###############################################################
#                                                             # 
#          	Build up Force Filed Martix (ff_M)            #
#                                                             #
###############################################################
#natom_build = len(itpProtein)
ntype_build = len(ff)
## dictionary for Charmm36  
ff_charmm36 = {"C" :  0 ,"CA" :  1 ,"CC" :  2 ,"CD" :  3 ,"CE1" :  4 ,"CE2" :  5 ,"CM" :  6 ,"CP1" :  7 ,"CP2" :  8 ,"CP3" :  9 ,
		"CPA" :  10 ,"CPB" :  11 ,"CPH1" :  12 ,"CPH2" :  13 ,"CPM" :  14 ,"CPT" :  15 ,"CS" :  16 ,"CST" :  17 ,
		"CT1" :  18 ,"CT2" :  19 ,"CT3" :  20 ,"CY" :  21 ,"CT" :  22 ,"CT1x" :  23 ,"CT2x" :  24 ,"CT3x" :  25 ,
		"H" :  26 ,"HA" :  27 ,"HE1" :  28 ,"HE2" :  29 ,"HB" :  30 ,"HC" :  31 ,"HP" :  32 ,"HR1" :  33 ,
		"HR2" :  34 ,"HR3" :  35 ,"HS" :  36 ,"HT" :  37 ,"HA1" :  38 ,"HA2" :  39 ,"HA3" :  40 ,"N" :  41 ,
		"NC2" :  42 ,"NH1" :  43 ,"NH2" :  44 ,"NH3" :  45 ,"NP" :  46 ,"NPH" :  47 ,"NR1" :  48 ,"NR2" :  49 ,
		"NR3" :  50 ,"NY" :  51 ,"O" :  52 ,"OB" :  53 ,"OC" :  54 ,"OH1" :  55 ,"OM" :  56 ,"OS" :  57 ,
		"OST" :  58 ,"OT" :  59 ,"S" :  60 ,"SM" :  61 ,"SS" :  62 ,"SOD" :  63 ,"POT" :  64 ,"CLA" :  65 ,
		"CAL" :  66 ,"MG" :  67 ,"CES" :  68 ,"ZN" :  69 ,"FE" :  70 ,"HE" :  71 ,"NE" :  72 ,"CLAL" :  73 ,
		"DUM" :  74 ,"CAP" :  75 ,"FA" :  76 ,"CN" :  77 ,"NC" :  78 ,"OCA" :  79 ,"COA" :  80 ,"CF1" :  81 ,
		"CF2" :  82 ,"CF3" :  83 ,"HF1" :  84 ,"HF2" :  85 ,"F1" :  86 ,"F2" :  87 ,"F3" :  88 ,"C3" :  89 ,
		"CC1A" :  90 ,"CC1B" :  91 ,"CC2" :  92 ,"NS1" :  93 ,"NS2" :  94 ,"HOL" :  95 ,"HAL1" :  96 ,
		"HAL2" :  97 ,"HAL3" :  98 ,"HCL" :  99 ,"HL" :  100 ,"HEL1" :  101 ,"HEL2" :  102 ,"HBL" :  103 ,
		"CCL" :  104 ,"CL" :  105 ,"CTL1" :  106 ,"CTL2" :  107 ,"CTL3" :  108 ,"CTL5" :  109 ,"CEL1" :  110 ,
		"CEL2" :  111 ,"OBL" :  112 ,"OCL" :  113 ,"O2L" :  114 ,"OHL" :  115 ,"OSL" :  116 ,"OSLP" :  117 ,
		"NH3L" :  118 ,"NTL" :  119 ,"SL" :  120 ,"PL" :  121 ,"OWT3" :  122 ,"HWT3" :  123 ,"OWT4" :  124 ,
		"HWT4" :  125 ,"MWT4" :  126 ,"OWT5" :  127 ,"HWT5" :  128 ,"MWT5" :  129 ,"OW" :  130 ,"HW" :  131 }

#### Check if all typr defined in dictionary
#for i in range(0, natom_build) :
#	print ff_charmm36[itpProtein.type[i]]

#### Calculation combination rule (type2)
ff_M = np.zeros([132,132,2]) # [i,j,0=c6, 1=c12] # There are 132 atom type in Charmm36 (exclude DNA)
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

###############################################################
#                                                             # 
#               	Load trajectory                       #
#                                                             #
###############################################################
u = MDAnalysis.Universe(gro,xtc)
Natoms = len(u.atoms) ### total number of atoms (or particles) 
nframes = len(u.trajectory)
Protein = u.select_atoms('protein')
part = ['Protein']*len(Protein.indices) 
#### System information ###
N_mol = {x:len(u.select_atoms('resname ' + x ).residues) for x in nonProtein_type}
N_mol['Protein'] = len(Protein.residues)
#### Add data in u.types and u.charges for all atomes in universe 
u.atoms[Protein.indices].types = itpProtein['type'].tolist()
u.atoms[Protein.indices].charges = itpProtein['charge'].tolist()
for sel in nonProtein_type:
	if N_mol[sel] != 0 :
		u.select_atoms('resname ' + sel).types = np.array(itpLip[itpLip['resname']==sel]['type'].tolist()*N_mol[sel])
		u.select_atoms('resname ' + sel).charges = np.array(itpLip[itpLip['resname']==sel]['charge'].tolist()*N_mol[sel])
		part = part + itpLip[itpLip['resname']==sel]['part'].tolist()*N_mol[sel]
del sel
part = np.array(part)

## To do check if the program built intranal topoloty correctly
'''
add check if all type and chard apply correctly
u.select_atoms('resname ' + 'POPC').types[:134]
u.atoms[Protein.indices].types
change method above then using :
_.names and use dictionary to define part: head tail or PO4  
'''

##### Build protein information ############
seq = Protein.residues.resnames
polar = [d_polar[i] for i in seq]
aminonum = [d_amino2num[i] for i in seq]
polar_text = np.array([d_num2polar[i] for i in polar])

print "number of frame is", nframes, "\n"
print "number of total atoms is ", Natoms , "\nThere are", N_mol," (molecules, residues). \n"
print "program is ruuning \n"

L = MDAnalysis.analysis.leaflet.LeafletFinder(u, 'name ' + P)
if (L.groups(0).center_of_mass()[2] > L.groups(1).center_of_mass()[2]):
	upper = L.groups(0) #### upper
	lower = L.groups(1) #### lower
else:
	upper = L.groups(1) #### upper
	lower = L.groups(0) #### lower

leaflet_profile = {}
sel_upper = 'resid'
for rid in upper.resids :
    sel_upper = sel_upper + ' ' +str(rid)
    leaflet_profile[rid] = 'u'

sel_lower = 'resid'
for rid in lower.resids :
    sel_lower = sel_lower + ' ' +str(rid)   
    leaflet_profile[rid] = 'l' 
    
##### Create list(s) to collect variable overtime 
#Lipid = u.select_atoms('resname POPC')
Fri_t, PF_t, ATPF_t = [], [], []

check1 = timeit.default_timer()
print "time to setup " + str(check1 - start )
######################################################
##              Loop over time                      ##
######################################################
# Run on each new frame:
### set time step parameter
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
### Lood over all atoms ###
	check4 = timeit.default_timer()
	pairs, Charge_ij, C6_ij, C12_ij = [], [], [], []
	RattimeT = np.empty((0,3))
	for k in Protein.residues.resids : ##change from index
		b_pro = u.select_atoms("resid " + str(k))
		a_lip = u.select_atoms("not protein" + " and around " + str(r_cutoff) + " resid " + str(k))
		if (len(a_lip) != 0) :
			for n in range(len(b_pro)) :
				r = a_lip.positions - b_pro.positions[n]
				RattimeT = np.vstack((RattimeT,r)) # for data less than 2 millions points
				for j in a_lip.indices:
					i = b_pro.indices[n]
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

############     Calculations   #################
	RattimeT_nm = RattimeT/float(10.)  ## A to nm
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
	
	df = pd.DataFrame({"P_index": pairs[:,0] , "L_index": pairs[:,1] , "P_resid": u.atoms[pairs[:,0]].resids, "L_resid": u.atoms[pairs[:,1]].resids , "P_atom": u.atoms[pairs[:,0]].names , "L_atom": u.atoms[pairs[:,1]].names , "P_resname" : u.atoms[pairs[:,0]].resnames, "F": F_norm, "F_direction" : F_direction,"r": r_norm, "frame": [int(ts.frame)]*len(pairs)})
	df['Polarity'] = polar_text[u.atoms[pairs[:,0]].resids-1]
	df['Part']= part[pairs[:,1]]
	df['Part2'] = [d_Lip_part[x] for x in df['Part']]
	df['leaflet'] = [leaflet_profile[x] for x in df['L_resid']]
	PF = df.groupby(['P_resid', 'L_resid','Part','frame','Polarity', 'leaflet', "P_resname"])['F'].sum().reset_index()
	Fri_tmp = df.groupby(['P_resid'])['F'].sum().reset_index()
	Fri = np.zeros(N_mol['Protein'])
	for i in Fri_tmp['P_resid'] :
		Fri[i-1] = Fri_tmp[Fri_tmp['P_resid']==i]['F']

### save data into files (for each time step) ###
	if save_all_data == 'Yes' :
		if t_step%save_step == 0 :
			df.to_csv('./FDA_Pro_Lip/t_step/Pairs_resid' + str(ts.frame) + '.dat')

### keep time series data in list ###
	Fri_t.append(Fri)	
	PF_t.append(PF)
	ATPF_t.append(df)
	print 'Frame = ', ts.frame
	t_step = t_step + 1 

######################################################
##                  Stat over time                  ##
######################################################
Fri_t = np.array(Fri_t)
Fri_av = np.average(Fri_t, axis = 0)
Fri_std = np.std(Fri_t, axis = 0)

PF_all = pd.concat(PF_t, ignore_index= True)
ATPF_t = pd.concat(ATPF_t, ignore_index= True)
PF_av = PF_all.groupby(['P_resid', 'L_resid','Part','Polarity', 'P_resname'])['F'].mean().reset_index()
MolPF_t = PF_all.groupby(['P_resid', 'L_resid','frame','Polarity', 'P_resname'])['F'].sum().reset_index()
MolPF_av = MolPF_t.groupby(['P_resid', 'L_resid','Polarity', 'P_resname'])['F'].mean().reset_index()
PF_AminoCG_av = PF_all.groupby(['P_resname', 'Part'])['F'].mean().reset_index()
PF_AminoCG_av['Part2'] = [d_Lip_part[x] for x in PF_AminoCG_av['Part']]
PF_AminoMol_av = PF_AminoCG_av.groupby(['P_resname', 'Part2'])['F'].sum()
PF_AminoMol_av.unstack()[['Head','GL','Tail']].plot.bar()
plt.savefig('./FDA_Pro_Lip/Amino_Mol_av' + '.svg', dpi=300)
plt.close()
PF_AminoMol_av = PF_AminoMol_av.reset_index()
PF_AminoMol_av['Polarity'] = [d_num2polar[d_polar[x]] for x in PF_AminoMol_av['P_resname']]
PF_Polar_av = PF_AminoMol_av.groupby(['Polarity','Part2'])['F'].mean()
PF_Polar_av = PF_Polar_av.unstack()[['Head','GL','Tail']]
PF_Polar_av.plot.bar(rot= 0)
plt.savefig('./FDA_Pro_Lip/Polarity_Mol_av' + '.svg', dpi=300)
plt.close()
######################################################
##                Save data to files                ##
######################################################
np.savetxt("./FDA_Pro_Lip/FriasfuntionofTime.csv", Fri_t)
PF_all.to_csv('./FDA_Pro_Lip/Pairs_CG_t' + '.csv')
PF_av.to_csv('./FDA_Pro_Lip/Pairs_CG_av' + '.csv')
MolPF_t.to_csv('./FDA_Pro_Lip/Pairs_Mol_t' + '.csv')
MolPF_av.to_csv('./FDA_Pro_Lip/Pairs_Mol_av' + '.csv')
PF_AminoCG_av.to_csv('./FDA_Pro_Lip/AminoAcid_CG_av' + '.csv')
PF_AminoMol_av.to_csv('./FDA_Pro_Lip/AminoAcid_Mol_av' + '.csv')
PF_Polar_av.to_csv('./FDA_Pro_Lip/Polarity_Mol_av' + '.csv')

################# Visualisation (VMD) ###################
#### make PDB file for average [Sum|Force|] ########## 
for i in range(N_mol['Protein']) :
    Protein.select_atoms('resid ' + str(i+1)).bfactors = [Fri_av[i]] * len(Protein.select_atoms('resid ' + str(i+1)))
#print BackBone.bfactors to PDB file
Protein.write('./FDA_Pro_Lip/Fri.pdb')

####  SELECTION for VMD ######
## To do : change to class for multiple selection
con = Fri_av > np.percentile(Fri_av, 95)
Im_resid = (np.array(np.where(con)))[0]+1
file1 = open('./FDA_Pro_Lip/Simid_vis' + '.vmd', 'w')
file1.write("mol representation Licorice 0.300000 10.000000 10.000000\n")
file1.write("mol selection {resid ")
for j in Im_resid :
	i = int(j)
	file1.write("%s " %i)
file1.write("}\nmol addrep top")
file1.close()

stop = timeit.default_timer()
t_run = stop - start
t_run = str(t_run)
print "      Program running time " + t_run + " s"
