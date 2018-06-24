#!/usr/bin/env python

import sys,os
import MDAnalysis
import numpy as np
import timeit, argparse
start = timeit.default_timer()

#         Set up parameters and locate all input files       #

parser = argparse.ArgumentParser(description='FDA--Protein-Lpid interaction use as : python Lipid_contacted1.py -f trj.gro -x trj.xtc -v ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
### required parameters
parser.add_argument('-f', '--gro', help='input gro file', required=True)
parser.add_argument('-x', '--xtc', help='input xtc files', required=True)
### optional parameters
parser.add_argument('-r', '--r_cutoff', help='distance cutoff in Angstrom', default = 6)
parser.add_argument('-b', help='First frame to read from trajectory', default = 0, nargs='?')
parser.add_argument('-e', help='Last frame to read from trajectory', default = -1, nargs='?')
parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")

args = parser.parse_args()
print "\n\n\n"
print args, "\n\n"

gro = args.gro
xtc = args.xtc
r_cutoff = args.r_cutoff

##### Create directories for files ####
if not os.path.exists('./Lipid_contact_VJ'):
	os.makedirs('Lipid_contact_VJ')

u = MDAnalysis.Universe(gro,xtc)
nframes = len(u.trajectory)
Protein = u.select_atoms('protein')

print "number of frame is", nframes, "\n"
print "program is ruuning \n"

check1 = timeit.default_timer()
print "time to setup " + str(check1 - start )
######################################################
##              Loop over time                      ##
######################################################
# Run on each new frame:
### set time step parameter
ContactC_t = []
ContactP_t = []
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
	num_C_around, num_P_around  = np.zeros(len(Protein.residues.resids)), np.zeros(len(Protein.residues.resids))
	for k in Protein.residues.resids : ##change from index
		a_lip = u.select_atoms("not protein" + " and around " + str(r_cutoff) + " resid " + str(k))
		foo = [name[0] for name in a_lip.names if name not in ['C1', 'C2', 'C3', 'C11','C12', 'C13', 'C14', 'C15']] ## do not count fist C of each chains and C1*
		num_C_around[k-1] = foo.count('C')
		num_P_around[k-1] = foo.count('P')
	check2 = timeit.default_timer()
	if args.verbose :
		print "time to loop all residues " + str(check2 - check4 )
	ContactC_t.append(num_C_around)
	ContactP_t.append(num_P_around)

ContactC_t = np.array(ContactC_t)
ContactP_t = np.array(ContactP_t)
ContactC_av, ContactP_av  = np.average(ContactC_t, axis = 0), np.average(ContactP_t, axis = 0)
np.savetxt("./Lipid_contact_VJ/Lipid_Ccontact.csv", ContactC_t)
np.savetxt("./Lipid_contact_VJ/Lipid_Pcontact.csv", ContactP_t)

for i in Protein.residues.resids :
    Protein.select_atoms('resid ' + str(i)).bfactors = [ContactC_av[i-1]] * len(Protein.select_atoms('resid ' + str(i)))
#print BackBone.bfactors to PDB file
Protein.write('./Lipid_contact_VJ/Lipid_contactC.pdb')

for i in Protein.residues.resids :
    Protein.select_atoms('resid ' + str(i)).bfactors = [ContactP_av[i-1]] * len(Protein.select_atoms('resid ' + str(i)))
#print BackBone.bfactors to PDB file
Protein.write('./Lipid_contact_VJ/Lipid_contactP.pdb')

f = open('./Lipid_contact_VJ/log' + '.dat', 'w')
f.write("Lipid-contact analyis was done with following parameter" + str(args) + "\n Output are Lipid_contact.pdb and Lipid_contact.csv")

