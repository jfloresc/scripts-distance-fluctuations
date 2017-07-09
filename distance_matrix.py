#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import distance 

t = md.load('stripped_gpa_wat.nc',top='step5_charmm2amber_stripped.pdb')
pdbref = md.load('step5_charmm2amber_stripped.pdb')
topology = pdbref.topology

#atoms_to_keep = [a.index for a in topology.atoms if a.name == 'CA']
# second selection method; VMD style selection
atoms_to_keep = topology.select('name CA')

t.restrict_atoms(atoms_to_keep)
pdbref.restrict_atoms(atoms_to_keep)

nume = distance.distance_fluc(t.xyz.reshape(t.n_frames,t.n_atoms*3)*10., t.n_frames, t.n_atoms)
np.savetxt("distance_matrix.dat",nume)

plt.figure(1)

nn = nume.shape[0]
plt.pcolor(nume, cmap=plt.cm.seismic, vmin=np.amin(nume), vmax=np.amax(nume))
axes = plt.gca()
axes.set_xlim([0,nn])
axes.set_ylim([0,nn])
axes.invert_yaxis()
axes.xaxis.tick_top()

# sets new tick labels
z = np.arange(70, 97, 2) # start residue index, #final residue index  + 1, stride
labels = np.concatenate((z,z),axis=0) 

# figure sets initial labels from 0 to nn-1
# sets ticks for the first monomer
z2 = np.arange(0, 27, 2) #  start, final + 1, stride
# sets ticks for the second monomer
z3 = np.arange(27, nn, 2) # start, final + 1, stride

x = np.concatenate((z2,z3),axis=0) 
y = np.concatenate((z2,z3),axis=0) 

plt.xticks(x, labels, fontsize=9)# rotation='vertical')
plt.yticks(y, labels, fontsize=9)# rotation='vertical')

plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
cbar = plt.colorbar()
cbar.set_label(r'Pairwise Distance Fluctuations ($\AA$)')
nameout='distance_matrix.png'
plt.savefig(nameout)
