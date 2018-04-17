#!/usr/bin/env python
# plots Calpha Atom correlation matrix
# by Jose Flores-Canales, 01/2017

from __future__ import print_function
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def get_cov(A,K):
  MA=np.mean(A.T,axis=1)
  M = (A-MA).T # subtract the mean (along columns)
  print(M.shape)
#1050, 14000
  covM=M.dot(M.T)
  covM=covM/K
  np.savetxt("covariance.dat",covM)
  return covM

def correlation_numpy(A, n_frames, n_atoms):
  #3000, 54, 3
  print("size of traj: ", A.shape)
  MA=np.mean(A, axis=0)
  print("mean(A): ", MA.shape)
  M = (A-MA) # subtract the mean (along columns)
#14000, 1050
  print("A-MA.T: ", M.shape)
#1050,
  cov_nn = np.tensordot(M,M.T, axes=([2,0],[0,2]))/n_frames
  diag = np.copy(cov_nn.diagonal())
  for i in xrange(n_atoms):
    for j in xrange(n_atoms):
      cov_nn[i,j] = cov_nn[i,j]/np.sqrt(diag[i]*diag[j])
  print(diag)
  return cov_nn   

def get_average_structure(traj):
  n_frames = traj.n_frames 
  n_atoms = traj.n_atoms
  average = np.zeros((n_atoms,3))
  for frame in xrange(n_frames):
    average += traj.xyz[frame] 
  return average/n_frames

t = md.load('stripped_gpa_G79L_wat.nc',top='step5_charmm2amber_stripped.pdb')
pdbref = md.load('step5_charmm2amber_stripped.pdb')
atom_to_keep = [a.index for a in t.topology.atoms if a.name == 'CA']
t.restrict_atoms(atom_to_keep)
pdbref.restrict_atoms(atom_to_keep)

pdbref.xyz = get_average_structure(t)
#print(pdbref.xyz)
rmsds = md.rmsd(t,pdbref,0)
index = np.argmin(rmsds)
pdbref.xyz = t.xyz[index] 
rmsds = md.rmsd(t,pdbref,0)
#for i in xrange(rmsds.shape[0]):
#  print(rmsds[i])
#index = np.argmin(rmsds)
t.superpose(pdbref)
#t.save_dcd('T0EG5_aligned.dcd')
pdbref.save_pdb('average.pdb')
nume = correlation_numpy(t.xyz*10, t.n_frames, t.n_atoms)
print(np.amax(nume), np.amin(nume))
np.savetxt("covariance.dat",nume)

plt.figure(1)

nn = nume.shape[0]
plt.pcolor(nume, cmap=plt.cm.seismic, vmin=-1, vmax=1)
axes = plt.gca()
axes.set_xlim([0,nn])
axes.set_ylim([0,nn])

# DEFINE LABELS RANGE 
# Residue 70 to 97, every 2 residues
z = np.arange(70, 97, 2)
labels = np.concatenate((z,z),axis=0) 
#z2 = np.arange(0, 27, 2)
#z3 = np.arange(27, nn, 2)
#x = np.concatenate((z2,z3),axis=0) 
#y = np.concatenate((z2,z3),axis=0) 
x=np.arange(0,nn,2)
y=np.arange(0,nn,2)
plt.xticks(x, labels, fontsize=9)# rotation='vertical')
plt.yticks(y, labels, fontsize=9)# rotation='vertical')
#plt.xticks(np.arange(70, 97, 2))
#plt.yticks(np.arange(70, 97, 2))
plt.xlabel('Residue')
plt.ylabel('Residue')
cbar = plt.colorbar()
cbar.set_label('Normalized Correlation')
nameout='correlation_np'
plt.savefig(nameout+'.png')
