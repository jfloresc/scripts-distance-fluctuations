#!/usr/bin/env python

# modified 25.07.18
# added user interface

from __future__ import print_function
import os
import sys
import optparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.decomposition import PCA


def parse_cmdline(cmdlineargs):
    """parsing user defined parameters"""
    parser = optparse.OptionParser("Usage: python plot_2drmsd.py -p pdbfile -t trajfile")

    parser.add_option("-p", "--pdbfile", action="store", dest="pdbfile")
    parser.add_option("-t", "--trajfile", action="store", dest="trajfile")
    
    opts, _ = parser.parse_args(cmdlineargs)
    pdb_file = opts.pdbfile
    traj_file = opts.trajfile
    if (pdb_file is None) or (traj_file is None): # or (pdb_ref is None):
        parser.print_help()
        exit()
return pdb_file, traj_file


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
    print("mean(A) shape: ", MA.shape)
    M = (A-MA) # subtract the mean (along columns)
#14000, 1050
    print("M=(A-MA.T) shape: ", M.shape)
#1050,
    cov_nn = np.tensordot(M,M.T, axes=([2,0],[0,2]))/n_frames
    diag = np.copy(cov_nn.diagonal())
    for i in xrange(n_atoms):
        for j in xrange(n_atoms):
            cov_nn[i,j] = cov_nn[i,j]/np.sqrt(diag[i]*diag[j])
    print("cov_nn shape:",cov_nn.shape)
    return cov_nn   


def get_average_structure(traj):
    n_frames = traj.n_frames 
    n_atoms = traj.n_atoms
    average = np.zeros((n_atoms,3))
    for frame in xrange(n_frames):
        average += traj.xyz[frame] 
    return average/n_frames


def rmsd2d(referencepdb, trajfile):
    """ Calculates 2drmsd plot """
    #t = md.load('All_T310_restart_16_25.dcd',top='ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
    #pdbref = md.load('ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
    #t = md.load('./ampar_active_cry_CA.mdcrd',top='./ampar_active_cry_CA.pdb')
    #pdbref = md.load('./ampar_active_cry_CA.pdb')
    t = md.load(trajfile, top=referencepdb)
    pdbref = md.load(referencepdb)
    topology = pdbref.topology
    atom_to_keep = [a.index for a in topology.atoms if a.name == 'CA']
    at_to_res_map = {a.index:a.residue.index for a in topology.atoms if a.name =='CA'}
    print(len(atom_to_keep))
    #print(at_to_res_map)
    t.restrict_atoms(atom_to_keep)
    pdbref.restrict_atoms(atom_to_keep)


##############################################################
    similarities = [md.rmsd(t, t, i) for i in range(t.n_frames)]
    similarities = np.array(similarities).reshape((t.n_frames, t.n_frames)) * 10. # converting nm to Angstroms
    for i in range(t.n_frames):
        for j in range(i + 1, t.n_frames):
            similarities[j, i] = similarities[i, j]
    print("dimensions of the 2 RMSD matrix:", similarities.shape)
    print("2D RMSD:", similarities)

#############################################################

    #print(np.amax(nume), np.amin(nume))
    np.savetxt("2drmsd.dat", similarities)
    rmsd_max = float(np.max(similarities))

    plt.figure(figsize=(9,7))

    nn = similarities.shape[0]
    cmap = plt.get_cmap('seismic', 10)
    #bounds = [-1,-0.75,-.6,-0.475,0.475,.6,.75,1.0]
    bounds = np.arange(0, rmsd_max + 0.5, .5)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    #cmap = plt.cm.seismic
    cmaplist = [cmap(i) for i in range(cmap.N)]
    print(cmaplist)
    #plt.pcolor(nume, cmap=plt.cm.seismic, vmin=-1, vmax=1)
    img = plt.pcolor(similarities, cmap=cmap, norm=norm, vmin=0, vmax= rmsd_max)
    axes = plt.gca()
    axes.set_xlim([0,nn])
    axes.set_ylim([0,nn])
    # THIS INVERTS THE PLOT
    #axes.invert_yaxis()
    #axes.xaxis.tick_top()
    # frame stride to be used for labeling on the X and Y axes
    step = 200
    z = np.arange(0, nn, step)
    #labels = np.concatenate((z,z),axis=0) 
    labels = z
    x = z 
    y = z
    #x=np.arange(0,nn,2)
    #y=np.arange(0,nn,2)
    plt.xticks(x, labels, fontsize=9)# rotation='vertical')
    plt.yticks(y, labels, fontsize=9)# rotation='vertical')
    #plt.xticks(np.arange(70, 97, 2))
    #plt.yticks(np.arange(70, 97, 2))
    plt.xlabel('Frames Index')
    plt.ylabel('Frames Index')
    # THIS is used when the plot is inverted
    #axes.xaxis.set_label_position('top')
    cbar = plt.colorbar(img, cmap=cmap, norm=norm,spacing='proportional',ticks=bounds,boundaries=bounds)
    cbar.set_label('2D RMSD')
    nameout='2d_rmsd.png'
    plt.savefig(nameout, format='png', dpi=500)


###########################################################
# __main__
###########################################################
if __name__ == "__main__":
    REFPDB, TRAJ_FILENAME, = parse_cmdline(sys.argv[1:]) 
    rmsd2d(REFPDB, TRAJ_FILENAME)
