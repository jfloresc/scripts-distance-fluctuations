#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv)<2:
  print("plot_distance_2d.py matrix_updiag matrix_lowdiag")
  sys.exit()
matrix_up = sys.argv[1]
matrix_low = sys.argv[2]

m_up = np.loadtxt(matrix_up)
m_low = np.loadtxt(matrix_low)

nx_up = m_up.shape[0] 
ny_up = m_up.shape[1] 
nx_low = m_low.shape[0] 
ny_low = m_low.shape[1]
print(nx_up,ny_up,nx_low,ny_low)

if (nx_up != ny_up) or (nx_up != nx_low) or (ny_up != ny_low):
  print("matrices hve different dimensions")
  exit()

data = np.zeros((nx_up, nx_up),dtype=float)
print(data.shape, m_up.shape, m_low.shape)

for i in xrange(nx_up):
  for j in xrange(nx_up):
    if i <= j:
      data[i,j] = m_up[i,j]
    else:
      data[i,j] = m_low[i,j]
    
plt.figure(1)

plt.pcolor(data, cmap=plt.cm.seismic, vmin=np.amin(data), vmax=np.amax(data))
#plt.pcolor(data, cmap=plt.cm.seismic, vmin=-1.0, vmax=1.0)
axes = plt.gca()
axes.set_xlim([0,nx_up])
axes.set_ylim([0,nx_up])
axes.invert_yaxis()
axes.xaxis.tick_top()

# sets new tick labels
z = np.arange(70, 97, 2) # start residue index, #final residue index  + 1, stride
labels = np.concatenate((z,z),axis=0) 

# figure sets initial labels from 0 to nx_up-1
# sets ticks for the first monomer
z2 = np.arange(0, 27, 2) #  start, final + 1, stride
# sets ticks for the second monomer
z3 = np.arange(27, nx_up, 2) # start, final + 1, stride

x = np.concatenate((z2,z3),axis=0) 
y = np.concatenate((z2,z3),axis=0) 
plt.xticks(x, labels, fontsize=9)# rotation='vertical')
plt.yticks(y, labels, fontsize=9)# rotation='vertical')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
cbar = plt.colorbar()
cbar.set_label(r'Pairwise Distance Fluctuations ($\AA$)')
nameout='distance_matrix_2d.png'
plt.savefig(nameout)
