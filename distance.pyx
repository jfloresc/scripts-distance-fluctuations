#!/usr/bin/env python

import numpy as np
cimport numpy as np
DTYPE = np.float
ctypedef np.float_t DTYPE_t

def euclidian2(DTYPE_t vx1, DTYPE_t vx2, DTYPE_t vx3, DTYPE_t vy1, DTYPE_t vy2, DTYPE_t vy3):
  cdef DTYPE_t d2
  d2 = (vx1-vy1)*(vx1-vy1) + (vx2-vy2)*(vx2-vy2) + (vx3-vy3)*(vx3-vy3)
  return d2

def distance_fluc(np.ndarray A, int n_frames, int n_atoms):
  cdef np.ndarray dist = np.zeros([n_atoms, n_atoms], dtype=DTYPE)
  cdef np.ndarray dist2 = np.zeros([n_atoms, n_atoms], dtype=DTYPE)
  cdef np.ndarray f = np.zeros([n_atoms, n_atoms], dtype=DTYPE)
  cdef DTYPE_t s, s2, temp, temp2, vx1, vx2, vx3, vy1, vy2, vy3
  cdef int i, j, iframe
  cdef n_i, n_x
  
  n_i = A.shape[0] 
  n_x = A.shape[1] 
  print(n_i, n_x)
  for i in xrange(0,n_atoms*3,3):
    for j in xrange(i,n_atoms*3,3):
      s = 0.0
      s2 = 0.0
      for iframe in xrange(n_frames):
        vx1 = A[iframe][i]
        vx2 = A[iframe][i+1]
        vx3 = A[iframe][i+2]
        vy1 = A[iframe][j]
        vy2 = A[iframe][j+1]
        vy3 = A[iframe][j+2]
        temp2 = euclidian2(vx1, vx2, vx3, vy1, vy2, vy3)
        temp = np.sqrt(temp2)
        s += temp
        s2 += temp2
 
      dist[i/3, j/3] = s*s 
      dist[j/3, i/3] = s*s
      dist2[i/3, j/3] = s2 
      dist2[j/3, i/3] = s2
      
  dist = dist / ((n_frames)*(n_frames))
  dist2 = dist2 /( n_frames)
  for i in xrange(n_atoms):
    for j in xrange(n_atoms):
      #if (i==j):
        #print(i, j, dist2[i,j], dist[i,j])
      #else:
      #print(i, j, dist2[i,j], dist[i,j])
      f[i,j] = np.sqrt(dist2[i,j]-dist[i,j])
      f[j,i] = f[i,j]
      #print(i, j, f[i,j])
  return f   
