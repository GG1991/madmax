#!/usr/bin/python
#
# Homogenization examples by FEM 
# of 2D structured cell 
#
# Guido Giuntoli
#

import numpy as np
import math as math
from scipy.sparse import csr_matrix

nx = 3
ny = 3
lx = 1.0
ly = 1.0

A = np.zeros( (nx*ny*2,nx*ny*2) )

elem = np.zeros( ((nx-1)*(ny-1),4) )
coor = np.zeros( (nx*ny,2) )

for i in range(0, ny-1):
  for j in range(0,nx-1):
    elem[i*(nx-1)+j,0] = j    + i*nx
    elem[i*(nx-1)+j,1] = j+1  + i*nx
    elem[i*(nx-1)+j,2] = j+1  + (i+1)*nx
    elem[i*(nx-1)+j,3] = j    + (i+1)*nx

for i in range(0, ny):
  for j in range(0,nx):
    coor[i*nx+j,0] = j*lx/(nx-1)
    coor[i*nx+j,1] = i*ly/(ny-1)


xp = np.zeros( (4,2) )

xp[0,:] = [(-1/math.sqrt(3.0)+1)/2 , (-1/math.sqrt(3.0)+1)/2] 
xp[1,:] = [(+1/math.sqrt(3.0)+1)/2 , (-1/math.sqrt(3.0)+1)/2] 
xp[2,:] = [(+1/math.sqrt(3.0)+1)/2 , (+1/math.sqrt(3.0)+1)/2] 
xp[3,:] = [(-1/math.sqrt(3.0)+1)/2 , (+1/math.sqrt(3.0)+1)/2] 

wp = np.zeros( (4) )

wp[0] = 0.25
wp[1] = 0.25
wp[2] = 0.25
wp[3] = 0.25


print elem, "\n"
print coor, "\n"
print xp, "\n"
print wp, "\n"
