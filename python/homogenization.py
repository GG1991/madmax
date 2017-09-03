#!/usr/bin/python
#
# Homogenization examples by FEM 
# of 2D structured cell 
#
# Guido Giuntoli
#

import numpy as np
import math as math

############################################################
def elem_matrix(ke):

  xp = np.zeros( (4,2) )
  xp[0,:] = [(-1/math.sqrt(3.0)+1)/2 , (-1/math.sqrt(3.0)+1)/2] 
  xp[1,:] = [(+1/math.sqrt(3.0)+1)/2 , (-1/math.sqrt(3.0)+1)/2] 
  xp[2,:] = [(+1/math.sqrt(3.0)+1)/2 , (+1/math.sqrt(3.0)+1)/2] 
  xp[3,:] = [(-1/math.sqrt(3.0)+1)/2 , (+1/math.sqrt(3.0)+1)/2] 

  wp = np.ones( (4) )
  wp = np.multiply(wp,0.25)

  dsh = np.zeros( (4,2,4) ) # num_sh, x_dir, num_gp
  for gp in range(0, xp.shape[0]):
    dsh[0,0,gp] = -1.0*(1-xp[gp,1]); dsh[0,1,gp] = (1-xp[gp,0])*-1.0
    dsh[1,0,gp] = +1.0*(1-xp[gp,1]); dsh[1,1,gp] = (0+xp[gp,0])*-1.0
    dsh[2,0,gp] = +1.0*(0+xp[gp,1]); dsh[2,1,gp] = (0+xp[gp,0])*+1.0
    dsh[3,0,gp] = -1.0*(0+xp[gp,1]); dsh[3,1,gp] = (1-xp[gp,0])*+1.0

  # define constitutive tensor
  nu = 0.3; E  = 1e6
  C = np.array([
  [1         ,nu/(1-nu) ,0                  ],
  [nu/(1-nu) ,1         ,0                  ],
  [0         ,0         ,(1-2*nu)/(2*(1-nu))]
  ])
  C = np.multiply(C,E*(1-nu)/((1+nu)*(1-2*nu)))

  for gp in range(0, xp.shape[0]):

    B = np.zeros( (3,4*2) )
    for sh in range(0, dsh.shape[0]):
      B[0,[sh*2+0, sh*2+1]] = [dsh[sh,0,gp], 0           ]
      B[1,[sh*2+0, sh*2+1]] = [0           , dsh[sh,0,gp]]
      B[2,[sh*2+0, sh*2+1]] = [dsh[sh,1,gp], dsh[sh,0,gp]]

    ke += reduce(np.dot,[B.transpose(),C,B])

  print ke, "\n"
  return;

############################################################

nx = 3
ny = 3
lx = 1.0
ly = 1.0

A  = np.zeros( (nx*ny*2,nx*ny*2) )
ke = np.zeros( (4*2,4*2) )

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




print elem, "\n"
print coor, "\n"

elem_matrix( ke )

