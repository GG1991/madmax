#!/usr/bin/python
#
# Homogenization examples by FEM 
# of 2D structured cell 
#
# Guido Giuntoli
#

import numpy as np
import matplotlib.pyplot as plt
import math as math

############################################################
#
# creates elemental matrix
#
def elem_matrix(ke):

  for gp in range(0, xp.shape[0]):

    B = np.zeros( (3,4*2) )
    for sh in range(0, dsh.shape[0]):
      B[0,[sh*2+0, sh*2+1]] = [dsh[sh,0,gp], 0           ]
      B[1,[sh*2+0, sh*2+1]] = [0           , dsh[sh,0,gp]]
      B[2,[sh*2+0, sh*2+1]] = [dsh[sh,1,gp], dsh[sh,0,gp]]

    ke += reduce(np.dot,[B.transpose(),C,B])

  return;
############################################################
#
# creates elemental residue
#
def elem_residue(e, elem, x, be):

  for n in range(0, 4):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  elem_disp = x[index]
  stress_gp = np.zeros( 3 )

  for gp in range(0, xp.shape[0]):

    B = np.zeros( (3,4*2) )
    for sh in range(0, dsh.shape[0]):
      B[0,[sh*2+0, sh*2+1]] = [dsh[sh,0,gp], 0           ]
      B[1,[sh*2+0, sh*2+1]] = [0           , dsh[sh,0,gp]]
      B[2,[sh*2+0, sh*2+1]] = [dsh[sh,1,gp], dsh[sh,0,gp]]

    stress_gp = np.dot(B,elem_disp)
    be += np.dot(B.transpose(), stress_gp)

  return;

############################################################
#
# main program
#
nx = 3
ny = 3
lx = 1.0
ly = 1.0

n_bc = ny*2 + (nx-2)*2

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

J  = np.zeros( (nx*ny*2 + n_bc*2,nx*ny*2 + n_bc*2) )
x  = np.zeros( nx*ny*2 + n_bc*2 )
dx = np.zeros( nx*ny*2 + n_bc*2)
b  = np.zeros( nx*ny*2 + n_bc*2 )
ke = np.zeros( (4*2,4*2) )
be = np.zeros( 4*2 )
D  = np.zeros( (3,n_bc*2) )

elem = np.zeros( ((nx-1)*(ny-1),4), dtype=np.int )
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

# Boundary conditions indeces
ux_x0_ind = np.arange(0          ,ny*nx*2,nx*2)
ux_x1_ind = np.arange((nx-1)*2   ,ny*nx*2,nx*2)
uy_x0_ind = np.arange(1          ,ny*nx*2,nx*2)
uy_x1_ind = np.arange((nx-1)*2+1 ,ny*nx*2,nx*2)

ux_y0_ind = np.arange(2          ,nx*2-2,2)
uy_y0_ind = np.arange(3          ,nx*2-2,2)
ux_y1_ind = np.arange(((ny-1)*nx+1)*2,ny*nx*2-2,2)
uy_y1_ind = np.arange(((ny-1)*nx+1)*2+1,ny*nx*2-2,2)
#print ux_x0_ind, "\n"
#print ux_x1_ind, "\n"
#print uy_x0_ind, "\n"
#print uy_x1_ind, "\n"
#print ux_y0_ind, "\n"
#print uy_y0_ind, "\n"
#print ux_y1_ind, "\n"
#print uy_y1_ind, "\n"

n_ind_bc = np.concatenate((ux_y0_ind, ux_x1_ind, ux_y1_ind, ux_x0_ind))/2
x_bc = np.concatenate((ux_y0_ind,uy_y0_ind,ux_x1_ind,uy_x1_ind,ux_y1_ind,uy_y1_ind,ux_x0_ind,uy_x0_ind))

print "n_ind_bc",n_ind_bc

for i in range(0, n_bc):
  D[0,i*2+0] = coor[n_ind_bc[i],0]   ; D[0,i*2+1] = 0
  D[1,i*2+0] = 0                     ; D[1,i*2+1] = coor[n_ind_bc[i],1]
  D[2,i*2+0] = coor[n_ind_bc[i],1]/2 ; D[2,i*2+1] = coor[n_ind_bc[i]/2,0]

#calculate elemental matrix
elem_matrix( ke )

# assembly J 
index = np.zeros(4*2, dtype=np.int)
for e in range(0, elem.shape[0]):
  for n in range( 0, elem.shape[1]):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  J[np.ix_(index,index)] += ke

# set BCs on J
J[x_bc[:], nx*ny*2 + np.arange(x_bc.size)] = 1e7;
J[nx*ny*2 + np.arange(x_bc.size), x_bc[:]] = 1e7;

e_mac = np.array([[1.0],[0.0],[0.0]])

# assembly b 
for e in range(0, elem.shape[0]):
  elem_residue(e, elem, x, be)
  for n in range( 0, elem.shape[1]):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  b[index] += be

ub = x[x_bc] 
print x_bc
print nx*ny*2 + x_bc
b[nx*ny*2 + np.arange(x_bc.size)] = ub - np.transpose(np.dot(D.transpose(),e_mac))
print b

# plot the matrix
plt.matshow(J)
plt.show()

print elem, "\n"
print coor, "\n"

############################################################
