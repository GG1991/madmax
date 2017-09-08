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

  ke.fill(0.0)
  for gp in range(0, xp.shape[0]):

    B = np.zeros( (3,4*2) )
    for sh in range(0, 4):
      B[0,[sh*2+0, sh*2+1]] = [dsh[sh,0,gp], 0           ]
      B[1,[sh*2+0, sh*2+1]] = [0           , dsh[sh,1,gp]]
      B[2,[sh*2+0, sh*2+1]] = [dsh[sh,1,gp], dsh[sh,0,gp]]

    ke += reduce(np.dot,[B.transpose(),C,B])*wp[gp]*area

  return;
############################################################
#
# creates elemental residue
#
def elem_residue(e, elem, x, be):

  be.fill(0.0)
  for n in range(0, 4):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  elem_disp = x[index]
  stress_gp = np.zeros( 3 )

  B = np.zeros( (3,4*2) )
  for gp in range(0,4):

    for sh in range(0, dsh.shape[0]):
      B[0,[sh*2+0, sh*2+1]] = [dsh[sh,0,gp], 0           ]
      B[1,[sh*2+0, sh*2+1]] = [0           , dsh[sh,1,gp]]
      B[2,[sh*2+0, sh*2+1]] = [dsh[sh,1,gp], dsh[sh,0,gp]]

    strain_gp = np.dot(B,elem_disp)
    stress_gp = np.dot(C,strain_gp)
    #print "strain_gp",strain_gp
    #print "stress_gp",stress_gp
    be += np.dot(B.transpose(), stress_gp)*wp[gp]*area

  return;

############################################################
#
# main program
#
nx = 2
ny = 2 
lx = 1.0
ly = 1.0
hx = lx/(nx-1)
hy = ly/(ny-1)
area = hx*hy

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
  dsh[0,0,gp] = -0.5*(1-xp[gp,1])/hx; dsh[0,1,gp] = (1-xp[gp,0])*-0.5/hy
  dsh[1,0,gp] = +0.5*(1-xp[gp,1])/hx; dsh[1,1,gp] = (0+xp[gp,0])*-0.5/hy
  dsh[2,0,gp] = +0.5*(0+xp[gp,1])/hx; dsh[2,1,gp] = (0+xp[gp,0])*+0.5/hy
  dsh[3,0,gp] = -0.5*(0+xp[gp,1])/hx; dsh[3,1,gp] = (1-xp[gp,0])*+0.5/hy

# define constitutive tensor
nu = 0.3; E  = 1e6
#C = np.array([
#    [1         ,nu/(1-nu) ,0                  ],
#    [nu/(1-nu) ,1         ,0                  ],
#    [0         ,0         ,(1-2*nu)/(2*(1-nu))]
#])
#C = np.multiply(C,E*(1-nu)/((1+nu)*(1-2*nu)))
C = np.array([
    [1         ,nu        ,0                  ],
    [nu        ,1         ,0                  ],
    [0         ,0         ,(1-nu)/2           ]
])
C = np.multiply(C,E/(1-nu*nu))

J  = np.zeros( (nx*ny*2 + n_bc*2,nx*ny*2 + n_bc*2) )
x  = np.zeros( nx*ny*2 + n_bc*2 )
b  = np.zeros( nx*ny*2 + n_bc*2 )
dx = np.zeros( nx*ny*2 + n_bc*2)
ke = np.zeros( (4*2,4*2) )
be = np.zeros( 4*2 )
D  = np.zeros( (3,n_bc*2) )

nelm = (nx-1)*(ny-1)
nnod = nx*ny
elem = np.zeros( (nelm,4), dtype=np.int)
coor = np.zeros( (nnod,2) )

for i in range(0, ny-1):
  for j in range(0,nx-1):
    elem[i*(nx-1)+j,0] = j    + i*nx
    elem[i*(nx-1)+j,1] = j+1  + i*nx
    elem[i*(nx-1)+j,2] = j+1  + (i+1)*nx
    elem[i*(nx-1)+j,3] = j    + (i+1)*nx

for i in range(0, ny):
  for j in range(0,nx):
    coor[i*nx+j,0] = j*hx
    coor[i*nx+j,1] = i*hy

# Boundary nodes index
y0_ind = np.arange(2          ,nx*2-2,2)
x1_ind = np.arange((nx-1)*2   ,ny*nx*2,nx*2)
y1_ind = np.arange(((ny-1)*nx+1)*2,ny*nx*2-2,2)
x0_ind = np.arange(0          ,ny*nx*2,nx*2)

bc_nods = np.sort(np.concatenate((y0_ind, x1_ind, y1_ind, x0_ind))/2)
bc_inds = np.zeros(bc_nods.size*2, dtype=np.int)
for i in range(0, bc_nods.size):
  bc_inds[i*2+0] = bc_nods[i]*2+0
  bc_inds[i*2+1] = bc_nods[i]*2+1

print
print "bc_nods",bc_nods
print "bc_inds",bc_inds
print

for i in range(0, n_bc):
  D[0,i*2+0] = coor[bc_nods[i],0]   ; D[0,i*2+1] = 0
  D[1,i*2+0] = 0                    ; D[1,i*2+1] = coor[bc_nods[i],1]
  D[2,i*2+0] = coor[bc_nods[i],1]/2 ; D[2,i*2+1] = coor[bc_nods[i],0]/2

e_mac = np.array([1.0,0.0,0.0])
e_mac = e_mac.transpose()
print "D^T*e",np.transpose(np.dot(D.transpose(),e_mac))

#calculate elemental matrix
elem_matrix( ke )

# assembly J 
index = np.zeros(4*2, dtype=np.int)
for e in range(0, elem.shape[0]):
  for n in range( 0, 4):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  J[np.ix_(index,index)] += ke

# set BCs on J
J[bc_inds[:], nx*ny*2 + np.arange(bc_inds.size)] = -1.0;
J[nx*ny*2 + np.arange(bc_inds.size), bc_inds[:]] = +1.0;
print "J ",J

# assembly b 
b.fill(0.0)
for e in range(0, elem.shape[0]):
  elem_residue(e, elem, x, be)
  for n in range( 0, 4):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  b[index] += be

u_bc = x[bc_inds] 
lamb = x[nx*ny*2:]
b[bc_inds] -= lamb
b[nx*ny*2 + np.arange(bc_inds.size)] = u_bc - np.transpose(np.dot(D.transpose(),e_mac))
b = -b

dx = np.linalg.solve(J, b)
x  = x + dx
print 
print "|b|",np.linalg.norm(b),"\n"

b.fill(0.0)
for e in range(0, elem.shape[0]):
  elem_residue(e, elem, x, be)
  for n in range(0, 4):
    index[[n*2+0, n*2+1]] = [elem[e,n]*2+0,elem[e,n]*2+1]
  b[index] += be

u_bc = x[bc_inds] 
lamb = x[nx*ny*2:]
b[bc_inds] -= lamb
b[nx*ny*2 + np.arange(bc_inds.size)] = u_bc - np.transpose(np.dot(D.transpose(),e_mac))
b = -b
dx = np.linalg.solve(J, b)
x = x + dx
print
print "|b|",np.linalg.norm(b),"\n"
print "u_bc", u_bc
print "b",b
print "forces",x[nnod*2:], "\n"

# plot the matrix
plt.matshow(J)
#plt.show()

#print "elem",elem, "\n"
#print "coor",coor, "\n"

############################################################
