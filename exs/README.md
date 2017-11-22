### macro_alone

Example of a boundary condition problem using the macro code. 

Needs : `PETSC`, `SLEPC`, `ParMETIS`, `MPI`, `Gmsh`

The "-eigen" or a "-normal" execution can be selected.

```bash
cd macro_alone
```

## Eigensystem

Lets solve the eigen-system associated with this problem. First we run `gmsh` and
generate `cube_2d.msh`. You can change the mesh parameters, for example, the number 
of elements in *x* and *y* directions ( `vim cube_2d.geo`):

```bash
N = 30;
```

then:

```bash
gmsh -2 cube_2d.geo
```

will generate `cube_2d.msh`.

time to run the eigensystem ! Only configure in `run.sh` the option `-eigen` to appear
as an argument, `-normal` should *not* appear for this case. You can select the number 
of eigenvalues to compute with the option `-eps_nev 2`. Remember to comfigure your `mpiexec` 
path correctly in the line `MPIEXEC="/home/guido/libs/openmpi-install/bin/mpiexec"`. Finally
do:

```bash
./run.sh
```
The results are stored in files `macro_eigen_*.pvtu` and you can see the solution in paraview 

<img src="../doc/sputnik-man/figures/macro_alone_a.jpg"/>
<img src="../doc/sputnik-man/figures/macro_alone_b.jpg"/>

## Normal execution


Dimension: 2

Geometry : Square, we can use triangles or quads.

### micro_hom

Example of a simple homogenization that can be perform using :

  *  uniform strains 
  *  parallel mixture theory 
  *  serial mixture theory 

Dimension: 2

Geometry : the microstructure is a simple circular fiber.

### fron_fib

Example where the heterogeneities are solved by the direct method :

<img src="../doc/sputnik-man/figures/front_fib_d.jpg" 
alt="example with fibers embedded in a matrix" />

[//]: # (add this to control image size width="400" height="400")

Using a mesh like:

![mesh in the *direct* simulation](../doc/sputnik-man/figures/front_fib_c.jpg?raw=true "Title")

If we made a zoom:

<img src="../doc/sputnik-man/figures/front_fib_b.jpg" 
alt="zoom to the mesh" width="400" height="400" />

Then the solution can be obtain by doing:

```bash
./direct
```
if you check inside you can notice you can vary the material properties:

```bash
-material "MATRIX TYPE_0 1.0e7 1.0e6 0.3","FIBER TYPE_0 1.0e7 1.0e6 0.3","MICRO TYPE_1" \
```
and the boundary conditions:

```bash
-boundary "X0 11 0 0","X1 11 1 0" \
```
These are two different boundaries that are present in the domain.
In the first one the first string `X0` is the name of the physical entity specified equally
in the gmsh file the second number `11` means that you will impose Dirichlet on `y` and `x` 
direction respectivelly and `0 0` means to that those boundary condition follows the function of time 
number `0` for `x` and `0` for `y`.

These functions are defined in the next line
```bash
-function "0 2 0.0 0.0 1.0 0.0","1 2 0.0 0.0 1.0 0.01" \
```
Here we have two functions `0` and `1` (are the first number that appear on each one).
Each function has `2` values and we give them as `x0 y0 x1 y1` then the interpolation is linear in time.

Then the solution is given in plot in `macro_t_1.pvtu`

The stress distribution is

![stress distribution](../doc/sputnik-man/figures/front_fib_a.jpg?raw=true "Title")

and by the homogenization using:

  *  uniform strains 
  *  parallel mixture theory 
  *  serial mixture theory 

Dimension: 2

Geometry : A square piece of composite with fiber perpendicular to the plane


### tran_fib

Example where the heterogeneities are solved by the direct method and by the homogenization using:

  *  uniform strains 
  *  parallel mixture theory 
  *  serial mixture theory 

Dimension: 2

Geometry : A square piece of composite with fiber transversal to the plane
