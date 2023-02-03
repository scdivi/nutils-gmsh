# nutils-gmsh
This is an example of using [gmsh](https://gmsh.info/) (with boundary tags) in [nutils](https://github.com/evalf/nutils).

In this script we solve the Laplace equation $$- k \quad \Delta T = 0$$ on a 
annulus ring domain $\Omega$ with Dirichlet boundary $\Gamma$,
subject to boundary conditions:

$$ T = T_1     \qquad {\rm on} \quad Γ_{\rm inner}, $$

$$ T = T_2     \qquad {\rm on} \quad Γ_{\rm outer}, $$

where k denotes the thermal conductivity with unit W/(m·K).
On the Dirichlet boundaries, the temperature is set to T1 
and T2, respectively. The exact solution is:

$$ T_{\rm exact} = \frac{ \ln{(r_1)} -\frac{1}{2} \ln{(x_0^2+x_1^2)} }{ \ln{(  r_1 )} - \ln{( r_2 )}}  ( T_2 - T_1 ) + T_1 .$$

## Instructions to create a `.geo` with boundary tags
We consider an annulus ring (see figure below) and label the inner (in blue) and outer (in red)
rings for applying the boundary conditions.

![Alt text](https://github.com/scdivi/nutils-gmsh/blob/main/assets/annulus_ring.png "Optional title")

Consider $r_1 = 0.1$ and $r_2 = 0.15$, construct points, lines, arcs, line loops and 
plane surfaces in a `.geo` file as following:

```
h  = 0.00125
r1 = 0.1;
r2 = 0.15;
Point(0) = {0,0,0,h};
Point(1) = {r1,0,0,h};
Point(2) = {r2,0,0,h};
Point(3) = {0,r2,0,h};
Point(4) = {0,r1,0,h};
Line(1) = {1,2};
Circle(2) = {2,0,3};
Line(3) = {3,4};
Circle(4) = {4,0,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
```

**Remark:** The order of point creation should be in anti-clockwise direction.
This makes sure that the normal is always pointing towards outward direction.

Note that mesh size can be identicataed using the $h$ parameter. Now tag the boundaries as following:

```
Physical Line("bottom") = {1};
Physical Line("outer")  = {2};
Physical Line("left")   = {3};
Physical Line("inner")  = {4};
Physical Surface("interior") = {1};
```

Use [gmsh](https://gmsh.info/) to generate `.msh` file using this `.geo` file.
For example, on a comman line run

```
gmsh <name_of_your_geo_file.geo> -2 -o <name_of_the_output_file>.msh
```
Here “-2” means it generates a 2D mesh.

## Importing `.msh` file into nutils

In a [nutils](https://github.com/evalf/nutils) script you can import the gmsh generated `.msh` and read the boundary tags as:

```
domain, geom = mesh.gmsh('<name_of_the_file>.msh')
boundary_integral = domain.boundary['inner'].integral('(T - T1)^2 dS' @ ps, degree=degree*2)
```

For a Laplace example see [script](https://github.com/scdivi/nutils-gmsh/blob/main/annulus.py). 
Run `python annulus.py main` to evaluate solution for $h=0.01$ (see figure below).

![Alt text](https://github.com/scdivi/nutils-gmsh/blob/main/assets/annulus_solution.png "Optional title")

Run `python annulus.py convergence` to perform convergence study of $L^2$ and  $H^1$ error norms.
For a first-order basis function this results in a 
slope of $2$ and $1$ respectively (see figure below).

![Alt text](https://github.com/scdivi/nutils-gmsh/blob/main/assets/annulus_convergence.png "Optional title")

