# Hydro1D #

Solver for 1D Hyperbolic Conservation Laws.  This is a testbed for ideas to be
implemented in larger codes, practice in building codes like this, and
hopefully a place to do some useful science like generate equilibrated initial
conditions for other projects.

-Geoff

## Installation ##

To install, just type `make`!

## Running ##

`Hydro1D` takes a single command line argument: a parameter file with the
desired setup specified. Sample parameter files are included in `parfiles/`.
To run the cartesian isentropic wave example, from the project directory run:

    $ bin/hydro1d parfiles/isentrope.par

### Visualization Requirements ###

Some simple python visualization scripts are included in `vis/`. They require
`numpy` and `matplotlib` python libraries to function correctly.

## Code Structure ##

The `Hydro1D` executable lives in `bin/`. Source code is contained in `src/`.
Sample parameter files are in `parfiles/` and visualization scripts in `vis/`.

Most of the code is modular. Boundary conditions, hydro systems, riemann solvers,
and more can all be varied at run time. Each module contains one or more
generic function pointers which are set to the appropriate specific functions
during initialization by a `set_...()` function.  For example `timestep.h`
declares the `timestep` function `void (*timestep)(struct grid *, double, struct parList *)`.  This pointer is initialized by the 
`set_timestep(struct parList *pars)` function.  Possible options are `step_fe`,
`step_rk2_mp`, `step_rk3_tvd`, etc. all defined in `timestep.c`.

The source code is organized as follows:

* `main.c`: Contains `main()`. Sets paramters and runs the code.
* `boundary.h/c`: Boundary conditions.  Currently implements fixed
    (e.g. Dirichlet), outflow (e.g. zero-gradient), and geodesic (relativistic
    only). 
* `geom.h/c`: Geometrical functions. Defines `CM` (center of mass of cell), `DA`
    (area of cell-cell interface), and `DV` (cell volume). Implements cartesian,
    cylindrical, and spherical geometry.
* `grid.h/c`: Defines `struct grid`, a container for the simulation domain and
    data. Sets the spatial reconstruction. Possible reconstruction schemes
    are piecewise constant and piecewise linear. Various grid utility functions.
* `hydro.h/c`: Sets hydro equations and defines data indices `RHO`, `PPP`, etc.
    * `hydro/newt_cart.c`: Compressible, adiabatic, euler equations in cartesian
        geometry.
    * `hydro/newt_cyl.c`: Compressible, adiabatic, axisymmetric euler equations
        in cylindrical geometry.
    * `hydro/newt_sph.c`: Compressible, adiabatic, spherically symmetric euler
        equations in spherical geometry.
    * `hydro/rel_cart.c`: Compressible, adiabatic, relativistic euler equations
        in cartesian geometry.
    * `hydro/rel_metric.c`: Compressible, adiabatic, relativistic euler equations
        in arbitrary geometry.
* `initial.h/c`: Initial conditions.







