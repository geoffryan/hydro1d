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






