#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "hydro.h"
#include "timestep.h"

void substep(struct grid *g, double dt, struct parList *pars)
{
    int nx = g->nx;
    int nq = g->nq;

    //Calculate new slopes.
    reconstruction(g);

    //Solve Riemann problems.
    
    //Add Sources.
    
    //Boundary Conditions.

    calc_prim(g, pars);
    calc_cons(g, pars);
}
