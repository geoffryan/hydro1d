#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "hydro.h"
#include "riemann.h"
#include "timestep.h"

//Local Functions

void add_fluxes(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars);
void add_sources(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars);

//Definitions

void substep(struct grid *g, double dt, struct parList *pars)
{
    int nx = g->nx;
    int nq = g->nq;

    //Calculate new slopes.
    reconstruction(g);

    //Solve Riemann problems.
    add_fluxes(g, 1, nx-1, nq, dt, pars);
    
    //Add Sources.
    add_sources(g, 0, nx-1, nq, dt, pars);
    
    //Boundary Conditions.

    //Update prims.
    calc_prim(g, pars);
    calc_cons(g, pars);
}

//Local Definitions.

void add_fluxes(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars)
{
    int i;

    for(i=a; i<=b; i++)
    {
        double xL = g->x[i-1]; 
        double xC = g->x[i]; 
        double xR = g->x[i+1];
        double x1 = 0.5*(xL+xC);
        double x2 = 0.5*(xC+xR);
        riemann_solve(&(g->prim[(i-1)*nq]), &(g->cons[(i-1)*nq]), x1,
                &(g->prim[i*nq]), &(g->prim[i*nq]), x2, xC, dt, pars);
    }
}

void add_sources(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars)
{
    int i;

    for(i=a; i<=b; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        double dV = xp-xm;

        add_source(&(g->prim[nq*i]), &(g->cons[nq*i]), x, dV*dt, pars);
    }
}
