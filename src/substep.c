#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "boundary.h"
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
    add_fluxes(g, 2, nx-2, nq, dt, pars);
    
    //Add Sources.
    add_sources(g, 1, nx-2, nq, dt, pars);
    
    //Update prims.
    calc_prim(g, pars);

    //Boundary Conditions.
    bc_inner(g, pars);
    bc_outer(g, pars);

    //Re-update cons.
    calc_cons(g, pars);
}

//Local Definitions.

void add_fluxes(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars)
{
    int i, q;

    for(i=a; i<=b; i++)
    {
        double xLL = g->x[i-1]; 
        double x = g->x[i]; 
        double xRR = g->x[i+1];
        double xL = 0.5*(xLL+x);
        double xR = 0.5*(x+xRR);
        double primL[nq], primR[nq], F[nq];
        for(q=0; q<nq; q++)
        {
            primL[q] = g->prim[(i-1)*nq+q] + g->grad[(i-1)*nq+q] * (x-xL);
            primR[q] = g->prim[  i  *nq+q] + g->grad[  i  *nq+q] * (x-xR);
        }
        riemann_flux(primL, primR, F, nq, x, dt, pars);
        double dA = 1.0;
        for(q=0; q<nq; q++)
        {
            g->cons[(i-1)*nq+q] -= F[q] * dA * dt;
            g->cons[  i  *nq+q] += F[q] * dA * dt;
        }
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
