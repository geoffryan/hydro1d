#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "geom.h"
#include "boundary.h"
#include "hydro.h"
#include "movement.h"
#include "riemann.h"
#include "timestep.h"

//Local Functions

void add_fluxes(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars);
void add_sources(struct grid *g, int a, int b, int nq, double dt, 
                struct parList *pars);
void move_grid(struct grid *g, double rkfac1, double rkfac2, double dt,
                struct parList *pars);

//Definitions

void substep(struct grid *g, double rkfac1, double rkfac2, double dt,
                struct parList *pars)
{
    int nx = g->nx;
    int nq = g->nq;

    //Set grid velocities.
    calc_grid_movement(g, pars);

    //Calculate new slopes.
    reconstruction(g);

    //Solve Riemann problems.
    add_fluxes(g, 2, nx-2, nq, dt, pars);
    
    //Add Sources.
    add_sources(g, 1, nx-2, nq, dt, pars);

    //Move Grid.
    move_grid(g, rkfac1, rkfac2, dt, pars);
    
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
        double xL = CM(xLL,x);
        double xR = CM(x,xRR);
        double primL[nq], primR[nq], F[nq];
        for(q=0; q<nq; q++)
        {
            primL[q] = g->prim[(i-1)*nq+q] + g->grad[(i-1)*nq+q] * (x-xL);
            primR[q] = g->prim[  i  *nq+q] + g->grad[  i  *nq+q] * (x-xR);
        }
        riemann_flux(primL, primR, F, nq, x, g->w[i], pars);
        double dA = DA(x);
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
        double x = CM(xm,xp);
        double dV = DV(xm,xp);

        add_source(&(g->prim[nq*i]), &(g->cons[nq*i]), x, dV*dt, pars);
    }
}

void move_grid(struct grid *g, double rkfac1, double rkfac2, double dt, 
                struct parList *pars)
{
    int i;
    int nx = g->nx;
    for(i=0; i<nx+1; i++)
        g->x[i] = rkfac1*g->x[i] + rkfac2*g->x_rk[i] + g->w[i] * dt;
}
