#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "hydro.h"
#include "par.h"
#include "timestep.h"

int set_timestep(struct parList *pars)
{
    int err = 0;
    int choice = pars->step;

    if(choice == 0)
        timestep = &step_fe;
    else if(choice == 1)
        timestep = &step_rk2;
    else
    {
        printf("ERROR - Invalid Timestep choice: %d\n", choice);
        err = 1;
    }

    return err;
}

double get_dt(struct grid *g, struct parList *pars)
{
    int i;
    int nx = g->nx;
    int nq = g->nq;

    double dtmin = 1.0e100;

    for(i=0; i<nx; i++)
    {
        double xp = g->x[i+1];
        double xm = g->x[i];
        double x = 0.5*(xm+xp);
        double dx = xp-xm;

        double dt = mindt(&(g->prim[i*nq]), x, dx, pars);

        dtmin = dt < dtmin ? dt : dtmin;
    }

    return dtmin;
}

void step_fe(struct grid *g, double dt, struct parList *pars)
{
    substep(g, dt, pars);
    g->t += dt;
}

//TODO: Implement this.
void step_rk2(struct grid *g, double dt, struct parList *pars)
{
    substep(g, dt, pars);
    g->t += dt;
}

void calc_cons(struct grid *g, struct parList *pars)
{
    int i;
    int nq = g->nq;
    int nx = g->nx;

    for(i=0; i<nx; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        double dV = xp-xm;

        prim2cons(&(g->prim[nq*i]), &(g->cons[nq*i]), x, dV, pars);
    }
}

void calc_prim(struct grid *g, struct parList *pars)
{
    int i;
    int nq = g->nq;
    int nx = g->nx;

    for(i=0; i<nx; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        double dV = xp-xm;

        cons2prim(&(g->cons[nq*i]), &(g->prim[nq*i]), x, dV, pars);
    }
}

