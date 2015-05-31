#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "geom.h"
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
        timestep = &step_rk2_mp;
    else if(choice == 2)
        timestep = &step_rk2_tvd;
    else if(choice == 3)
        timestep = &step_rk3_tvd;
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
        double x = CM(xm,xp);
        double dx = xp-xm;

        double dt = mindt(&(g->prim[i*nq]), x, dx, pars);

        dtmin = dt < dtmin ? dt : dtmin;
    }

    return pars->cfl * dtmin;
}

void step_fe(struct grid *g, double dt, struct parList *pars)
{
    substep(g, dt, pars);
    g->t += dt;
}

void step_rk2_mp(struct grid *g, double dt, struct parList *pars)
{
    // The Midpoint RK2 Method.

    copy_to_rk(g);
    substep(g, 0.5*dt, pars);
    g->t += 0.5*dt;

    update_cons(g, 0.0, 1.0);
    substep(g, dt, pars);
    g->t += 0.5*dt;
}

void step_rk2_tvd(struct grid *g, double dt, struct parList *pars)
{
    // The TVD RK2 Method by Gottlied & Shu

    copy_to_rk(g);
    substep(g, dt, pars);

    update_cons(g, 0.5, 0.5);
    substep(g, 0.5*dt, pars);
    g->t += dt;
}

void step_rk3_tvd(struct grid *g, double dt, struct parList *pars)
{
    // The TVD RK3 Method by Shu & Osher

    copy_to_rk(g);
    substep(g, dt, pars);

    update_cons(g, 0.25, 0.75);
    substep(g, 0.25*dt, pars);

    update_cons(g, 2.0/3.0, 1.0/3.0);
    substep(g, 2.0/3.0*dt, pars);

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
        double x = CM(xm,xp);
        double dV = DV(xm,xp);

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
        double x = CM(xm,xp);
        double dV = DV(xm,xp);

        cons2prim(&(g->cons[nq*i]), &(g->prim[nq*i]), x, dV, pars);
    }
}

