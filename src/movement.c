#include <stdio.h>
#include <stdlib.h>
#include "geom.h"
#include "grid.h"
#include "hydro.h"
#include "par.h"
#include "movement.h"

int set_movement(struct parList *pars)
{
    int err = 0;
    int choice = pars->movement;

    if(choice == 0)
        calc_grid_movement = &move_none;
    else if(choice == 1)
        calc_grid_movement = &move_local;
    else
    {
        printf("ERROR = Invalid Movement choice: %d\n", choice);
        err = 1;
    }

    return err;
}

void move_none(struct grid *g, struct parList *pars)
{
    int nx = g->nx;
    int i;

    for(i=0; i<nx+1; i++)
        g->w[i] = 0.0;
}

void move_local(struct grid *g, struct parList *pars)
{
    int nx = g->nx;
    int nq = g->nq;

    int i;
    double vL, vR;
    vR = grid_V(&(g->prim[0]), CM(g->x[0],g->x[1]), pars);
    for(i=1; i<nx; i++)
    {
        vL = vR;
        vR = grid_V(&(g->prim[nq*i]), CM(g->x[i],g->x[i+1]), pars);
        g->w[i] = 0.5*(vL+vR);
    }
    g->w[0] = g->w[1];
    g->w[nx] = g->w[nx-1];
}
