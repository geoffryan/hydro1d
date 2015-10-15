#include <stdio.h>
#include <stdlib.h>
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
    for(i=1; i<nx; i++)
        g->w[i] = 0.5*(g->prim[nq*(i-1)+VX1] + g->prim[nq*i+VX1]);
    g->w[0] = g->w[1];
    g->w[nx] = g->w[nx-1];
}
