#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "par.h"
#include "initial.h"

int set_initial(struct parList *pars)
{
    int err = 0;
    int choice = pars->init;

    if(choice == 0)
    {
        initial_value = &initial_shocktube;
    }
    else
    {
        printf("ERROR - Invalid Initial Condition choice: %d\n", pars->init);
        err = 1;
    }

    return err;
}

void initialize_grid(struct grid *g, struct parList *pars)
{
    int i;
    int nx = g->nx;
    int nq = g->nq;

    for(i=0; i<nx; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}
