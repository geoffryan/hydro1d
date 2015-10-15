#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "par.h"
#include "geom.h"
#include "initial.h"

int set_initial(struct parList *pars)
{
    int err = 0;
    int choice = pars->init;

    if(choice == 0)
    {
        initial_value = &initial_uniform;
    }
    else if(choice == 1)
    {
        initial_value = &initial_shocktube;
    }
    else if(choice == 2)
    {
        initial_value = &initial_shocktube_transverse;
    }
    else if(choice == 3)
    {
        initial_value = &initial_isentrope;
    }
    else if(choice == 4)
    {
        initial_value = &initial_disc;
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
        double x = CM(xm,xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}
