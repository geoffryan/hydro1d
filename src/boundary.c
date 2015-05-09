#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "par.h"
#include "initial.h"
#include "boundary.h"

int set_bc(struct parList *pars)
{
    int err = 0;
    int choiceInner = pars->bcInner;
    int choiceOuter = pars->bcOuter;

    if(choiceInner == 0)
        bc_inner = &bc_fixed_inner;
    else if(choiceInner == 1)
        bc_inner = &bc_outflow_inner;
    else
    {
        err++;
        printf("ERROR - Invalid choice for inner BC: %d\n", choiceInner);
    }

    if(choiceOuter == 0)
        bc_outer = &bc_fixed_outer;
    else if(choiceOuter == 1)
        bc_outer = &bc_outflow_outer;
    else
    {
        err++;
        printf("ERROR - Invalid choice for outer BC: %d\n", choiceOuter);
    }

    return err;
}

void bc_fixed_inner(struct grid *g, struct parList *pars)
{
    int i;
    int nq = g->nq;
    int ng = g->ng;

    for(i=0; i<ng; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}

void bc_fixed_outer(struct grid *g, struct parList *pars)
{
    int i;
    int nx = g->nx;
    int nq = g->nq;
    int ng = g->ng;

    for(i=nx-ng; i<nx; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}

void bc_outflow_inner(struct grid *g, struct parList *pars)
{
    int q;
    int nq = g->nq;

    for(q=0; q<nq; q++)
        g->prim[0*nq+q] = g->prim[1*nq+q];
}

void bc_outflow_outer(struct grid *g, struct parList *pars)
{
    int q;
    int nq = g->nq;
    int nx = g->nx;

    for(q=0; q<nq; q++)
        g->prim[(nx-1)*nq+q] = g->prim[(nx-2)*nq+q];
}
