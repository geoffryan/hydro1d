#include <stdlib.h>
#include <stdio.h>
#include "par.h"
#include "grid.h"

void make_grid(struct grid *g, struct parList *pars)
{
    //Set grid variables from parameter list and allocate
    //space for data arrays.
    g->nx_int = pars->nx;
    g->nx = pars->nx + 2*pars->nghost;
    g->ng = pars->nghost;
    g->xmin = pars->xmin;
    g->xmax = pars->xmax;
    g->nc = pars->nc;
    g->nq = pars->nc + pars->np;

    g->x = (double *) malloc((g->nx+1) * sizeof(double));
    g->prim = (double *) malloc(g->nx * g->nq * sizeof(double));
    g->cons = (double *) malloc(g->nx * g->nq * sizeof(double));
    g->cons_rk = (double *) malloc(g->nx * g->nq * sizeof(double));
    g->grad = (double *) malloc(g->nx * g->nq * sizeof(double));

    int i;
    double dx = (g->xmax - g->xmin) / g->nx_int;
    for(i=0; i<g->nx+1; i++)
        g->x[i] = g->xmin + (i - g->ng) * dx;
}

void print_grid(struct grid *g, char *filename)
{
    FILE *f;
    if(filename == NULL)
        f == stdout;
    else
        f = fopen(filename, "w");

    int i,q;

    int NQ = g->nq;

    fprintf(f, "%.12g", g->x[0]);
    for(i=1; i<g->nx+1; i++)
        fprintf(f, " %.12g", g->x[i]);
    fprintf(f, "\n");

    for(i=0; i<g->nx; i++)
    {
        fprintf(f, "%d", i);
        for(q=0; q<NQ; q++)
            fprintf(f, " %.12g", g->prim[NQ*i+q]);
        fprintf(f, "\n");
    }

    if(filename != NULL)
        fclose(f);
}

void free_grid(struct grid *g)
{
    free(g->x);
    free(g->prim);
    free(g->cons);
    free(g->cons_rk);
    free(g->grad);
}

void interpolate_constant(struct grid *g)
{}
