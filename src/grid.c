#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "par.h"
#include "grid.h"

//Local Functions
double minmod(double a, double b, double c);

int set_reconstruction(struct parList *pars)
{
    int err = 0;
    int choice = pars->recon;

    if(choice == 0)
        reconstruction = &interpolate_constant;
    else if(choice == 1)
        reconstruction = &interpolate_plm;
    else
    {
        err = 1;
        printf("ERROR - Invalid Reconstruction choice: %d\n", choice);
    }

    return err;
}

void make_grid(struct grid *g, struct parList *pars)
{
    //Set grid variables from parameter list and allocate
    //space for data arrays.
    g->nx_int = pars->nx;
    g->nx = pars->nx + 2*pars->nghost;
    g->ng = pars->nghost;
    g->xmin = pars->xmin;
    g->xmax = pars->xmax;
    g->t = pars->tmin;
    g->nc = pars->nc;
    g->nq = pars->nc + pars->np;
    g->PLM = pars->plm;

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
        f = stdout;
    else
        f = fopen(filename, "w");

    int i,q;

    int NQ = g->nq;

    fprintf(f, "%s\n", VERSION);

    fprintf(f, "t %.12g\n", g->t);    
    fprintf(f, "x");
    for(i=0; i<g->nx+1; i++)
        fprintf(f, " %.12g", g->x[i]);
    fprintf(f, "\n");

    for(i=0; i<g->nx; i++)
    {
        fprintf(f, "%d", i);
        for(q=0; q<NQ; q++)
            fprintf(f, " %.12g", g->prim[NQ*i+q]);
        for(q=0; q<NQ; q++)
            fprintf(f, " %.12g", g->cons[NQ*i+q]);
        for(q=0; q<NQ; q++)
            fprintf(f, " %.12g", g->grad[NQ*i+q]);
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
{
    int i,j;
    int nx = g->nx;
    int nq = g->nq;

    for(i=0; i<nx; i++)
        for(j=0; j<nq; j++)
            g->grad[nq*i+j] = 0.0;
}

void interpolate_plm(struct grid *g)
{
    int i,j;
    int nx = g->nx;
    int nq = g->nq;
    double PLM = g->PLM;

    for(j=0; j<nq; j++)
    {
        double dR = 2*(g->prim[1*nq+j] - g->prim[0*nq+j]) / (g->x[2]-g->x[0]);
        g->grad[j] = dR;
    }
    for(i=1; i<nx-1; i++)
    {
        double idxL = 2.0 / (g->x[i+1] - g->x[i-1]);
        double idxC = 2.0 / (g->x[i+2] + g->x[i+1] - g->x[i] - g->x[i-1]);
        double idxR = 2.0 / (g->x[i+2] - g->x[i]);
        int id = nq*i;
        
        for(j=0; j<nq; j++)
        {
            double SL = (g->prim[id   +j] - g->prim[id-nq+j]) * idxL;
            double SC = (g->prim[id+nq+j] - g->prim[id-nq+j]) * idxC;
            double SR = (g->prim[id+nq+j] - g->prim[id   +j]) * idxR;
            
            g->grad[id+j] = minmod(PLM*SL, SC, PLM*SR);
        }
    }
    for(j=0; j<nq; j++)
    {
        double dR = 2*(g->prim[(nx-1)*nq+j] - g->prim[(nx-2)*nq+j])
                        / (g->x[nx]-g->x[nx-2]);
        g->grad[j] = dR;
    }
}

//Local definitions

double minmod(double a, double b, double c)
{
    double m;

    if(a*b < 0)
        m = 0;
    else if(fabs(a) < fabs(b))
        m = a;
    else
        m = b;

    if(b*c < 0)
        m = 0;
    else if(fabs(c) < fabs(m))
        m = c;

    return m;
}
