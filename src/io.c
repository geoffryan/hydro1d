#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "io.h"

void io_init(struct io *iom, struct parList *par)
{
    iom->nTot = par->nChkpt;
    iom->current = 0;
    iom->tnext = par->tmin;
    iom->tmin = par->tmin;
    iom->tmax = par->tmax;
}

void io_out(struct io *iom, struct grid *g)
{
    if(g->t >= iom->tnext && iom->current <= iom->nTot)
    {
        char filename[128];

        sprintf(filename, "grid_%05d.txt", iom->current);

        io_print_grid_ascii(g, filename);

        iom->current += 1;
        iom->tnext = iom->current * (iom->tmax-iom->tmin)/(iom->nTot)
                        + iom->tmin;
        if(iom->current == iom->nTot)
            iom->tnext = iom->tmax;
    }
}

void io_print_grid_ascii(struct grid *g, char *filename)
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
