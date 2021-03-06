#include <stdlib.h>
#include "../par.h"
#include "../grid.h"
#include "../hydro.h"
#include "../initial.h"

void initial_shocktube(double *prim, double x, struct parList *pars)
{
    double x0 = pars->initPar1;
    double rhoL = pars->initPar2;
    double PL = pars->initPar3;
    double vL = pars->initPar4;
    double rhoR = pars->initPar5;
    double PR = pars->initPar6;
    double vR = pars->initPar7;
    double x1 = pars->initPar8;

    if(x < x0)
    {
        prim[RHO] = rhoL;
        prim[PPP] = PL;
        prim[VX1] = vL;
        prim[VX2] = 0.0;
    }
    else
    {
        prim[RHO] = rhoR;
        prim[PPP] = PR;
        prim[VX1] = vR;
        prim[VX2] = 0.0;
    }

    if(pars->np > 0)
    {
        prim[pars->nc] = x < x1 ? 1.0 : 0.0;

        int q;
        for(q = pars->nc+1; q < pars->nc+pars->np; q++)
            prim[q] = 0.0;
    }
}
