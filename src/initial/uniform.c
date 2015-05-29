#include <stdlib.h>
#include "../par.h"
#include "../grid.h"
#include "../hydro.h"
#include "../initial.h"

void initial_uniform(double *prim, double x, struct parList *pars)
{
    int q;
    double rho = pars->initPar1;
    double P = pars->initPar2;
    double vx = pars->initPar3;
    double vy = pars->initPar4;

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = vx;
    prim[VX2] = vy;

    for(q = pars->nc; q < pars->nc+pars->np; q++)
        prim[q] = 0.0;
}
