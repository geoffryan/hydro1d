#include <stdlib.h>
#include <math.h>
#include "../par.h"
#include "../grid.h"
#include "../hydro.h"
#include "../initial.h"

void initial_disc(double *prim, double x, struct parList *pars)
{
    double x0 = 1.0;

    double r0 = pars->initPar1;
    double rho0 = pars->initPar2;
    double rhoa = pars->initPar3;
    double P0 = pars->initPar4;
    double Pa = pars->initPar5;
    double M = pars->M;

    prim[RHO] = rho0 * pow(x/r0, rhoa);
    prim[PPP] = P0 * pow(x/r0, Pa);

    double X = 6*M/x - 1.0;

    if(x < 6*M)
    {
        prim[VX1] = (4*sqrt(2.0)*M/x - sqrt(X*X*X)) / (3.0 * (1.0-2.0*M/x));
        prim[VX2] = 2.0*sqrt(3.0)*M;
    }
    else
    {
        prim[VX1] = 2*M/(x*sqrt(1-3*M/x));
        prim[VX2] = sqrt(M*x / (1.0 - 3*M/x));
    }

    if(pars->np > 0)
    {
        prim[pars->nc] = x < x0 ? 1.0 : 0.0;

        int q;
        for(q = pars->nc+1; q < pars->nc+pars->np; q++)
            prim[q] = 0.0;
    }
}
