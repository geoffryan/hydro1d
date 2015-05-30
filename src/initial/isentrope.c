#include <stdlib.h>
#include <math.h>
#include "../par.h"
#include "../grid.h"
#include "../hydro.h"
#include "../initial.h"

// Initial condition for an isentropic wave.  See Section 4.6
// of the RAM paper (Zhang & MacFadyen 2006) for details.

void initial_isentrope(double *prim, double x, struct parList *pars)
{
    double rho0 = pars->initPar1; // Reference density
    double P0 = pars->initPar2;   // Reference pressure
    double x0 = pars->initPar3;   // Pulse location
    double L = pars->initPar4;    // Pulse width
    double a = pars->initPar5;    // Pulse strength

    double rho, P, vx, dx, K, GAM, cs, J;

    GAM = pars->gammalaw;
    K = P0 / pow(rho0, GAM);
    dx = (x-x0)/L;
    cs = sqrt(GAM*P0/rho0);

    J = -2.0 * cs / (GAM-1.0); 

    rho = rho0;
    P = P0;
    vx = 0.0;

    if(fabs(dx) < 1.0)
    {
        rho += a * rho0 * pow(dx*dx-1.0, 4);
        P = K * pow(rho, GAM);
        cs = sqrt(GAM*P/rho); 
        vx = J + 2.0 * cs / (GAM-1.0);
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = vx;
    prim[VX2] = 0.0;

    if(pars->np > 0)
    {
        prim[pars->nc] = x < x0 ? 1.0 : 0.0;

        int q;
        for(q = pars->nc+1; q < pars->nc+pars->np; q++)
            prim[q] = 0.0;
    }
}
