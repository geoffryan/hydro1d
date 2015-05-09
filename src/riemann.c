#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "par.h"
#include "hydro.h"
#include "riemann.h"

int set_riemann_solver(struct parList *pars)
{
    int choice = pars->riemann;
    int err = 0;

    if(choice == 0)
        riemann_flux = &lax_friedrichs_flux;
    else if(choice == 1)
        riemann_flux = &hll_flux;
    else
    {
        err = 1;
        printf("ERROR - Invalid choice for Riemann Solver: %d\n", choice);
    }

    return err;
}

void lax_friedrichs_flux(double primL[], double primR[], double F[], int nq,
                         double x, double dt, struct parList *pars)
{
    double sL, sR, sC, s;
    double UL[nq], UR[nq], FL[nq], FR[nq];

    prim2cons(primL, UL, x, 1.0, pars);
    prim2cons(primR, UR, x, 1.0, pars);
    flux(primL, FL, x, pars);
    flux(primR, FR, x, pars);
    wave_speeds(primL, primR, &sL, &sR, &sC, x, pars);

    s = fabs(sL)>fabs(sR) ? fabs(sL) : fabs(sR);

    int q;
    for(q=0; q<nq; q++)
        F[q] = 0.5*(FL[q] + FR[q] - s*(UR[q] - UL[q]));
}

void hll_flux(double primL[], double primR[], double F[], int nq,
                double x, double dt, struct parList *pars)
{
    double sL, sR, sC;
    double UL[nq], UR[nq], FL[nq], FR[nq];

    prim2cons(primL, UL, x, 1.0, pars);
    prim2cons(primR, UR, x, 1.0, pars);
    flux(primL, FL, x, pars);
    flux(primR, FR, x, pars);
    wave_speeds(primL, primR, &sL, &sR, &sC, x, pars);

    int q;
    if(sL > 0.0)
        for(q=0; q<nq; q++)
            F[q] = FL[q];
    else if(sR < 0.0)
        for(q=0; q<nq; q++)
            F[q] = FR[q];
    else
        for(q=0; q<nq; q++)
            F[q] = (sR*FL[q] - sL*FR[q] + sL*sR*(UR[q] - UL[q])) / (sR - sL);
}
