#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../hydro.h"

enum{RHO, PPP, VXX};
enum{DDD, TAU, SXX};

void prim2cons(double *prim, double *cons, double x, double dV,
                struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double v = prim[VXX];
    double gam = pars->gammalaw;

    cons[DDD] = rho * dV;
    cons[TAU] = (0.5*rho*v*v + P/(gam-1)) * dV;
    cons[SXX] = rho * v * dV;
}

void cons2prim(double *cons, double *prim, double x, double dV,
                struct parList *pars)
{
    double rho = cons[DDD] / dV;
    double en = cons[TAU] / dV;
    double m = cons[SXX] / dV;
    double gam = pars->gammalaw;

    prim[RHO] = rho;
    prim[PPP] = (gam-1.0)*(en - 0.5*m*m/rho);
    prim[VXX] = m / rho;
}

void flux(double *prim, double *F, double x, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double v = prim[VXX];
    double gam = pars->gammalaw;

    F[DDD] = rho * v;
    F[TAU] = (0.5*rho*v*v + gam/(gam-1.0)*P) * v;
    F[SXX] = rho*v*v + P;
}

void add_source(double *prim, double *cons, double x, double dV, 
                struct parList *pars)
{
}

void wave_speeds(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double x, struct parList *pars)
{
    double rho1 = prim1[RHO];
    double P1 = prim1[PPP];
    double v1 = prim1[VXX];
    double rho2 = prim2[RHO];
    double P2 = prim2[PPP];
    double v2 = prim2[VXX];
    double gam = pars->gammalaw;

    double cs1 = sqrt(gam*P1/rho1);
    double cs2 = sqrt(gam*P2/rho2);

    *sL = v1 - cs1;
    if(v2 - cs2 < *sL)
        *sL = v2 - cs2;

    *sR = v1 + cs1;
    if(v2 + cs2 > *sR)
        *sR = v2 + cs2;

    //TODO: Use proper HLLC wavespeed.
    *sC = 0.0;
}

double mindt(double *prim, double x, double dx, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double v = prim[VXX];
    double gam = pars->gammalaw;

    double cs = sqrt(gam*P/rho);
    double dt = dx/(cs + fabs(v));

    return dt;
}
