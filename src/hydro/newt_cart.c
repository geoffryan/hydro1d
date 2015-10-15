#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../hydro.h"

void prim2cons_newt_cart(double *prim, double *cons, double x, double dV,
                            struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vx = prim[VX1];
    double vy = prim[VX2];
    double gam = pars->gammalaw;

    cons[DDD] = rho * dV;
    cons[TAU] = (0.5*rho*(vx*vx+vy*vy) + P/(gam-1)) * dV;
    cons[SX1] = rho * vx * dV;
    cons[SX2] = rho * vy * dV;
}

void cons2prim_newt_cart(double *cons, double *prim, double x, double dV,
                            struct parList *pars)
{
    double rho = cons[DDD] / dV;
    double en = cons[TAU] / dV;
    double mx = cons[SX1] / dV;
    double my = cons[SX2] / dV;
    double gam = pars->gammalaw;

    prim[RHO] = rho;
    prim[PPP] = (gam-1.0)*(en - 0.5*(mx*mx+my*my)/rho);
    prim[VX1] = mx / rho;
    prim[VX2] = my / rho;
}

void flux_newt_cart(double *prim, double *F, double x, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vx = prim[VX1];
    double vy = prim[VX2];
    double gam = pars->gammalaw;

    F[DDD] = rho * vx;
    F[TAU] = (0.5*rho*(vx*vx+vy*vy) + gam/(gam-1.0)*P) * vx;
    F[SX1] = rho*vx*vx + P;
    F[SX2] = rho*vy*vx;
}

void add_source_newt_cart(double *prim, double *cons, double x, double dVdt, 
                struct parList *pars)
{
}

void wave_speeds_newt_cart(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double x, struct parList *pars)
{
    double rho1 = prim1[RHO];
    double P1 = prim1[PPP];
    double v1 = prim1[VX1];
    double rho2 = prim2[RHO];
    double P2 = prim2[PPP];
    double v2 = prim2[VX1];
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

double mindt_newt_cart(double *prim, double x, double dx, double cw, 
                        struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double v = prim[VX1];
    double gam = pars->gammalaw;

    double cs = sqrt(gam*P/rho);
    
    double sl = fabs(v + cs - cw);
    double sr = fabs(v - cs - cw);
    double s = sl > sr ? sl : sr;

    double dt = dx / s;

    return dt;
}
