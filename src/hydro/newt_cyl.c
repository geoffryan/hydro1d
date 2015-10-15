#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../hydro.h"

void prim2cons_newt_cyl(double *prim, double *cons, double r, double dV,
                        struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vr = prim[VX1];
    double vp = prim[VX2];
    double gam = pars->gammalaw;

    cons[DDD] = rho * dV;
    cons[TAU] = (0.5*rho*(vr*vr+r*r*vp*vp) + P/(gam-1)) * dV;
    cons[SX1] = rho * vr * dV;
    cons[SX2] = rho * r*r*vp * dV;
}

void cons2prim_newt_cyl(double *cons, double *prim, double r, double dV,
                        struct parList *pars)
{
    double rho = cons[DDD] / dV;
    double en = cons[TAU] / dV;
    double mr = cons[SX2] / dV;
    double mp = cons[SX1] / dV;
    double gam = pars->gammalaw;

    prim[RHO] = rho;
    prim[PPP] = (gam-1.0)*(en - 0.5*(mr*mr+mp*mp/(r*r))/rho);
    prim[VX1] = mr / rho;
    prim[VX2] = mp / (r*r*rho);
}

void flux_newt_cyl(double *prim, double *F, double r, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vr = prim[VX1];
    double vp = prim[VX2];
    double gam = pars->gammalaw;

    F[DDD] = rho * vr;
    F[TAU] = (0.5*rho*(vr*vr+r*r*vp*vp) + gam/(gam-1.0)*P) * vr;
    F[SX1] = rho*vr*vr + P;
    F[SX2] = rho*r*r*vp*vr;
}

void add_source_newt_cyl(double *prim, double *cons, double r, double dVdt, 
                            struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vp = prim[VX2];

    cons[SX1] += (P/r + rho*r*vp*vp) * dVdt;
}

void wave_speeds_newt_cyl(double *prim1, double *prim2, double *sL, double *sR,
                            double *sC, double r, struct parList *pars)
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

double mindt_newt_cyl(double *prim, double x, double dx, double cw, 
                        struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double vr = prim[VX1];
    double gam = pars->gammalaw;

    double cs = sqrt(gam*P/rho);
    
    double sl = fabs(vr + cs - cw);
    double sr = fabs(vr - cs - cw);
    double s = sl > sr ? sl : sr;

    double dt = dx / s;

    return dt;
}

double grid_V_newt_cyl(double *prim, double x, struct parList *pars)
{
    return prim[VX1];
}
