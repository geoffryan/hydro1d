#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../hydro.h"

void prim2cons_rel_cart(double *prim, double *cons, double x, double dV,
                            struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ux = prim[VX1];
    double uy = prim[VX2];
    double gam = pars->gammalaw;

    u2 = ux*ux + uy*uy
    double u0 = math.sqrt(1.0 + u2);
    double rhoh = rho + gam/(gam-1.0)*P;

    cons[DDD] = rho * u0 * dV;
    cons[TAU] = (rho * u0*u2/(1+u0) + P * (u0*u0/(gam-1.0)+u2)) * dV;
    cons[SX1] = rhoh * u0 * ux * dV;
    cons[SX2] = rhoh * u0 * uy * dV;
}

void cons2prim_rel_cart(double *cons, double *prim, double x, double dV,
                            struct parList *pars)
{
    double ERR = 1.0e-3;
    double D = cons[DDD] / dV;
    double tau = cons[TAU] / dV;
    double Sx = cons[SX1] / dV;
    double Sy = cons[SX2] / dV;
    double gam = pars->gammalaw;

    double S2 = Sx*Sx + Sy*Sy;
    double e = tau/D + 1.0;

    double c0 = e*e - s2;
    double c1 = (gam-1.0)*(gam-1.0)/(gam*gam) - 2*S2/gam + e*e;
    double c2 = -S2/(gam*gam);
    double c3 = -2*e*(gam-1.0)/gam;

    double u2 = prim[VX1]*prim[VX1] + prim[VX2]*prim[VX2];
    double du2;

    do
    {
        double f = c0*u2*u2 + c1*u2 + c2 + c3*u2*sqrt(1.0+u2);
        double df = 2*c0*u2 + c1 + c3*(1.0+1.5*u2)/sqrt(1.0+u2);

        if(f == 0.0)
            break;

        du2 = -df/f;
        u2 += du2;
    }
    while(fabs(du2) > ERR);

    double w = sqrt(1.0 + u2);
    double rho = D/w;
    double hmo = (w*e - w*w) / (u2 + 1.0/gam);

    double P = (gam-1.0)/gam * rho * hmo;
    double ux = Sx/(D*(hmo+1.0));
    double uy = Sy/(D*(hmo+1.0));

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = ux;
    prim[VX2] = uy;
}

void flux_rel_cart(double *prim, double *F, double x, struct parList *pars)
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

void add_source_rel_cart(double *prim, double *cons, double x, double dVdt, 
                struct parList *pars)
{
}

void wave_speeds_rel_cart(double *prim1, double *prim2, double *sL, double *sR,
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

double mindt_rel_cart(double *prim, double x, double dx, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double v = prim[VX1];
    double gam = pars->gammalaw;

    double cs = sqrt(gam*P/rho);
    double dt = dx/(cs + fabs(v));

    return dt;
}
