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

    double u2 = ux*ux + uy*uy;
    double u0 = sqrt(1.0 + u2);
    double rhoh = rho + gam/(gam-1.0)*P;

    cons[DDD] = rho * u0 * dV;
    cons[TAU] = (rho * u0*u2/(1.0+u0) + P * (u0*u0/(gam-1.0)+u2)) * dV;
    cons[SX1] = rhoh * u0 * ux * dV;
    cons[SX2] = rhoh * u0 * uy * dV;
}

void cons2prim_rel_cart(double *cons, double *prim, double x, double dV,
                            struct parList *pars)
{
    double ERR = 1.0e-12;
    double D = cons[DDD] / dV;
    double tau = cons[TAU] / dV;
    double Sx = cons[SX1] / dV;
    double Sy = cons[SX2] / dV;
    double gam = pars->gammalaw;

    double s2 = (Sx*Sx + Sy*Sy)/(D*D);
    double e = tau/D + 1.0;

    double c0 = e*e - s2;
    double c1 = (gam-1.0)*(gam-1.0)/(gam*gam) - 2*s2/gam + e*e;
    double c2 = -s2/(gam*gam);
    double c3 = -2*e*(gam-1.0)/gam;

    double u2 = prim[VX1]*prim[VX1] + prim[VX2]*prim[VX2];
    double du2;

    int i = 0;
    do
    {
        double f = c0*u2*u2 + c1*u2 + c2 + c3*u2*sqrt(1.0+u2);
        double df = 2*c0*u2 + c1 + c3*(1.0+1.5*u2)/sqrt(1.0+u2);

        if(f == 0.0)
            break;

        du2 = -f/df;
        u2 += du2;

        //printf("%d %.6g %.6g %.6g %.6g\n", i, u2, f, df, du2);
        i++;
    }
    while(fabs(du2) > ERR && i < 1000);

    if(i == 1000)
        printf("cons2prim failed. ERR = %.12lg\n", du2);

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
    double ux = prim[VX1];
    double uy = prim[VX2];
    double gam = pars->gammalaw;

    double u0 = sqrt(1.0 + ux*ux + uy*uy);
    double rhoh = rho + gam/(gam-1.0) * P;

    F[DDD] = rho*ux;
    F[TAU] = rhoh*u0*ux - rho*ux;
    F[SX1] = rhoh*ux*ux + P;
    F[SX2] = rhoh*uy*ux;
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
    double ux1 = prim1[VX1];
    double uy1 = prim1[VX2];
    double rho2 = prim2[RHO];
    double P2 = prim2[PPP];
    double ux2 = prim2[VX1];
    double uy2 = prim2[VX2];
    double gam = pars->gammalaw;

    double rhoh1 = rho1 + gam/(gam-1.0)*P1;
    double rhoh2 = rho2 + gam/(gam-1.0)*P2;

    double cs21 = gam*P1/rhoh1;
    double cs22 = gam*P2/rhoh2;

    double w1 = sqrt(1.0 + ux1*ux1 + uy1*uy1);
    double w2 = sqrt(1.0 + ux2*ux2 + uy2*uy2);

    double vx1 = ux1/w1;
    double vy1 = uy1/w1;
    double vx2 = ux2/w2;
    double vy2 = uy2/w2;

    double v21 = vx1*vx1 + vy1*vy1;
    double v22 = vx2*vx2 + vy2*vy2;

    double dv1 = sqrt(cs21 * (1.0 - (vx1*vx1 + vy1*vy1*cs21)) / (w1*w1));
    double dv2 = sqrt(cs22 * (1.0 - (vx2*vx2 + vy2*vy2*cs22)) / (w2*w2));

    double sL1 = (vx1*(1-cs21) - dv1) / (1.0 - v21*cs21);
    double sL2 = (vx2*(1-cs22) - dv2) / (1.0 - v22*cs22);
    double sR1 = (vx1*(1-cs21) + dv1) / (1.0 - v21*cs21);
    double sR2 = (vx2*(1-cs22) + dv2) / (1.0 - v22*cs22);

    *sL = sL1<sL2 ? sL1 : sL2;
    *sR = sR1>sR2 ? sR1 : sR2;

    //TODO: Use proper HLLC wavespeed.
    *sC = 0.0;
}

double mindt_rel_cart(double *prim, double x, double dx, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ux = prim[VX1];
    double uy = prim[VX2];
    double gam = pars->gammalaw;

    double rhoh = rho + gam/(gam-1.0)*P;
    double cs2 = gam*P/rhoh;

    double w = sqrt(1.0 + ux*ux + uy*uy);
    double vx = ux/w;
    double vy = uy/w;
    double v2=  vx*vx + vy*vy;

    double dv = sqrt(cs2 * (1.0 - (vx*vx + vy*vy*cs2)) / (w*w));

    double s;
    if(vx > 0)
        s = (vx*(1-cs2) + dv) / (1.0 - v2*cs2);
    else
        s = (vx*(1-cs2) - dv) / (1.0 - v2*cs2);

    double dt = dx / fabs(s);

    return dt;
}
