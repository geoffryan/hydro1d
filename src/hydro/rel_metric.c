#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../hydro.h"
#include "metric/schw_ks.h"

void prim2cons_rel_metric(double *prim, double *cons, double r, double dV,
                            struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ur = prim[VX1];
    double up = prim[VX2];
    double U[3];
    double gam = pars->gammalaw;
    double M = pars->M;

    frame_U(r, M, U);
    double w = sqrt(1.0 + IG3RR*ur*ur + IG3PP*up*up + 2*IG3RP*ur*up);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double W = AL * U[0];
    double rhoh = rho + gam/(gam-1.0)*P;
    double uU = u0*U[0] + ur*U[1] + up*U[2];

    cons[DDD] = rho * w * J3*dV;
    cons[TAU] = (-rhoh*uU*w - W*P - rho*w) * J3*dV;
    cons[SX1] = rhoh * w * ur * J3*dV;
    cons[SX2] = rhoh * w * up * J3*dV;
}

void cons2prim_rel_metric(double *cons, double *prim, double r, double dV,
                            struct parList *pars)
{
    double ERR = 1.0e-12;
    double D = cons[DDD] / dV;
    double tau = cons[TAU] / dV;
    double Sr = cons[SX1] / dV;
    double Sp = cons[SX2] / dV;
    double gam = pars->gammalaw;
    double M = pars->M;
    double U[3];
    frame_U(r, M, U);

    double W = AL * U[0];
    double Urd = G40R*U[0] + G4RR*U[1] + G4RP*U[2];
    double Upd = G40P*U[0] + G4RP*U[1] + G4PP*U[2];

    double s2 = (IG3RR*Sr*Sr + 2*IG3RP*Sr*Sp + IG3PP*Sp*Sp)/(D*D);
    double su = (IG3RR*Sr*Urd + IG3RP*(Sr*Upd+Sp*Urd) + IG3PP*Sp*Upd) / D;
    double e = (tau/D + su + 1.0)/W;

    double c0 = e*e - s2;
    double c1 = (gam-1.0)*(gam-1.0)/(gam*gam) - 2*s2/gam + e*e;
    double c2 = -s2/(gam*gam);
    double c3 = -2*e*(gam-1.0)/gam;

    double u2 = IG3RR*prim[VX1]*prim[VX1] + IG3PP*prim[VX2]*prim[VX2]
                + 2*IG3RP*prim[VX1]*prim[VX2];
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
    double rho = D/(J3 * w);
    double hmo = (w*e - w*w) / (u2 + 1.0/gam);

    double P = (gam-1.0)/gam * rho * hmo;
    double ur = Sr/(D*(hmo+1.0));
    double up = Sp/(D*(hmo+1.0));

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = ur;
    prim[VX2] = up;
}

void flux_rel_metric(double *prim, double *F, double r, struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ur = prim[VX1];
    double up = prim[VX2];
    double gam = pars->gammalaw;
    double M = pars->M;
    double U[3];

    frame_U(r, M, U);
    double w = sqrt(1.0 + IG3RR*ur*ur + IG3PP*up*up + 2*IG3RP*ur*up);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double rhoh = rho + gam/(gam-1.0)*P;
    double uU = u0*U[0] + ur*U[1] + up*U[2];

    double urr = IG40R*u0 + IG4RR*ur + IG4RP*up;

    F[DDD] = J4 * rho*urr;
    F[TAU] = J4 * (-rhoh*uU*urr - U[1]*P - rho*urr);
    F[SX1] = J4 * (rhoh*urr*ur + P);
    F[SX2] = J4 * rhoh*urr*up;
}

void add_source_rel_metric(double *prim, double *cons, double r, double dVdt, 
                struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ur = prim[VX1];
    double up = prim[VX2];
    double gam = pars->gammalaw;
    double M = pars->M;

    double U[3], dU[3];
    frame_U(r, M, U);
    frame_dU(r, M, dU);
    double w = sqrt(1.0 + IG3RR*ur*ur + IG3PP*up*up + 2*IG3RP*ur*up);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double rhoh = rho + gam/(gam-1.0)*P;

    double u00 = IG400*u0 + IG40R*ur + IG40P*up;
    double urr = IG40R*u0 + IG4RR*ur + IG4RP*up;
    double upp = IG40P*u0 + IG4RP*ur + IG4PP*up;

    double Sr = 0.5 * ((rhoh*u00*u00 + IG400*P)*DG400
                     + (rhoh*urr*urr + IG4RR*P)*DG4RR
                     + (rhoh*upp*upp + IG4PP*P)*DG4PP
                     + 2*(rhoh*u00*urr + IG40R*P)*DG40R
                     + 2*(rhoh*u00*upp + IG40P*P)*DG40P
                     + 2*(rhoh*urr*upp + IG4RP*P)*DG4RP);

    double Stau = -(rhoh * urr * (u0*dU[0] + ur*dU[1] + up*dU[2]) + P * dU[1]);

    cons[SX1] += Sr * J4 * dVdt;
    cons[TAU] += (Stau - U[1]*Sr) * J4 * dVdt;
}

void wave_speeds_rel_metric(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double r, struct parList *pars)
{
    double rho1 = prim1[RHO];
    double P1 = prim1[PPP];
    double ur1 = prim1[VX1];
    double up1 = prim1[VX2];
    double rho2 = prim2[RHO];
    double P2 = prim2[PPP];
    double ur2 = prim2[VX1];
    double up2 = prim2[VX2];
    double gam = pars->gammalaw;
    double M = pars->M;

    double rhoh1 = rho1 + gam/(gam-1.0)*P1;
    double rhoh2 = rho2 + gam/(gam-1.0)*P2;

    double cs21 = gam*P1/rhoh1;
    double cs22 = gam*P2/rhoh2;

    double w1 = sqrt(1.0 + IG3RR*ur1*ur1 + 2*IG3RP*ur1*up1 + IG3PP*up1*up1);
    double w2 = sqrt(1.0 + IG3RR*ur2*ur2 + 2*IG3RP*ur2*up2 + IG3PP*up2*up2);
    double u01 = (w1/AL - IG40R*ur1 - IG40P*up1) / IG400;
    double u02 = (w2/AL - IG40R*ur2 - IG40P*up2) / IG400;

    double urr1 = IG40R*u01 + IG4RR*ur1 + IG4RP*up1;
    double urr2 = IG40R*u02 + IG4RR*ur2 + IG4RP*up2;
    double vr1 = urr1*AL/w1;
    double vr2 = urr2*AL/w2;

    double br = - IG40R / IG400;

    double sig1 = cs21 / (w1*w1*(1.0-cs21));
    double sig2 = cs22 / (w2*w2*(1.0-cs22));

    double dv1 = sqrt(sig1*(1.0+sig1)*AL*AL*IG3RR - sig1*(vr1+br)*(vr1+br));
    double dv2 = sqrt(sig2*(1.0+sig2)*AL*AL*IG3RR - sig2*(vr2+br)*(vr2+br));

    double sL1 = (vr1 - sig1*br - dv1) / (1.0 + sig1);
    double sL2 = (vr2 - sig2*br - dv2) / (1.0 + sig2);
    double sR1 = (vr1 - sig1*br + dv1) / (1.0 + sig1);
    double sR2 = (vr2 - sig2*br + dv2) / (1.0 + sig2);

    *sL = sL1<sL2 ? sL1 : sL2;
    *sR = sR1>sR2 ? sR1 : sR2;

    //TODO: Use proper HLLC wavespeed.
    *sC = 0.0;
}

double mindt_rel_metric(double *prim, double r, double dx, 
                        struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double ur = prim[VX1];
    double up = prim[VX2];
    double gam = pars->gammalaw;
    double M = pars->M;

    double rhoh = rho + gam/(gam-1.0)*P;
    double cs2 = gam*P/rhoh;

    double u2 = IG3RR*ur*ur + 2*IG3RP*ur*up + IG3PP*up*up;
    double w = sqrt(1.0 + u2);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double urr = IG40R*u0 + IG4RR*ur + IG4RP*up;
    double vr = urr*AL/w;

    double br = - IG40R / IG400;
    double sig = cs2 / (w*w * (1.0-cs2));
    double dv = sqrt(sig*(1.0+sig)*AL*AL*IG3RR - sig*(vr+br)*(vr+br));

    double s;
    if(vr > 0)
        s = (vr - sig*br + dv) / (1.0 + sig);
    else
        s = (vr - sig*br - dv) / (1.0 + sig);

    double dt = dx / fabs(s);

    return dt;
}

void frame_normal_U(double r, double M, double U[])
{
    U[0] = sqrt(1.0+2.0*M/r);
    U[1] = -2.0*M / sqrt(r*(2.0*M+r));
    //U[0] = 1.0/sqrt(1.0-2*M/r);
    //U[1] = 0.0;
    //U[0] = 1.0;
    //U[1] = 0.0;
    U[2] = 0.0;
}

void frame_normal_dU(double r, double M, double dU[])
{
    double x = r*r + 2.0*M*r;
    dU[0] = -M/(r*sqrt(x));
    dU[1] = 2*M*(M+r) / sqrt(x*x*x);
    //dU[0] = -M / (r*r*sqrt((1-2*M/r)*(1-2*M/r)*(1-2*M/r)));
    //dU[1] = 0.0;
    //dU[0] = 0.0;
    //dU[1] = 0.0;
    dU[2] = 0.0;
}

void frame_geo_schw_U(double r, double M, double U[])
{
    if(r > 6*M)
    {
        U[0] = sqrt(r/(r-3.0*M));
        U[1] = 0.0;
        U[2] = sqrt(M / (r*r*(r-3.0*M)));
    }
    else
    {
        double x = 6*M/r - 1.0;
        U[0] = 2.0 * (sqrt(2.0)*r - M*sqrt(x*x*x)) / (3.0*(r-2.0*M));
        U[1] = -sqrt(x*x*x) / 3.0;
        U[2] = 2.0*sqrt(3.0) * M/(r*r);
    }
}

void frame_geo_schw_dU(double r, double M, double dU[])
{
    if(r > 6*M)
    {
        double x = 1.0 - 3.0*M/r;
        dU[0] = -1.5*M / (r*r * sqrt(x*x*x));
        dU[1] = 0.0;
        dU[2] = -1.5 * (r-2.0*M) * sqrt(M/(r*x*x*x)) / (r*r*r);
    }
    else
    {
        double x = 6*M/r - 1.0;
        double y = M/r;
        dU[0] = -2.0*M * (sqrt(x)*(18.0*y*y-15.0*y+1.0) + 2.0*sqrt(2.0))
                / (3.0*(r-2.0*M)*(r-2.0*M));
        dU[1] = 3.0 * sqrt(x) * M/(r*r);
        dU[2] = -4.0 * sqrt(3.0) * M / (r*r*r);
    }
}
