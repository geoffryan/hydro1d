#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grid.h"
#include "geom.h"
#include "par.h"
#include "initial.h"
#include "boundary.h"
#include "hydro.h"
#include "hydro/metric/schw_ks.h"

int set_bc(struct parList *pars)
{
    int err = 0;
    int choiceInner = pars->bcInner;
    int choiceOuter = pars->bcOuter;

    if(choiceInner == 0)
        bc_inner = &bc_fixed_inner;
    else if(choiceInner == 1)
        bc_inner = &bc_outflow_inner;
    else if(choiceInner == 2)
        bc_inner = &bc_geodesic_inner;
    else
    {
        err++;
        printf("ERROR - Invalid choice for inner BC: %d\n", choiceInner);
    }

    if(choiceOuter == 0)
        bc_outer = &bc_fixed_outer;
    else if(choiceOuter == 1)
        bc_outer = &bc_outflow_outer;
    else if(choiceOuter == 2)
        bc_outer = &bc_geodesic_outer;
    else
    {
        err++;
        printf("ERROR - Invalid choice for outer BC: %d\n", choiceOuter);
    }

    return err;
}

void bc_fixed_inner(struct grid *g, struct parList *pars)
{
    int i;
    int nq = g->nq;
    int ng = g->ng;

    for(i=0; i<ng; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}

void bc_fixed_outer(struct grid *g, struct parList *pars)
{
    int i;
    int nx = g->nx;
    int nq = g->nq;
    int ng = g->ng;

    for(i=nx-ng; i<nx; i++)
    {
        double xm = g->x[i];
        double xp = g->x[i+1];
        double x = 0.5*(xm+xp);
        initial_value(&(g->prim[nq*i]), x, pars);
    }
}

void bc_outflow_inner(struct grid *g, struct parList *pars)
{
    int q,i;
    int nq = g->nq;
    int ng = g->ng;

    for(i=ng-1; i>=0; i--)
        for(q=0; q<nq; q++)
            g->prim[i*nq+q] = g->prim[(i+1)*nq+q];
}

void bc_outflow_outer(struct grid *g, struct parList *pars)
{
    int q, i;
    int nq = g->nq;
    int nx = g->nx;
    int ng = g->ng;

    for(i=nx-ng; i<nx; i++)
        for(q=0; q<nq; q++)
            g->prim[i*nq+q] = g->prim[(i-1)*nq+q];
}

void bc_geodesic_inner(struct grid *g, struct parList *pars)
{
    int q, i;
    int nq = g->nq;
    int ng = g->ng;

    double M = pars->M;
    double gam = pars->gammalaw;
    double r = CM(g->x[ng], g->x[ng+1]);

    double rho = g->prim[nq*ng+RHO];
    double p = g->prim[nq*ng+PPP];
    double ur = g->prim[nq*ng+VX1];
    double up = g->prim[nq*ng+VX2];
    double w = sqrt(1.0 + IG3RR*ur*ur + 2*IG3RP*ur*up + IG3PP*up*up);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double urr = IG40R*u0 + IG4RR*ur + IG4RP*up;

    double Mdot = - r * rho * urr;
    double K = p / pow(rho, gam); 
    double l = up;
    double e = u0;

    for(i=0; i<ng; i++)
    {
        r = CM(g->x[i], g->x[i+1]);
        //double U[3];
        //frame_U(r, M, U);
        //double Ur = G40R * U[0] + G4RR * U[1] + G4RP * U[2];
        //double Up = G40P * U[0] + G4RP * U[1] + G4PP * U[2];
        //rho = -Mdot / (r * U[1]);
        
        double hb = IG40R*e + IG4RP*l;
        double a = IG4RR;
        double c = 1.0 + IG400*e*e + IG4PP*l*l + 2*IG40P*e*l;
        ur = (-hb - sqrt(hb*hb-a*c)) / a;
        urr = IG40R*e + IG4RR*ur + IG4RP*l;
        rho = -Mdot / (r * urr);
        
        p = K * pow(rho, gam);

        g->prim[i*nq + RHO] = rho;
        g->prim[i*nq + PPP] = p;
        g->prim[i*nq + VX1] = ur;
        g->prim[i*nq + VX2] = l; //Up

        for(q=4; q<nq; q++)
            g->prim[i*nq+q] = g->prim[ng*nq+q];
    }
}

void bc_geodesic_outer(struct grid *g, struct parList *pars)
{
    int q, i;
    int nq = g->nq;
    int nx = g->nx;
    int ng = g->ng;

    double M = pars->M;
    double gam = pars->gammalaw;
    double r = CM(g->x[nx-ng-1], g->x[nx-ng]);

    double rho = g->prim[nq*(nx-ng-1)+RHO];
    double p = g->prim[nq*(nx-ng-1)+PPP];
    double ur = g->prim[nq*(nx-ng-1)+VX1];
    double up = g->prim[nq*(nx-ng-1)+VX2];
    double w = sqrt(1.0 + IG3RR*ur*ur + 2*IG3RP*ur*up + IG3PP*up*up);
    double u0 = (w/AL - IG40R*ur - IG40P*up) / IG400;
    double urr = IG40R*u0 + IG4RR*ur + IG4RP*up;

    double Mdot = - r * rho * urr;
    double K = p / pow(rho, gam); 

    for(i=nx-ng; i<nx; i++)
    {
        double U[3];
        r = CM(g->x[i], g->x[i+1]);
        frame_U(r, M, U);

        double Ur = G40R * U[0] + G4RR * U[1] + G4RP * U[2];
        double Up = G40P * U[0] + G4RP * U[1] + G4PP * U[2];
        rho = -Mdot / (r * U[1]);
        p = K * pow(rho, gam);

        g->prim[i*nq + RHO] = rho;
        g->prim[i*nq + PPP] = p;
        g->prim[i*nq + VX1] = Ur;
        g->prim[i*nq + VX2] = Up;

        for(q=4; q<nq; q++)
            g->prim[i*nq+q] = g->prim[(i-1)*nq+q];
    }
}
