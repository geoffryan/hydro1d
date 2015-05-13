#include <stdlib.h>
#include <stdio.h>
#include "par.h"
#include "hydro.h"

int set_hydro(struct parList *pars)
{
    int err = 0;
    int choice = pars->hydro;

    if(choice == 0)
    {
        prim2cons = &prim2cons_newt_cart;
        cons2prim = &cons2prim_newt_cart;
        flux = &flux_newt_cart;
        add_source = &add_source_newt_cart;
        wave_speeds = &wave_speeds_newt_cart;
        mindt = &mindt_newt_cart;
    }
    else if(choice == 1)
    {
        prim2cons = &prim2cons_newt_cyl;
        cons2prim = &cons2prim_newt_cyl;
        flux = &flux_newt_cyl;
        add_source = &add_source_newt_cyl;
        wave_speeds = &wave_speeds_newt_cyl;
        mindt = &mindt_newt_cyl;
    }
    else
    {
        err = 1;
        printf("ERROR - Invalid choice for hydro: %d\n", choice);
    }

    return err;
}