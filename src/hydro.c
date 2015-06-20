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
    else if(choice == 2)
    {
        prim2cons = &prim2cons_newt_sph;
        cons2prim = &cons2prim_newt_sph;
        flux = &flux_newt_sph;
        add_source = &add_source_newt_sph;
        wave_speeds = &wave_speeds_newt_sph;
        mindt = &mindt_newt_sph;
    }
    else if(choice == 3)
    {
        prim2cons = &prim2cons_rel_cart;
        cons2prim = &cons2prim_rel_cart;
        flux = &flux_rel_cart;
        add_source = &add_source_rel_cart;
        wave_speeds = &wave_speeds_rel_cart;
        mindt = &mindt_rel_cart;
    }
    else
    {
        err = 1;
        printf("ERROR - Invalid choice for hydro: %d\n", choice);
    }

    return err;
}
