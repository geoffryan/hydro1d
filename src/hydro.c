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
        grid_V = &grid_V_newt_cart;
    }
    else if(choice == 1)
    {
        prim2cons = &prim2cons_newt_cyl;
        cons2prim = &cons2prim_newt_cyl;
        flux = &flux_newt_cyl;
        add_source = &add_source_newt_cyl;
        wave_speeds = &wave_speeds_newt_cyl;
        mindt = &mindt_newt_cyl;
        grid_V = &grid_V_newt_cyl;
    }
    else if(choice == 2)
    {
        prim2cons = &prim2cons_newt_sph;
        cons2prim = &cons2prim_newt_sph;
        flux = &flux_newt_sph;
        add_source = &add_source_newt_sph;
        wave_speeds = &wave_speeds_newt_sph;
        mindt = &mindt_newt_sph;
        grid_V = &grid_V_newt_sph;
    }
    else if(choice == 3)
    {
        prim2cons = &prim2cons_rel_cart;
        cons2prim = &cons2prim_rel_cart;
        flux = &flux_rel_cart;
        add_source = &add_source_rel_cart;
        wave_speeds = &wave_speeds_rel_cart;
        mindt = &mindt_rel_cart;
        grid_V = &grid_V_rel_cart;
    }
    else if(choice == 4)
    {
        prim2cons = &prim2cons_rel_metric;
        cons2prim = &cons2prim_rel_metric;
        flux = &flux_rel_metric;
        add_source = &add_source_rel_metric;
        wave_speeds = &wave_speeds_rel_metric;
        mindt = &mindt_rel_metric;
        grid_V = &grid_V_rel_metric;

        int frame = pars->frame;
        if(frame == 0)
        {
            frame_U = &frame_normal_U;
            frame_dU = &frame_normal_dU;
        }
        else if(frame == 1)
        {
            frame_U = &frame_geo_schw_U;
            frame_dU = &frame_geo_schw_dU;
        }

        else
        {
            err++;
            printf("ERROR - Invalid choice for frame: %d\n", frame);
        }
    }
    else
    {
        err++;
        printf("ERROR - Invalid choice for hydro: %d\n", choice);
    }

    return err;
}
