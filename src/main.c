#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "par.h"
#include "boundary.h"
#include "hydro.h"
#include "initial.h"
#include "riemann.h"
#include "timestep.h"

void getTheHellOuttaHere(struct grid *g);

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        printf("usage: hydro1c <parameter file>\n");
        return 0;
    }

    printf("\nWelcome to Hydro-1D!\n%s\n\n", VERSION);

    int err = 0;
    struct parList pars = PAR_DEFAULT;
    struct grid grid = GRID_DEFAULT;

    read_pars(&pars, argv[1]);
    print_pars(&pars);

    err += set_initial(&pars);
    err += set_reconstruction(&pars);
    err += set_riemann_solver(&pars);
    err += set_timestep(&pars);
    err += set_bc(&pars);
    err += set_hydro(&pars);

    if(err)
    {
        printf("Error during setup.\n");
        exit(0);
    }

    make_grid(&grid, &pars);
    initialize_grid(&grid, &pars);
    calc_cons(&grid, &pars);

    char filename[128];
    sprintf(filename, "grid_%04d.txt", 0);
    print_grid(&grid, filename);

    int i = 1;
    while(grid.t < pars.tmax)
    {
        double dt = get_dt(&grid, &pars);
        dt = grid.t+dt > pars.tmax ? pars.tmax-grid.t : dt;

        printf("t: %.6e dt: %.6e\n", grid.t, dt);

        step_fe(&grid, dt, &pars);

        sprintf(filename, "grid_%04d.txt", i);
        print_grid(&grid, filename);
        i++;
    }

    free_grid(&grid);

    return 0;
}

void getTheHellOuttaHere(struct grid *g)
{
    printf("Aborting...\n");
    print_grid(g, "abort_dump.txt");
    free_grid(g);
    exit(0);
}
